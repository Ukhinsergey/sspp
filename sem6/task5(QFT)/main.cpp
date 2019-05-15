#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "./logic_gate.h"



typedef std::complex<double> complexd;
using namespace std;

int main(int argc, char **argv) {
    if (argc != 6) {
        cout << "input: n, mode(1-file \"in.bin\" , 2-random)";
        cout << "numtreads, in.bin out.bin" << endl;
        return 0;
    }

    srand(time(NULL));
    MPI_Init(&argc, &argv);
    int n, numtreads;
    int mode;
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &mode);
    sscanf(argv[3], "%d", &numtreads);
    char * outfilename = argv[5];
    char * infilename = argv[4];

    omp_set_num_threads(numtreads);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    MPI_Datatype filetype;
    MPI_File file;

    int powproc = 0;

    for (int i = 0 ; i < 10; ++i) {
        if (1 << powproc == size) {
            break;
        }
        ++powproc;
    }

    long vec_length = 1 << (n - powproc);  // = 2^n/size (size = 2 ^ powproc)
    complexd *a = new complexd[vec_length];
    complexd *b = new complexd[vec_length];

    double timetemp;
    if (rank == 0) {
        timetemp = MPI_Wtime();
    }
    MPI_Bcast(&timetemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mode == 2) {
        double dlina = 0;
        long start = vec_length * rank;
#pragma omp parallel
        {
            unsigned int seed = timetemp + omp_get_num_threads() * rank +
            omp_get_thread_num();
#pragma omp for reduction(+ : dlina)
            for (int i = 0 ; i < vec_length; ++i) {
                // a[i] = complexd((double)rand()/RAND_MAX * MAXD,
                // (double) rand()/RAND_MAX * MAXD);
                a[i] = complexd(i + start, i + start);
                dlina += norm(a[i]);
            }
        }

        double temp;
        MPI_Allreduce(&dlina, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dlina = sqrt(temp);
#pragma omp parallel for
        for (int i = 0; i < vec_length; ++i) {
            a[i] = a[i] / dlina;
        }
    } else {
        int s = 1;
        int p = 0;
        MPI_Type_create_subarray(1,  &s, &s, &p, MPI_ORDER_C,
        MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);
        int offset = sizeof(int) + vec_length * 2 * rank * sizeof(double);
        MPI_File_open(MPI_COMM_WORLD, infilename, MPI_MODE_RDONLY,
        MPI_INFO_NULL, &file);
        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype,
        "native", MPI_INFO_NULL);
        MPI_File_read_all(file, a, vec_length * 2, MPI_DOUBLE, &status);
        MPI_File_close(&file);
        double dlina = 0;
    }

    double start_time, end_time;


    start_time = MPI_Wtime();


    for(int i = 1; i <= n; ++i) {
        Hadamar(a, b, n, i, powproc, rank);
        complexd * swap = a;
        a = b;
        b = swap;
        int m = 2;
        for(int j = i + 1; j <= n; ++j) {
            double phi = 2 * M_PI / (1 << m);
            C_rw(a, b, n, j, i, powproc, rank, phi);
            complexd * swap = a;
        	a = b;
        	b = swap;
            m++;
        }
    }
    end_time = MPI_Wtime();
    double time = end_time - start_time;
    double maxtime;
    MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        cout << " maxtime: " << maxtime << endl;
    }

    if ( rank == 0 ) {
        MPI_File_open(MPI_COMM_SELF, outfilename,
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        MPI_File_write(file, &n, 1, MPI_INT, &status);
        MPI_File_close(&file);
    }
    

    vector<int> bit_masks(n);
    bit_masks[n - 1] = 1;
    for (int i = n - 2; i >= 0; --i) {
        bit_masks[i] = bit_masks[i + 1] << 1;
    }
    std::vector<int> index_vec(vec_length);
    for (int i = 0 ; i < vec_length; ++i) {
        int new_index = 0;
        for (int j = 0; j < n; ++j) {
            new_index |= (((i + rank * vec_length) & bit_masks[j]) && 1) << j;
        }
        index_vec[i] = new_index;
    }


    int s = 1;
    int p = 0;
    MPI_Type_create_subarray(1,  &s, &s, &p, MPI_ORDER_C,
    MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    for(int i = 0; i < vec_length; ++i) {
        int offset = sizeof(int) + index_vec[i] * 2 * sizeof(double);
        MPI_File_open(MPI_COMM_WORLD, outfilename,
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        MPI_File_set_view(file, offset, MPI_DOUBLE,
        filetype, "native", MPI_INFO_NULL);
        MPI_File_write(file, a + i, 2, MPI_DOUBLE, &status);
        MPI_File_close(&file);
    }

    delete[] a;
    delete[] b;

    MPI_Finalize();
    return 0;
}