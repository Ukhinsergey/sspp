#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "./logic_gate.h"

typedef std::complex<double> complexd;
const double MAXD = 232323;
using namespace std;







int main(int argc, char **argv) {
    if (argc != 7) {
        cout << "input: n q1 q2 mode(1-file \"in.bin\" , 2-random)";
        cout << " numtreads, out.bin" << endl;
        return 0;
    }
    srand(time(NULL));
    MPI_Init(&argc, &argv);
    int n, q1, q2, numtreads;
    int mode;
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &q1);
    sscanf(argv[3], "%d", &q2);
    sscanf(argv[4], "%d", &mode);
    sscanf(argv[5], "%d", &numtreads);
    omp_set_num_threads(numtreads);

    complexd **u = new complexd *[2];
    for (int i = 0; i < 2; ++i) {
        u[i] = new complexd[2];
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            if ( i != 1 || j != 1 ) {
                u[i][j] = 1.0 / sqrt(2);
            } else {
                u[i][j] = - 1.0/sqrt(2);
            }
        }
    }


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
                a[i] = i + start;
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
        MPI_File_open(MPI_COMM_WORLD, "in.bin", MPI_MODE_RDONLY,
        MPI_INFO_NULL, &file);
        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype,
        "native", MPI_INFO_NULL);
        MPI_File_read_all(file, a, vec_length * 2, MPI_DOUBLE, &status);
        MPI_File_close(&file);
        double dlina = 0;
        #pragma omp parallel for reduction(+ : dlina)
        for (int i = 0 ; i < vec_length; ++i) {
            dlina += norm(a[i]);
        }
        double temp;
        MPI_Allreduce(&dlina, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dlina = sqrt(temp);
        #pragma omp parallel for
        for (int i = 0; i < vec_length; ++i) {
            a[i] = a[i] / dlina;
        }
    }


    /*
    for(int i = 0; i < vec_length; ++i) {
            cout << a[i] << endl;
    }
    */

    /*
    double dlina = 0;
    for(int i = 0 ; i < vec_length; ++i) {
                //a[i] = complexd((double)rand()/RAND_MAX * MAXD,
                //(double) rand()/RAND_MAX * MAXD);   
        dlina += norm(b[i]);
    } 
    cout << dlina << endl;
    */
    complexd *b = N_hadamar(a, n, powproc, rank);


    int s = 1;
    int p = 0;
    MPI_Type_create_subarray(1,  &s, &s, &p, MPI_ORDER_C,
    MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    int offset = sizeof(int) + vec_length * 2 * rank * sizeof(double);
    if ( rank == 0 ) {
        MPI_File_open(MPI_COMM_SELF, argv[6],
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        MPI_File_write(file, &n, 1, MPI_INT, &status);
        MPI_File_close(&file);
    }
    MPI_File_open(MPI_COMM_WORLD, argv[6],
    MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
    MPI_File_set_view(file, offset, MPI_DOUBLE,
    filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(file, b, vec_length * 2, MPI_DOUBLE, &status);
    MPI_File_close(&file);

    for (int i = 0 ; i < 2; ++i) {
        delete [] u[i];
    }
    delete []b;
    delete []a;
    delete []u;
    MPI_Finalize();
    return 0;
}
