#include <iostream>
#include <complex>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "mpi.h"
#include <omp.h>

typedef std::complex<double> complexd;
const double MAXD = 232323;
using namespace std;






complexd *qubit_transform(complexd *a, int  n, complexd u[2][2], int  k, long powproc,int rank) {
    long vec_length = 1 << (n - powproc);
    long dist = 1 << (n - k);
    long start = vec_length * rank; // 2*n*rank/size (part number of the whole vector)
    complexd *b = new complexd[vec_length];
    if (dist < vec_length) {
        //#pragma omp parallel for
        for(int i = 0; i < vec_length; ++i) {
            b[i] = u[((i + start) & dist) >> (n - k)][0] * a[(((i + start) | dist) ^ dist) - start] + u[((i + start) & dist) >> (n - k)][1] * a[ ((i + start) | dist) - start];
        }
    } else {
        int needrank;
        if ((start & dist) == 0) {                           //look at the bit that needs to be changed
            needrank = (start | dist) / vec_length;
        } else {
            needrank = (start & ~dist) / vec_length;
        }
        complexd *temp = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, temp, vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        complexd *vec0, *vec1;
        if (rank < needrank) { // bit needs to be changed in a == 0
            vec0 = a;
            vec1 = temp;
        } else {               // bit needs to be changed in temp == 0
            vec0 = temp;
            vec1 = a;
        }
        //#pragma omp parallel for
        for(int i = 0; i < vec_length; ++i) {
            b[i] = u[((i + start) & dist) >> (n - k)][0] * vec0[i] + u[((i + start) & dist) >> (n - k)][1] * vec1[i]; // 
        }
        delete []temp;
    }
    return b;
}

complexd *Hadamar(complexd *a, int  n, int  k, long powproc, int rank) {
    complexd u[2][2];
    u[0][0] = u[0][1] = u[1][0] = 1.0 / sqrt(2.0);
    u[1][1] = -u[0][0];
    return qubit_transform(a, n, u, k, powproc, rank);
}

complexd *N_hadamar(complexd *a, int  n, long powproc, int rank) {
    complexd u[2][2];
    u[0][0] = u[0][1] = u[1][0] = 1.0 / sqrt(2.0);
    u[1][1] = -u[0][0];
    complexd *temp = qubit_transform(a, n, u, 1, powproc, rank);
    complexd *b;
    for(int i = 2 ; i <= n; ++i) {
        b = qubit_transform(temp, n, u, i, powproc, rank);
        delete[] temp;
        temp = b;
    }
    return b;
}

complexd *Not(complexd *a, int  n, int  k, long powproc, int rank) {
    complexd u[2][2];
    u[0][0] = u[1][1] = 0;
    u[0][1] = u[1][0] = 1;
    return qubit_transform(a, n, u, k, powproc, rank);
}

complexd *Rw(complexd *a, int  n, int  k, long powproc, int rank, double phi) {
    complexd u[2][2];
    u[0][0] = 1;
    u[0][1] = u[1][0] = 0;
    u[1][1] = exp(complexd(0,1) * phi);
    return qubit_transform(a, n, u, k, powproc, rank);
}


int main(int argc, char **argv) {
    if( argc != 6) {
        cout << "input: n q1 q2 mode(1-file \"in.bin\" , 2-random) numtreads" << endl;
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

    complexd **u= new complexd *[2];
    for(int i = 0; i < 2; ++i) {
        u[i] = new complexd[2];
    }

    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < 2; ++j) {
            if ( i != 1 || j != 1) {
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

    for(int i = 0 ; i < 10; ++i) {
        if (1 << powproc == size) {
            break;
        }
        ++powproc;
    }

    long vec_length = 1 << (n - powproc); // = 2^n/size (size = 2 ^ powproc)
    complexd *a = new complexd[vec_length]; 
    
    double timetemp;
    if(rank == 0) {
        timetemp = MPI_Wtime ();
    }
    MPI_Bcast(&timetemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (mode == 2) {
        double dlina = 0;
        long start = vec_length * rank;
        //#pragma omp parallel
        {
            unsigned int seed = timetemp + omp_get_num_threads() * rank + omp_get_thread_num();
            //#pragma omp for reduction(+ : dlina)
            for(int i = 0 ; i < vec_length; ++i) {
                //a[i] = complexd((double)rand()/RAND_MAX * MAXD, (double) rand()/RAND_MAX * MAXD);   
                a[i] = i + start;
                dlina += norm(a[i]);
            } 
        }
        
        double temp;
        MPI_Allreduce(&dlina, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dlina = sqrt(temp); 
        //#pragma omp parallel for
        for(int i = 0; i < vec_length; ++i) {
            a[i] = a[i] / dlina;
        }

    } else {
        int s = 1;
        int p = 0;
        MPI_Type_create_subarray(1,  &s, &s, &p, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);
        int offset = sizeof(int) + vec_length * 2 * rank * sizeof(double);
        MPI_File_open(MPI_COMM_WORLD, "in.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
        MPI_File_read_all(file, a, vec_length * 2, MPI_DOUBLE, &status);
        MPI_File_close(&file);
        double dlina = 0;
        //#pragma omp parallel for reduction(+ : dlina)
        for(int i = 0 ; i < vec_length; ++i) {
            dlina += norm(a[i]);
        } 
        double temp;
        MPI_Allreduce(&dlina, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dlina = sqrt(temp); 
        //#pragma omp parallel for
        for(int i = 0; i < vec_length; ++i) {
            a[i] = a[i] / dlina;
        }
    }
    /*
    complexd d = 0;
    for(int i = 0; i < vec_length; ++i) {
            d += abs(a[i] * a[i]);
    }
    cout << d << endl;
    */
    //end_time1 = MPI_Wtime();

    /*double sum = 0;
    for(int i = 0; i < vec_length; ++i) {
            cout << a[i] << endl;
    }
    */

    /*
    double dlina = 0;
    for(int i = 0 ; i < vec_length; ++i) {
                //a[i] = complexd((double)rand()/RAND_MAX * MAXD, (double) rand()/RAND_MAX * MAXD);   
        dlina += norm(b[i]);
    } 
    cout << dlina << endl;
    */
    complexd *b = Rw(a, n, q1, powproc, rank, 3.1415);


    int s = 1;
    int p = 0;
    MPI_Type_create_subarray(1,  &s, &s, &p, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    int offset = sizeof(int) + vec_length * 2 * rank * sizeof(double);
    if( rank == 0) {
        MPI_File_open(MPI_COMM_SELF, "out.bin", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        MPI_File_write(file, &n, 1, MPI_INT, &status);
        MPI_File_close(&file);
    }
    MPI_File_open(MPI_COMM_WORLD, "out.bin", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
    MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(file, b, vec_length * 2, MPI_DOUBLE, &status);
    MPI_File_close(&file);
    


    
    

    for(int i = 0 ; i < 2; ++i) {
        delete [] u[i];
    }
    delete []b;
    delete []a;
    delete []u;
    MPI_Finalize();
    return 0;
}