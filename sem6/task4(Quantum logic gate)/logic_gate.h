#include <complex>
#include "mpi.h"
#include <omp.h>


typedef std::complex<double> complexd;

complexd *qubit_transform(complexd *a, int  n, complexd u[2][2], int  k, long powproc, int rank) {
    long vec_length = 1 << (n - powproc);
    long dist = 1 << (n - k);
    long start = vec_length * rank; // 2*n*rank/size (part number of the whole vector)
    complexd *b = new complexd[vec_length];
    if (dist < vec_length) {
        #pragma omp parallel for
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
        #pragma omp parallel for
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




complexd *TwoQubitsEvolution(complexd *a, int n, int q1, int q2, complexd u[4][4], long powproc, int rank)
{
    //n - число кубитов
    //q1, q2 - номера кубитов, над которыми производится преобразование
    int shift1 = n - q1;
    int shift2 = n - q2;

    //Все биты нулевые, за исключением соответсвующего номеру первого изменяемого кубита
    int pow2q1 = 1 << (shift1);
    //Все биты нулевые, за исключением соответсвующего номеру второгоизменяемого кубита
    int pow2q2 = 1 << (shift2);

    long vec_length = 1 << (n - powproc);
    long start = vec_length * rank; 
    complexd *b = new complexd[vec_length];

    for (int i = 0; i < vec_length; i++)
    {
        //Установка изменяемых битов во все возможные позиции
        int i00 = i & ~pow2q1 & ~pow2q2;
        int i01 = i & ~pow2q1 | pow2q2;
        int i10 = (i | pow2q1) & ~pow2q2;
        int i11 = i | pow2q1 | pow2q2;
        

        
        //Получение значений изменяемых битов
        int iq1 = (i & pow2q1) >> shift1;
        int iq2 = (i & pow2q2) >> shift2;

        //Индекс в матрице
        int iq = (iq1 << 1 ) + iq2;
        b[i] = u[iq][(0 << 1) + 0] * a[i00] + u[iq][( 0 << 1 ) + 1] * a[i01] + u[iq][( 1 << 1 ) + 0] * a[i10] + u[iq][( 1 << 1 ) + 1] * a[i11];
    }
    return b;
}

complexd *C_not(complexd *a, int  n, int  q1, int q2, long powproc, int rank) {
    complexd u[4][4];
    u[0][0] = u[1][1] = u[2][3] = u[3][2] = 1;
    u[0][1] = u[0][2] = u[0][3] = u[1][0] = u[1][2] = u[1][3] = u[2][0] = u[2][1] = u[2][2] = u[3][0] = u[3][1] = u[3][3] = 0;
    return TwoQubitsEvolution(a, n, q1, q2,u, powproc, rank);
}