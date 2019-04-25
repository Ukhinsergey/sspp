#include <omp.h>
#include <mpi.h>
#include <complex>
#ifndef SEM6_TASK4_QUANTUM_LOGIC_GATE__LOGIC_GATE_H_
#define SEM6_TASK4_QUANTUM_LOGIC_GATE__LOGIC_GATE_H_
#endif   //  SEM6_TASK4_QUANTUM_LOGIC_GATE__LOGIC_GATE_H_


typedef std::complex<double> complexd;

complexd *qubit_transform(complexd *a, int  n, complexd u[2][2], int  k,
    int powproc, int rank) {

    long vec_length = 1 << (n - powproc);
    long dist = 1 << (n - k);
    long start = vec_length * rank;
    complexd *b = new complexd[vec_length];
    if (dist < vec_length) {
        #pragma omp parallel for
        for (int i = 0; i < vec_length; ++i) {
            b[i] = u[((i + start) & dist) >> (n - k)][0] *
            a[(((i + start) | dist) ^ dist) - start] +
            u[((i + start) & dist) >> (n - k)][1] * a[ ((i + start) | dist) -
                start];
        }
    } else {
        int needrank;
        if ((start & dist) == 0) {  //  look at the bit that needs to be changed
            needrank = (start | dist) / vec_length;
        } else {
            needrank = (start & ~dist) / vec_length;
        }
        complexd *temp = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, temp,
            vec_length * 2, MPI_DOUBLE, needrank,
            0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        complexd *vec0, *vec1;
        if (rank < needrank) {  //   bit needs to be changed in a == 0
            vec0 = a;
            vec1 = temp;
        } else {               //   bit needs to be changed in temp == 0
            vec0 = temp;
            vec1 = a;
        }
        #pragma omp parallel for
        for (int i = 0; i < vec_length; ++i) {
            b[i] = u[((i + start) & dist) >> (n - k)][1] * vec0[i] +
            u[((i + start) & dist) >> (n - k)][0] * vec1[i];
        }
        delete []temp;
    }
    return b;
}


complexd *TwoQubitsEvolution(complexd *a, int n, int q1, int q2,
    complexd u[4][4], int powproc, int rank) {
    //  n - число кубитов
    //  q1, q2 - номера кубитов, над которыми производится преобразование
    int shift1 = n - q1;
    int shift2 = n - q2;


    //  Все биты нулевые, за исключением
    //  соответсвующего номеру первого изменяемого кубита
    int pow2q1 = 1 << (shift1);


    //  Все биты нулевые, за исключением
    //  соответсвующего номеру второгоизменяемого кубита
    int pow2q2 = 1 << (shift2);

    long vec_length = 1 << (n - powproc);
    long start = vec_length * rank;
    complexd *b = new complexd[vec_length];

    if (pow2q1 < vec_length && pow2q2 < vec_length) {
        #pragma omp parallel for
        for (int i1 = 0; i1 < vec_length; i1++) {
            //  Установка изменяемых битов во все возможные позиции
            int i = i1 + start;
            int i00 = (i & ~pow2q1 & ~pow2q2) - start;
            int i01 = (i & ~pow2q1 | pow2q2) - start;
            int i10 = ((i | pow2q1) & ~pow2q2) - start;
            int i11 = (i | pow2q1 | pow2q2) - start;


            //  Получение значений изменяемых битов
            int iq1 = (i & pow2q1) >> shift1;
            int iq2 = (i & pow2q2) >> shift2;


            //  Индекс в матрице
            int iq = (iq1 << 1) + iq2;
            b[i1] = u[iq][(0 << 1) + 0] * a[i00] + u[iq][(0 << 1) + 1] *
            a[i01] + u[iq][(1 << 1) + 0] * a[i10] + u[iq][(1 << 1) + 1] *
            a[i11];
        }
        return b;
    } else if (pow2q1 < vec_length && pow2q2 >= vec_length) {
        int needrank;
        if ((start & pow2q2) == 0) {  // look at the bit that need to change
            needrank = (start | pow2q2) / vec_length;
        } else {
            needrank = (start & ~pow2q2) / vec_length;
        }
        complexd *temp = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, temp,
            vec_length * 2, MPI_DOUBLE, needrank, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        complexd *vecq2_0, *vecq2_1;
        if (rank < needrank) {  //   bit needs to be changed in a == 0
            vecq2_0 = a;
            vecq2_1 = temp;
        } else {               //   bit needs to be changed in temp == 0
            vecq2_0 = temp;
            vecq2_1 = a;
        }
        #pragma omp parallel for
        for (int i1 = 0; i1 < vec_length; i1++) {
            //  Установка изменяемых битов во все возможные позиции
            int i = i1 + start;
            int i00 = (i & ~pow2q1) - start;
            int i01 = (i & ~pow2q1) - start;
            int i10 = (i | pow2q1) - start;
            int i11 = (i | pow2q1) - start;


            //  Получение значений изменяемых битов
            int iq1 = (i & pow2q1) >> shift1;
            int iq2 = (i & pow2q2) >> shift2;


            //  Индекс в матрице
            int iq = (iq1 << 1) + iq2;
            b[i1] = u[iq][(0 << 1) + 0] * vecq2_0[i00] + u[iq][(0 << 1) + 1] *
            vecq2_1[i01] + u[iq][(1 << 1) + 0] * vecq2_0[i10] +
            u[iq][(1 << 1) + 1] * vecq2_1[i11];
        }
        delete []temp;
        return b;
    } else if (pow2q2 < vec_length && pow2q1 >= vec_length) {
        int needrank;
        if ((start & pow2q1) == 0) {  // look at the bit that need to change
            needrank = (start | pow2q1) / vec_length;
        } else {
            needrank = (start & ~pow2q1) / vec_length;
        }
        complexd *temp = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, temp,
            vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
        complexd *vecq1_0, *vecq1_1;
        if (rank < needrank) {  // bit needs to be changed in a == 0
            vecq1_0 = a;
            vecq1_1 = temp;
        } else {               //   bit needs to be changed in temp == 0
            vecq1_0 = temp;
            vecq1_1 = a;
        }
        #pragma omp parallel for
        for (int i1 = 0; i1 < vec_length; i1++) {
            //  Установка изменяемых битов во все возможные позиции
            int i = i1 + start;
            int i00 = (i & ~pow2q2) - start;
            int i01 = (i | pow2q2) - start;
            int i10 = (i & ~pow2q2) - start;
            int i11 = (i | pow2q2) - start;


            //  Получение значений изменяемых битов
            int iq1 = (i & pow2q1) >> shift1;
            int iq2 = (i & pow2q2) >> shift2;


            //  Индекс в матрице
            int iq = (iq1 << 1) + iq2;
            b[i1] = u[iq][(0 << 1) + 0] * vecq1_0[i00] + u[iq][(0 << 1) + 1] *
            vecq1_0[i01] + u[iq][(1 << 1) + 0] * vecq1_1[i10] +
            u[iq][(1 << 1) + 1] * vecq1_1[i11];
        }
        delete []temp;
        return b;
    } else if (pow2q1 >= vec_length && pow2q2 >= vec_length) {
        int needrank;
        int value00 = start & ~pow2q1 & ~pow2q2;
        needrank = value00 / vec_length;
        complexd *vec00 = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, vec00,
            vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        int value01 = start & ~pow2q1 | pow2q2;
        needrank = value01 / vec_length;
        complexd *vec01 = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, vec01,
            vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        int value10 = (start | pow2q1) & ~pow2q2;
        needrank = value10 / vec_length;
        complexd *vec10 = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, vec10,
            vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        int value11 = start | pow2q1 | pow2q2;
        needrank = value11 / vec_length;
        complexd *vec11 = new complexd[vec_length];
        MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, vec11,
            vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        #pragma omp parallel for
        for (int i1 = 0; i1 < vec_length; i1++) {
            int i = i1 + start;
            //  Получение значений изменяемых битов
            int iq1 = (i & pow2q1) >> shift1;
            int iq2 = (i & pow2q2) >> shift2;

            //  Индекс в матрице
            int iq = (iq1 << 1) + iq2;
            b[i1] = u[iq][(0 << 1) + 0] * vec00[i1] + u[iq][(0 << 1) + 1] *
            vec01[i1] + u[iq][(1 << 1) + 0] *
            vec10[i1] + u[iq][(1 << 1) + 1] * vec11[i1];
        }
        delete []vec00;
        delete []vec01;
        delete []vec10;
        delete []vec11;
        return b;
    }
}

complexd *Hadamar(complexd *a, int  n, int  k, int powproc, int rank) {
    complexd u[2][2];
    u[0][0] = u[0][1] = u[1][0] = 1.0 / sqrt(2.0);
    u[1][1] = -u[0][0];
    return qubit_transform(a, n, u, k, powproc, rank);
}

complexd *N_hadamar(complexd *a, int  n, int powproc, int rank) {
    complexd u[2][2];
    u[0][0] = u[0][1] = u[1][0] = 1.0 / sqrt(2.0);
    u[1][1] = -u[0][0];
    complexd *temp = qubit_transform(a, n, u, 1, powproc, rank);
    complexd *b;
    for (int i = 2 ; i <= n; ++i) {
        b = qubit_transform(temp, n, u, i, powproc, rank);
        delete[] temp;
        temp = b;
    }
    return b;
}

complexd *Not(complexd *a, int  n, int  k, int powproc, int rank) {
    complexd u[2][2];
    u[0][0] = u[1][1] = 0;
    u[0][1] = u[1][0] = 1;
    return qubit_transform(a, n, u, k, powproc, rank);
}

complexd *Rw(complexd *a, int  n, int  k, int powproc, int rank, double phi) {
    complexd u[2][2];
    u[0][0] = 1;
    u[0][1] = u[1][0] = 0;
    u[1][1] = exp(complexd(0, 1) * phi);
    return qubit_transform(a, n, u, k, powproc, rank);
}





complexd *C_not(complexd *a, int  n, int  q1, int q2, long powproc, int rank) {
    complexd u[4][4];
    u[0][0] = u[1][1] = u[2][3] = u[3][2] = 1;
    u[0][1] = u[0][2] = u[0][3] = u[1][0] = u[1][2] = u[1][3] = u[2][0] =
    u[2][1] = u[2][2] = u[3][0] = u[3][1] = u[3][3] = 0;
    return TwoQubitsEvolution(a, n, q1, q2, u, powproc, rank);
}


complexd *C_rw(complexd *a, int  n, int  q1, int q2, long powproc, int rank,
    double phi) {

    complexd u[4][4];
    u[0][0] = u[1][1] = u[2][2] = 1;
    u[0][1] = u[0][2] = u[0][3] = u[1][0] = u[1][2] = u[1][3] = u[2][0] =
    u[2][1] = u[2][3]= u[3][0] = u[3][1] = u[3][2] = 0;
    u[3][3] = exp(complexd(0, 1) * phi);
    return TwoQubitsEvolution(a, n, q1, q2, u, powproc, rank);
}
