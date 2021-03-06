#include <iostream>
#include <complex>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "mpi.h"

typedef std::complex<double> complexd;
const double MAXD = 232323;
using namespace std;




complexd *qubit_transform(complexd *a, int  n, complexd **u, int  k, long powproc,int rank) {
	long vec_length = 1 << (n - powproc);
	long dist = 1 << (n - k);
	long start = vec_length * rank; // 2*n*rank/size (part number of the whole vector)
	complexd *b = new complexd[vec_length];
	//cout << "vec_length: " << vec_length << " dist " << dist << " start " << start << endl;
	if (dist < vec_length) {
		for(int i = 0; i < vec_length; ++i) {
			b[i] = u[(i + start) & dist >> (n - k)][0] * a[(((i + start) | dist) ^ dist) - start] + u[((i + start) & dist) >> (n - k)][1] * a[ ((i + start) | dist) - start];
		}
	} else {
		int needrank;
		if ((start & dist) == 0) {                           //look at the bit that needs to be changed
			needrank = (start | dist) / vec_length;
		} else {
			needrank = (start & ~dist) / vec_length;
		}
		//cout << "rank "<< rank << " needrank" << needrank << endl<<endl;
		complexd *temp = new complexd[vec_length];
		//MPI_Send(a, vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD);
		//MPI_Recv(temp, vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(a, vec_length * 2, MPI_DOUBLE, needrank, 0, temp, vec_length * 2, MPI_DOUBLE, needrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		complexd *vec0, *vec1;
		if (rank < needrank) { // bit needs to be changed in a == 0
			vec0 = a;
			vec1 = temp;
		} else {			   // bit needs to be changed in temp == 0
			vec0 = temp;
			vec1 = a;
		}
		for(int i = 0; i < vec_length; ++i) {
			b[i] = u[(i + start) & dist >> (n - k)][0] * vec0[i] + u[((i + start) & dist) >> (n - k)][1] * vec1[i]; // 
		}
		delete []temp;
	}
	
	return b;
}


int main(int argc, char **argv) {
	if( argc != 4) {
		cout << "input: n k mode (1-file \"in.bin\" , 2-random)" << endl;
		return 0;
	}
	double start_time1, start_time2, end_time1, end_time2;
	int n, k;
	int mode;
	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%d", &k);
	sscanf(argv[3], "%d", &mode);

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
	MPI_Init(&argc, &argv);
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
	start_time1 = MPI_Wtime();
	
	if (mode == 2) {
		double dlina = 0;
		srand(timetemp * (rank + 1));
		long start = vec_length * rank;
		for(int i = 0 ; i < vec_length; ++i) {
			a[i] = complexd((double)rand()/RAND_MAX * MAXD, (double) rand()/RAND_MAX * MAXD);	
			//a[i] = i + start;
			dlina += norm(a[i]);
		} 
		double temp;
		MPI_Allreduce(&dlina, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		dlina = sqrt(temp); 
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
		for(int i = 0 ; i < vec_length; ++i) {
			dlina += norm(a[i]);
		} 
		double temp;
		MPI_Allreduce(&dlina, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		dlina = sqrt(temp); 
		for(int i = 0; i < vec_length; ++i) {
			a[i] = a[i] / dlina;
		}
	}
	
	end_time1 = MPI_Wtime();
	



	start_time2 = MPI_Wtime();

	complexd *b = qubit_transform(a,n,u,k, powproc, rank);
	
	end_time2 = MPI_Wtime();


	double time1 = end_time1 - start_time1;
	double time2 = end_time2 - start_time2;
	double timelocal = time1 + time2;
	double sumtime1, sumtime2, maxtime;
	MPI_Reduce(&time1, &sumtime1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time2, &sumtime2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&timelocal, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		cout << "sumtime1: " << sumtime1 <<" sumtime2: " << sumtime2 << " maxtime: " << maxtime << endl;
	}
	
	//file output
	
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
	delete []u;
	delete []a;
	delete []b;
	MPI_Finalize();
	return 0;
}
