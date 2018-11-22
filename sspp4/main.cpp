#include "mpi.h"
#include <iostream>
#include <cstdint>
#include <fstream>

using namespace std;

int tagdata = 1;
int tagb = 2;
int tagc = 3;

int main(int argc, char **argv)
{
	if (argc != 4) {
		cout << "A.dat b.dat c.dat" << endl;
		return 0;
	}
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	double time_start, time_finish, time;
	fstream finb(argv[2]);
	uint64_t mb, nb;
	char type;
	finb.read(&type, sizeof(type));
	finb.read((char *) &mb, sizeof(mb));
	finb.read((char *) &nb, sizeof(nb));
	double sumtime = 0;

	uint64_t n, m;
	if (rank == 0) {
		fstream fin(argv[1]);
		fin.read(&type, sizeof(type));
		fin.read((char *) &m, sizeof(m));
		fin.read((char *) &n, sizeof(n));
		int sizes[2];
		sizes[1] = n;
		sizes[0] = m; 

		if (m >= n) {

			long long begin = (long long) rank * m / size ;
			long long end = (long long) (rank + 1) * m /  size - 1;

			float a[end - begin + 1][n];
			for(int i = begin ; i <= end; ++i) {
				for (int j = 0 ; j < n; ++j) {
					fin.read((char *) &a[i][j], sizeof(a[i][j]));
				}
			}
    		for(int k = 1; k < size; ++k) { // k = rank 
    			MPI_Send(&sizes[0], 2, MPI_INT, k, tagdata, MPI_COMM_WORLD);
    			begin = (long long) k * m / size ;
    			end = (long long) (k + 1) * m /  size - 1;
    			float temp [end - begin + 1][n];
    			for(int i = begin ; i <= end; ++i) {
    				for (int j = 0 ; j < n; ++j) {
    					fin.read((char *) &temp[i-begin][j], sizeof(temp[i-begin][j]));
    				}
    				
    				MPI_Send(&temp[i-begin][0], n, MPI_FLOAT, k, tagdata, MPI_COMM_WORLD);
    			}
    			
    		}
    		float b[mb];
    		for (int i = 0 ; i < mb; ++i) {
    			finb.read((char *) &b[i], sizeof(b[i]));
    		}
    		float c[m];
    		begin = (long long) rank * m / size ;
    		end = (long long) (rank + 1) * m /  size - 1;
    		for (int i = begin; i <=end; ++i) {
    			c[i - begin] = 0;
    			for(int j = 0 ; j < n ; ++j) {
    				c[i - begin] += a[i - begin][j] * b[j];
    			}
    		}
    		for ( int i = 1; i < size; ++i) {
    			begin = (long long) i * m / size ;
    			end = (long long) (i + 1) * m /  size - 1;
    			MPI_Recv(&c[begin], end - begin  + 1, MPI_FLOAT, i, tagc, MPI_COMM_WORLD, &status);
    		}
    		fstream foutc;
    		foutc.open(argv[3], ios::out | ios::binary);
    		type = 'f';
    		foutc.write(&type, sizeof(type));
    		foutc.write((char *) &m, sizeof(m));
    		n = 1;
    		foutc.write((char *) &n, sizeof(n));
    		for(int i = 0 ; i < m ; ++i) {
    			foutc.write((char *) &c[i], sizeof(c[i]));
    		}
    	} else {
    		float a[m][n];
    		for(int i = 0 ; i < m; ++i) {
    			for (int j = 0 ; j < n; ++j) {
    				fin.read((char *) &a[i][j], sizeof(a[i][j]));
    			}
    		}
    		float b[mb];
    		for (int i = 0 ; i < mb; ++i) {
    			finb.read((char *) &b[i], sizeof(b[i]));
    		}
    		for(int k = 1; k < size; ++k) {
    			MPI_Send(&sizes[0], 2, MPI_INT, k, tagdata, MPI_COMM_WORLD);
    			long long begin = (long long) k * n / size;
    			long long end = (long long) (k + 1) * n /  size - 1;
    			float temp[end-begin +1][m];
    			for(int j = begin; j <= end; ++j) {
    				for(int i = 0 ; i < m; ++i) {
    					temp[j - begin][i] = a[i][j];
    				}
    				MPI_Send(&temp[j - begin][0], n, MPI_FLOAT, k, tagdata, MPI_COMM_WORLD);
    			}
    			MPI_Send(&b[begin], end - begin + 1, MPI_FLOAT, k, tagb, MPI_COMM_WORLD);
    		}
    		float c[m];
    		long long begin = (long long) rank * m / size ;
    		long long end = (long long) (rank + 1) * m /  size - 1;
    		for(int i = 0;  i < m; ++i) {
    			c[i] = 0;
    			for(int j = begin; j <= end; ++j) {
    				c[i] += a[i][j - begin] * b[j - begin];
    			}
    		}

    		//MPI_Reduce

    		fstream foutc;
    		foutc.open(argv[3], ios::out | ios::binary);
    		type = 'f';
    		foutc.write(&type, sizeof(type));
    		foutc.write((char *) &m, sizeof(m));
    		n = 1;
    		foutc.write((char *) &n, sizeof(n));
    		for(int i = 0 ; i < m ; ++i) {
    			foutc.write((char *) &c[i], sizeof(c[i]));
    		}
    	}
    	

    } else {
    	int sizes[2];
    	MPI_Recv(&sizes[0], 2, MPI_INT, 0, tagdata, MPI_COMM_WORLD, &status);
    	m = sizes[0];
    	n = sizes[1];
    	if (m >= n) {
    		float b[mb];
    		for (int i = 0 ; i < mb; ++i) {
    			finb.read((char *) &b[i], sizeof(b[i]));
    		}
    		long long begin = (long long) rank * m / size ;
    		long long end = (long long) (rank + 1) * m /  size - 1;
    		float a[end - begin + 1][n];
    		for(int i = begin; i <= end; ++i) {
    			MPI_Recv(&a[i - begin][0], n, MPI_FLOAT, 0, tagdata, MPI_COMM_WORLD, &status);
    		}
    		float c[end - begin + 1];
    		time_start = MPI_Wtime();
    		for(int i = begin ; i <= end; ++i) {
    			c[i - begin] = 0;
    			for(int j = 0 ; j < n ; ++j) {
    				c[i - begin] += a[i - begin][j] * b[j];
    			}
    		}
    		time_finish = MPI_Wtime();
    		MPI_Send(&c[0], end - begin  + 1, MPI_FLOAT, 0, tagc, MPI_COMM_WORLD);
    		time = time_finish - time_start;
    	} else {
    		long long begin = (long long) rank * n / size ;
    		long long end = (long long) (rank + 1) * n /  size - 1;
    		float a[end - begin + 1][m];
    		for(int i = begin; i <= end; ++i) {
    			MPI_Recv(&a[i - begin][0], n, MPI_FLOAT, 0, tagdata, MPI_COMM_WORLD, &status);
    		}
    		float b[end - begin + 1];
    		MPI_Recv(&b[0], end - begin + 1, MPI_FLOAT, 0, tagb, MPI_COMM_WORLD, &status);
    		float c[m];
    		time_start = MPI_Wtime();
    		for(int i = 0 ; i < m; ++i) {
    			c[i] = 0;
    			for (int j = begin ; j <= end; ++j) {
    				c[i] += b[j - begin] * a[j-begin][i];
    			}
    		}
    		time_finish = MPI_Wtime();
    		
    		time = time_finish - time_start;
    	}
    }
    MPI_Reduce(&time, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0 ) {
    	cout << "all time :" << sumtime << endl;
    }
    MPI_Finalize();
    return 0;	
}