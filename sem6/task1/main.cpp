#include <omp.h>
#include <iostream>
#include <complex>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cstdlib>

typedef std::complex<double> complexd;
const double MAXD = 232323;
using namespace std;




complexd *qubit_transform(complexd *a, int n, complexd **u, int k, double dlina) {
	int num_qubits = 1 << n;
	int temp = 1 << (n - k); //(int) pow(2, n - k);
	complexd *b = new complexd[num_qubits];
	#pragma omp parallel
	{	
		#pragma omp for 
		for(int i = 0 ; i < num_qubits; ++i) {
		b[i] = (u[(i & temp) >> (n - k)][1] * a[ i | temp] + u[(i & temp) >> (n - k)][0] * a[(i | temp) ^ temp]) / dlina; // 
		}
	}
	
	return b;
}


int main(int argc, char **argv) {
	if( argc != 4) {
		cout << "input: n, k, numtreads" << endl;
		return 0;
	}
	double start_time1, start_time2, end_time1, end_time2;
	int n, k, numtreads;
	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%d", &k);
	sscanf(argv[3], "%d", &numtreads);
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
	int num_qubits = 1 << n; //(int) pow(2, n);
	complexd *a = new complexd[num_qubits];
	start_time1 = omp_get_wtime();
	double dlina = 0;
	#pragma omp parallel
	{
		unsigned int seed = omp_get_thread_num();
		#pragma omp for reduction(+ : dlina)
		for(int i = 0 ; i < num_qubits; ++i) {
			a[i] = complexd(((double)rand_r(&seed))/RAND_MAX * MAXD, ((double)rand_r(&seed))/RAND_MAX * MAXD);	
			//a[i] = i;
			dlina += norm(a[i]);
		} 
	}
	end_time1 = omp_get_wtime();
	start_time2 = omp_get_wtime();

	complexd *b = qubit_transform(a,n,u,k, sqrt(dlina));
	end_time2 = omp_get_wtime(); 
	cout << "work time: " << end_time1 - start_time1 <<' ' << end_time2 - start_time2 << endl;

	
	/*ofstream out("out.txt");
	for(int i = 0 ; i < num_qubits; ++i) {
		out << i << ' ' << b[i] << ' ' << endl;
	}
	out.close();
	*/

	for(int i = 0 ; i < 2; ++i) {
		delete [] u[i];
	}
	delete []u;
	delete []a;
	delete []b;
	
	return 0;
}
