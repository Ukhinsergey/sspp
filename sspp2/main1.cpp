#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <papi.h>

float ptime;
float rtime;
long long flpops;
float mflops;
uint64_t blocksize = 32;
uint64_t myblock = 72;//myblock = sqrt(ml/3);
using namespace std;


void mulmatr(fstream &a, fstream &b, fstream &c, int r){
	uint64_t na, ma, nb, mb, i, j, k;
	uint64_t  n;
	a.read((char *) &na, sizeof(na));
	a.read((char *) &ma, sizeof(ma));
	b.read((char *) &nb, sizeof(nb));
	b.read((char *) &mb, sizeof(mb));
	if (ma != nb) {
		cout << "cant multiply matrs because of size\n";
		return;
	}
	float **matra = new float*[na];
	float **matrb = new float*[nb];
	float **matrc = new float*[na];
	for (i = 0; i < na; ++i) {
		matra[i] = new float[ma];
		matrc[i] = new float[mb];
	}
	for (i = 0; i < nb; ++i) {
		matrb[i] = new float[mb];
	}

	for (i = 0 ; i < na; ++i) {
		for(j = 0; j < ma; ++j) {
			a.read((char *) &matra[i][j], sizeof(matra[i][j]));
		}
	}
	for (i = 0 ; i < nb; ++i){
		for (j = 0; j < mb; ++j) {
			b.read((char *) &matrb[i][j], sizeof(matrb[i][j]));
		}
	}
	for( i = 0 ; i < na; ++i) {
		for(j = 0 ; j < mb; ++j) {
			matrc[i][j] = 0;
		}
	}

	ofstream time, flops;
	time.open("time.txt", ios::out | ios::app);
	flops.open("flops.txt", ios::out | ios::app);
	//tlb.open("tlb.txt", ios::out | ios::app);
	n = na;
	uint64_t i1, j1, k1;
	switch (r) {
		case 0:
			if (PAPI_OK != PAPI_flops(&rtime, &ptime, &flpops, &mflops)) {
				cout << "error papi start flops" << endl;
				exit(1);
			}
			for(i = 0 ; i < n; i+=blocksize) {
				for( j = 0 ; j < n; j+=blocksize) {
					for ( k = 0; k < n; k+=blocksize) {
						//block mul
						uint64_t imin = min(n, i + blocksize);
						uint64_t jmin = min(n, j + blocksize);
						uint64_t kmin = min(n, k + blocksize);
						for(i1 = i; i1 < imin; ++i1) {
							for (j1 = j; j1 < jmin ; ++j1) {
								for (k1 = k; k1 < kmin; ++k1) {
									matrc[i1][j1] += matra[i1][k1] * matrb[k1][j1];
								}
							}
						}
					}
				}
			}
			if (PAPI_OK != PAPI_flops(&rtime, &ptime, &flpops, &mflops)) {
				cout << "error papi stop flops" << endl;
				exit(1);
			}
			break;

		case 1:
			if (PAPI_OK != PAPI_flops(&rtime, &ptime, &flpops, &mflops)) {
				cout << "error papi start" << endl;
				exit(1);
			}
			for(i = 0 ; i < n; i+=blocksize) {
				for( j = 0 ; j < n; j+=blocksize) {
					for ( k = 0; k < n; k+=blocksize) {
						//block mul
						uint64_t imin = min(n, i + blocksize);
						uint64_t jmin = min(n, j + blocksize);
						uint64_t kmin = min(n, k + blocksize);
						for(i1 = i; i1 < imin; ++i1) {
							for (k1 = k; k1 < kmin; ++k1) { 
								for (j1 = j; j1 < jmin; ++j1) {
									matrc[i1][j1] += matra[i1][k1] * matrb[k1][j1];
								}
							}
						}
					}
				}
			}

			if (PAPI_OK != PAPI_flops(&rtime, &ptime, &flpops, &mflops)) {
				cout << "error papi stop" << endl;
				exit(1);
			}
			break;

		case 2:
			blocksize = myblock;
			if (PAPI_OK != PAPI_flops(&rtime, &ptime, &flpops, &mflops)) {
				cout << "error papi start" << endl;
				exit(1);
			}
			for(i = 0 ; i < n; i+=blocksize) {
				for( j = 0 ; j < n; j+=blocksize) {
					for ( k = 0; k < n; k+=blocksize) {
							//block mul
						uint64_t imin = min(n, i + blocksize);
						uint64_t jmin = min(n, j + blocksize);
						uint64_t kmin = min(n, k + blocksize);
						for(i1 = i; i1 < imin; ++i1) {
							for (k1 = k; k1 < kmin; ++k1) { 
								for (j1 = j; j1 < jmin; ++j1) {
									matrc[i1][j1] += matra[i1][k1] * matrb[k1][j1];
								}
							}
						}
					}
				}
			}
			if (PAPI_OK != PAPI_flops(&rtime, &ptime, &flpops, &mflops)) {
				cout << "error papi stop" << endl;
				exit(1);
			}
			break;
		default : 
			cout << "error: mode" << r << endl;
			return;
	}
	
	//tlb << r  << ' ' << values[3] << endl; 
	time << n << ' ' << r << ' ' << ptime << endl;
	flops << n << ' ' << r << ' ' << flpops << endl;
	
	char type ='f';

	c.write((char*) &type, sizeof(type));
	c.write((char*) &na, sizeof(na));
	c.write((char*) &mb, sizeof(mb));
	for (i = 0 ; i < na; ++i) {
		for ( j = 0; j < mb; ++j) {
			c.write((char*) &matrc[i][j], sizeof(matrc[i][j]));
		}
	}


	for (i = 0; i < na; ++i) {
		delete []matra[i];
		delete []matrc[i];
	}
	for (i = 0; i < nb; ++i) {
		delete []matrb[i];
	}
	delete []matra;
	delete []matrb;
	delete []matrc;
}

int main(int argc, char **argv) {
	//ijk == 0 , ikj == 1 
	
	if (argc != 5) {
		cout << "format: A.dat B.dat input.dat  r( 0/1/2 )" << endl;
		return 0;
	}

	fstream a,b,c;
	a.open(argv[1], ios::in | ios::binary);
	b.open(argv[2], ios::in | ios::binary);
	char ta, tb;

	c.open(argv[3], ios::out | ios::binary);
	int r;
	sscanf(argv[4], "%d", &r);
	a.read(&ta, sizeof(ta));
	b.read(&tb, sizeof(tb));
	mulmatr(a, b, c, r);
	a.close();
	b.close();
	c.close();

	return 0;
}