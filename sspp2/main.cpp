#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <papi.h>

const size_t EVENT_COUNT = 3;
int events[] = {PAPI_TOT_CYC, PAPI_L1_DCM, PAPI_L2_DCM};//{PAPI_TOT_INS, PAPI_L1_DCM, PAPI_L2_DCM, PAPI_TLB_DM};
long long values[EVENT_COUNT];
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

	ofstream dat, ins, l1, l2, tlb;
	dat.open("dat.txt", ios::out | ios::app);
	ins.open("ins.txt", ios::out | ios::app);
	l1.open("l1.txt", ios::out | ios::app);
	l2.open("l2.txt", ios::out | ios::app);
	tlb.open("tlb.txt", ios::out | ios::app);
	n = na;
	uint64_t i1, j1, k1;
	clock_t time;
	switch (r) {
		case 0:
			time = clock();
			if (int retval = PAPI_start_counters(events, EVENT_COUNT) != PAPI_OK) {
				cout << "error papi start" << endl;
				exit(1);
			}
			for(i = 0 ; i < n; i+=blocksize) {
				for( j = 0 ; j < n; j+=blocksize) {
					for ( k = 0; k < n; k+=blocksize) {
						//block mul
						for(i1 = i; i1 < min(n, i + blocksize); ++i1) {
							for (j1 = j; j1 < min(n, j + blocksize); ++j1) {
								for (k1 = k; k1 < min(n, k + blocksize); ++k1) {
									matrc[i1][j1] += matra[i1][k1] * matrb[k1][j1];
								}
							}
						}
					}
				}
			}
			if (int retval = PAPI_stop_counters(values, EVENT_COUNT) != PAPI_OK) {
				cout << "error papi stop" << endl;
				exit(1);
			}
			time -= clock();
			break;

		case 1:
			time = clock();
			if (int retval = PAPI_start_counters(events, EVENT_COUNT) != PAPI_OK) {
				cout << "error papi start" << endl;
				exit(1);
			}
			for(i = 0 ; i < n; i+=blocksize) {
				for( j = 0 ; j < n; j+=blocksize) {
					for ( k = 0; k < n; k+=blocksize) {
						//block mul
						for(i1 = i; i1 < min(n, i + blocksize); ++i1) {
							for (k1 = k; k1 < min(n, k + blocksize); ++k1) { 
								for (j1 = j; j1 < min(n, j + blocksize); ++j1) {
									matrc[i1][j1] += matra[i1][k1] * matrb[k1][j1];
								}
							}
						}
					}
				}
			}

			if (int retval = PAPI_stop_counters(values, EVENT_COUNT) != PAPI_OK) {
				cout << "error papi stop" << endl;
				exit(1);
			}
			time -= clock();
			break;

		case 2:
			blocksize = myblock;
			time = clock();
			if (int retval = PAPI_start_counters(events, EVENT_COUNT) != PAPI_OK) {
				cout << "error papi start" << endl;
				exit(1);
			}
			for(i = 0 ; i < n; i+=blocksize) {
				for( j = 0 ; j < n; j+=blocksize) {
					for ( k = 0; k < n; k+=blocksize) {
							//block mul
						for(i1 = i; i1 < min(n, i + blocksize); ++i1) {
							for (k1 = k; k1 < min(n, k + blocksize); ++k1) { 
								for (j1 = j; j1 < min(n, j + blocksize); ++j1) {
									matrc[i1][j1] += matra[i1][k1] * matrb[k1][j1];
								}
							}
						}
					}
				}
			}
			if (int retval = PAPI_stop_counters(values, EVENT_COUNT) != PAPI_OK) {
			cout << "error papi stop" << endl;
				exit(1);
			}
			time -= clock();
			break;
		default : 
			cout << "error: mode" << r << endl;
			return;
	}
	dat << n << ' ' << r << ' ' << fixed << setprecision(6) << ((double) -time)/CLOCKS_PER_SEC << ' ' << endl;
	ins << n << ' ' << r << ' ' << values[0] << endl;
	l1 << n << ' ' << r  << ' ' << values[1] << endl;
	l2 << n << ' ' << r  << ' ' << values[2] << endl;
	//tlb << n << ' ' << r  << ' ' << values[3] << endl; */ l1<< n << ' ' << r  << ' ' << values[0] << endl; l2<< n << ' ' << r  << ' ' << values[1] << endl;
	dat.close();
	ins.close();
	l1.close();
	l2.close();
	tlb.close();

	
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