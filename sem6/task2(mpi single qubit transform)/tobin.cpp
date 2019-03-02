#include <iostream>
#include <fstream>
#include <complex>

typedef std::complex<double> complexd;

using namespace std;

int main(int argc, char **argv) 
{
	if (argc != 2) {
		cout << " format: in.txt in.bin" ;
		return 0;
	}
	fstream filein(argv[1]);
	int n;
	filein >> n;
	fstream fout(argv[2], ios::out | ios::binary | ios::trunc);
	fout.write((char *) &n, sizeof(n));
	long vec_length = 1 << n;
	for(int i = 0 ; i < vec_length; ++i) {
		complexd temp;
		filein >> temp;
		fout.write((char *)	&temp, sizeof(double));
		fout.write((char *) &temp + sizeof(double), sizeof(double));
	}
	fout.close();
	return 0;


}