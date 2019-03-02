#include <fstream>
#include <complex>
#include <iostream>


typedef std::complex<double> complexd;

using namespace std;

int main(int argc, char **argv)
{
	if (argc !=2) {
		cout << " format : input.bin" << endl;
		return 0;
	}

	fstream fin(argv[1], ios::in | ios::binary);
	fstream fout("input.txt", ios::out | ios::trunc);
	int n;
	fin.read((char *) &n, sizeof(n));
	fout << n << endl;

	long vec_length = 1 << n;
	for(int i = 0 ; i < vec_length; ++i) {
		double re, im;
		fin.read((char *) &re, sizeof(re));
		fin.read((char *) &im, sizeof(im));
		fout << "(" << re <<"," << im <<")" << endl;
	}
	fout.close();
	fin.close();
	return 0;
}