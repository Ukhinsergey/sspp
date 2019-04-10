#include <fstream>
#include <complex>
#include <iostream>


typedef std::complex<double> complexd;

using namespace std;

int main(int argc, char **argv)
{
	if (argc != 3) {
		cout << " format : in1.bin in2.txt " << endl;
		return 0;
	}

	fstream fin1(argv[1], ios::in | ios::binary);
	fstream fin2(argv[2], ios::in | ios::binary);
	int n1;
	int n2;
	fin1.read((char *) &n1, sizeof(n1));
	fin2.read((char *) &n2, sizeof(n2));
	if (n1 != n2) {
		cout << "false 0" << endl;
		fin1.close();
		fin2.close();
		return 0;
	}
	int n = n1;
	long vec_length = 1 << n;
	for(int i = 0 ; i < vec_length; ++i) {
		double re1, im1, re2, im2;
		fin1.read((char *) &re1, sizeof(re1));
		fin1.read((char *) &im1, sizeof(im1));
		fin2.read((char *) &re2, sizeof(re2));
		fin2.read((char *) &im2, sizeof(im2));
		if ( re1 != re2 || im1 != im2) {
			cout << "false 1" << endl;
			fin1.close();
			fin2.close();
			return 0;
		}
	}
	cout << "True" << endl;
	fin1.close();
	fin2.close();
	return 0;
}