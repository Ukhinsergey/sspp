#include <fstream>
#include <complex>
#include <iostream>


typedef std::complex<double> complexd;

using namespace std;

int main(int argc, char **argv) {
    if (argc != 3) {
        cout << " format : input.bin output.txt " << endl;
        return 0;
    }

    fstream fin(argv[1], ios::in | ios::binary);
    fstream fout(argv[2], ios::out | ios::trunc);
    int n;
    fin.read(reinterpret_cast<char *> &n, sizeof(n));
    fout << n << endl;

    long vec_length = 1 << n;
    for (int i = 0 ; i < vec_length; ++i) {
        double re, im;
        fin.read(reinterpret_cast<char *> &re, sizeof(re));
        fin.read(reinterpret_cast<char *> &im, sizeof(im));
        fout << "(" << re <<"," << im <<")" << endl;
    }
    fout.close();
    fin.close();
    return 0;
}
