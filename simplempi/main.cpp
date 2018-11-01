#include "mpi.h"
#include <stdio.h>
#include <vector>
#include <fstream>


using namespace std;

int main(int argc, char **argv)
{   

    double time_start, time_finish;
    time_start = MPI_Wtime();
    long long  begin, end;
    int tagdata = 1;
    int tagtime = 2;
    sscanf(argv[1], "%llu", &begin);
    sscanf(argv[2], "%llu", &end);
    vector<bool> prime(end + 1, true);
    int rank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    prime[0] = prime[1] = false; 
    long long  temp = 0;
    for(temp = 2; temp * temp <= end; ++temp) {
        if (prime[temp]) {
            for(long long  j = temp * temp; j * j <= end; j += temp ) {
                prime[j] = false;
            }
        }

    }
    long long  ibegin = (rank - 1) * (end - max(begin,temp) + 1) / (size - 1) + max(begin,temp);
    long long  iend = rank * (end - max(begin,temp) + 1) / (size - 1) + max(begin,temp) - 1;  
    //cout << rank << ' ' << ibegin << ' ' << iend << endl;
    if (rank != 0) {
        for(long long  j = 2; j * j <= end; ++j) {
            if(prime[j]) {
                for(long long  i = (ibegin/j + 1 * (ibegin % j != 0)) * j; i <= min(iend, end); i += j) {
                    prime[i] = false;
                }
            }
            
        }
        /*for(long long  i = ibegin; i <= min(iend,end); ++i) {
            cout << rank << ' ' << i << ' ' << prime[i] << endl;
        }*/

        int i = -1;
        MPI_Send(&i, 1, MPI_LONG_LONG, 0, tagdata, MPI_COMM_WORLD);
        time_finish = MPI_Wtime();
        double time = time_finish - time_start;
        MPI_Send(&time, 1, MPI_DOUBLE, 0, tagtime, MPI_COMM_WORLD);
    } else {
        std::vector<long long > simple;
        long long  count = 0;
        std::vector<double> time;
        int k = 0;
       
        for(temp = 2; temp * temp <= end; ++temp) {
            if (prime[temp]) {
                for(long long  j = temp * temp; j * j <= end; j += temp ) {
                    prime[j] = false;
                }
            }

        }
        for (k = 1; k < size; ++k) {
             {
                long long  tmp;
                MPI_Recv(&temp, 1, MPI_LONG_LONG,k, tagdata, MPI_COMM_WORLD, &status);
                if (tmp != -1) {
                    simple.push_back(tmp);
                    count++;
                } else {
                    double tmptime;
                    MPI_Recv(&tmptime, 1, MPI_DOUBLE ,k, tagdata, MPI_COMM_WORLD, &status);
                    time.push_back(tmptime);
                    break;
                }
            }
        }
        ofstream out;
        for(long long  i = ibegin; i <= min(iend, end); ++i) {
            out << simple[i] << ' ';
            count++;
        }
        out.open(argv[3]);
        for(long long i = 0; i < count; ++i) {
            out << simple[i] << ' ';
        }
        cout << count << endl;
        for (k = 1; k < size ; ++k) {
           // cout << time[k-1] << ' ' ;
        }
        
    }
    MPI_Finalize();
}