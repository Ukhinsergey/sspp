
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <time.h>


using namespace std;



struct smpl{
    long long  start;
    long long end;
    int size;
    int rank;
    int count;
    double time;
    std::vector<bool> * prime;
    smpl(int s = 1 , int r = 1, int c = 0, long long st = 0, long long e = 1, double t = 0, std::vector<bool> *p = NULL){
        size = s;
        rank = r;
        count = c;
        start = st;
        end = e;
        time = t;   
        prime = p;
    }
};

void *worker(void *atr) {

    smpl *a = (smpl*) atr;
    //cout << "worke # " << a->rank << " count: " << a->count << "size: " << a->size << "start: " << a->start << endl;
    long long  ibegin = (long long) (a->rank - 1)  * (a->end - a->start + 1)/(a->size - 1)  + a->start;
    long long  iend = (long long ) a->rank  * (a->end - a->start + 1)/(a->size - 1)  + a->start- 1; 
    std::vector<bool> prime1(iend - ibegin + 1, true);
    timespec time_start;
    timespec time_finish;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_start);

    for(long long  j = 2; j * j <= a->end; ++j) {
        if((*(*a).prime)[j]) {
            for(long long  i = (ibegin/j + 1 * (ibegin % j != 0)) * j; i <= min(iend, a->end); i += j) {
                prime1[i - ibegin] = false;
            }
        }

    }

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_finish);
    for(long long  i = ibegin; i <= min(iend, a->end); ++i) {

        if (prime1[i - ibegin]) {   
            a->count = a->count + 1;
        }
    }
    a->time = double (time_finish.tv_sec - time_start.tv_sec + 1e-9 * (time_finish.tv_nsec - time_start.tv_nsec));
    cout << "worke # " << a->rank << " count: " << a->count << " time: " << (double)a->time << endl;
    return a;
}


int main(int argc, char **argv)
{   

    long long  begin, end;
    sscanf(argv[1], "%llu", &begin);
    sscanf(argv[2], "%llu", &end);
    int numthread;
    sscanf(argv[4], "%d", &numthread);
    //vector<bool> prime(end + 1, true);
    std::vector<bool> prime((int) sqrt(end + 1.0), true);
    prime[0] = prime[1] = false; 
    long long  temp;
    int count = 0;
    for(temp = 2; temp * temp <= end; ++temp) {
        if (prime[temp]) {
            count++;
            for(long long  j = temp * temp; j * j <= end; j += temp ) {
                prime[j] = false;
            }
        }

    }

    pthread_t tid[numthread];
    smpl param[numthread];
    for(int i = 0 ; i < numthread; ++i) {
        param[i] = smpl(numthread + 1, i + 1, 0, max(begin, temp), end, 0, &prime);
        pthread_create(&tid[i], NULL, worker, &param[i]);
    }
    double maxtime = 0;
    double sumtime = 0;
    for(int i = 0; i < numthread; ++i) {
        smpl *smtemp;
        pthread_join(tid[i],(void **) &smtemp);
        count += smtemp->count;
        sumtime += smtemp->time;
        if (smtemp->time > maxtime) {
        	maxtime = smtemp->time;
        }
    }
    cout << count << endl;
    cout << "maxtime: " << maxtime << endl << "sumtime: " << sumtime << endl;
    return 0;
}
