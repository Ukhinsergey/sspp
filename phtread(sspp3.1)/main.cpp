
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <pthread.h>

using namespace std;



struct smpl{
    long long  start;
    long long end;
    int size;
    int rank;
    int count;
    clock_t time;
    std::vector<bool> * prime;
    smpl(int s = 1 , int r = 1, int c = 0, long long st = 0, long long e = 1, clock_t t = 0, std::vector<bool> *p = NULL){
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
    clock_t time_start = clock();

    for(long long  j = 2; j * j <= a->end; ++j) {
        if((*(*a).prime)[j]) {
            for(long long  i = (ibegin/j + 1 * (ibegin % j != 0)) * j; i <= min(iend, a->end); i += j) {
                prime1[i - ibegin] = false;
            }
        }

    }
    clock_t time_finish = clock();
    for(long long  i = ibegin; i <= min(iend, a->end); ++i) {

        if (prime1[i - ibegin]) {   
            a->count = a->count + 1;
        }
    }
    a->time = time_finish - time_start;
    cout << "worke # " << a->rank << " count: " << a->count << " time: " << (double)a->time/CLOCKS_PER_SEC << endl;
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
    }
    for(int i = 0 ; i < numthread; ++i) {
        pthread_create(&tid[i], NULL, worker, &param[i]);
    }
    for(int i = 0; i < numthread; ++i) {
        smpl *smtemp;
        pthread_join(tid[i],(void **) &smtemp);
        count += smtemp->count;
    }
    cout << count << endl;
    return 0;
}
