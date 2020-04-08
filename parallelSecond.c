#include <stdio.h>
#include <omp.h>
static long num_steps = 100000;
double step;
#define NUM_THREADS 4
void main()
{
    int i, nthreads;
    double pi, sum[NUM_THREADS];
    double time_start, time_finish;
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    time_start = omp_get_wtime();
    #pragma omp parallel
    {
        int i, nthrds;
        int ID = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if(ID==0) nthreads = nthrds;
        double x;
        for (i=ID, sum[ID]=0.0;i<num_steps;i+=nthrds){
            x = (i+0.5)*step;
            sum[ID] += 4.0/(1.0+x*x);
        }
    }
    for (i=0, pi=0.0; i<nthreads; i++)pi += step * sum[i];   
    time_finish = omp_get_wtime();
    printf("Area is: %f", pi);
    printf("Time required: %f", time_finish-time_start);
}