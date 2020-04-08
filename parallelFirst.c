#include <stdio.h>
#include <omp.h>

int main()
{
    double A[1000];
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf("hello(%d)", ID);
        printf(" world(%d) \n",ID);
    }

}