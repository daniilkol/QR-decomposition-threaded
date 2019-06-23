#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "qr_razlog.h"
#include <pthread.h>
#include "get_time.h"
#include <time.h>

#define N 100


int main(int argc, char* argv[])
{
    int n,nthreads,indicator,i,j,t;
    double *a,*copy_of_a,*x,t_full,time,*c,*s;
    pthread_t * threads;
    ARGS * args;
    struct timespec t_full_begin, t_full_end;
    char* name = new char[1024];
    if (argc == 2 && strcmp(argv[1], "help") == 0)
    {
        printf("First argument: the number of threads\n");
        printf("Second argument: 1 for file, 0 for formula\n");
        printf("Third argument: if from file, then file name; if from formula, then matrix size\n");
        return 0;
    } else if (argc != 4)
            {
                fprintf(stderr, "The number of arguments is wrong\n");
                return -1;
            }

    if (sscanf(argv[1], "%d", &nthreads) != 1 || nthreads <= 0)
    {
        fprintf(stderr, "Failure to read nthreads\n");
        return -2;
    }

    if (sscanf(argv[2], "%d", &indicator) != 1 || !(indicator == 1 || indicator == 0))
    {
        fprintf(stderr, "Failure to read from_indicator\n");
        return -3;
    }

    if (indicator == 1 && sscanf(argv[3], "%s", name) != 1)
    {
        fprintf(stderr, "Failure to read the file\n");
        return -4;
    }  else if (indicator == 0 && (sscanf(argv[3], "%d", &n) != 1 || n <= 1))
            {
                fprintf(stderr, "Failure to read size\n");
                return -5;
            }

    if (indicator)
    {
        FILE*fp;
        if(!(fp=fopen(name,"r")))
        {
            printf("Check your filename,please\n");
            return -7;
        }
        if(fscanf(fp,"%d",&n)!=1)
        {
            printf("Cannot read your n from file,sorry\n");
            return -8;
        }
        if(n<=0)
        {
            printf("Your n is negative, sorry\n");
            return -9;
        }
        a=new double[n*n];
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
            {
                if(fscanf(fp,"%lf",a+i*n+j)!=1)
                {
                    printf("Cannot read %s\n",name);
                    fclose(fp);
                    delete [] a;
                    return -2;
                }
            }
        fclose(fp);
        } else
            {
                a=new double[n*n];
                init_matrix(a,n);
            }
    printf("MATRIX A BEFORE:\n");
    print_matrix(a,n);
    copy_of_a=new double [n*n];
    x=new double [n*n];
    c=new double [n];
    s=new double [n];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            copy_of_a[i*n+j]=a[i*n+j];
        }
    }
    matrix_E(x, n);
    threads=new pthread_t [nthreads];
    args=new ARGS [nthreads];
    int error=0;
    for (t = 0; t < nthreads; t++)
    {
        args[t].a = a;
        args[t].copy_of_a = copy_of_a;
        args[t].x = x;
        args[t].n = n;
        args[t].thread_num = t;
        args[t].total_threads = nthreads;
        args[t].error = &error;
        args[t].sin=s;
        args[t].cos=c;
    }
    clock_gettime(CLOCK_MONOTONIC,&t_full_begin);
    for (i = 0; i < nthreads; i++)
    {
      if (pthread_create (threads + i, 0,qr_razlog_threaded,args + i))
        {
          fprintf (stderr, "cannot create thread #%d!\n",i);
          delete [] a;
          delete [] copy_of_a;
          delete [] x;
          delete [] threads;
          delete [] args;
          delete [] c;
          delete [] s;
          return -10;
        }
    }
    for (i = 0; i < nthreads; i++)
    {
        if (pthread_join (threads[i], 0))
            fprintf (stderr, "cannot wait thread #%d!\n", i);
    }
    clock_gettime(CLOCK_MONOTONIC,&t_full_end);
    t_full = t_full_end.tv_sec-t_full_begin.tv_sec+(double)(t_full_end.tv_nsec-t_full_begin.tv_nsec)/1000000000;
    printf("\n\n");
    printf ("Total full time = %f\n",t_full);
    if(error==0)
    {
        printf("\nREVERSE MATRIX:\n");
        print_matrix(x,n);
        time=clock();
        printf("\nNEVYAZKA=%e\n",nevyazka(copy_of_a,x,a,n));
        time=(clock()-time)/CLOCKS_PER_SEC;
        printf("Elapse time for counting discrepancy=%.6f\n",time);
    }
    delete [] a;
    delete [] copy_of_a;
    delete [] x;
    delete [] threads;
    delete [] args;
    delete [] name;
    delete [] c;
    delete [] s;
    return 1;

}
