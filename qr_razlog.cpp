#include "qr_razlog.h"
#include "synchronize.h"
#include "get_time.h"
#include <pthread.h>
#include <time.h>

double threads_total_time;
static pthread_mutex_t threads_total_time_mutex
    = PTHREAD_MUTEX_INITIALIZER;
int key=0;

void matrix_E(double* a, int n)
{
    int i, j, p;
    for(i=0;i<n;i++)
    {
        p=i*n;
        for(j=0;j<n;j++)
        {
            if(i==j)
                a[p+j]=1;
            else a[p+j]=0;
        }
    }
}
void transponirovanie(double* a,int n)
{
    int i,j;
    double p;
    for (i=0;i<n;i++)
    {
      for(j=i;j<n;j++)
      {
          p=a[i*n+j];
          a[i*n+j]=a[j*n+i];
          a[j*n+i]=p;
      }
    }
}
#define MAX_N 8
void print_qr (double* a,double* x,int n)
{
    int i,j,p;
    int n_max=(n>MAX_N?MAX_N:n);
    transponirovanie(x,n);
    for(i=0;i<n_max;i++)
    {
        p=i*n;
        printf("|");
        for(j=0;j<n_max;j++)
            printf("  %.2f",x[p+j]);
        printf("\t|");
        printf("         ");
        printf("|");
        for(j=0;j<n_max;j++)
            printf("  %.2f",a[p+j]);
        printf("|");
        printf("\n");
    }
    transponirovanie(x,n);
}
void obratnii_xod (double* a,double* x,int n, int thread_num, int nthreads)
{
    int i,j,l;
    double tmp,b;
    double* t,*p1,*p2,h;
    int first_column, last_column;
    first_column = n * thread_num;
    first_column /= nthreads;
    last_column = n * (thread_num + 1);
    last_column = last_column / nthreads - 1;
    for(i=n-1;i>=0;i--)
    {
      tmp=a[i*n+i];
      if((fabs(tmp)> 1e+308) || (fabs(tmp)< 1e-16))
      {
          //printf("\n\n\n i = %d\n a[i*n+i] = % f\n", i, a[i*n+i]);
          printf("Matrica virozgdena\n");
          key=1;
          return ;
      }
      b=1/tmp;
      t=x+i*n;
if( (fabs(b)> 1e+200) || ((fabs(b)< 1e-100) && (fabs(b)>1e-100)) )
    {
        key = 1;
        return;
    }
      for (l=first_column;l<=last_column;l++)
      {
          t[l]=b*t[l];
      }
      p2=x+i*n;
      for(j=0;j<i;j++)
      {
         b=-a[j*n+i];
         p1=x+j*n;
         for (l=first_column;l<=last_column;l++)
         {
             h=p2[l];
             p1[l]=b*h+p1[l];
         }
       }
     }
}
void change (double* x,double* a,int n,int i,double *c,double *s, int thread_num, int nthreads)
{
    int k,p=i*n,t,j;
    double tmp1,tmp2;
    int first_column, last_column;
    first_column = n * thread_num;
    first_column /= nthreads;
    last_column = n * (thread_num + 1);
    last_column = last_column / nthreads - 1;
    for (k=first_column;k<=last_column;k++)
    {
       for(j=i+1;j<n;j++)
       {
	  t=j*n;
          tmp1=a[p+k];
          tmp2=x[p+k];
          a[p+k]=c[j]*a[p+k]+s[j]*a[t+k];
          a[t+k]=-s[j]*tmp1+c[j]*a[t+k];
	  x[p+k]=c[j]*x[p+k]+s[j]*x[t+k];
          x[t+k]=-s[j]*tmp2+c[j]*x[t+k];
       }
    }
}
void * qr_razlog_threaded(void *pa)
{
    ARGS *pargs = (ARGS*)pa;
    struct timespec begin, end;
    printf("Thread %d started\n", pargs->thread_num);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&begin);
    if(!(qr_razlog(pargs->a, pargs->x, pargs->n, pargs->thread_num, pargs->total_threads,pargs->error,pargs->cos,pargs->sin)))
    {
        printf("Unable to search the inverse matrix\n");
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&end);
    pthread_mutex_lock (&threads_total_time_mutex);
    threads_total_time += (end.tv_sec-begin.tv_sec)+(double)(end.tv_nsec-begin.tv_nsec)/1000000000;
    pthread_mutex_unlock (&threads_total_time_mutex);
    printf("Thread %d finished, time = %f\n", pargs->thread_num, ((end.tv_sec-begin.tv_sec)+(double)(end.tv_nsec-begin.tv_nsec)/1000000000));
    return 0;
}
int qr_razlog (double* a, double* x, int n, int thread_num, int nthreads,int *error,double* cosi,double* sini)
{
    int i,j;
    double p,t=0,kor;
    key = 0;
    for(i=0;i<n;i++)
    {
       p=a[i*n+i];
       cosi[i]=1;
       sini[i]=0;
       for(j=i+1;j<n;j++)
       {
            p=p*cosi[j-1]+t*sini[j-1];
	    t=a[j*n+i];
	    if(fabs(p)<1e-16&&fabs(t)<1e-16)
            {
                sini[j]=0;
                cosi[j]=1;

            }else
            {
                kor=sqrt(p*p+t*t);
                cosi[j]=p/kor;
                sini[j]=t/kor;
            }
        }
        synchronize(nthreads);
        change(x,a,n,i,cosi,sini,thread_num,nthreads);
        synchronize(nthreads);
    }
    if(thread_num==0){
    printf("QR decomposition of matrix A:\n");
    print_qr(a,x,n);}
    synchronize(nthreads);
    obratnii_xod(a,x,n,thread_num,nthreads);
    synchronize(nthreads);
    if(key)
    {
      *error=-1;
      return 0;
    }
    return 1;
}
#define MAX_N 8
void print_matrix(double* a, int n)
{
  if(MAX_N > n-2)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
               printf("%.2f ",a[i*n+j]);
            }
            printf("\n");
        }
    }
    else
    {
        for(int i = 0;i<MAX_N;i++)
        {
            for(int j = 0; j<MAX_N;j++)
            {
               printf("%.2f ",a[i*n+j]);
            }
           printf(".. %.2f\n",a[i*n + n -1]);
        }
        printf("..\n");
        for(int j = 0; j<MAX_N;j++)
        {
            printf("%.2f ",a[(n-1)*n+j]);
        }
        printf(".. %.2f\n",a[(n-1)*n + n -1]);

    }

}
double f(int i,int j,int n)
{

    /*if(i<n/2&&j<n/2)
    {
        if(j-i==1)
            return 2;
        if(i==j)
            return 1;
    }
    if(i+j==(3*n)/2-1)
      return 1;
    return 0;*/
(void) n;
return fabs(i-j);

}
void init_matrix (double* a,int n)
{
    int i,j,p;
    for(i=0;i<n;i++)
    {
       p=i*n;
       for(j=0;j<n;j++)
          a[p+j]=f(i,j,n);
    }
}
/*void fill_matrix(enum FUNC f, int n, double* a)
{
    int i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            switch(f) {
            case symmetric:
                a[i*n+j] = abs(i-j);
                break;
            case hilbert:
                a[i*n+j] = 1./(double)(i+j+1);
                break;
            case triangle:
		        if(i==j) a[i*n+j]=1;
                if(i<j) a[i*n+j]=-1;
                if(i>j) a[i*n+j]=0;
                break;
            case jordan:
                if(i>j) a[i*n+j]=0;
                if(i==j) a[i*n+j]=1;
                if(j==i+1 && j<n) a[i*n+j]= B;
                if(j>i+1) a[i*n+j]=0;
                break;
            case angle:
                if(i>j) a[i*n+j]=n-i;
                if(i<=j) a[i*n+j]=n-j;
                break;
			case errf:
				break;
              }
        }
    }
}*/
void multiplication_of_matrix (double* a,double* b,double* c,int n)
{
    int i,j,t,p;
    double sum = 0;
    for (i=0;i<n;i++)
    {
        p = i * n;
        for(j=0;j<n;j++)
        {
            for(t=0;t<n;t++)
            {
                sum += a[p+t]*b[t*n+j];
            }
            c[p+j] = sum;
            sum = 0;
        }
    }
    return;
}
double norma (double* a,int n)
{
    int i,j;
    double p,max;
    max = 0;
    for (i=0;i<n;i++)
    {
        p = 0;
        for(j=0;j<n;j++)
        {
            p += fabs(a[i*n+j]);
        }
        if (p > max){
            max = p;
        }
    }
    return max;
}
double nevyazka (double* a,double* x,double* c,int n)
{
    int i;
    double p;
    multiplication_of_matrix(a,x,c,n);
    for(i=0;i<n;i++)
    {
        c[i*n+i] = c[i*n+i]-1;
    }

    p = norma(c,n);
    return p;
}
