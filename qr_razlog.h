#ifndef QR_RAZLOG_H
#define QR_RAZLOG_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


typedef struct _ARGS
{
    double* a;
    double* copy_of_a;
    double* x;
    double* sin;
    double* cos;
    int n;
    int thread_num;
    int total_threads;
    int *error;
} ARGS;

void init_matrix (double* a,int n);
double f(int i,int j,int n);
void matrix_E(double* a, int n);
int qr_razlog (double* a,double* x,int n,int thread_num, int nthreads,int *error,double* cosi,double* sini);
void obratnii_xod (double* a,double* x,int n,int thread_num, int nthreads);
void print_matrix(double* a,int n);
void change (double* x,double* a,int n,int i,double *c,double *s, int thread_num, int nthreads);
void multiplication_of_matrix (double* a,double* b,double* c,int n);
double norma (double* a,int n);
double nevyazka (double* a,double* x,double* c,int n);
void transponirovanie(double* a,int n);
void * qr_razlog_threaded(void *pa);
void print_qr (double* a,double* x,int n);


#endif
