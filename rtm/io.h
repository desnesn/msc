#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "segy.h"

/**********************************************************/
/*********************** MACROS ***************************/
/**********************************************************/

#ifndef NULL
#define NULL	((void *)0)
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif

#ifndef SEEK_SET
#define SEEK_SET (0)
#endif

#ifndef SEEK_CUR
#define SEEK_CUR (1)
#endif

#ifndef SEEK_END
#define SEEK_END (2)
#endif

#ifndef PI
#define PI (3.141592653589793)
#endif

#ifndef D_PI             
#define D_PI 3.14159265358979323846264338328
#endif

#ifndef GOLDEN_RATIO 
#define GOLDEN_RATIO (1.618034)   /* the golden ratio */
#endif

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

#ifndef YES
#define YES (1)
#endif

#ifndef NO
#define NO (0)
#endif

#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif

/*
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
*/
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#define CLOSETO(x, y) ((ABS((x) - (y)) <= FLT_EPSILON*ABS(y))?cwp_true:cwp_false)

#define ISODD(n) ((n) & 01)

#define ISIZE sizeof(int)

#define FSIZE sizeof(float)

#define DSIZE sizeof(double)

#define	STREQ(s,t) (strcmp(s,t) == 0)

#define	STRLT(s,t) (strcmp(s,t) < 0)

#define	STRGT(s,t) (strcmp(s,t) > 0)

#define	DIM(a) (sizeof(a)/sizeof(a[0]))

#define SQR(x) ((x))*((x))

float sinc(float);

double interp_trace(float, int, float, float, float *);

void get_vel_model(double**, FILE*, size_t, size_t, size_t, size_t, size_t);

int max(int, int);
int min(int, int);

float max_f(float, float);
float min_f(float, float);

//int ABS(int);

int maxval(int *, int);
int minval(int *, int);

/**********************************************************/
/********************* TRACES IO's ************************/
/**********************************************************/

int get_tr(int, su_header*, float*, int, FILE*);

void put_tr(int, su_header*, float*, int, FILE*);

void get_sx_gx(su_header*, int*, int*);

//void get_vel_model(float**, FILE *, size_t, size_t);

void get_min_max(double**,float*,float*, size_t, size_t, size_t, size_t);

/**********************************************************/
/**************COPY & PRINT FUNCTIONS *********************/
/**********************************************************/

void copy2fdata(float **, size_t, float *, size_t);

void copy1float(float *, float *, size_t);

void print1int(int*, int);

void print1float(float*, int);

void print2int(int**, int, int);

void print2float(float**, int, int);

/**********************************************************/
/********************ALOCATE POINTERS**********************/
/**********************************************************/

/* allocate a 1-d array of ints - zerado */
int *alloc1int(size_t);

/* allocate a 2-d array of ints - zerado*/
int **alloc2int(size_t n_lin, size_t n_col);

/* allocate a 3-d array of ints - zerado*/
int ***alloc3int(size_t, size_t, size_t);

/* allocate a 1-d array of floats - zerado */
float *alloc1float(size_t);

/* allocate a 2-d array of floats - zerado*/
float **alloc2float(size_t, size_t);

/* allocate a 3-d array of floats - zerado*/
float ***alloc3float(size_t, size_t, size_t);

/* allocate a 1-d array of double - zerado */
double *alloc1double(size_t);

/* allocate a 2-d array of double - zerado*/
double **alloc2double(size_t, size_t);

/* allocate a 3-d array of double - zerado*/
double ***alloc3double(size_t, size_t, size_t);

/* allocate a 1-d array of su_trace - zerado*/
su_trace* alloc1su_trace(int);

/**********************************************************/
/**********************FREE POINTERS***********************/
/**********************************************************/

/* free a 1-d array of ints */
void free1int(int *);

/* free a 2-d array of ints */
void free2int(int **, size_t, size_t);

/* free a 3-d array of int */
void free3int(int ***, size_t, size_t, size_t);

/* free a 1-d array of floats */
void free1float(float *);

/* free a 2-d array of floats */
void free2float(float **, size_t, size_t);

/* free a 3-d array of floats */
void free3float(float ***, size_t, size_t, size_t);

/* free a 1-d array of double */
void free1double(double *vet);

/* free a 2-d array of double */
void free2double(double **mat, size_t, size_t);

/* free a 3-d array of double */
void free3double(double ***mat, size_t, size_t, size_t);

/* free a 1-d array of su_trace */
void free1su_trace(su_trace*);

#endif
