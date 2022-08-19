#ifndef IO_H
#define IO_H

#include "alloc.h"

//////////////////////////////////////////////////////////
//			MACROS				//
//////////////////////////////////////////////////////////

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

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/*
#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif

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

#define ISIZE sizeof(int)

#define FSIZE sizeof(float)

#define DSIZE sizeof(double)

#define SQR(x) ((x))*((x))

//////////////////////////////////////////////////////////////////
//			SEISMIC FUNCTIONS			//
//////////////////////////////////////////////////////////////////

double interp_trace(float, int, float, float, float *);

void get_vel_model(float*, FILE*, size_t, size_t, size_t);

void get_tr(int, su_trace*, int, FILE*);
void put_tr(int, su_trace*, int, FILE*);

//////////////////////////////////////////////////////////////////
//		MATH FUNCTIONS AND SEARCHES			//
//////////////////////////////////////////////////////////////////

void get_min_max(float*,float*,float*, size_t, size_t, size_t, size_t);

float sinc(float);

int max(int, int);
int min(int, int);

float max_f(float, float);
float min_f(float, float);

int maxval(int *, int);
int minval(int *, int);

//////////////////////////////////////////////////////////////////
//			COPY FUNCTION				//
//////////////////////////////////////////////////////////////////

void copy1float(float *, float *, size_t);

#endif
