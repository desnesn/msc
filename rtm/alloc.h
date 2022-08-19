#ifndef ALLOC_H
#define ALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "segy.h"

#include "memcpy.h"

//////////////////////////////////////////////////
//		ALLOC - INT			//
//////////////////////////////////////////////////

//	allocate a 1-d array of ints		//
int *alloc1int(size_t);

//	allocate a 2-d array of ints		//
int *alloc2int(size_t n_lin, size_t n_col);

//	allocate a 3-d array of ints		//
int *alloc3int(size_t, size_t, size_t);

//////////////////////////////////////////////////
//		ALLOC - FLOAT			//
//////////////////////////////////////////////////

//	allocate a 1-d array of float		//
float *alloc1float(size_t);

// 	allocate a 2-d array of float		//
float *alloc2float(size_t, size_t);

//	allocate a 3-d array of ints		//
float *alloc3float(size_t, size_t, size_t);

//////////////////////////////////////////////////
//		ALLOC - SU_TRACE		//
//////////////////////////////////////////////////

//	allocate a 1-d array of SU_TRACE	//
su_trace* alloc1su_trace(int);

//////////////////////////////////////////////////
//		DEALLOC - INT			//
//////////////////////////////////////////////////

//	free a 1-d array of ints 		//
void free1int(int *);

// 	free a 2-d array of ints 		//
void free2int(int *);

// 	free a 3-d array of ints 		//
void free3int(int *);

//////////////////////////////////////////////////
//		DEALLOC - FLOAT			//
//////////////////////////////////////////////////

// 	free a 1-d array of floats 		//
void free1float(float *);

// 	free a 2-d array of floats 		//
void free2float(float *);

// 	free a 3-d array of floats 		//
void free3float(float *);

//////////////////////////////////////////////////
//		DEALLOC - SU_TRACE		//
//////////////////////////////////////////////////

// 	free a 1-d array of su_trace 		//
void free1su_trace(su_trace*);

#endif

