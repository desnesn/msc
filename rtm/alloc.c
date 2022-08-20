#include "alloc.h"

//////////////////////////////////////////////////
//		ALLOC - INT			//
//////////////////////////////////////////////////

//	allocate a 1-d array of ints		//
int* alloc1int(size_t length)
{
    	return  (int*) malloc(sizeof(int)*length);
}

//	allocate a 2-d array of ints		//
int* alloc2int(size_t width, size_t length)
{
    	return  (int*) malloc(sizeof(int)*width*length);
}

//////////////////////////////////////////////////
//		ALLOC - FLOAT			//
//////////////////////////////////////////////////

//	allocate a 1-d array of float		//
float* alloc1float(size_t length)
{
    	return  (float*) malloc(sizeof(float)*length);
}

// 	allocate a 2-d array of float		//
float* alloc2float(size_t width, size_t length)
{
    	return  (float*) malloc(sizeof(float)*width*length);
}

// 	allocate a 2-d array of float		//
float* alloc3float(size_t width, size_t length, size_t dim)
{
    	return  (float*) malloc(sizeof(float)*width*length*dim);
}

//////////////////////////////////////////////////
//		ALLOC - SU_TRACE		//
//////////////////////////////////////////////////

//	allocate a 1-d array of SU_TRACE	//
su_trace* alloc1su_trace(int length)
{
	return (su_trace*) malloc(sizeof(su_trace) * length);
}

//////////////////////////////////////////////////
//		DEALLOC - INT			//
//////////////////////////////////////////////////

//	free a 1-d array of ints 		//
void free1int(int *vet) 
{
    free(vet);
}

// 	free a 2-d array of ints 		//
void free2int(int *vet)
{
    free(vet);
}

//////////////////////////////////////////////////
//		DEALLOC - FLOAT			//
//////////////////////////////////////////////////

// 	free a 1-d array of floats 		//
void free1float(float *vet)
{
	free(vet);
}

// 	free a 2-d array of floats 		//
void free2float(float *vet)
{
	free(vet);
}

// 	free a 3-d array of floats 		//
void free3float(float *vet)
{
	free(vet);
}

//////////////////////////////////////////////////
//		DEALLOC - SU_TRACE		//
//////////////////////////////////////////////////

// 	free a 1-d array of su_trace 		//
void free1su_trace(su_trace *vet)
{
	free(vet);
}

