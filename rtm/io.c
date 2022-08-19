#include "io.h"

//////////////////////////////////////////////////////////////////
//			SEISMIC FUNCTIONS			//
//////////////////////////////////////////////////////////////////

//nt= number of samples ! t = target sample; ! t0= amostra inicial ! dt= intervalo de amostragem original //
double interp_trace(float t, int nt, float t0, float dt, float *trace)
{
	float interp_trace = 0.0;
	//float tol = 1.0e-3;
	int nlag = 10;
	int it,itb,ite;
	float aux1, aux0;

	aux0 = (t-t0)/dt;

	it  = ( (int)floor(aux0) ) + 1;

	if ( (t0 <= t) && (t <= (t0 + nt*dt) ) )
	{
		itb = max(it-nlag,0); 
		ite = min(it+nlag,nt);

		for(it=itb-1;it<ite;it++)
		{
			aux1 = aux0 - ((float)it);
			interp_trace = interp_trace + trace[it] * sinc(aux1);
   		}
	}

	return ((double)interp_trace);
}

void get_vel_model(float* vel, FILE *vel_file, size_t nz, size_t nx, size_t nborda)
{
	unsigned int ix;
	unsigned int nzz = nz + 2 * nborda;

	float *vel_rd = (float*) malloc(nz * nx * sizeof(float));

       	fread((float*) vel_rd, sizeof(float), nz * nx, vel_file);

	//TODO: calcular min e max aqui

	for(ix=0;ix<nx;++ix)
		memcpy(vel + (ix*nzz+nzz*nborda) + nborda, vel_rd + (ix*nz), nz*sizeof(float) );
}

void get_tr(int tr_index, su_trace *tr, int ns, FILE *su_file)
{
	fseek(su_file, tr_index * ( (sizeof(float)*ns) + HDRBYTES ), SEEK_SET);

	fread(tr->tr_header, HDRBYTES, 1, su_file);	
	fread(tr->tr_data, sizeof(float), ns, su_file);	
}

void put_tr(int tr_index, su_trace *tr, int ns, FILE *su_file)
{
	fseek(su_file, tr_index * ( (sizeof(float)*ns) + HDRBYTES ), SEEK_SET);
	
	fwrite(tr->tr_header, HDRBYTES, 1, su_file);
	fwrite(tr->tr_data, sizeof(float), ns, su_file);
}

//////////////////////////////////////////////////////////////////
//		MATH FUNCTIONS AND SEARCHES			//
//////////////////////////////////////////////////////////////////

void get_min_max(float* mat,float *vmin, float *vmax, size_t ixb, size_t ixe, size_t izb, size_t ize)
{
	unsigned ix,iz;

	unsigned int nzz = ize + izb - 1;

	*vmin = mat[(izb-1) + ((ixb-1)*(nzz))];
    	*vmax = mat[(izb-1) + ((ixb-1)*(nzz))];

	for(ix=ixb-1;ix<ixe;++ix)
    	{
        	for(iz=izb-1;iz<ize;++iz)
        	{
        	    	*vmax = max_f(mat[(iz) + (ix*nzz)], *vmax);    
        	    	*vmin = min_f(mat[(iz) + (ix*nzz)], *vmin);
        	}
   	}
}

float sinc(float x)
{
	return ( fabs(x) > 0.0 ) ? sin(D_PI*x) / (D_PI*x) : 1.0;
}

int max(int a, int b) 
{ 
	return (a>b)? a:b; 
}

int min(int a, int b) 
{ 
	return (a<b)? a:b; 
}

float max_f(float a, float b) 
{	
	return (a>b)? a:b;
}

float min_f(float a, float b) 
{ 
	return (a<b)? a:b; 
}

int minval(int *v, int size)
{
	int i, min;

	min=v[0];
	for(i=1;i<size;i++)
		if (v[i]<min) min=v[i];

	return min;
}

int maxval(int *v, int size)
{
	int i, max;

	max=v[0];
	for(i=1;i<size;i++)
		if (v[i]>max) max=v[i];

	return max;
}

//////////////////////////////////////////////////////////////////
//			COPY FUNCTION				//
//////////////////////////////////////////////////////////////////

void copy1float(float *v1, float *v2, size_t size)
{
    memcpy(v1,v2,sizeof(float)*size);
}

