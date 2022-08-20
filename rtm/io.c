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

		//omp_set_num_threads(24);
		//#pragma omp parallel for default(none) private(it,aux1) shared(itb,ite,aux0,trace) reduction(+: interp_trace)
		for(it=itb-1;it<ite;it++)
		{
			aux1 = aux0 - ((float)it);
			interp_trace = interp_trace + trace[it] * sinc(aux1);
			//interp_trace = trace[it] * sinc(aux1);
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

	//omp_set_num_threads(24);
	//#pragma omp parallel for default(none) private(ix) shared(vel,nz,nx,nzz,nborda,vel_rd)
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

	//condicao de corrida de dados
	//omp_set_num_threads(24);
	//#pragma omp parallel for default(none) private(ix,iz) shared(ixb,ixe,izb,ize,mat,vmax,vmin,nzz)
	for(ix=ixb-1;ix<ixe;++ix)
    	{
        	for(iz=izb-1;iz<ize;++iz)
        	{
			//#pragma omp critical
			//{
				*vmax = (mat[(iz) + (ix*nzz)]>*vmax) ? mat[(iz) + (ix*nzz)]:*vmax;    	
			//}

			//#pragma omp critical
			//{			
				*vmin = (mat[(iz) + (ix*nzz)]<*vmin) ? mat[(iz) + (ix*nzz)]:*vmin;
			//}

			//*vmax = max_f(mat[(iz) + (ix*nzz)], *vmax);
			//*vmin = min_f(mat[(iz) + (ix*nzz)], *vmin);
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
/*
void copy1float(float *v1, float *v2, size_t size)
{
    memcpy(v1,v2,sizeof(float)*size);
}
*/

//////////////////////////////////////////////////////////////////
//			TIME FUNCTION				//
//////////////////////////////////////////////////////////////////

/* Return 1 if the difference is negative, otherwise 0.  */
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

/*
void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;

    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}

int main()
{
    struct timeval tvBegin, tvEnd, tvDiff;

    // begin
    gettimeofday(&tvBegin, NULL);
    timeval_print(&tvBegin);

    // lengthy operation
    int i,j;
    for(i=0;i<999999L;++i) {
        j=sqrt(i);
    }

    //end
    gettimeofday(&tvEnd, NULL);
    timeval_print(&tvEnd);

    // diff
    timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
    printf("%ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

    return 0;
}

*/


/*
int get_exec_time(struct timeval start,struct timeval finish)
{
	int msec;
	msec = finish.tv_sec * 1000 + finish.tv_usec / 1000;
	msec-= start.tv_sec * 1000 + start.tv_usec / 1000;
	return msec;
}

int get_exec_time_sec(struct timeval start,struct timeval finish)
{
	int msec;
	msec = finish.tv_sec * 1000 + finish.tv_usec / 1000;
	msec-= start.tv_sec * 1000 + start.tv_usec / 1000;
	return msec;
}
*/


