#include "io.h"

//const double _pi = 3.14159265358979323846264338328;

float sinc(float x)
{
	float sinc_v;
	sinc_v = ( fabs(x) > 0.0 ) ? sin(D_PI*x) / (D_PI*x) : 1.0;
/*
	if ( fabs(x) > 0.0 ) 
		sinc_v = ( sin( x*_pi ) / (_pi*x) );
	else
		sinc_v = 1.0;
*/
	return sinc_v;
}

//nt= number of samples ! t = target sample; ! t0= amostra inicial ! dt= intervalo de amostragem original //
// TODO - verificar esta funcao //
double interp_trace(float t, int nt, float t0, float dt, float *trace)
{		//ttotal-t,ns,0.0,dtrec,shotgather[itrc].tr_data
	float interp_trace;
	float tol = 1.0e-3;
	int nlag = 10;	// ???
	int it,itb,ite;
	float aux1, aux0;

	aux0 = (t-t0)/dt;
	//printf("dt: %lf\n",dt);
	//printf("aux0: %lf\n",aux0);				//724

	it  = ( (int)floor(aux0) ) + 1;
//	printf("it: %d\t",it);				//725

	if ( fabs(aux0-(it-1)) < tol )
     		interp_trace = trace[it];
	// ???? else ????
	interp_trace = 0.0;
	
	//printf("(t0+(nt-1)*dt): %lf\n",(t0+(nt-1)*dt));	//2.896

	if ( (t0 <= t) && (t <= (t0 + nt*dt) ) )//-1)*dt) ) )
	{
		itb = max(it-nlag,0); //   itb = max(it-nlag,1)
		ite = min(it+nlag,nt);//   ite = min(it+nlag,nt)

		//printf("itb[%d]\tite[%d]\n",itb,ite);

		for(it=itb-1;it<ite;it++)
		{
			aux1 = aux0 - ((float)it);//      aux1 = aux0 - real(it-1)
			//printf("aux0(%lf) - it(%d) - aux1(%lf) - sinc(aux1): %.25lf\n", aux0,it,aux1,sinc(aux1));

			interp_trace = interp_trace + trace[it] * sinc(aux1); // melhor analise pendente
		//	printf("sinc(aux1): %lf\n", sinc(aux1));
			
   		}
	}

	//printf("interp_trace: %lf\n",interp_trace);

	return ((double)interp_trace);
}

//////////////////////////
//	nx: 2301	//
//	nz: 751		//
//	nxx: 2401	//
//	nzz: 801	//
//	nborda: 50	//
//	ixb: 51		//
//	ixe: 2351	//
//	izb: 1		//
//	ize: 751	//
//////////////////////////

void get_vel_model(double** vel, FILE *vel_file, size_t ixb, size_t ixe, size_t izb, size_t ize, size_t nborda)
{
	unsigned ix,iz,nz,ntrace=1;
    
	float *ptr;

    	nz=ize-nborda;
    	ptr = (float*) calloc(sizeof(float),nz);

	for(ix=ixb-1;ix<ixe;ix++)
    	{
        	fread((float*) ptr, sizeof(float), nz, vel_file);
        	for(iz=izb-1;iz<ize;iz++)
        	{
            		vel[iz][ix]=((double)ptr[iz-nborda]);		// otimizar com memcpy
        	}
        	fseek ( vel_file , nz * sizeof(float) * ntrace , SEEK_SET);
        	ntrace++;
    	}
}

void get_min_max(double** mat,float *vmin,float *vmax, size_t ixb, size_t ixe, size_t izb, size_t ize)
{
    unsigned ix,iz;

    *vmin=(float)(mat[izb-1][ixb-1]);
    *vmax=(float)(mat[izb-1][ixb-1]);

    for(ix=ixb-1;ix<ixe;ix++)
    {
        for(iz=izb-1;iz<ize;iz++)
        {
            *vmax = max_f((float)(mat[iz][ix]), *vmax);    
            *vmin = min_f((float)(mat[iz][ix]), *vmin);
        }
   } 
}

int max(int a, int b) { return (a>b)? a:b; }
int min(int a, int b) { return (a<b)? a:b; }

float max_f(float a, float b) {	return (a>b)? a:b; }
float min_f(float a, float b) { return (a<b)? a:b; }

/*
int ABS(int x)
{
	return (x<0)?(-1)*(x):(x);
}
*/

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

/**********************************************************/
/********************* TRACES IO's ************************/
/**********************************************************/

int get_tr(int tr_index, su_header *tr_hdr, float *tr_data, int ns, FILE *su_file)
{

	fseek(su_file, tr_index * ( (sizeof(float)*ns) + HDRBYTES ), SEEK_SET);
	
	int nread;

	void *ptr;
	ptr = (su_header*) malloc( sizeof(su_header) );

	nread = fread((su_header*)ptr, HDRBYTES, 1, su_file);

	switch(nread)
    	{
        	case 0:
            		free(ptr);
            		return 0;
            	break;

		default:
            		memcpy(tr_hdr, ptr, sizeof(su_header));    
        
			free(ptr);
            		ptr = (float*) malloc ( sizeof(float) * ns );
 
			fread((float*) ptr, sizeof(float) , ns, su_file);
            		memcpy(tr_data, ptr, sizeof(float) * ns);
   
            free(ptr);
            return nread;
    }
}

void put_tr(int tr_index, su_header *tr_hdr, float *tr_data, int ns, FILE *su_file)
{
	fseek(su_file, tr_index * ( (sizeof(float)*ns) + HDRBYTES ), SEEK_SET);
	
	fwrite(tr_hdr, HDRBYTES, 1, su_file);
	fwrite(tr_data, sizeof(float), ns, su_file);
}

void get_sx_gx(su_header *tr, int *sx, int *gx)
{
    *sx = tr->sx;
    *gx = tr->gx;
}

/*
void get_vel_model(float** vel, FILE *vel_file, size_t n_lin, size_t n_col)
{
    unsigned i,j,ntrace=1;
    
    float *ptr;
    ptr = (float*) calloc(sizeof(float),n_col);

    for(j=0;j<n_col;j++)
    {
        fread((float*) ptr, sizeof(float), n_col, vel_file);
        for(i=0;i<n_lin;i++)
        {
            vel[i][j]=ptr[i];
        }
        fseek ( vel_file , 751 * sizeof(float) * ntrace , SEEK_SET );

        ntrace++;
    }
    fseek ( vel_file , 0 , SEEK_SET );
}


/**********************************************************/
/********************* TRACES IO's ************************/
/**********************************************************/

/**********************************************************/
/**************COPY & PRINT FUNCTIONS *********************/
/**********************************************************/

// copia um vetor em uma coluna de uma matrix //
void copy2fdata(float **shot, size_t column, float *tr_data, size_t ns)
{
    unsigned i;
    for(i=0;i<ns;i++)
        shot[i][column]=tr_data[i];
}

void copy1float(float *v1, float *v2, size_t size)
{
    unsigned i;
    for(i=0;i<size;i++) v1[i] = v2[i];
}

void print1int(int* vet, int length)
{
    unsigned j;
    for(j=0;j<length;j++) printf("%d\n",vet[j]);
}

void print2int(int** mat, int n_lin, int n_col)
{
    unsigned i,j;
    for(i=0;i<n_lin;i++)
    {
        for(j=0;j<n_col;j++) 
            printf("%d\t",mat[i][j]);
        printf("\n");
    }
}

void print1float(float* vet, int length)
{
    unsigned j;
    for(j=0;j<length;j++) printf("%.5f\n",vet[j]);
}

void print2float(float** mat, int n_lin, int n_col)
{
    unsigned i,j;
    for(i=0;i<n_lin;i++)
    {
        for(j=0;j<n_col;j++) 
            printf("%.5f\t",mat[i][j]);
        printf("\n");
    }
}

/**********************************************************/
/**************COPY & PRINT FUNCTIONS *********************/
/**********************************************************/

/**********************************************************/
/********************ALOCATE POINTERS**********************/
/**********************************************************/

/* allocate a 1-d array of ints - zerado */
int *alloc1int(size_t length)
{
	int *vet;
    	vet = (int*) calloc(sizeof(int), length);
    	return vet;
}

/* allocate a 2-d array of ints - zerado*/
int **alloc2int(size_t n_lin, size_t n_col)
{
    	unsigned int i;
    	int **mat = (int**) malloc(sizeof(int*)*n_lin);
    	for(i=0;i<n_lin;i++)
       		mat[i] = (int*) calloc(sizeof(int),n_col);
	return mat;
}

/* allocate a 3-d array of ints - zerado*/
int ***alloc3int(size_t x, size_t y, size_t t)
{
	unsigned int i,j;
    	int ***mat = (int***) malloc(sizeof(int**) * x);
    	for(i=0;i<x;i++)
    	{
        	mat[i] = (int**) malloc(sizeof(int*) * y);
        	for(j=0;j<y;j++)
        	{
            		mat[i][j] = (int*) calloc(sizeof(int), t);
        	}
    	}
    	return mat;
}

/* allocate a 1-d array of floats - zerado */
float *alloc1float(size_t length)
{
	float *vet;
    	vet = (float*) calloc(sizeof(float), length);
	return vet;
}

/* allocate a 2-d array of floats - zerado*/
float **alloc2float(size_t n_lin, size_t n_col)
{
	unsigned int i,j;
    	float **mat = (float**) malloc(sizeof(float*) * n_lin);
    	for(i=0;i<n_lin;i++)
       		mat[i] = (float*) calloc(sizeof(float), n_col);
    	return mat;
}

/* allocate a 3-d array of floats - zerado*/
float ***alloc3float(size_t x, size_t y, size_t t)
{
	unsigned int i,j;
    	float ***mat = (float***) malloc(sizeof(float**) * x);
    	for(i=0;i<x;i++)
    	{
        	mat[i] = (float**) malloc(sizeof(float*) * y);
        	for(j=0;j<y;j++)
        	{
            		mat[i][j] = (float*) calloc(sizeof(float), t);
        	}
    	}
    	return mat;
}

/* allocate a 1-d array of double - zerado */
double *alloc1double(size_t length)
{
	double *vet;
    	vet = (double*) calloc(sizeof(double), length);
	return vet;
}

/* allocate a 2-d array of double - zerado*/
double **alloc2double(size_t n_lin, size_t n_col)
{
	unsigned int i,j;
    	double **mat = (double**) malloc(sizeof(double*) * n_lin);
    	for(i=0;i<n_lin;i++)
       		mat[i] = (double*) calloc(sizeof(double), n_col);
    	return mat;
}

/* allocate a 3-d array of double - zerado*/
double ***alloc3double(size_t x, size_t y, size_t t)
{
	unsigned int i,j;
    	double ***mat = (double***) malloc(sizeof(double**) * x);
    	for(i=0;i<x;i++)
    	{
        	mat[i] = (double**) malloc(sizeof(double*) * y);
        	for(j=0;j<y;j++)
        	{
            		mat[i][j] = (double*) calloc(sizeof(double), t);
        	}
    	}
    	return mat;
}

su_trace* alloc1su_trace(int length)
{
	su_trace *vet;
    	vet = (su_trace*) calloc(sizeof(su_trace), length);
	return vet;
}


/**********************************************************/
/********************ALOCATE POINTERS**********************/
/**********************************************************/

/**********************************************************/
/**********************FREE POINTERS***********************/
/**********************************************************/

/* free a 1-d array of ints */
void free1int(int *vet)
{
    free(vet);
}

/* free a 2-d array of ints */
void free2int(int **mat, size_t n_lin,size_t n_col)
{
    int i;
    for(i=(n_lin-1); i>=0 ; i--)
        free(mat[i]);
    free(mat);
}

/* free a 3-d array of int */
void free3int(int ***mat, size_t x, size_t y, size_t t)
{
    int i,j;
    for(i=0;i<x;i++)
    {
        for(j=0;j<y;j++)
	    free(mat[i][j]);
	free(mat[i]);
    }
    free(mat);
}

/* free a 1-d array of floats */
void free1float(float *vet)
{
	free(vet);
}

/* free a 2-d array of floats */
void free2float(float **mat, size_t n_lin, size_t n_col)
{
    int i;
    for(i=(n_lin-1); i>=0 ; i--)
        free(mat[i]);
    free(mat);
}

/* free a 3-d array of floats */
void free3float(float ***mat, size_t x, size_t y, size_t t)
{

    int i,j;
    for(i=0;i<x;i++)
    {
        for(j=0;j<y;j++)
	    free(mat[i][j]);
	free(mat[i]);
    }
    free(mat);
}

/* free a 1-d array of double */
void free1double(double *vet)
{
	free(vet);
}

/* free a 2-d array of double */
void free2double(double **mat, size_t n_lin, size_t n_col)
{
    int i;
    for(i=(n_lin-1); i>=0 ; i--)
        free(mat[i]);
    free(mat);
}

/* free a 3-d array of double */
void free3double(double ***mat, size_t x, size_t y, size_t t)
{
    int i,j;
    for(i=0;i<x;i++)
    {
        for(j=0;j<y;j++)
	    free(mat[i][j]);
	free(mat[i]);
    }
    free(mat);
}

/* free a 1-d array of su_trace */
void free1su_trace(su_trace *vet)
{
	free(vet);
}

/**********************************************************/
/**********************FREE POINTERS***********************/
/**********************************************************/
