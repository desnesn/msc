/********************************************************************************/
/*	RTM Migration 2.0 - Baseado no codigo em Fortran do Prof. Jesse Costa  	*/
/*	Autor: Desnes A N do Rosário			(jesse@ufpa.br)		*/
/*	e-mail: username@gmail.com						*/
/*										*/
/*	FD segunda ordem no tempo-oitava ordem no espaco			*/
/*										*/
/*	pulso fonte : Ricker 							*/
/*	f(t) = exp(-(pi*freq*(t-t0)^2)*(1-2.0*(pi*freq*(t-t0)^2)		*/
/*										*/
/*  	interpolacao de tracos usando funcao sinc 				*/
/*										*/
/*  	condicoes PML (Perfect Matched Layer) para absorcao na fronteira	*/
/*  	Nborda=50 largura da borda de absorcao 					*/
/*  										*/
/*  	condicao de imagem: cross-correlacao					*/
/*										*/
/*  	filtro laplaciano aplicado a imagem migrada				*/
/*										*/
/*  	campo retropropagado armazenado para cada tiro 				*/
/*  	armazenado no arquivo teporario fbackward.bin				*/
/*										*/
/*	campo propagado forward no tempo pode ser armazenado 			*/
/*  	no arquivo forward.bin fazero o parametro MVFWRD /= 0			*/
/*	observe que isso aumenta o tempo e custo de armazenamento		*/	
/*										*/
/********************************************************************************/
/*	dtrtm  = 0.008 ms  intervalo de gravacao dos frames para rtm		*/
/********************************************************************************/
#include <string.h>


#include <stdio.h>
#include <stdlib.h>
//#include <sys/time.h>
#include <omp.h>

#include "io.h"

struct timeval start, finish, diff;
struct timeval start_princ, finish_princ, diff_princ;

//int get_exec_time(struct timeval,struct timeval);
//int get_exec_time_sec(struct timeval,struct timeval);

#define UNUSED __attribute__ ((unused))

int main(int argc, char *argv[])
{
	//////////////////////////////////////////////////////////////////////////////////
	//				PRE-INIT					//
	//////////////////////////////////////////////////////////////////////////////////    

	if (argc!=3) { printf("Erro nas opcoes de entrada\n"); return -1; }

	int n_threads = atoi(argv[2]);

	//	g - modelo grande		//
	//	p - modelo pequeno		//
	char modelo = (argv[1])[1];

	//	0 - programa nao paralelo	//
	//	1 - programa paralelo		//
	int UNUSED is_parallel;
	is_parallel = (n_threads==1) ? 0 : 1;

	//////////////////////////////////////////////////////////////////////////////////
	//				RTM PRE-INIT					//
	//////////////////////////////////////////////////////////////////////////////////    
	
	printf("|---------------------------------------------------------------------------|\n");
	printf("|                     Reverse Time Migration - RTM 4.0                      |\n");
	printf("|---------------------------------------------------------------------------|\n");

	//	largura da borda da ordem de 3xlambda/2, lambda=Vmax/freq		//
	const unsigned int UNUSED nborda = 50;
	
	//	LFRAC: {0,1,2,3,4}							//
	//	0 - migra cada tiro usando toda extensao do modelo de velocidade	//
	//	caso contrario a janela de migracao limitada ao intervalo do modelo de	//
	//	velocidade: lanco/NFRAC+lanco+lanco/NFRAC				//
	const unsigned int UNUSED LFRAC = 3;

	// intervalo para gravacao dos frames na migracao RTM //
	const float UNUSED dtrtm = 0.008;// seconds

	float UNUSED *gama_x;
	float UNUSED *gama_z;
	float UNUSED *window_x;

	float UNUSED *vel;

	float UNUSED *p_1;
	float UNUSED *p_2;
	float UNUSED *p_aux;

	float UNUSED dswap;

	float UNUSED *psrc;
	float UNUSED *prcv;
	float UNUSED *imag;
	float UNUSED *imag_filter;

	//	Operador de segunda derivada de oitava ordem				//
	//	serve para resolver numericamente a equacao da onda (EDP) 		//
	const unsigned int UNUSED DRVLEN = 5;
	float UNUSED deriv2[5] = {-2.84722E+00,1.60000E+00,-2.00000E-01,2.53968E-02,-1.785714E-03};

	//	Blackmann window 							//
	//	utiliza um taper nas bordas de cada imagem suavizando as derivadas	//
	// 	janela Blackmann harris calculada para decaimento suave 		//
	//	nas bordas das derivadas						//
	const unsigned int UNUSED WDWLEN = 5;
	float UNUSED blackmann[5] = {0.893012701892219E+00,0.63000E+00,0.34000E+00,0.130000E+00,0.026987298107781E+00};
	
	su_trace UNUSED *shotgather;

	//	0 - nao grava o campo emitido	//
	//	1 - grava o campo emitido
	int UNUSED write_forward_file = 1;//0;

	int UNUSED reclen;
	int UNUSED nrec, ilanco, nmigtrc, ixmig1, ixmig0, ixx0, ixx1;

	const float UNUSED sec2mcs = 1.0e6;

	float UNUSED beta, delt, gama, invpgama, mgama, laplacian, fonte, source;
	float UNUSED ttotal, t, vmin, vmax, dxmax;
	float UNUSED dx, dt, dz, freq;
	int UNUSED nx, nz, nt, z0, x0, it;
	float UNUSED dtrec;
	int UNUSED ndtrec, ndtrtm;
	int UNUSED ix, iz, it0;
	int UNUSED isx, isz;
	int UNUSED *igx , *igz;
	int UNUSED *shotidx;
	int UNUSED iconv, ixb, ixe, izb, ize, nxx, nzz;
	int UNUSED iframe, nframes, irec, ircv, itrc, ns;
	int UNUSED ntrc, ngeophones, ishot, nshots;
	float UNUSED factor;

	if(modelo=='p')
	{
		//	Dados do modelo de velocidade-PEQUENO					//TODO: argc e argv
    		nz = 243;     //n1
    		nx = 767;     //n2
    		z0 = 0.0;     //??
    		x0 = 0.0;     //??
    		dz = 12.0;    //d1
    		dx = 12.0;    //d2
		//	Dados do modelo de velocidade-PEQUENO 						//TODO: argc e argv
	}
	else if(modelo=='g')
	{
		//	Dados do modelo de velocidade-GRANDE					//TODO: argc e argv
    		nz = 751;    //n1
    		nx = 2301;   //n2
    		z0 = 0.0;    //??
    		x0 = 0.0;    //??
    		dz = 4.0;    //d1
    		dx = 4.0;    //d2
		//	Dados do modelo de velocidade-GRANDE 						//TODO: argc e argv
	}
	else{ printf("Erro na opcao do modelo de velocidades\n"); return -1; }

	// 	Frequencia da fonte 							//
	freq = 12.0;

	// 	PML nos 4 lados do grid - Atenuar a energia das bordas com uma 		//
	//	exponencial quadratica da velocidade  - Bordas para x e z no modelo	//
	nxx = nx+2*nborda;
	nzz = nz+2*nborda;
	ixb = nborda+1;
	ixe = nborda+nx;
	izb = nborda+1;
	ize = nborda+nz;

	printf("\tInitial parameters\n");
	printf("\tfreq: %.3f   ",freq);    //freq: 12
	printf("nz: %d   ",nz);          //nz: 751
	printf("nx: %d   ",nx);          //nx: 2301
	printf("nzz: %d   ",nzz);        //nzz: 801
	printf("nxx: %d\n",nxx);        //nxx: 2401
	printf("\tnborda: %d   ",nborda);  //nborda: 50
	printf("ixb: %d   ",ixb);        //ixb: 51
	printf("ixe: %d   ",ixe);        //ixe: 2351
	printf("izb: %d   ",izb);        //izb: 1
	printf("ize: %d\n",ize);        //ize: 751

	printf("|---------------------------------------------------------------------------|\n");

	//	Alocacao de vetores e matrizes 						//

	//p = alloc3float(nzz,nxx,2);    // Campos de pressão: p[851][2401][2]
	//if(p==NULL) { printf("Allocation of p[%d,%d,%d] failed!!!\n",nzz,nxx,2); return -1;}

	p_2 = alloc2float(nzz,nxx);    // Campos de pressão: p[851][2401][2]
	if(p_2==NULL) { printf("Allocation of p_2[%d,%d] failed!!!\n",nzz,nxx); return -1;}

	p_1 = alloc2float(nzz,nxx);    // Campos de pressão: p[851][2401][2]
	if(p_1==NULL) { printf("Allocation of p_1[%d,%d] failed!!!\n",nzz,nxx); return -1;}

	vel = alloc2float(nzz,nxx);      // Matriz de velocidade: vel[851,2401]
	if(vel==NULL) { printf("Allocation of vel[%d,%d] failed!!!\n",nzz,nxx); return -1;}
    	
    	gama_z = alloc1float(nzz);          
	if(gama_z==NULL) { printf("Allocation of gama_z[%d] failed!!!\n",nzz); return -1;}

    	gama_x = alloc1float(nxx);
	if(gama_x==NULL) { printf("Allocation of gama_x[%d] failed!!!\n",nxx); return -1;}

	window_x = alloc1float(nx);
	if(window_x==NULL) { printf("Allocation of window_x[%d] failed!!!\n",nx); return -1;}

	imag = alloc2float(nz,nx);
	if(imag==NULL) { printf("Allocation of imag[%d,%d] failed!!!\n",nz,nx); return -1;}

	imag_filter = alloc2float(nz,nx);
	if(imag_filter==NULL) { printf("Allocation of imag_filter[%d,%d] failed!!!\n",nz,nx); return -1;}

	//////////////////////////////////////////////////////////////////////////////////    
	//				RTM INIT					//
	//////////////////////////////////////////////////////////////////////////////////    
	char *vel_name;

	if(modelo=='p')
	{
		//	Leitura e criação de bordas no modelo de velocidades - PEQUENO		//
		vel_name = "marm_vel_smooth_05_05.bin";// vel[243][767]
	}
	else if(modelo=='g')
	{
		//	Leitura e criação de bordas no modelo de velocidades - GRANDE		//
		//char *vel_name = "velocity.h@";// vel[751][2301]
		vel_name = "velocity_smooth_05_05.bin";
	}

	FILE *vel_file;
	
	if( (vel_file = fopen(vel_name, "rb")) == NULL ) { printf("Error opening vel. model file\n"); return -1; }	

//	gettimeofday(&start,NULL);

	//	Leitura do modelo de velocidades 					//
	get_vel_model(vel, vel_file, nz, nx, nborda);	

	fclose(vel_file);

	get_min_max(vel, &vmin, &vmax, ixb, ixe, izb, ize);

	// 	Preenchimento das bordas laterais e superiores				//
	float *flt_vet_left;
	flt_vet_left = (float*) malloc(sizeof(float)*nz);
	memcpy(flt_vet_left, vel + (nborda*nzz + nborda), nz*sizeof(float));

	float *flt_vet_right;
	flt_vet_right = (float*) malloc(sizeof(float)*nz);
	memcpy(flt_vet_right, vel + ((ixe-1)*nzz + nborda), nz*sizeof(float));

	//omp_set_num_threads(n_threads);
	//#pragma omp parallel for default(none) private(ix) shared(vel,flt_vet_left,flt_vet_right,nz,nzz,ixe)
 	for(ix=0;ix<nborda;++ix) // collum order ... super otimizado
	{
		memcpy(vel + (ix*nzz)     + nborda, flt_vet_left, nz*sizeof(float));
		memcpy(vel + (ixe+ix)*nzz + nborda, flt_vet_right, nz*sizeof(float));
	}

	free(flt_vet_left);
	free(flt_vet_right);

	//omp_set_num_threads(n_threads);
	//#pragma omp parallel for default(none) private(ix,iz) shared(vel,nzz,nxx,ixe,izb,ize)
	for(ix=0;ix<nxx;++ix) // line order ... nao tem como super otimizar
    	{
        	for(iz=0;iz<(izb-1);++iz)
		{
			vel[(iz) + (ix*nzz)] =  vel[(izb-1) + (ix*nzz)];
			vel[(iz + ize) + (ix*nzz)] =  vel[(ize-1) + (ix*nzz)];
		}
   	}

//	gettimeofday(&finish,NULL);
//	printf("Tempo de criacao do modelo de velocidades: %d ms\n",get_exec_time(start,finish));

	printf("\tModelo de velocidades\n");
	printf("\tnz: %d   ",nz);
	printf("nx: %d   ",nx);
	printf("z0: %d   ",z0);
	printf("x0: %d   ",x0);
	printf("dz: %.2f   ",dz);
	printf("dx: %.2f\n",dx);
	printf("\tvmin: [%.2f]   vmax: [%.2f]\n",vmin, vmax);

	printf("|---------------------------------------------------------------------------|\n");

	//////////////////////////////////////////////////////////////////////////////////
        //				GETTING SEISMIC DATA	 			//
	//////////////////////////////////////////////////////////////////////////////////

	//	Dados dos dados sismicos						//TODO: argc e argv
    	nshots 		= 1;//240;
    	ngeophones 	= 96;
	ns		= 725;
	ntrc   		= nshots*ngeophones;
	char *su_name = "MARMOUSI_xtLINUX_CORRECT2.su";

	FILE *su_file;
	if( (su_file = fopen(su_name, "rb")) == NULL ) { printf("Error opening seismic data file\n"); return -1; }

	su_trace trace;
	trace.tr_header = (su_header*) malloc( sizeof(su_header) ); 
	memset(trace.tr_header, 0, sizeof(su_header));
	trace.tr_data = alloc1float(ns);
	if(trace.tr_data==NULL) { printf("Allocation of trace.tr_data[%d] failed!!!\n",ns); return -1; }
	
	//	Dados do primeiro traco							//
	get_tr(0, &trace, ns, su_file);

	isx = trace.tr_header->sx;

	ns = trace.tr_header->ns;
	dtrec  = (trace.tr_header->dt)/sec2mcs; // Intervalo de amostragem dos dados //
	ttotal = (ns-1)*dtrec;
	
	printf("\tAtributos do modelo de dados\n");
	printf("\tns: %d   ", ns);
	printf("dt: %.3lf   ", (float)(trace.tr_header->dt));
	printf("sec2mcs: %.3lf   ", sec2mcs);
	printf("dtrec: %lf\n", dtrec);
	printf("\tttotal: %lf   ", ttotal);

	// 	fator de escala para coordenadas scalco 				//
	printf("scalco: %d   ", (trace.tr_header->scalco));

	factor = 1.0;
	if ( (trace.tr_header->scalco) < 0 )
		factor = 1.0 / fabs((float)(trace.tr_header->scalco)); //!SEG-Y standard
	else if( (trace.tr_header->scalco) > 0 )
		factor = ((float)(trace.tr_header->scalco));// ! SEG-Y standard
	
	printf("factor: %.3f\n",factor);

	//////////////////////////////////////////////////////////////////////////////////
        //			ALLOCATING DYNAMIC DATA 				//
	//////////////////////////////////////////////////////////////////////////////////

	printf("\tNumber of traces: %d\n",ntrc);

	shotidx = alloc1int(nshots+1);
    	if(shotidx==NULL) { printf("Allocation of shotidx[%d] failed!!!\n",nshots); return -1;}

	//	First trace from each shot position					//
	for(ix=0;ix<(nshots+1);++ix)
		shotidx[ix]=ngeophones*ix;

	printf("\tNumber of shots: %d\n",nshots);
	printf("\tNumber of geophones: %d\n",ngeophones);

	printf("|---------------------------------------------------------------------------|\n");

	//TODO: transformar shotgather em uma matriz e criar um vetor shotgather_header

    	shotgather = alloc1su_trace(ntrc);
	if(shotgather==NULL) { printf("Allocation of shotgather[%d] failed!!!\n",ntrc); return -1;}

	//TODO: otimizar isto
	for(itrc=0;itrc<ntrc;++itrc)
	{
		shotgather[itrc].tr_header = (su_header*) malloc( sizeof(su_header) );
		shotgather[itrc].tr_data = alloc1float(ns);

		if(shotgather[itrc].tr_data==NULL) 
		{ printf("Allocation of shotgather[%d].tr_data[%d] failed!!!\n",itrc,ns); return -1; }
	}
	

	igz = alloc1int(ngeophones); // correcao devida a profundidade da fonde e do receptor numa levantamento marinho //
	if(igz==NULL) { printf("Allocation of igz[%d] failed!!!\n",ngeophones); return -1;}

	igx = alloc1int(ngeophones); // correcao devida a profundidade da fonde e do receptor numa levantamento marinho //
	if(igx==NULL) { printf("Allocation of igx[%d] failed!!!\n",ngeophones); return -1;}

	//////////////////////////////////////////////////////////////////////////////////
        //				SETTING PARAMETERS 				//
	//////////////////////////////////////////////////////////////////////////////////

	// check grid interval com o objetivo de evitar dispersao numerica, e garantir acuidade!! //
	// limitar valores referentes a amostragem do campo em relacao ao dx e dt max //
	dxmax    = vmin / ( ( 6.0 * freq ) );
	dt       = dx   / ( ( 4.0 * vmax ) );
	ndtrec   = (int) ceil( dtrec/dt );
	ndtrtm   = (int) ceil( dtrtm/dt );

	dt = (ndtrec > 1) ? dtrec/ndtrec : dtrec; 

	nt = ceil( ((double)ttotal)/((double)dt) );// - 1 ;// POG

	if(dx > dxmax) { printf("Freq muito alta!!!"); return -1; }

	printf("\tAtributos do RTM\n");
	printf("\tdxmax: %lf   ",dxmax);
	printf("dt: %lf   ",dt);
	printf("dtrec: %lf   ",dtrec);
	printf("dtrtm: %lf\n",dtrtm);
	printf("\tndtrec: %d   ",ndtrec);
	printf("ndtrtm: %d   ",ndtrtm);
	printf("nt: %d   ",nt);

	//	pulso fonte - ricker 						//
	beta = D_PI*freq*dt;
    	it0 = (int) ceil( 1.0/(freq*dt) );

	printf("beta: %lf\t",beta);
	printf("it0: %d\n",it0);

	#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) private(ix,iz) shared(vel,nxx,nzz,dt,dx)
    	for(ix=0;ix<nxx;++ix) // (v[z][x])^2
        	for(iz=0;iz<nzz;++iz)
            		vel[(iz) + (ix*nzz)] = pow( (vel[(iz) + (ix*nzz)]*dt/dx), 2.0);

	memset(gama_z, 0, nzz * sizeof(float) );

	for(ix=1;ix<(nborda+1);++ix)
    	{
        	gama = beta * pow( ( ((float)ix)/((float)nborda) ), 2);
        	gama_z[izb-ix-1] = gama;
        	gama_z[ize+ix-1] = gama;
    	}	

	su_trace imagtrace;
	imagtrace.tr_header = (su_header*) malloc( sizeof(su_header) );
	memset(imagtrace.tr_header, 0, sizeof(su_header));
	imagtrace.tr_data = alloc1float(nz);
	if(imagtrace.tr_data==NULL) { printf("Allocation of imagtrace.tr_data[%d] failed!!!\n",nz); return -1;}

    	// 	header do traco imagem 						//
    	imagtrace.tr_header->ns = (unsigned short) nz;
    	imagtrace.tr_header->dt = (unsigned short) (sec2mcs*dtrtm);
    	imagtrace.tr_header->f1 = (float) z0;
    	imagtrace.tr_header->f2 = (float) x0;
    	imagtrace.tr_header->d1 = (float) dz;
    	imagtrace.tr_header->d2 = (float) dx;

	//////////////////////////////////////////////////////////////////////////////////
        //				RTM MIGRATION	                     	 	//
	//////////////////////////////////////////////////////////////////////////////////

	memset(imag, 0, sizeof(float) * nz * nx);
	memset(imag_filter, 0, sizeof(float) * nz * nx);

	//-FILE *fdbackward_file;	// campo retropropagado dos dados //
	//-char *fdbackward_name = "fdbackward.bin";

	//-FILE *fdforward_file;	// campo da fonte propagado //
	//-char *fdforward_name;

	//-if(write_forward_file)
	//-	fdforward_name = "fdforward.bin";

	FILE *out_file;		// saida - dados migrados
	char *out_name = "rtm_migrated_otimizado.su";

    	ilanco=0;

	printf("|---------------------------------------------------------------------------|\n");

	nrec = ntrc/nshots; // 96

	//TODO: testar depois sem fseek

	// lendo os dados do marmousi //	
	for(itrc=0;itrc<ntrc;++itrc)
     	{
		get_tr(itrc, &shotgather[itrc], ns, su_file);
	}

	gettimeofday(&start_princ,NULL);

    	for(ishot=0;ishot<nshots;++ishot)
 	{
		printf("\n\n");
		printf("|###########################################################################|\n");
		printf("\tMIGRATING shot %d out of %d...\n",ishot+1,nshots);
		printf("|###########################################################################|\n");

		ircv = 0; // numero de tracos para common-shot //
	
		//#pragma omp for
		for(itrc=shotidx[ishot];itrc<shotidx[ishot+1];++itrc)
		{
			//	armazena posicao do traco ma malha do modelo 		//
			//	converte unidades usando a palavra scalco do header 	//
			// 	=== ATENCAO === 					//
			// 	as palavras gelev e sdepth devem ser setadas no header 	//

			igx[ircv] = (int) floor( ( ((float)(shotgather[itrc].tr_header->gx))-((float)x0) )*\
												(factor/dx) )+ixb;

			igz[ircv] = (int) floor( (-1.0)*((float)(shotgather[itrc].tr_header->gelev)) *\
												(factor/dz) )+izb;

			isx 	  = (int) floor( ( ((float)(shotgather[itrc].tr_header->sx))-((float)x0) )*\
												(factor/dx) )+ixb;

			isz       = (int) floor( ((float)(shotgather[itrc].tr_header->sdepth)) * (factor/dz) ) + izb;

			//#pragma omp critical
			//{
				ilanco = max(ilanco, fabs(isx-igx[ircv]));
 				++ircv;
			//}
  		}


		// 	DEFINIR JANELA DE MIGRACAO 					//
		if( LFRAC == 0 )
		{
			ixmig0 = 1;
			ixmig1 = nx;
		}
		else
		{
			ixmig0 = max( min(minval(igx,nrec),isx)-ilanco/LFRAC, ixb);
			ixmig1 = min( max(maxval(igx,nrec),isx)+ilanco/LFRAC, ixe);
		}

		printf("\tJANELA DE MIGRACAO do tiro\n");
		printf("\tJanela s/ borda - ixmig0: %d   ixmig1: %d\n",ixmig0,ixmig1);

		nmigtrc = ixmig1 - ixmig0 + 1;

		// 	DEFINIR JANELA DE PROPAGACAO 					//
		ixx0 = ixmig0 - nborda - 1;
		ixx1 = ixmig1 + nborda;
	
		printf("\tJanela c/ borda - ixx0: %d   ixx1: %d\n",ixx0, ixx1);
	

		//memset(gama_x,0,sizeof(float)*nxx);
		//#pragma omp for
		for(ix=0;ix<nxx;++ix) gama_x[ix]=0;

		//#pragma omp for
		for(ix=1;ix<nborda+1;++ix)
		{
			gama = beta * pow( ( ((float)ix)/((float)nborda) ) , 2 );
			gama_x[ixmig0-ix-1] = gama;
			gama_x[ixmig1+ix-1] = gama;
		}

		// JANELA BLACKMANN para suavizar borda da imagem //
		//#pragma omp for nowait
		for(ix=0;ix<nx;++ix) window_x[ix] = 1.0; // window_x(:) = 1.0

		//#pragma omp for
		for(ix=0;ix<WDWLEN;++ix)
		{
			window_x[ixmig0-1-nborda+WDWLEN-1-ix] = blackmann[ix];
			window_x[ixmig1-nborda-WDWLEN+ix] = blackmann[ix];
		}

		//////////////////////////////////////////////////////////////////////////////////
		//	       	RETROPROPAGACAO DO CAMPO DOS RECEPTORES				//
		//////////////////////////////////////////////////////////////////////////////////

		// 	condiçao inicial: repouso 		//
		// 	p(:,:,1) pressure field at time t-dt 	//
		// 	p(:,:,2) pressure field at time t 	//

		//memset(p_2, 0, sizeof(float)*nzz*nxx);
		//memset(p_1, 0, sizeof(float)*nzz*nxx);
		#pragma omp parallel for if(is_parallel) num_threads(n_threads) default(none) \
			shared(nzz,nxx,p_2,p_1) private(iz,ix)
		for(ix=0;ix<nxx;++ix)
		{
			for(iz=0;iz<nzz;++iz)
			{
				p_2[ (iz) + (ix*nzz)]  = 0.0;
				p_1[ (iz) + (ix*nzz)]  = 0.0;
			}
		}

		reclen = nz*nmigtrc*sizeof(float);

		printf("\tnmigtrc: %d   ",nmigtrc);
		printf("reclen: %d\n",reclen);

		printf("|###########################################################################|\n");

		printf("\tRETROPROPAGATION START\n");

		//-if( (fdbackward_file = fopen(fdbackward_name, "wb")) == NULL ) 
		//-{ printf("Error fdbackward file\n"); }//return -1; }

		nframes = (nt%ndtrtm==0) ? nt/ndtrtm : (int) ceil( ((float)nt)/((float)ndtrtm) );

		prcv = alloc3float(nframes,nz,nmigtrc);
		if(prcv==NULL) { printf("Allocation of prcv[%d,%d,%d] failed!!!\n",nframes,nz,nmigtrc); }//return -1; }
		//-prcv = alloc2float(nz,nmigtrc);
		//-if(prcv==NULL) { printf("Allocation of prcv[%d,%d] failed!!!\n",nz,nmigtrc); }//return -1; }

		irec=0;
		iframe=0;

		// medindo tempo de retropropagacao
		gettimeofday(&start,NULL);
	

		for(it=0;it<nt;++it) // modelando a evolucao do campo de pressao //
		{
			// FD scheme: forward time //
			t = ((float)it) * dt;//((float)(it-1)) * dt;

  			#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) \
			shared(ixx0,ixx1,nzz,nxx,nrec,igx,igz,dx,ttotal,t,ns,dtrec,shotgather,gama_x,  \
				gama_z,deriv2,p_1,p_2,vel,ishot)				       \
			private(source,laplacian,ix,iz,itrc,gama,invpgama,mgama,iconv)
			for(ix=(ixx0+DRVLEN);ix<(ixx1-DRVLEN);++ix)
			{
				for(iz=(DRVLEN-1);iz<(nzz-DRVLEN);++iz)
				{
					source = 0.0;

					for(itrc=0;itrc<nrec;++itrc)
					{
						if( (ix == (igx[itrc]-1)) && (iz == (igz[itrc]-1)) )
						{
							source = (dx*dx) * \
							interp_trace(ttotal-t,ns,0.0,dtrec,		\
									shotgather[ itrc + (ishot*nrec) ].tr_data);
									//shotgather[itrc].tr_data);
						}
					}

					gama = gama_x[ix] + gama_z[iz]; 
					invpgama  = ( 1.0 / (1.0 + gama) );
					mgama     = 1.0 - gama;

					//	aumentar a ordem do operador de diferencas finitas: 		//
					laplacian = 2.0 * deriv2[0] * p_2[ (iz) + (ix*nzz)];

					for(iconv=1;iconv<DRVLEN;++iconv)
					{
						laplacian = laplacian + deriv2[iconv] 				* \
			      				(							  \
								p_2[ (iz-iconv) + (ix*nzz) ] 	+ \
								p_2[ (iz+iconv) + (ix*nzz) ] 	+ \
								p_2[ (iz) + ((ix-iconv)*nzz) ]	+ \
								p_2[ (iz) + ((ix+iconv)*nzz) ] 	  \
							);
					}

					p_1[ (iz) + (ix*nzz) ] = invpgama * ( 	 2.0 * \
							    p_2[ (iz) + (ix*nzz) ] - mgama * \
							    p_1[ (iz) + (ix*nzz) ] + 	       \
						            vel[(iz) + (ix*nzz)] * ( laplacian - source ) );
				}
			}

			// mega - new swap fields
			p_aux  = p_2;
			p_2    = p_1;
			p_1    = p_aux;

			if( ( it % ndtrtm) == 0)
			{
				if( (iframe%100==0) )
					printf("\tit: %d - backward frames %d completed\n", it, iframe);

				//-#pragma omp for private(ix)
				for(ix=0;ix<nmigtrc;++ix)
					memcpy(prcv + (iframe*nz*nmigtrc) + (ix*nz),p_2 + ((ixmig0-1+ix)*nzz) + (nborda), \
													nz * sizeof(float) );

				//-for(ix=0;ix<nmigtrc;++ix)
				//-	memcpy(prcv+(ix*nz), p_2 + ((ixmig0-1+ix)*nzz)+(nborda), nz*sizeof(float));

				//-fwrite(prcv, sizeof(float), nz * nmigtrc, fdbackward_file);

				++iframe;
			}

		} // fim de modelagem

		gettimeofday(&finish,NULL);

		timeval_subtract(&diff, &finish, &start);
		printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);

		//-free2float(prcv);	//-
		//-fclose(fdbackward_file);//-

		nframes = iframe;

		printf("\tRETROPROPAGATION OF %d FRAMES END\n", nframes);

		printf("|###########################################################################|\n");
	
		//////////////////////////////////////////////////////////////////////////////////
		//		       	PROPAGACAO DO CAMPO DA FONTE				//
		//////////////////////////////////////////////////////////////////////////////////

		printf("\tSOURCE PROPAGATION START\n");

		//-if( (fdbackward_file = fopen(fdbackward_name, "rb")) == NULL )
		//-{printf("Error opening model file\n"); }//return -1;}

		//-if(write_forward_file)
		//-	if((fdforward_file=fopen(fdforward_name,"wb"))==NULL)
		//-		{printf("Error opening forward source file\n");}//return -1;}
		
		// 	condiçao inicial: repouso 		//
		// 	p(:,:,1) pressure field at time t-dt 	//
		// 	p(:,:,2) pressure field at time t 	//

		//memset(p_2, 0, sizeof(float)*nzz*nxx);
		//memset(p_1, 0, sizeof(float)*nzz*nxx);
		#pragma omp parallel for if(is_parallel) num_threads(n_threads) default(none) \
			shared(nzz,nxx,p_2,p_1) private(iz,ix)
		for(ix=0;ix<nxx;++ix)
		{
			for(iz=0;iz<nzz;++iz)
			{
				p_2[ (iz) + (ix*nzz)] = 0.0;
				p_1[ (iz) + (ix*nzz)] = 0.0;
			}
		}	

		//-prcv = alloc2float(nz,nmigtrc);
		//-if(prcv==NULL) { printf("Allocation of prcv[%d,%d] failed!!!\n",nz,nmigtrc); }//return -1;}

		psrc = alloc3float(nframes,nz,nmigtrc);
		if(psrc==NULL) { printf("Allocation of psrc[%d,%d,%d] failed!!!\n",nframes,nz,nmigtrc); }//return -1;}
		//-psrc = alloc2float(nz,nmigtrc);
		//-if(psrc==NULL) { printf("Allocation of psrc[%d,%d] failed!!!\n",nz,nmigtrc); }//return -1;}		

		irec   = 0;
		iframe = 0;

		gettimeofday(&start,NULL);

		// modelando a evolucao do campo de pressao //
		for(it=0;it<nt;it++)
		{
			delt  = pow( beta * ((float)(1+it-it0)), 2.0);

			// Pulso fonte Ricker frequancia pico freq //
		        fonte = (dx*dx) * exp( -delt ) * ( 1.0 - (2.0 * delt) );   

  			#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) \
			shared(ixx0,ixx1,nzz,nxx,nrec,isx,isz,dx,ttotal,t,ns,dtrec,shotgather,gama_x,  \
				gama_z,deriv2,p_1,p_2,vel,fonte)				       \
			private(source,laplacian,ix,iz,itrc,gama,invpgama,mgama,iconv)
			for(ix=(ixx0+DRVLEN);ix<(ixx1-DRVLEN);ix++)
			{
				for(iz=(DRVLEN-1);iz<(nzz-DRVLEN);iz++)
				{
					source = 0.0;

					if ( (ix == (isx-1)) && (iz == (isz-1)) )
	   				{
						source = fonte;
					}

					gama = gama_x[ix] + gama_z[iz]; 
					invpgama  = 1.0 / ( 1.0 + gama);
					mgama     = 1.0 - gama;

					//	aumentar a ordem do operador de diferencas finitas: 		//
					laplacian = 2.0 * deriv2[0] * p_2[ (iz) + (ix*nzz)];

					for(iconv=1;iconv<DRVLEN;++iconv)
					{
						laplacian = laplacian + deriv2[iconv] 		* \
			      				(					  \
								p_2[ (iz-iconv) + (ix*nzz) ] 	+ \
								p_2[ (iz+iconv) + (ix*nzz) ] 	+ \
								p_2[ (iz) + ((ix-iconv)*nzz) ]	+ \
								p_2[ (iz) + ((ix+iconv)*nzz) ] 	  \
							);
					}

					p_1[ (iz) + (ix*nzz) ] = invpgama * ( 	 2.0 * \
							    p_2[ (iz) + (ix*nzz) ] - mgama * \
							    p_1[ (iz) + (ix*nzz) ] + 	       \
						            vel[(iz) + (ix*nzz)] * ( laplacian - source ) );
				}
			}
		
			// mega - new swap fields
			p_aux      = p_2;
			p_2    = p_1;
			p_1 = p_aux;

			if( ( it % ndtrtm) == 0)
			{
				if( (iframe%100==0) ) 
					printf("\tit: %d - forward frames %d completed\n", it, iframe);

				//%irec = nframes-iframe-1;

				//%fseek(fdbackward_file, irec * reclen , SEEK_SET);
				//%fread(prcv, sizeof(float), nz * nmigtrc, fdbackward_file);

				//#pragma omp for private(ix)
				for(ix=0;ix<nmigtrc;++ix)
					memcpy(psrc + (iframe*nz*nmigtrc) + (ix*nz), p_2 + ((ixmig0-1+ix)*nzz) + (nborda), \
											nz * sizeof(float) );

				//-if(write_forward_file)
				//-	fwrite(psrc, sizeof(float), nz * nmigtrc, fdforward_file);
			
			/*%
			%	//////////////////////////////////////////////////////////////////
			%	//	 	CONDICAO DE IMAGEM CROSS-CORRELACAO 		//
			%	//////////////////////////////////////////////////////////////////
			%
			%	itrc=0;

			%	//--#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) \
			%	//-		private(ix,iz,itrc) \
			%	//-		shared(window_x,psrc,prcv,ixmig0,ixmig1,nz,nmigtrc,imag)
			%	//-for(ix=(ixmig0-1-nborda);ix<(ixmig1-nborda);++ix)
			%	//-{
			%	//-//	itrc = ix % ( (ixmig1-nborda) - (ixmig0-1-nborda) ); //TODO: rever isto
			%	//-
			%	//-	for(iz=0;iz<nz;++iz)
	  		%	//-		imag[(iz) + (ix*nz)]  = imag[(iz) + (ix*nz)] + ( window_x[ix] * \
			%	//-				psrc[ (iframe*nz*nmigtrc) + (iz) + (itrc*nz)] * \
			%	//-				prcv[ (irec*nz*nmigtrc)   + (iz) + (itrc*nz)] );
			%	//-	++itrc;
			%	//-}
%
			%	//#pragma omp for
			%	for(ix=(ixmig0-1-nborda);ix<(ixmig1-nborda);++ix)
			%	{
			%		itrc = ix % ( (ixmig1-nborda) - (ixmig0-1-nborda) ); //TODO: rever isto
			%		//printf("tn:%d|itrc: %d\n",omp_get_thread_num(),itrc);
%
			%		for(iz=0;iz<nz;++iz)
			%			imag[(iz) + (ix*nz)]  = imag[(iz) + (ix*nz)] + ( window_x[ix] * \
			%					psrc[ (iz) + (itrc*nz)] * prcv[ (iz) + (itrc*nz)] );
			%		
			%		//++itrc;
			%	}
			%*/
				++iframe;
			}
		} // fim de modelagem //

		nframes = iframe;
		//-fclose(fdforward_file);

		printf("\tSOURCE PROPAGATION OF %d FRAMES END\n",nframes);

		gettimeofday(&finish,NULL);

		timeval_subtract(&diff, &finish, &start);
		printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);

		printf("|###########################################################################|\n");
		printf("\tCROSS CORRELATION START\n");

		//-if(write_forward_file)
		//-	if((fdforward_file=fopen(fdforward_name,"rb"))==NULL)
		//-		{printf("Error opening forward source file\n");}//return -1;}

		irec=0;
		itrc=0;

		gettimeofday(&start,NULL);

		for(iframe=0;iframe<nframes;++iframe)
		{
		
			//-fseek(fdforward_file, iframe * reclen , SEEK_SET);
			//-fread(psrc, sizeof(float), nz * nmigtrc, fdforward_file);

			irec = nframes-iframe-1;
			//-fseek(fdbackward_file, irec * reclen , SEEK_SET);
			//-fread(prcv, sizeof(float), nz * nmigtrc, fdbackward_file);

			//////////////////////////////////////////////////////////////////
			//			NOVA CROSS-CORRELACAO 			//
			//////////////////////////////////////////////////////////////////

			itrc = 0;

  			#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none)	\
				shared(ixx0,nz,ixmig0,ixmig1,iframe,imag,psrc,prcv,nmigtrc,irec,window_x)\
				private(ix,iz,itrc)
			for(ix=(ixmig0-1-nborda);ix<(ixmig1-nborda);++ix)
			{
				itrc = (ix-ixx0) % (nmigtrc);

				for(iz=0;iz<nz;++iz)
					imag[(iz) + (ix*nz)]  = imag[(iz) + (ix*nz)] + ( window_x[ix] * \
							psrc[ (iframe*nz*nmigtrc) + (iz) + (itrc*nz)] * \
							prcv[ (irec*nz*nmigtrc)   + (iz) + (itrc*nz)] );

				//++itrc;
			}

			//-#pragma omp for
			//-for(ix=(ixmig0-1-nborda);ix<(ixmig1-nborda);++ix)
			//-{
			//-	//itrc = ix % ( (ixmig1-nborda) - (ixmig0-1-nborda) ); //TODO: rever isto
			//-					
			//-	itrc = (ix-ixx0) % (nmigtrc);//ixx0=ixmig0-1-borda
			//-	
			//-	//printf("ix:%d|itrc: %d\n",ix,itrc);
			//-
			//-	for(iz=0;iz<nz;++iz)
			//-		imag[(iz) + (ix*nz)]  = imag[(iz) + (ix*nz)] + ( window_x[ix] * \
			//-				psrc[ (iz) + (itrc*nz)] * prcv[ (iz) + (itrc*nz)] );
			//-
			//-	//++itrc;
			//-}

		}

		printf("\tCROSS CORRELATION END\n");
		gettimeofday(&finish,NULL);

		timeval_subtract(&diff, &finish, &start);
		printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);
	
		//-free2float(prcv);
		//-free2float(psrc);

		free3float(prcv);
		free3float(psrc);

		//-fclose(fdbackward_file);

		//-if(write_forward_file)
		//-	fclose(fdforward_file);

		//////////////////////////////////////////////////////////////////////////////////
        	//				FILTRO LAPLACIANO				//
		//////////////////////////////////////////////////////////////////////////////////

		#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) \
		shared(nz,nx,deriv2,vel,imag,imag_filter)				       \
		private(laplacian,ix,iz,iconv)
		for(ix=(DRVLEN-1);ix<(nx-DRVLEN);ix++)
		{
			for(iz=(DRVLEN-1);iz<(nz-DRVLEN);iz++)
			{
				laplacian = 2.0 * deriv2[0] * imag[(iz) + (ix*nz)];

				for(iconv=1;iconv<DRVLEN;++iconv)
				{
					laplacian = laplacian + deriv2[iconv] 			* \
						(						  \
							imag[(iz-iconv) + (ix*nz)        ]	+ \
							imag[(iz+iconv) + (ix*nz)        ]	+ \
							imag[(iz)       + ((ix-iconv)*nz)]	+ \
							imag[(iz)       + ((ix+iconv)*nz)]  	  \
						);
				}

				imag_filter[(iz) + (ix*nz)] = laplacian;
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
        	//				ESCRITA FINAL DOS DADOS				//
		//////////////////////////////////////////////////////////////////////////////////

		if(ishot==nshots-1)
		{

			if( (out_file = fopen(out_name, "wb")) == NULL ) { printf("Error opening out file\n"); }//return -1; }

			//#pragma omp for 
			for(ix=0;ix<nx;ix++)
			{
				imagtrace.tr_header->tracl = ix + 1;

				//memcpy(imagtrace.tr_data, imag_filter + (ix*nz), nz * sizeof(float) );
				//put_tr(ix, &imagtrace, nz, out_file);

				fwrite(imagtrace.tr_header, HDRBYTES, 1, out_file);
				fwrite(imag_filter + (ix*nz), sizeof(float), nz, out_file);
			}

			fclose(out_file);
		}

		printf("|###########################################################################|\n");

	} // END LOOP OF SHOTS

	gettimeofday(&finish_princ,NULL);

	timeval_subtract(&diff_princ, &finish_princ, &start_princ);
	printf("%ld.%06ld\n", diff_princ.tv_sec, diff_princ.tv_usec);

	fprintf(stderr,"%ld.%06ld\n", diff_princ.tv_sec, diff_princ.tv_usec);

	printf("\n\n|---------------------------------------------------------------------------|\n");
	printf("\tend of migration\n");
	printf("|---------------------------------------------------------------------------|\n");

	//////////////////////////////////////////////////////////////////////////////////
        //		DEALOCATING VECTORS & CLOSING FILES				//
	//////////////////////////////////////////////////////////////////////////////////

	//TODO: checar se todo mundo foi desalocado
	
	fclose(su_file);

	free1int(shotidx);
	free1int(igz);
	free1int(igx);

	free1float(gama_x);
	free1float(gama_z);
	free1float(window_x);
	free1float(trace.tr_data);
	free1float(imagtrace.tr_data);

	free2float(imag);
	free2float(imag_filter);
	free2float(vel);

	free2float(p_1);
	free2float(p_2);

	//TODO: dealocate all shotgather
	free1su_trace(shotgather);

    	return 0;
}

