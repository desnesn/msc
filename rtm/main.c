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
	printf("|                     Reverse Time Migration - RTM 6.0                      |\n");
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

	//&float UNUSED *p_1;
	//&float UNUSED *p_2;
	float UNUSED *p_aux;
	float UNUSED p_swap;

	wave_fields *fields;

	float UNUSED dswap;

	//&float UNUSED *psrc;
	//&float UNUSED *prcv;
	p_frames *frames;

	//&float UNUSED *imag;
	//&float UNUSED *imag_filter;
	mig_section *sections;

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
	int UNUSED *shotidx;
	int UNUSED iconv, ixb, ixe, izb, ize, nxx, nzz;
	int UNUSED iframe, nframes, irec, ircv, itrc, ns;
	int UNUSED ntrc, ngeophones, ishot, nshots;
	float UNUSED factor;

	//&int UNUSED *igx , *igz;
	trace_pos *pos;

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

	fields = (wave_fields*) malloc(sizeof(wave_fields)*nzz*nxx);
	//&p_2 = alloc2float(nzz,nxx);    // Campos de pressão: p[851][2401][2]
	//&if(p_2==NULL) { printf("Allocation of p_2[%d,%d] failed!!!\n",nzz,nxx); return -1;}

	//&p_1 = alloc2float(nzz,nxx);    // Campos de pressão: p[851][2401][2]
	//&if(p_1==NULL) { printf("Allocation of p_1[%d,%d] failed!!!\n",nzz,nxx); return -1;}

	vel = alloc2float(nzz,nxx);      // Matriz de velocidade: vel[851,2401]
	if(vel==NULL) { printf("Allocation of vel[%d,%d] failed!!!\n",nzz,nxx); return -1;}
    	
    	gama_z = alloc1float(nzz);          
	if(gama_z==NULL) { printf("Allocation of gama_z[%d] failed!!!\n",nzz); return -1;}

    	gama_x = alloc1float(nxx);
	if(gama_x==NULL) { printf("Allocation of gama_x[%d] failed!!!\n",nxx); return -1;}

	window_x = alloc1float(nx);
	if(window_x==NULL) { printf("Allocation of window_x[%d] failed!!!\n",nx); return -1;}

	/*&
	imag = alloc2float(nz,nx);
	if(imag==NULL) { printf("Allocation of imag[%d,%d] failed!!!\n",nz,nx); return -1;}

	imag_filter = alloc2float(nz,nx);
	if(imag_filter==NULL) { printf("Allocation of imag_filter[%d,%d] failed!!!\n",nz,nx); return -1;}
	&*/

	sections = (mig_section*) malloc(sizeof(mig_section)*nz*nx);	

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

 	for(ix=0;ix<nborda;++ix) // collum order ... super otimizado
	{
		memcpy(vel + (ix*nzz)     + nborda, flt_vet_left, nz*sizeof(float));
		memcpy(vel + (ixe+ix)*nzz + nborda, flt_vet_right, nz*sizeof(float));
	}

	free(flt_vet_left);
	free(flt_vet_right);

	for(ix=0;ix<nxx;++ix) // line order ... nao tem como super otimizar
    	{
        	for(iz=0;iz<(izb-1);++iz)
		{
			vel[(iz) + (ix*nzz)] =  vel[(izb-1) + (ix*nzz)];
			vel[(iz + ize) + (ix*nzz)] =  vel[(ize-1) + (ix*nzz)];
		}
   	}

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
	
	/*&
	igz = alloc1int(ngeophones); // correcao devida a profundidade da fonde e do receptor numa levantamento marinho //
	if(igz==NULL) { printf("Allocation of igz[%d] failed!!!\n",ngeophones); return -1;}

	igx = alloc1int(ngeophones); // correcao devida a profundidade da fonde e do receptor numa levantamento marinho //
	if(igx==NULL) { printf("Allocation of igx[%d] failed!!!\n",ngeophones); return -1;}
	&*/

	pos = (trace_pos*) malloc(sizeof(trace_pos)*ngeophones);

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

	//&memset(imag, 0, sizeof(float) * nz * nx);
	//&memset(imag_filter, 0, sizeof(float) * nz * nx);

	memset(sections, 0, sizeof(float) * nz * nx * 2);

	//-FILE *fdbackward_file;	// campo retropropagado dos dados //
	//-char *fdbackward_name = "fdbackward.bin";

	//-FILE *fdforward_file;	// campo da fonte propagado //
	//-char *fdforward_name;

	//-if(write_forward_file)
	//-	fdforward_name = "fdforward.bin";

	FILE *out_file;		// saida - dados migrados
	char *out_name = "rtm_migrated_6.0.su";

	int iframe_cross;
	int chunk=0;

    	ilanco=0;

	printf("|---------------------------------------------------------------------------|\n");

	nrec = ntrc/nshots; // 96

	//TODO: testar depois sem fseek
	// lendo os dados do marmousi //	
	for(itrc=0;itrc<ntrc;++itrc)
     	{
		get_tr(itrc, &shotgather[itrc], ns, su_file);
	}

	//TEMPORARIO
	for(ix=0;ix<nxx;++ix)
	{
		for(iz=0;iz<nzz;++iz)
		{
			fields[(iz)+(nzz*ix)].vel = vel[(iz)+(nzz*ix)];
		}
	}

	gettimeofday(&start_princ,NULL);

    	for(ishot=0;ishot<nshots;++ishot)
 	{
	//&prcv,psrc,igx,igz,imag,imag_filter,p_1,p_2,p_aux(shared),p_swap
	#pragma omp parallel if (is_parallel) num_threads(n_threads) default(none) \
		shared(beta,blackmann,diff,dx,dtrec,deriv2,dt,dz,factor,gama_x,gama_z,ixx0,ixx1,pos,		\
			ixmig1,ishot,ixmig0,iframe,ixb,ixe,isx,ilanco,izb,isz,finish,it0,sections,fields,	\
			imagtrace,ndtrtm,nshots,nzz,nxx,nrec,ns,nx,nz,nmigtrc,nt,nframes,shotidx,out_file,	\
			out_name,reclen,start,ttotal,shotgather,vel,window_x,x0,chunk,n_threads,frames)\
		private(delt,fonte,gama,ix,iz,itrc,invpgama,iconv,it,ircv,irec,iframe_cross,laplacian,mgama,source,t,p_aux,p_swap)
	{
		
		#pragma omp single nowait
		{

			printf("\n\n");
			printf("|###########################################################################|\n");
			printf("\tMIGRATING shot %d out of %d...\n",ishot+1,nshots);
			printf("|###########################################################################|\n");
		
			//ircv = 0; // numero de tracos para common-shot //

			printf("\tINIT START\n");

			// medindo tempo de retropropagacao
			gettimeofday(&start,NULL);

			//chunk = nrec/n_threads;
		}

		#pragma omp for schedule(static,10)//schedule(static,100)
		for(itrc=shotidx[ishot];itrc<shotidx[ishot+1];++itrc)
		{
			ircv = (itrc-shotidx[ishot]) % nrec;
			
			//	armazena posicao do traco ma malha do modelo 		//
			//	converte unidades usando a palavra scalco do header 	//
			// 	=== ATENCAO === 					//
			// 	as palavras gelev e sdepth devem ser setadas no header 	//
			
			/*&
			igx[ircv] = (int) floor( ( ((float)(shotgather[itrc].tr_header->gx))-((float)x0) )*\
												(factor/dx) )+ixb;

			igz[ircv] = (int) floor( (-1.0)*((float)(shotgather[itrc].tr_header->gelev)) *\
												(factor/dz) )+izb;
			&*/

			pos[ircv].igx = (int) floor( ( ((float)(shotgather[itrc].tr_header->gx))-((float)x0) )*\
												(factor/dx) )+ixb;

			pos[ircv].igz = (int) floor( (-1.0)*((float)(shotgather[itrc].tr_header->gelev)) *\
												(factor/dz) )+izb;

			#pragma omp critical
			{
				isx 	  = (int) floor( ( ((float)(shotgather[itrc].tr_header->sx))-((float)x0) )*\
												(factor/dx) )+ixb;

				isz       = (int) floor( ((float)(shotgather[itrc].tr_header->sdepth)) * (factor/dz) ) + izb;

				ilanco = max(ilanco, fabs(isx-pos[ircv].igx));
 				//++ircv;
			}
  		}

		#pragma omp single nowait
		{
			gettimeofday(&finish,NULL);

			timeval_subtract(&diff, &finish, &start);

			//-free2float(prcv);
			//-fclose(fdbackward_file);

			printf("\tINIT END\n");
			printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);

			printf("|###########################################################################|\n");

			// 	DEFINIR JANELA DE MIGRACAO 					//
			if( LFRAC == 0 )
			{
				ixmig0 = 1;
				ixmig1 = nx;
			}
			else
			{
				ixmig0 = max( min(minval_igx(pos,nrec),isx)-ilanco/LFRAC, ixb);
				ixmig1 = min( max(maxval_igx(pos,nrec),isx)+ilanco/LFRAC, ixe);
			}

			printf("\tJANELA DE MIGRACAO do tiro\n");
			printf("\tJanela s/ borda - ixmig0: %d   ixmig1: %d\n",ixmig0,ixmig1);

			nmigtrc = ixmig1 - ixmig0 + 1;

			// 	DEFINIR JANELA DE PROPAGACAO 					//
			ixx0 = ixmig0 - nborda - 1;
			ixx1 = ixmig1 + nborda;
	
			printf("\tJanela c/ borda - ixx0: %d   ixx1: %d\n",ixx0, ixx1);
		}

		//#pragma omp single nowait
		//{
		//	printf("|###########################################################################|\n");
		//	printf("\tTESTE START\n");
		//	chunk = nxx/n_threads;
		//	gettimeofday(&start,NULL);
		//}

		#pragma omp for //schedule(static,10)
		for(ix=0;ix<nxx;++ix) gama_x[ix]=0;

		//#pragma omp single nowait
		//{
		//	printf("\tTESTE END\n");
		//	gettimeofday(&finish,NULL);
		//	timeval_subtract(&diff, &finish, &start);
		
		//	printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);
		//	printf("|###########################################################################|\n");
		//}

		//#pragma omp single nowait
		//{
		//	printf("|###########################################################################|\n");
		//	printf("\tTESTE START\n");
		//	chunk = (nborda+1)/n_threads;
		//	gettimeofday(&start,NULL);
		//}

		#pragma omp for	//schedule(static,1)
		for(ix=1;ix<nborda+1;++ix)
		{
			gama = beta * pow( ( ((float)ix)/((float)nborda) ) , 2 );
			gama_x[ixmig0-ix-1] = gama;
			gama_x[ixmig1+ix-1] = gama;
		}

		//#pragma omp single nowait
		//{
		//	printf("\tTESTE END\n");
		//	gettimeofday(&finish,NULL);
		//	timeval_subtract(&diff, &finish, &start);
		//
		//	printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);
		//	printf("|###########################################################################|\n");
		//}

		//#pragma omp single nowait
		//{
		//	printf("|###########################################################################|\n");
		//	printf("\tTESTE START\n");
		//	chunk = nx/n_threads;
		//	gettimeofday(&start,NULL);
		//}

		// JANELA BLACKMANN para suavizar borda da imagem //
		#pragma omp for schedule(guided)
		for(ix=0;ix<nx;++ix) window_x[ix] = 1.0; // window_x(:) = 1.0

		//#pragma omp single nowait
		//{
		//	printf("\tTESTE END\n");
		//	gettimeofday(&finish,NULL);
		//	timeval_subtract(&diff, &finish, &start);
		//
		//	printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);
		//	printf("|###########################################################################|\n");
		//}

		#pragma omp single nowait
		{
			for(ix=0;ix<WDWLEN;++ix)
			{
				window_x[ixmig0-1-nborda+WDWLEN-1-ix] = blackmann[ix];
				window_x[ixmig1-nborda-WDWLEN+ix] = blackmann[ix];
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
		//	       	RETROPROPAGACAO DO CAMPO DOS RECEPTORES				//
		//////////////////////////////////////////////////////////////////////////////////

		// 	condiçao inicial: repouso 		//
		// 	p(:,:,1) pressure field at time t-dt 	//
		// 	p(:,:,2) pressure field at time t 	//

		//#pragma omp single nowait
		//{
		//	printf("|###########################################################################|\n");
		//	printf("\tTESTE START\n");
		//	chunk = nxx/n_threads;
		//	gettimeofday(&start,NULL);
		//}

		//memset(p_2, 0, sizeof(float)*nzz*nxx);
		//memset(p_1, 0, sizeof(float)*nzz*nxx);
		#pragma omp for schedule(static,100)
		for(ix=0;ix<nxx;++ix)
		{
			for(iz=0;iz<nzz;++iz)
			{
				fields[(iz) + (ix*nzz)].p_1  = 0.0;
				fields[(iz) + (ix*nzz)].p_2  = 0.0;
			}
		}

		//#pragma omp single nowait
		//{
		//	printf("\tTESTE END\n");
		//	gettimeofday(&finish,NULL);
		//	timeval_subtract(&diff, &finish, &start);
		//
		//	printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);
		//	printf("|###########################################################################|\n");
		//}

		#pragma omp single nowait
		{
			reclen = nz*nmigtrc*sizeof(float);

			//-if( (fdbackward_file = fopen(fdbackward_name, "wb")) == NULL ) 
			//-{ printf("Error fdbackward file\n"); }//return -1; }

			nframes = (nt%ndtrtm==0) ? nt/ndtrtm : (int) ceil( ((float)nt)/((float)ndtrtm) );

			//prcv = alloc3float(nframes,nz,nmigtrc);
			//if(prcv==NULL) { printf("Allocation of prcv[%d,%d,%d] failed!!!\n",nframes,nz,nmigtrc); }//return -1; }

			frames = (p_frames*) malloc(sizeof(p_frames)*nframes*nz*nmigtrc);

			//-prcv = alloc2float(nz,nmigtrc);
			//-if(prcv==NULL) { printf("Allocation of prcv[%d,%d] failed!!!\n",nz,nmigtrc); }//return -1; }

			iframe=0;
		}
		
		#pragma omp single nowait
		{
			printf("\tnmigtrc: %d   ",nmigtrc);
			printf("reclen: %d\n",reclen);

			printf("|###########################################################################|\n");

			printf("\tRETROPROPAGATION START\n");

			// medindo tempo de retropropagacao
			gettimeofday(&start,NULL);
		}


		for(it=0;it<nt;++it) // modelando a evolucao do campo de pressao //
		{
			// FD scheme: forward time //
				t = ((float)it) * dt;//((float)(it-1)) * dt;ndtrtm

			#pragma omp for schedule(guided)
			for(ix=(ixx0+DRVLEN);ix<(ixx1-DRVLEN);++ix)
			{
				for(iz=(DRVLEN-1);iz<(nzz-DRVLEN);++iz)
				{
					source = 0.0;

					for(itrc=0;itrc<nrec;++itrc)
					{
						if( (ix == (pos[itrc].igx-1)) && (iz == (pos[itrc].igz-1)) )
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
					laplacian = 2.0 * deriv2[0] * fields[(iz) + (ix*nzz)].p_2;

					for(iconv=1;iconv<DRVLEN;++iconv)
					{
						laplacian = laplacian + deriv2[iconv] 				* \
			      				(							  \
								fields[ (iz-iconv) + (ix*nzz) ].p_2 	+ \
								fields[ (iz+iconv) + (ix*nzz) ].p_2 	+ \
								fields[ (iz) + ((ix-iconv)*nzz) ].p_2	+ \
								fields[ (iz) + ((ix+iconv)*nzz) ].p_2 	  \
							);
					}

					fields[ (iz) + (ix*nzz) ].p_1 = invpgama * ( 	 2.0 * \
							    fields[ (iz) + (ix*nzz) ].p_2 - mgama * \
							    fields[ (iz) + (ix*nzz) ].p_1 + 	       \
						            fields[ (iz) + (ix*nzz) ].vel * ( laplacian - source ) );
				}
			}

			/*
			#pragma omp single
			{
				// mega - new swap fields
				p_aux  = p_2;
				p_2    = p_1;
				p_1    = p_aux;
			}
			*/

			#pragma omp for schedule(static,100)
			for(p_aux=&(fields[ (ixx0*nzz) ].p_1);p_aux<=&(fields[ (ixx1*nzz) ].p_1);p_aux+=3)
			{
				p_swap     = *(p_aux);
				*(p_aux)   = *(p_aux+1);
				*(p_aux+1) = p_swap;
			}
			
			
			/*
			#pragma omp for
			for(ix=ixx0;ix<ixx1;++ix)
			{
				for(iz=0;iz<nzz;++iz)
				{
					p_swap = fields[ (iz) + (ix*nzz) ].p_2;
					fields[ (iz) + (ix*nzz) ].p_2 = fields[ (iz) + (ix*nzz) ].p_1;
					fields[ (iz) + (ix*nzz) ].p_1 = p_swap;
				}
			}
			*/

		//#pragma omp single nowait
		//{
			if( ( it % ndtrtm) == 0)
			{
				//if( (iframe%100==0) )
				//	printf("\tit: %d - backward frames %d completed\n", it, iframe);
	  			
				//#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) \
					shared(nmigtrc,prcv,iframe,nz,p_2,ixmig0,nzz) private(ix)					
				//for(ix=0;ix<nmigtrc;++ix)
				//	memcpy(prcv + (iframe*nz*nmigtrc) + (ix*nz),
				//		p_2 + ((ixmig0-1+ix)*nzz) + (nborda), nz * sizeof(float) );

				#pragma omp for
				for(ix=0;ix<nmigtrc;++ix)
				{
					for(iz=0;iz<nz;++iz)
					{
						frames[(iframe*nz*nmigtrc) + (iz) + (ix*nz)].prcv = \
							fields[((ixmig0-1+ix)*nzz) + (iz) + (nborda)].p_2;
							//*(p_2 + ((ixmig0-1+ix)*nzz) + (iz) + (nborda));
					}
				}

				//-for(ix=0;ix<nmigtrc;++ix)
				//-	memcpy(prcv+(ix*nz), p_2 + ((ixmig0-1+ix)*nzz)+(nborda), nz*sizeof(float));

				//-fwrite(prcv, sizeof(float), nz * nmigtrc, fdbackward_file);
				#pragma omp single nowait
				{
					if( (iframe%100==0) )
						printf("\tit: %d - backward frames %d completed\n", it, iframe);
					++iframe;
				}
			}
		//}
		} // fim de modelagem

		#pragma omp single nowait
		{
			nframes = iframe;

			gettimeofday(&finish,NULL);

			timeval_subtract(&diff, &finish, &start);

			//-free2float(prcv);
			//-fclose(fdbackward_file);

			printf("\tRETROPROPAGATION OF %d FRAMES END\n", nframes);
			printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);

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
		}
		
		// 	condiçao inicial: repouso 		//
		// 	p(:,:,1) pressure field at time t-dt 	//
		// 	p(:,:,2) pressure field at time t 	//

		//memset(p_2, 0, sizeof(float)*nzz*nxx);
		//memset(p_1, 0, sizeof(float)*nzz*nxx);
		#pragma omp for schedule(static,100)
		for(ix=0;ix<nxx;++ix)
		{
			for(iz=0;iz<nzz;++iz)
			{
				fields[ (iz) + (ix*nzz)].p_1 = 0.0;
				fields[ (iz) + (ix*nzz)].p_2 = 0.0;
			}
		}	

		#pragma omp single nowait
		{
			//-prcv = alloc2float(nz,nmigtrc);
			//-if(prcv==NULL) { printf("Allocation of prcv[%d,%d] failed!!!\n",nz,nmigtrc); }//return -1;}

			// cache
			//psrc = alloc3float(nframes,nz,nmigtrc);
			//if(psrc==NULL) { printf("Allocation of psrc[%d,%d,%d] failed!!!\n",nframes,nz,nmigtrc); }//return -1;}
			
			//-psrc = alloc2float(nz,nmigtrc);
			//-if(psrc==NULL) { printf("Allocation of psrc[%d,%d] failed!!!\n",nz,nmigtrc); }//return -1;}		

			iframe = 0;

			gettimeofday(&start,NULL);
			//chunk = (ixx0-ixx1)/n_threads;;
		}

		// modelando a evolucao do campo de pressao //
		for(it=0;it<nt;it++)
		{
			delt  = pow( beta * ((float)(1+it-it0)), 2.0);

			// Pulso fonte Ricker frequancia pico freq //
			fonte = (dx*dx) * exp( -delt ) * ( 1.0 - (2.0 * delt) );   
		
			#pragma omp for schedule(static,1)
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
					laplacian = 2.0 * deriv2[0] * fields[ (iz) + (ix*nzz)].p_2;

					for(iconv=1;iconv<DRVLEN;++iconv)
					{
						laplacian = laplacian + deriv2[iconv] 			* \
			      				(					  	  \
								fields[ (iz-iconv) + (ix*nzz) ].p_2 	+ \
								fields[ (iz+iconv) + (ix*nzz) ].p_2 	+ \
								fields[ (iz) + ((ix-iconv)*nzz) ].p_2	+ \
								fields[ (iz) + ((ix+iconv)*nzz) ].p_2 	  \
							);
					}

					fields[ (iz) + (ix*nzz) ].p_1 = invpgama * ( 	 2.0 * \
							    fields[ (iz) + (ix*nzz) ].p_2 - mgama * \
							    fields[ (iz) + (ix*nzz) ].p_1 + 	       \
						            fields[(iz) + (ix*nzz)].vel * ( laplacian - source ) );
				}
			}
			
			/*
			#pragma omp single
			{
				// mega - new swap fields
				p_aux  = p_2;
				p_2    = p_1;
				p_1    = p_aux;
			}
			*/

			
			#pragma omp for schedule(static,100)
			for(p_aux=&(fields[ (ixx0*nzz) ].p_1);p_aux<=&(fields[ (ixx1*nzz) ].p_1);p_aux+=3)
			{
				p_swap     = *(p_aux);
				*(p_aux)   = *(p_aux+1);
				*(p_aux+1) = p_swap;
			}
			
			/*
			#pragma omp for
			for(ix=ixx0;ix<ixx1;++ix)
			{
				for(iz=0;iz<nzz;++iz)
				{
					p_swap = fields[ (iz) + (ix*nzz) ].p_2;
					fields[ (iz) + (ix*nzz) ].p_2 = fields[ (iz) + (ix*nzz) ].p_1;
					fields[ (iz) + (ix*nzz) ].p_1 = p_swap;
				}
			}
			*/
		//#pragma omp single nowait
		//{
			if( ( it % ndtrtm) == 0)
			{
				//if( (iframe%100==0) ) 
				//	printf("\tit: %d - forward frames %d completed\n", it, iframe);

				//#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none) \
					shared(nmigtrc,psrc,iframe,nz,p_2,ixmig0,nzz) private(ix)
				//for(ix=0;ix<nmigtrc;++ix)
				//	memcpy(psrc + (iframe*nz*nmigtrc) + (ix*nz), 
				//		p_2 + ((ixmig0-1+ix)*nzz) + (nborda), nz * sizeof(float) );

//PAREI AQUI
				/*
				#pragma omp for //fields[ (ixx1*nzz) ].p_1
				for(p_aux=&(frames[(iframe*nz*nmigtrc)].psrc);p_aux<=&(frames[(iframe*nz*nmigtrc)+(nzz)].psrc);p_aux+=2)
				{
					*(p_aux2) = fields[((ixmig0-1+ix)*nzz) + (nborda)].p_2;
					*(p_aux) = *(p_aux2);
					p_aux2++;
				}
				*/
				
				#pragma omp for
				for(ix=0;ix<nmigtrc;++ix)
				{
					for(iz=0;iz<nz;++iz)
					{
						frames[(iframe*nz*nmigtrc) + (iz) + (ix*nz)].psrc = \
							fields[((ixmig0-1+ix)*nzz) + (iz) + (nborda)].p_2;
							//*(p_2 + ((ixmig0-1+ix)*nzz) + (iz) + (nborda));
					}
				}
				

				//-if(write_forward_file)
				//-	fwrite(psrc, sizeof(float), nz * nmigtrc, fdforward_file);
				
				#pragma omp single nowait
				{
					if( (iframe%100==0) ) 
						printf("\tit: %d - forward frames %d completed\n", it, iframe);
					++iframe;
				}
			}
		//}
		} // fim de modelagem //

		#pragma omp single nowait
		{
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

			gettimeofday(&start,NULL);

			//chunk=nmigtrc/n_threads;
		}

		for(iframe_cross=0;iframe_cross<nframes;++iframe_cross)
		{

			//-fseek(fdforward_file, iframe_cross * reclen , SEEK_SET);
			//-fread(psrc, sizeof(float), nz * nmigtrc, fdforward_file);

			irec = nframes-iframe_cross-1;

			//-fseek(fdbackward_file, irec * reclen , SEEK_SET);
			//-fread(prcv, sizeof(float), nz * nmigtrc, fdbackward_file);

			//////////////////////////////////////////////////////////////////
			//			NOVA CROSS-CORRELACAO 			//
			//////////////////////////////////////////////////////////////////

			#pragma omp for schedule(static,1)
			for(ix=(ixmig0-1-nborda);ix<(ixmig1-nborda);++ix)
			{
				itrc = (ix-ixx0) % (nmigtrc);

				//for(iz=0;iz<nz;++iz)
				//	imag[(iz) + (ix*nz)]  = imag[(iz) + (ix*nz)]    + ( window_x[ix] * \
							psrc[ (iframe_cross*nz*nmigtrc) + (iz) + (itrc*nz)] * \
							prcv[ (irec*nz*nmigtrc)         + (iz) + (itrc*nz)] );

				for(iz=0;iz<nz;++iz)
					sections[(iz) + (ix*nz)].imag  = sections[(iz) + (ix*nz)].imag + ( window_x[ix] * \
							frames[(iframe_cross*nz*nmigtrc) + (iz) + (itrc*nz)].psrc * \
							frames[(irec*nz*nmigtrc) + (iz) + (itrc*nz)].prcv );
				//++itrc;
			}


			//-itrc = 0;
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
		
		#pragma omp single nowait
		{
			printf("\tCROSS CORRELATION END\n");
			gettimeofday(&finish,NULL);

			timeval_subtract(&diff, &finish, &start);
			printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);
	
			//-free2float(prcv);
			//-free2float(psrc);

			//&free3float(prcv);
			//&free3float(psrc);
			free(frames);

			//-fclose(fdbackward_file);

			//-if(write_forward_file)
			//-	fclose(fdforward_file);
		

		//////////////////////////////////////////////////////////////////////////////////
        	//				FILTRO LAPLACIANO				//
		//////////////////////////////////////////////////////////////////////////////////
			printf("|###########################################################################|\n");
			printf("\tLAPLACIAN FILTER START\n");
			gettimeofday(&start,NULL);		
		}

		#pragma omp for	schedule(guided)
		for(ix=(DRVLEN-1);ix<(nx-DRVLEN);ix++)
		{
			for(iz=(DRVLEN-1);iz<(nz-DRVLEN);iz++)
			{
				laplacian = 2.0 * deriv2[0] * sections[(iz) + (ix*nz)].imag;

				for(iconv=1;iconv<DRVLEN;++iconv)
				{
					laplacian = laplacian + deriv2[iconv] 			* \
						(						  \
							sections[(iz-iconv) + (ix*nz)        ].imag	+ \
							sections[(iz+iconv) + (ix*nz)        ].imag	+ \
							sections[(iz)       + ((ix-iconv)*nz)].imag	+ \
							sections[(iz)       + ((ix+iconv)*nz)].imag  	  \
						);
				}

				sections[(iz) + (ix*nz)].imag_filter = laplacian;
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
        	//				ESCRITA FINAL DOS DADOS				//
		//////////////////////////////////////////////////////////////////////////////////
		#pragma omp single nowait
		{
			printf("\tLAPLACIAN FILTER END\n");
			gettimeofday(&finish,NULL);

			timeval_subtract(&diff, &finish, &start);
			printf("%ld.%06ld\n", diff.tv_sec, diff.tv_usec);

			if(ishot==nshots-1)
			{

				if( (out_file = fopen(out_name, "wb")) == NULL ) { printf("Error opening out file\n"); }//return -1; }
				//#pragma omp parallel for if (is_parallel) num_threads(n_threads) default(none)	\
					shared(nz,nx,out_file,imag_filter) private(ix) firstprivate(imagtrace)
				for(ix=0;ix<nx;++ix)
				{
					imagtrace.tr_header->tracl = ix + 1;

					for(iz=0;iz<nz;++iz)
					{
						imagtrace.tr_data[iz] = sections[(iz)+(ix*nz)].imag_filter;
					}

					//memcpy(imagtrace.tr_data, imag_filter + (ix*nz), nz * sizeof(float) );
					//put_tr(ix, &imagtrace, nz, out_file);
				
					//&fwrite(imagtrace.tr_header, HDRBYTES, 1, out_file);
					//&fwrite(imag_filter + (ix*nz), sizeof(float), nz, out_file);

					fwrite(imagtrace.tr_header, HDRBYTES, 1, out_file);
					fwrite(imagtrace.tr_data, sizeof(float), nz, out_file);
				}

				fclose(out_file);
			}

			printf("|###########################################################################|\n");
		}

		} // end of parallel section

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
	
	//&free1int(igz);
	//&free1int(igx);
	free(pos);

	free1float(gama_x);
	free1float(gama_z);
	free1float(window_x);
	free1float(trace.tr_data);
	free1float(imagtrace.tr_data);

	//&free2float(imag);
	//&free2float(imag_filter);
	free(sections);

	//free2float(vel);
	//free2float(p_1);
	//free2float(p_2);
	free(fields);

	//TODO: dealocate all shotgather
	free1su_trace(shotgather);

    	return 0;
}

