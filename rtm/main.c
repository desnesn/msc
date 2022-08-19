/********************************************************************************/
/*	RTM Migration - Baseado no codigo em Fortran do Prof. Jesse Costa  	*/
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

#include <stdio.h>
#include <stdlib.h>

#include "io.h"

#define UNUSED __attribute__ ((unused))

int main(int argc, char *argv[])
{
	//////////////////////////////////////////////////////////////////////////////////
	//				RTM INIT					//
	//////////////////////////////////////////////////////////////////////////////////    

	// variaveis nao usadas //
	// 1 caso queira gravar campo propagado forward //
	// const unsigned int MVFWRD = 0; // 1 caso queira gravar campo propagado forward //
	// int ierr,izf, ixf, swap, itrcmv idx;
	// float swap;
	// double xf, zf;
	// Unidades de arquivo //
	//
	// int MDID=72;
	// int MVID=61;
	// int MVFWID=63;
	// int SUID=90;
	// int SUCSID=81;
	// int SUIMGID=69;
	//
	// variaveis nao usadas //

	// largura da borda da ordem de 3xlambda/2, lambda=Vmax/freq //
	const unsigned int UNUSED nborda = 50;
	
	// LFRAC: {0,1,2,3,4}
	// 0 - migra cada tiro usando toda extensao do modelo de velocidade
	// caso contrario a janela de migracao limitada ao intervalo do modelo de
	// velocidade: lanco/NFRAC+lanco+lanco/NFRAC
	const unsigned int UNUSED LFRAC = 3;

	// intervalo para gravacao dos frames na migracao RTM //
	const float UNUSED dtrtm = 0.008;// seconds

	double UNUSED *gama_x;
	double UNUSED *gama_z;
	double UNUSED *window_x;

	double UNUSED **vel;
//	double **sourcemask;

	double UNUSED ***p;

	double UNUSED dswap;

	float UNUSED **psrc;
	float UNUSED **prcv;
	float UNUSED **imag;
	float UNUSED **imag_filter;

	// operador de segunda derivada de oitava ordem //
	// serve para resolver numericamente a equacao da onda (EDP) //
	const unsigned int UNUSED DRVLEN = 5;
	double UNUSED deriv2[5] = {-2.84722E+00,1.60000E+00,-2.00000E-01,2.53968E-02,-1.785714E-03};

	// blackmann window 
	// utiliza um taper nas bordas de cada imagem suavizando as derivadas
	// janela Blackmann harris calculada para decaimento suave nas bordas das derivadas
	const unsigned int UNUSED WDWLEN = 5;
	double UNUSED blackmann[5] = {0.893012701892219E+00,0.63000E+00,0.34000E+00,0.130000E+00,0.026987298107781E+00};
	
	//su_trace imagtrace = {0};
	//su_trace trace = {0};
	su_trace trace,imagtrace;
	memset(&imagtrace, 0, sizeof(su_trace));
	memset(&trace, 0, sizeof(su_trace));

	su_trace UNUSED *shotgather;

	int UNUSED reclen;
	int UNUSED nrec, ilanco, nmigtrc, ixmig1, ixmig0, ixx0, ixx1;

	const float UNUSED sec2mcs = 1.0e6;

	double UNUSED beta, delt, gama, invpgama, mgama, laplacian, fonte, source;
	float UNUSED ttotal, t, vmin, vmax, dxmax;
	double UNUSED dx, dt, dz, freq;
	int UNUSED nx, nz, nt, z0, x0, it;
	float UNUSED dtrec;
	int UNUSED ndtrec, ndtrtm;
	int UNUSED ix, iz, it0;
	int UNUSED isx, isz;
	int UNUSED *igx , *igz;
	int UNUSED *shotidx;
	int UNUSED iconv, ixb, ixe, izb, ize, nxx, nzz;
	int UNUSED iframe, nframes, irec, ircv, itrc, ns;
	int UNUSED ntrc, ntrcmax, ishot, nshots;
	float UNUSED factor;

	FILE *su_file;		// complete su data - marmousi - 240 tiros - 23040 traços
    	char *f_name = "MARMOUSI_xtLINUX_CORRECT2.su";

	if( (su_file = fopen(f_name, "rb")) == NULL )
    	{
        	printf("Error opening seismic data file\n");
	        return -1;
    	}

	FILE *vel_file;
	char *vel_name = "velocity.h@";// vel[751][2301]
//	char *vel_name = "velocity_smooth_05_05.bin";

	if( (vel_file = fopen(vel_name, "rb")) == NULL )
    	{
        	printf("Error opening vel. model file\n");
        	return -1;
    	}

	// Dados do modelo de velocidade //
    	nz = 751;    //n1
    	nx = 2301;   //n2
    	z0 = 0.0;    //??
    	x0 = 0.0;    //??
    	dz = 4.0;    //d1
    	dx = 4.0;    //d2
    	// Velocity model constants //

	// Frequencia da fonte //
	freq = 12.0;

	// PML nos 4 lados do grid // - Atenuar a energia das bordas com uma exponencial quadratica da velocidade
	nxx=nx+2*nborda;  // Incliur o modelo incluindo as bordas para x e z
	nzz=nz+2*nborda;
	ixb=nborda+1;
	ixe=nborda+nx;
	izb=nborda+1;
	ize=nborda+nz;

	printf("\nInitial parameters\n");
	printf("freq: %.3f\n",freq);    //freq: 12
	printf("nz: %d\n",nz);          //nz: 751
	printf("nx: %d\n",nx);          //nx: 2301
	printf("nzz: %d\n",nzz);        //nzz: 801
	printf("nxx: %d\n",nxx);        //nxx: 2401
	printf("nborda: %d\n",nborda);  //nborda: 50
	printf("ixb: %d\n",ixb);        //ixb: 51
	printf("ixe: %d\n",ixe);        //ixe: 2351
	printf("izb: %d\n",izb);        //izb: 1
	printf("ize: %d\n",ize);        //ize: 751

	// Alocacao de vetores e matrizes //

	p = alloc3double(nzz,nxx,2);    // Campos de pressão: p[851][2401][2]
	if(p==NULL) { printf("Allocation of p[%d,%d,%d] failed!!!\n",nzz,nxx,2); return -1;}

    	vel = alloc2double(nzz,nxx);      // Matriz de velocidade: vel[851,2401]
	if(vel==NULL) { printf("Allocation of vel[%d,%d] failed!!!\n",nzz,nxx); return -1;}

    	gama_z = alloc1double(nzz);          
	if(gama_z==NULL) { printf("Allocation of gama_z[%d] failed!!!\n",nzz); return -1;}

    	gama_x = alloc1double(nxx);
	if(gama_x==NULL) { printf("Allocation of gama_x[%d] failed!!!\n",nxx); return -1;}

	window_x = alloc1double(nx);
	if(window_x==NULL) { printf("Allocation of window_x[%d] failed!!!\n",nx); return -1;}

	psrc = alloc2float(nz,nx);
	if(psrc==NULL) { printf("Allocation of psrc[%d,%d] failed!!!\n",nz,nx); return -1;}

	prcv = alloc2float(nz,nx);
	if(prcv==NULL) { printf("Allocation of prcv[%d,%d] failed!!!\n",nz,nx); return -1;}

	imag = alloc2float(nz,nx);
	if(imag==NULL) { printf("Allocation of imag[%d,%d] failed!!!\n",nz,nx); return -1;}

	imag_filter = alloc2float(nz,nx);
	if(imag_filter==NULL) { printf("Allocation of imag_filter[%d,%d] failed!!!\n",nz,nx); return -1;}

	imagtrace.tr_data = alloc1float(nz);
	if(imagtrace.tr_data==NULL) { printf("Allocation of imagtrace.tr_data[%d] failed!!!\n",nz); return -1;}

//    	sourcemask = alloc2double(nzz,nxx);      // Injeção da fonte (sourcemask)
//	if(sourcemask==NULL) { printf("Allocation of sourcemask[%d,%d] failed!!!\n",nzz,nxx); return -1;}

	//////////////////////////////////////////////////////////////////////////////////    
	//		Leitura e criação de bordas no modelo de velocidades		//
	//////////////////////////////////////////////////////////////////////////////////    

	// Leitura do modelo de velocidades //
	get_vel_model(vel, vel_file, ixb, ixe, izb, ize, nborda);
    	fclose(vel_file);
	// Obtendo o valor maximo e minimo do modelo de velocidades //
	get_min_max(vel,&vmin,&vmax, ixb, ixe, izb, ize);
	
	// preenchendo  o modelo de velocidades nas bordas //
	for(ix=0;ix<nborda;ix++)
    	{
        	for(iz=izb-1;iz<ize;iz++)
        	{ 
            		vel[iz][ix]     =  vel[iz][ixb-1]; 
            		vel[iz][ixe+ix] =  vel[iz][ixe-1];
        	}
   	}

	for(ix=0;ix<nxx;ix++)
    	{
        	for(iz=0;iz<(izb-1);iz++)
			vel[iz][ix]     =  vel[izb-1][ix]; 
		for(iz=ize;iz<nzz;iz++)
            		vel[iz][ix] 	=  vel[ize-1][ix];
        	
   	}
	
	printf("\nModelo de velocidades\n");
	printf("nz: %d\n",nz);
	printf("nx: %d\n",nx);
	printf("z0: %d\n",z0);
	printf("x0: %d\n",x0);
	printf("dz: %.2f\n",dz);
	printf("dx: %.2f\n",dx);
	printf("vmin: [%.2f]\nvmax: [%.2f]\n",vmin, vmax);

	//////////////////VERIFICAÇÃO DO MODELO DE VELOCIDADE COM BORDA///////////////////
	//	 	valores: ix[0;8], ix[46;54], ix[2340;2354], ix[2392,2401]	//
	// 			valores: iz[0;70], iz[730;nzz]				//
	//////////////////////////////////////////////////////////////////////////////////    

	//
	//for(ix=2392;ix<2401;ix++)
	//	printf("\t[%d]",ix);
	//printf("\n");
	//
	//for(iz=730;iz<nzz;iz++)//Lado inferior esquerdo do modelo
    	//{
        //	printf("[%d]: ",iz);
        //	for(ix=2392;ix<2401;ix++)
        //  		printf("%.2f\t",vel[iz][ix]);
        //	printf("\n");
    	//}
	//

	//////////////////////////////////////////////////////////////////////////////////    

	//////////////////////////////////////////////////////////////////////////////////
        //			Getting relevant data from first trace 			//
	//////////////////////////////////////////////////////////////////////////////////

	//int min_isx_igx;	// min(sx,gx)			//

  	void *ptr;
	ptr = (su_header*) malloc( sizeof(su_header) );
	fread(ptr, HDRBYTES, 1, su_file);

	isx = ((su_header*)ptr)->sx;

	//
	//get_sx_gx((su_header*)ptr, &isx, &igx);
	//min_isx_igx = MIN(isx,igx);
	//isx = isx - min_isx_igx;
	//igx = igx - min_isx_igx;
	//nxshot=240;//number of shots

	//ng=96;//number  of geophones
    	//nt = ((segy*) ptr)->ns;

	//    dx = ( ((segy*) ptr)->d2 ) ? ((segy*) ptr)->d2 : 1.0;// sample spacing between traces
	//

	ns = ((su_header*)ptr)->ns;
	dtrec  = ((float)(((su_header*)ptr)->dt))/sec2mcs; // Intervalo de amostragem dos dados //
	ttotal = ((float)ns-1)*dtrec; // POG: ((float)(ns-1))*dtrec; //
	
	printf("\nData from first trace\n");
	printf("ns: %d\n",ns);
	printf("dt: %lf\n",(float)(((su_header*)ptr)->dt));
	printf("sec2mcs: %lf\n", sec2mcs);
	printf("dtrec: %lf\n", dtrec);
	printf("ttotal: %lf\n", ttotal); 

	// fator de escala para coordenadas scalco
	factor = 1.0;

	trace.tr_data = alloc1float(ns);
	if(trace.tr_data==NULL) { printf("Allocation of trace.tr_data[%d] failed!!!\n",ns); return -1;}

	printf("(((su_header*)ptr)->scalco): %d\n",(((su_header*)ptr)->scalco));

	if ( (((su_header*)ptr)->scalco) < 0 )
		//!factor = 1.0 / (10.0**abs(trace%trcheader%scalco)) ! SU manual
		factor = 1.0 / fabs((float)(((su_header*)ptr)->scalco)); //!SEG-Y standard
	else if( (((su_header*)ptr)->scalco) > 0 )
		//factor = 10.0**abs( trace%trcheader%scalco) ! SU manual
		factor = (float) (((su_header*)ptr)->scalco);// ! SEG-Y standard

	printf("factor: %.3f\n",factor);

	free(ptr);
    	fseek(su_file,0,SEEK_SET);

	//////////////////////////////////////////////////////////////////////////////////
        //		Reading dynamic data and info from su file	 		//
	//////////////////////////////////////////////////////////////////////////////////

	int oldisx=0;		// 	old sx position		//
	//int oldigx=0;		// 	old sx position		//

	ntrc=0;
	nshots=1;

	get_tr(ntrc,&(trace.tr_header), trace.tr_data, ns, su_file);

       	oldisx=trace.tr_header.sx;
        //oldigx=igx;

	// Testing data read //
	//
	//printf("\nTesting reading data from the nth trace\n");
	//printf("hdr.sx=%d\n",trace.tr_header.sx);
	//printf("hdr.sy=%d\n",trace.tr_header.sy);
        //printf("hdr.gx=%d\n",trace.tr_header.gx);
        //printf("hdr.gy=%d\n\n",trace.tr_header.gy);

	//for(iz=700;iz<ns;iz++)
    	//{
      	//	printf("[%d]: %.3f\n",iz,trace.tr_data[iz]);
    	//}
	//

        // Finding out the number of traces in the data //
        do{
            	//get_sx_gx(trace.tr_header, &sx,&gx);
 
            	//isx = (isx - min_sx_gx);
		//igx = (igx - min_sx_gx);
		isx = trace.tr_header.sx;

		ntrc++;

		if(isx!=oldisx) 
		{
			nshots++;
			oldisx=isx;
		}
        }while( get_tr(ntrc, &(trace.tr_header),trace.tr_data, ns, su_file) != 0 );

	fseek(su_file,0,SEEK_SET);

	printf("\nData from gather\n");
	printf("Number of traces: %d\n",ntrc);
	
	//////////////////////////////////////////////////////////////////

	shotidx = alloc1int(nshots+1);
    	if(shotidx==NULL) { printf("Allocation of shotidx[%d] failed!!!\n",nshots); return -1;}

	ntrc=0;	
	nshots=1;

	get_tr(ntrc, &(trace.tr_header), trace.tr_data, ns, su_file);
	
	oldisx=trace.tr_header.sx;
        //oldigx=igx;

	shotidx[0]=0;

        // Finding out the number of geofones and shots in the data
        do{
            	//get_sx_gx(trace.tr_header, &sx,&gx);
 
            	//isx = (isx - min_sx_gx);
		//igx = (igx - min_sx_gx);
		isx = trace.tr_header.sx;

		if(isx!=oldisx) 
		{	
//			printf("ntrc: %d\n",ntrc);
			shotidx[nshots] = ntrc;
			nshots++;
			oldisx=isx;
		}
		ntrc++;
        }while( get_tr(ntrc, &(trace.tr_header),trace.tr_data,ns,su_file) != 0 );
	
	shotidx[nshots] = ntrc;

	// First trace from each shot position
	//
	//for(iz=0;iz<nshots;iz++)
	//	printf("shotidx[%d]: %d\n",iz,shotidx[iz]);
	//

	fseek(su_file,0,SEEK_SET);

	ntrcmax=ntrc/nshots;

	//						// DESNECESSARIO AO MEU VER - ntrcmax //
	//printf("meu ntrcmax: %d\n\n",ntrcmax);
	//ntrcmax=0;
	//for(ishot=0;ishot<nshots;ishot++)
	//{
	//	ntrcmax = MAX(ntrcmax,shotidx[ishot+1]-shotidx[ishot]);
	//	printf("[%d] - ntrcmax: %d\n",ishot,ntrcmax);
	//}
	//

	printf("Number of shots: %d\n",nshots);
	printf("Number of geophones: %d\n",ntrcmax);

    	shotgather = alloc1su_trace(ntrcmax);          
	if(shotgather==NULL) { printf("Allocation of shotgather[%d] failed!!!\n",ntrcmax); return -1;}

	for(itrc=0;itrc<ntrcmax;itrc++)
	{	
		shotgather[itrc].tr_data = alloc1float(ns);
		if(shotgather[itrc].tr_data==NULL) 
		{ printf("Allocation of shotgather[itrc].tr_data[%d] failed!!!\n",ntrcmax); return -1; }
	}

	igz = alloc1int(ntrcmax); // correcao devida a profundidade da fonde e do receptor numa levantamento marinho //
	if(igz==NULL) { printf("Allocation of igz[%d] failed!!!\n",ntrcmax); return -1;}

	igx = alloc1int(ntrcmax); // correcao devida a profundidade da fonde e do receptor numa levantamento marinho //
	if(igx==NULL) { printf("Allocation of igx[%d] failed!!!\n",ntrcmax); return -1;}

	// check grid interval com o objetivo de evitar dispersao numerica, e garantir acuidade!! //
	// limitar valores referentes a amostragem do campo em relacao ao dx e dt max //
	dxmax    = vmin / ( (float) ( ((double)6.0) * freq) );
	dt       = dx   / ( ((double) 4.0) * ((double) vmax) );
	ndtrec   = (int) ceil( ((double)dtrec)/dt );
	ndtrtm   = (int) ceil( ((double)dtrtm)/dt );

	printf("\ndxmax: %lf\n",dxmax);
	printf("dt: %lf\n",dt);
	printf("dtrec: %lf\n",dtrec);
	printf("dtrtm: %lf\n",dtrtm);
	printf("ndtrec: %d\n",ndtrec);
	printf("ndtrtm: %d\n",ndtrtm);

	// escolha dt tal que o tempo de registro seja multiplo intereiro de dt //
	if(ndtrec > 1)
   		dt = ((double)dtrec)/((double)(ndtrec));
	else
   		dt = ((double)dtrec);
	
	printf("novo dt: %.12lf\n",dt);

	nt = (int) ceil( ((double)ttotal)/dt ) - 1;// -1 == POG	

	if(dx > dxmax) {printf("Freq muito alta!!!"); return -1;}

	printf("\nMigration Parameters\n");
	printf("ttotal: %f\n", ttotal); 
	printf("nt: %d\n",nt);

    	for(ix=0;ix<nxx;ix++) // (v[z][x])^2
    	{
        	for(iz=0;iz<nzz;iz++)
        	{
            		vel[iz][ix] = pow( (vel[iz][ix]*dt/dx), (double) 2);
        	}
    	}

	// pulso fonte - ricker //

	beta = D_PI*freq*dt;
    	it0 = (int) ceil( ((double)1.0)/(freq*dt) );

	printf("beta: %lf\t",beta);	// ???
	printf("it0: %d\n",it0);	// ???

//
	//! mascara de absorcao:
    	// gama_x(:) = 0.0_dp    
    	// gama_z(:) = 0.0_dp
	//
	//for(iz=0;iz<nzz;iz++)
	//	gama_z[iz] = ((double)0.0);
	//
	//for(ix=0;ix<nxx;ix++)
	//	gama_x[ix] = ((double)0.0);
	//

	for(ix=1;ix<(nborda+1);ix++)
    	{
        	gama = beta * pow( ( ((double)ix)/((double)nborda) ), 2);
	       	//gama_x[ixb-ix-1] = gama;
        	//gama_x[ixe+ix-1] = gama;
        	gama_z[izb-ix-1] = gama;
        	gama_z[ize+ix-1] = gama;
    	}	

// 
//    // Testing gamas //
//    for(ix=1;ix<(nborda+1);ix++)
//        printf("gama_x[%d]: %lf\tgama_x[%d]: %lf\n",ixb-ix-1,gama_x[ixb-ix-1],ixe+ix-1,gama_x[ixe+ix-1]);
//    printf("\n");
//    for(iz=1;iz<(nborda+1);iz++)
//       printf("gama_z[%d]: %lf\tgama_z[%d]: %lf\n",izb-iz-1,gama_z[izb-iz-1],ize+iz-1,gama_z[ize+iz-1]);
//
//

    	// header do traco imagem //
    	imagtrace.tr_header.ns = (unsigned short) nz;
    	imagtrace.tr_header.dt = (unsigned short) (sec2mcs*dtrtm); //TODO: =int(sec2mcs*dtrtm,2)
    	imagtrace.tr_header.f1 = (float) z0;
    	imagtrace.tr_header.f2 = (float) x0;
    	imagtrace.tr_header.d1 = (float) dz;
    	imagtrace.tr_header.d2 = (float) dx;

//
//    printf("\nDados do Imagtrace\n");
//    printf("imagtrace.tr_header.ns: %d\n",imagtrace.tr_header.ns);
//    printf("imagtrace.tr_header.dt: %d\n",imagtrace.tr_header.dt);
//    printf("imagtrace.tr_header.f1: %f\n",imagtrace.tr_header.f1);
//    printf("imagtrace.tr_header.f2: %f\n",imagtrace.tr_header.f2);
//    printf("imagtrace.tr_header.d1: %f\n",imagtrace.tr_header.d1);
//    printf("imagtrace.tr_header.d2: %f\n",imagtrace.tr_header.d2);
//

	//////////////////////////////////////////////////////////////////////////////////
        //		Inicializacao da imagem RTM                     	 	//
	//////////////////////////////////////////////////////////////////////////////////

    	// imag(:,:) = 0.0 //
    	// imag_filter(:,:) = 0.0 //
	//
	//for(iz=0;iz<nz;iz++)
	//{
	//	for(ix=0;ix<nx;ix++)
	//	{
	//		imag[iz][ix]=((float)0.0);
	//		imag_filter[iz][ix]=((float)0.0);
	//	}
	//}
	//

//	int ntrace;

	//float pt_aux;
	float *flt_vet;
    	//flt_vet = (float*) malloc(sizeof(float)*nz);
	flt_vet = (float*) calloc(sizeof(float),nz);

	FILE *cshot_file;
	char *cshot_name = "cshot.su";

	FILE *fdbackward_file;
	char *fdbackward_name = "fdbackward.bin";

	FILE *fdforward_file;
	char *fdforward_name = "fdforward.bin";

	FILE *out_file;		// saida - dados migrados
	char *out_name = "rtm_migrated.su";

	if( (out_file = fopen(out_name, "ab")) == NULL )
    	{
        	printf("Error opening rtm_migrated file\n");
        	return -1;
    	}

    	ilanco=0;

	int ntrmig=0; // minha variavel para escrita

	// Loop sobre shots //
    	for(ishot=0;ishot<1;ishot++) // nshots
    	{
        	// numero de tracos para common-shot //
		ircv = 0;
		nrec = shotidx[ishot+1]-shotidx[ishot]; // ntrc/nshots; // 96

		if( (cshot_file = fopen(cshot_name, "ab")) == NULL )
    		{
        		printf("Error opening shot file\n");
	        	return -1;
    		}


        	for(itrc=shotidx[ishot];itrc<shotidx[ishot+1];itrc++) // ??? NAO MIGRA O ULTIMO TIRO ??? //
        	{
			//printf("ircv: %d\titrc: %d\n",ircv,itrc);
			
			// lendo do dado e escrevendo o commom-shot //
            		get_tr(itrc, &(trace.tr_header), trace.tr_data, ns, su_file);
			put_tr(ircv, &(trace.tr_header), trace.tr_data, ns, cshot_file);
			
		      	shotgather[ircv].tr_header     = trace.tr_header;
			//printf("[%d]shotgather[%d].tr_header.gx: %d\t",itrc,ircv,shotgather[ircv].tr_header.gx);
			//printf("[%d]shotgather[%d].tr_header.sx: %d\t",itrc,ircv,shotgather[ircv].tr_header.sx);

			copy1float(shotgather[ircv].tr_data, trace.tr_data, ns);
			//printf("copy1float\n");

			// armazena posicao do traco ma malha do modelo //
			// converte unidades usando a palavra scalco do header //
			// ===============ATENCAO
			// as palavras gelev e sdepth devem ser setadas no header // isso eu vejo aonde ???
			// =========================================================
			igx[ircv] = (int) floor( ( ((float)(trace.tr_header.gx)) - ((float)x0) ) * factor/((float)dx) ) + ixb;

			//
			//printf("igx[%d]: %d\n", ircv, igx[ircv]);
			//printf("((float)(trace.tr_header.gx)): %f\t", (float)(trace.tr_header.gx) );
			//printf("x0: %f\t",((float)x0));
			//printf("factor: %f\t",factor);
			//printf("dx: %f\t",dx);
			//printf("ixb: %d\n",ixb);
			//

			igz[ircv] = (int) floor( ((float)(-1.0)) * ((float)(trace.tr_header.gelev)) * factor/((float)dz) ) + izb;
			//printf("igz[%d]: %d\n",ircv,igz[ircv]);
		
			isx 	  = (int) floor( ( ((float)(trace.tr_header.sx)) - ((float)x0) ) * factor/((float)dx) ) + ixb;
			//printf("isx: %d\t",isx);

			isz       = (int) floor( ((float)(trace.tr_header.sdepth)) * factor/((float)dz) ) + izb;
			//printf("isz: %d\n",isz);

			ilanco    = max(ilanco,fabs(isx-igx[ircv]));
			
			//printf("ilanco: %d\n", ilanco);

	 		ircv++;
			// TESTAR TUDO ISSO: OK
    		} // fim do loop dos tracos

		fclose(cshot_file);

		// DEFINIR JANELA DE MIGRACAO //

		if( LFRAC == 0 )
		{
			ixmig0 = 1;
			ixmig1 = nx;
		}
		else
		{
			ixmig0 = max( min(minval(igx,nrec),isx)-ilanco/LFRAC, ixb);
			//for(ix=0;ix<96;ix++)	printf("igx[%d]: %d\t",ix,igx[ix]);printf("\n");
			//printf("\nMIN(minval(igx,nrec),isx): %d\n",MIN(minval(igx,nrec),isx));
			//printf("minval(igx,nrec): %d\n",minval(igx,nrec));
			//printf("isx: %d\n",isx);
			//printf("ilanco: %d\tLFRAC: %d\tilanco/LFRAC: %d\n",ilanco, LFRAC, ilanco/LFRAC);
			//printf("MIN(minval(igx,nrec),isx)-ilanco/LFRAC: %d\n",MIN(minval(igx,nrec),isx)-ilanco/LFRAC);
			//printf("ixb: %d\n",ixb);
			//printf("MAX(-57,51): %d\n",MAX(-57,51));
			//printf("MAX( MIN(minval(igx,nrec),isx)-ilanco/LFRAC, ixb): %d\n",MAX( (int)(MIN(minval(igx,nrec),isx)-ilanco/LFRAC), ixb));
			ixmig1 = min( max(maxval(igx,nrec),isx)+ilanco/LFRAC, ixe);
			//printf("\nMAX(maxval(igx,nrec),isx): %d\n",max(maxval(igx,nrec),isx) );
			//printf("maxval(igx,nrec): %d\n",maxval(igx,nrec));
			//printf("isx: %d\n",isx);
			//printf("ilanco/LFRAC: %d\n",ilanco/LFRAC);
			//printf("MAX(maxval(igx,nrec),isx)+ilanco/LFRAC: %d\n",max(maxval(igx,nrec),isx)+ilanco/LFRAC);
			//printf("ixe: %d\n",ixe);
			//printf("MIN( MAX(maxval(igx,nrec),isx)+ilanco/LFRAC, ixe): %d\n",min( max(maxval(igx,nrec),isx)+ilanco/LFRAC, ixe));
		}

		printf("[%d][%d]ixmig0: %d\tixmig1: %d\t\n",itrc,ircv,ixmig0,ixmig1); // TODO: REVER ISTO... oscilação do tamanho da janela de mig

   		nmigtrc = ixmig1 - ixmig0 + 1;

		//printf("nmigtrc: %d\n",nmigtrc);

		// DEFINIR JANELA DE PROPAGACAO //
		memset(gama_x,0,sizeof(double)*nxx); // zerar gama_x????    gama_x(:) = 0.0 

		//for(ix=0;ix<nxx;ix++)
		//	gama_x[ix] = ((double)0.0);

		ixx0 = ixmig0 - nborda - 1;
		ixx1 = ixmig1 + nborda;
		
		printf("ixx0: %d\t",ixx0);
		printf("ixx1: %d\n",ixx1);

		for(ix=1;ix<nborda+1;ix++)
		{
			gama = beta * pow( ( ((double)ix)/((double)nborda) ) , 2 );
			gama_x[ixmig0-ix-1] = gama;
			//printf("gama_x[%d]=%f\t",ixmig0-ix-1,gama);
			gama_x[ixmig1+ix-1] = gama;
			//printf("gama_x[%d]=%f\n",ixmig1+ix-1,gama);
		}

		// JANELA BLACKMANN para suavizar borda da imagem //
		for(ix=0;ix<nx;ix++) window_x[ix] = ((double)1.0); // window_x(:) = 1.0

		//printf("\n");
		for(ix=0;ix<WDWLEN;ix++)
		{
			window_x[ixmig0-1-nborda+WDWLEN-1-ix] = blackmann[ix];
			//printf("window_x[%d]=%lf\t",ixmig0-nborda+WDWLEN-ix, blackmann[ix]);
			window_x[ixmig1-nborda-WDWLEN+ix] = blackmann[ix];
			//printf("window_x[%d]=%lf\n",ixmig1-nborda-WDWLEN+ix, blackmann[ix]);
		}

		//for(ix=0;ix<nx;ix++) printf("window_x[%d]: %lf\n",ix,window_x[ix]);

		
		//printf("ishot-base: [%d;%d]\n",ishot, ishot+1);
		//printf("OK: trace[%d;%d]\n",(shotidx[ishot]-1),(shotidx[ishot+1]-1));
		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		/////	       RETROPROPAGACAO DO CAMPO DOS RECEPTORES: INICIO		/////
		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////

		printf("\nRETROPROPAGATION\n");

		reclen = nz*nmigtrc*sizeof(float);
		printf("nz: %d\n",nz);
		printf("nmigtrc: %d\n",nmigtrc);
		printf("reclen: %d\n",reclen);

		printf("Migrando tiro %d de %d\n",ishot+1, nshots);
		printf("Janela de Migracao: nz[%d] e nx[%d]\n",nz, nmigtrc);
		printf("Numero de Receptores: %d\n", nrec);

		if( (fdbackward_file = fopen(fdbackward_name, "wb")) == NULL )
    		{
        		printf("Error opening shot file\n");
	        	return -1;
    		}

		// condiçao inicial: repouso //
		// p(:,:,1) pressure field at time t-dt //
		// p(:,:,2) pressure field at time t //

		// p(:,:,:)= 0.0_dp // zerar p
		for(it=0;it<2;it++)
			for(iz=0;iz<nzz;iz++)
				for(ix=0;ix<nxx;ix++)
					p[iz][ix][it] = ((double)(0.0));

		irec=0;
		//itrcmv=0;
		iframe=0;

		//printf("nt: %d\n",nt);
		
		for(it=0;it<nt;it++) // modelando a evolucao do campo de pressao //
		{
			// FD scheme:   //
			// forward time //
			t = ((float)it)*((float)dt);//((float)(it-1)) * dt;

			//printf("it[%d] e t:%lf\n",it,t); // t e it OK

			//printf("ix: [%d]-[%d]\n",(ixx0+DRVLEN-1),(ixx1-DRVLEN-1)); // OK
			//printf("iz: [%d]-[%d]\n\n",(DRVLEN-1),(nzz-DRVLEN-1));     // OK

			for(ix=(ixx0+DRVLEN);ix<(ixx1-DRVLEN);ix++) // nova sementinha do mal ... ixx0
			{
				for(iz=(DRVLEN-1);iz<(nzz-DRVLEN);iz++)
				{
					source = 0.0;

					for(itrc=0;itrc<nrec;itrc++)
					{
						if( (ix == (igx[itrc]-1)) && (iz == (igz[itrc]-1)) )
						{
							source = (dx*dx) * interp_trace(ttotal-t,ns,0.0,dtrec,shotgather[itrc].tr_data);
							//printf("ix: %d\tiz: %d\titrc: %d\n",ix,iz,itrc);
						}
					}

					gama = gama_x[ix] + gama_z[iz]; 
					invpgama  = ( ((double)1.0) / ( ((double)1.0) + gama) );
					mgama     = ((double)1.0)-gama;

					//if(source!=0.0) printf("[%d]source[%d][%d]: %lf   |ttotal: %lf   |t:%lf\n",it,ix,iz,source,ttotal,t); //OK
					//if(gama!=0.0) printf("[%d]gama[%d][%d]: %lf\n",it,ix,iz,gama); //OK
					//if(invpgama!=0.0) printf("[%d]invpgama[%d][%d]: %lf\n",it,ix,iz,invpgama); //OK
					//if(mgama!=0.0) printf("[%d]mgama[%d][%d]: %lf\n",it,ix,iz,mgama); //OK

					// aumentar a ordem do operador de diferencas finitas: //
      					
					laplacian = ((double)2.0) * deriv2[0] * p[iz][ix][1];
////					printf("deriv2[0]: %lf\n",deriv2[0]);
					//if( (laplacian!=0.0) ) printf("it-[%d]\tp[%d][%d]: %lf\tinit-laplacian: %lf\n",it,iz,ix,p[iz][ix][1],laplacian);	

					for(iconv=1;iconv<DRVLEN;iconv++)
					{
						laplacian = laplacian + deriv2[iconv] * \
									( p[iz-iconv][ix][1] + p[iz+iconv][ix][1] + \
									  p[iz][ix-iconv][1] + p[iz][ix+iconv][1] ) ;
//						printf("p[%d][%d]\n",(iz-iconv),ix);
//						printf("p[%d][%d]\n",(iz+iconv),ix);
//						printf("p[%d][%d]\n",iz,ix-iconv);
//						printf("p[%d][%d]\n",iz,ix+iconv);
//						printf("deriv2[%d]\n",iconv);
//						printf("lap: %lf\n\n",laplacian);
					} 

					p[iz][ix][0] = invpgama * ( ((double)2.0) * p[iz][ix][1] - \
							mgama * p[iz][ix][0] + vel[iz][ix] * ( laplacian - source ) );

			//printf("vel[%d][%d]: %lf\n",iz,ix,vel[iz][ix]);

			//if(p[iz][ix][1]!=((double)0.0)) printf("it: %d\tix: %d\tiz: %d\tp[iz][ix][1]: %lf\n",it,ix,iz,p[iz][ix][1]);

//			if( laplacian!=0.0 ) printf("it: %d\tiz: %d\tix: %d\tlaplacian: %lf\n",it,iz,ix,laplacian);
////					printf("******************************************\n");

				}
			}

			//printf("it[%d] iz[%d] ix[%d]\n",it,iz,ix); //

			// swap fields //
			//printf("ix[%d][%d]\n",ixx0-1,ixx1-1);     //OK
			//printf("iz[%d][%d]\n",0,nzz-1);      	    //OK				
			
			for(ix=ixx0; ix<ixx1; ix++)
			{
				for(iz=0; iz<nzz; iz++)
				{
					dswap = p[iz][ix][1];
					p[iz][ix][1] = p[iz][ix][0];
					p[iz][ix][0] = dswap;
				}
			}

			if( ( it % ndtrtm) == 0) // if( ((it-1) % ndtrtm) == 0)
			{
				//printf("it: %d\tndtrtm: %d\n",it,ndtrtm);
				if( (iframe%50==0) ) 
					printf("ishot[%d]\tit: %d\tbackward frames %d completed\n",ishot, it, iframe);
				
				fseek(fdbackward_file, iframe * reclen , SEEK_SET);
    				iframe++;
////				printf("izb[%d]\tize[%d]\n",izb,ize);
////				printf("ixmig0[%d]\tixmig1[%d]\n",ixmig0,ixmig1);

				for(ix=ixmig0-1;ix<ixmig1;ix++)
				{
					for(iz=izb-1;iz<ize;iz++)
					{
						//printf("lendo p[%d][%d]: %lf\n",iz,ix,((float)p[iz][ix][1]));
						//pt_aux = ((float)p[iz][ix][1]);
						//fwrite(&pt_aux, sizeof(float), nz, fdforward_file);
						//printf("writing p[%d][%d]: %f\n",iz,ix,pt_aux);
						flt_vet[iz-nborda] = ((float)p[iz][ix][1]);
					}
					fwrite(flt_vet, sizeof(float), nz, fdbackward_file);
					//
					//for(iz=0;iz<nz;iz++)
					//{
					//	if(flt_vet[iz]!=0.0) printf("[%d][%d]: %.15lf\n",iz,ix,flt_vet[iz]);
					//	//if(p[iz][ix][1]!=0.0) printf("[%d][%d]: %.15lf\n",iz,ix,p[iz][ix][1]);
					//}
					//
				}

				//for(ix=ixmig0-1;ix<ixmig1;ix++) // escriva velha
				//{
				//	for(iz=izb-1;iz<ize;iz++)
				//	{
						//printf("lendo p[%d][%d]: %lf\n",iz,ix,((float)p[iz][ix][1]));
				//		pt_aux = ((float)p[iz][ix][1]);
////				//		printf("writing p[%d][%d]: %f\n",iz,ix,pt_aux);
				//		fwrite( &pt_aux, sizeof(float), 1, fdbackward_file);
				//	}
				//}
			}

		}// fim de modelagem

		fclose(fdbackward_file);

		nframes = iframe;
		
		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		/////	       RETROPROPAGACAO DO CAMPO DOS RECEPTORES: FIM		/////
		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		printf("(FIM RETROPROPAGACAO nframes = %d)\n", nframes);

		reclen = nz*nmigtrc*sizeof(float);

		if( (fdbackward_file = fopen(fdbackward_name, "rb")) == NULL )
   		{
        		printf("Error opening model file\n");
	        	return -1;
   		}
		
		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		/////	       	PROPAGACAO DO CAMPO DA FONTE: INICIO			/////
		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////

		printf("\nSOURCE PROPAGATION\n");

		if( (fdforward_file = fopen(fdforward_name, "wb")) == NULL )
    		{
        		printf("Error opening forward source file\n");
	        	return -1;
    		}


		if( (out_file = fopen(out_name, "wb")) == NULL )
    		{
        		printf("Error opening rtm_migrated file\n");
        		return -1;
    		}


		// condiçao inicial: repouso //
		// p(:,:,1) pressure field at time t-dt //
		// p(:,:,2) pressure field at time t //

		// p(:,:,:)= 0.0_dp // nzz,nxx,2
		for(it=0;it<2;it++)
			for(ix=0;ix<nxx;ix++)
				for(iz=0;iz<nzz;iz++)
					p[iz][ix][it] = ((double)0.0);
		
		//free3double(p,nxx,nzz,2);
		//p = alloc3double(nzz,nxx,2);    // Campos de pressão: p[851][2401][2]
		//if(p==NULL) { printf("Reallocation of p[%d,%d,%d] failed!!!\n",nzz,nxx,2); return -1;}

		irec   = 0;
		//itrcmv = 0;
		iframe = 0;

		for(it=0;it<nt;it++) // modelando a evolucao do campo de pressao //
		{
			delt  = pow( beta * ((double)(1+it-it0)), 2);
			//printf("\nit: %d\t",it);
			//printf("it0: %d\t",it0);
			//printf("delt: %lf\n",delt);

		        fonte = (dx*dx) * exp( -delt ) * ( ((double)1.0) - ((double)2.0)*delt );   // Pulso fonte Ricker frequancia pico freq //

			//printf("ix: [%d]-[%d]\n",(ixx0+DRVLEN-1),(ixx1-DRVLEN-1)); // OK
			//printf("iz: [%d]-[%d]\n",(DRVLEN-1),(nzz-DRVLEN-1));     // OK

			// testar fonte //
			for(ix=(ixx0+DRVLEN);ix<(ixx1-DRVLEN);ix++)
			{
				for(iz=(DRVLEN-1);iz<(nzz-DRVLEN);iz++)
				{
					source = 0.0;

					if ( (ix == (isx-1)) && (iz == (isz-1)) )
           				{
						source = fonte;
						//printf("ix: %d\tisx: %d\t",ix,isx-1);
						//printf("iz: %d\tisz: %d\n",iz,isz-1);
					}

					//if(source!=0.0) printf("[%d]source[%d][%d]: %lf\n",it,ix,iz,source); //OK
					
					gama = gama_x[ix] + gama_z[iz]; 
					invpgama  = ( ((double)1.0) / ( ((double)1.0) + gama) );
					mgama     = ((double)1.0)-gama;

					laplacian = ((double)2.0) * deriv2[0] * p[iz][ix][1];

					for(iconv=1;iconv<DRVLEN;iconv++)
					{
						laplacian = laplacian + deriv2[iconv] * \
										( p[iz-iconv][ix][1] + p[iz+iconv][ix][1] + \
										  p[iz][ix-iconv][1] + p[iz][ix+iconv][1] );
					}
					
					p[iz][ix][0] = invpgama * ( ((double)2.0) * p[iz][ix][1] - \
								mgama * p[iz][ix][0] + vel[iz][ix] * ( laplacian - source ) );

					//if(p[iz][ix][0]!=((double)0.0)) printf("it: %d\tix: %d\tiz: %d\tp[iz][ix][1]: %.12lf\n",it,ix,iz,p[iz][ix][0]);
				}
			}
			
			// swap fields //
			for(ix=ixx0; ix<ixx1; ix++)
			{
				for(iz=0; iz<nzz; iz++)
				{
					dswap = p[iz][ix][1];
					p[iz][ix][1] = p[iz][ix][0];
					p[iz][ix][0] = dswap;
				}
			}

			if( ( it % ndtrtm) == 0) // if( ((it-1) % ndtrtm) == 0)
			{
				if( (iframe%50==0) ) printf("ishot[%d]\tit: %d\tforward frame %d completed\n",ishot, it, iframe);

				//read(MVID,rec=irec) prcv(1:nz,1:nmigtrc)
				irec = nframes-iframe-1;
				fseek(fdbackward_file, irec * reclen , SEEK_SET);
				for(ix=0;ix<nmigtrc;ix++)
    				{
        				fread(flt_vet, sizeof(float), nz, fdbackward_file);
        				for(iz=0;iz<nz;iz++)
            					prcv[iz][ix] = flt_vet[iz];
    				}
				//read(MVID,rec=irec) prcv(1:nz,1:nmigtrc)

				//psrc(1:nz,1:nmigtrc) = real(p(izb:ize,ixmig0:ixmig1,2),sp) //
				fseek(fdforward_file, iframe * reclen , SEEK_SET);
    				iframe++;
				for(ix=0;ix<nmigtrc;ix++)
				{
					for(iz=0;iz<nz;iz++)
					{
						psrc[iz][ix]=((float)p[izb-1+iz][ixmig0-1+ix][1]);
						flt_vet[iz] = psrc[iz][ix];
					}
					fwrite(flt_vet, sizeof(float), nz, fdforward_file);
				}
				//psrc(1:nz,1:nmigtrc) = real(p(izb:ize,ixmig0:ixmig1,2),sp) //

//
//				//write(MVFWID,rec=iframe) real(p(izb:ize,ixmig0:ixmig1,2),sp)
//				fseek(fdforward_file, iframe * reclen , SEEK_SET);
//  				iframe++;
//				for(ix=ixmig0-1;ix<ixmig1;ix++)
//				{
//					for(iz=izb-1;iz<ize;iz++)
//					{
//						//printf("lendo p[%d][%d]: %lf\n",iz,ix,((float)p[iz][ix][1]));
//						//pt_aux = ((float)p[iz][ix][1]);
//						//fwrite(&pt_aux, sizeof(float), nz, fdforward_file);
//						//printf("writing p[%d][%d]: %f\n",iz,ix,pt_aux);
//						flt_vet[iz-nborda] = ((float)p[iz][ix][1]);
//					}
//					fwrite(flt_vet, sizeof(float), nz, fdforward_file);
//				}
//				//write(MVFWID,rec=iframe) real(p(izb:ize,ixmig0:ixmig1,2),sp)
//	

				// CONDICAO DE IMAGEM CROSS-CORRELACAO //
    				itrc=0;
    				for(ix=(ixmig0-1-nborda);ix<(ixmig1-nborda);ix++)
				{
       					for(iz=0;iz<nz;iz++)
          					imag[iz][ix] = imag[iz][ix] + \
							( ((float)window_x[ix]) * psrc[iz][itrc] * prcv[iz][itrc] );
    				       	itrc++;
				}

			} // fim do if //
		} // fim de modelagem //
		printf("(FIM SOURCE MODELING nframes = %d)\n", nframes);

		fclose(fdbackward_file);
		fclose(fdforward_file);

		for(ix=(DRVLEN-1);ix<(nx-DRVLEN);ix++)
		{
			for(iz=(DRVLEN-1);iz<(nz-DRVLEN);iz++)
			{

				laplacian = ((double)2.0) * deriv2[0] * imag[iz][ix];

				for(iconv=1;iconv<DRVLEN;iconv++)
				{
					laplacian = laplacian + deriv2[iconv] * \
							( ((double)imag[iz-iconv][ix]) + ((double)imag[iz+iconv][ix]) + \
							  ((double)imag[iz][ix-iconv]) + ((double)imag[iz][ix+iconv]) );
				}

				imag_filter[iz][ix] = ((float)laplacian);
			}
		}
	

		//////////////////////////////////////////////////////////////////////////////////
        	//			WRITE DATA TO OUTPUT SU FILE				//
		//////////////////////////////////////////////////////////////////////////////////


		// Arquivo para armazenar imagem //
		for(ix=0;ix<nx;ix++)
		{
			imagtrace.tr_header.tracl = ix+1;
			//imagtrace.tr_header.tracl = ntrmig+1;
			//imagtrace.tr_header.tracf = ix+1;

			for(iz=0;iz<nz;iz++)
				imagtrace.tr_data[iz] = imag_filter[iz][ix];

			//put_tr(shotidx[ishot]+ix, &(imagtrace.tr_header) , imagtrace.tr_data, nz, out_file);

			//fseek(out_file, ntrmig * ( (sizeof(float)*nz) + HDRBYTES ), SEEK_SET);
			ntrmig++;

			fwrite(&(imagtrace.tr_header), HDRBYTES, 1, out_file);
			fwrite(imagtrace.tr_data, sizeof(float), nz, out_file);
	
		}
	
	} // END LOOP OF SHOTS

	printf("end of migration\n");

	fclose(out_file);

	fclose(su_file);

//
//	int tr_number=0;
//
//    	segy *tr_hdr;
//    	float *tr_data;
//    	float **shotgather;
//
//	tr_hdr      = (segy*) malloc(sizeof(segy));
//    	tr_data     = alloc1float(nt);
//    	shotgather  = alloc2float(nt, ng);
//
//    	do{
//    
//        	oldsx=sx;
//        	oldgx=gx;
//
//        	// Looping over traces to get data //
//        
//        	do{
//            		
//            		//printf("hdr.sx=%d\n",tr_hdr->sx);
//            		//printf("hdr.sy=%d\n",tr_hdr->sy);
//            		//printf("hdr.gx=%d\n",tr_hdr->gx);
//            		//printf("hdr.gy=%d\n\n",tr_hdr->gy);
//            		//for(i=720;i<725;i++)
//            		//	printf("[%i]: %f\n",i,tr_data[i]);
//            		
//            		copy2fdata(shotgather, (tr_number%ng), tr_data, nt);
//            
//			//        print2float(shot,nt,2);
//            
//            		get_sx_gx(tr_hdr, &sx,&gx);
//            
//            		sx = (sx - min_sx_gx);
//            		gx = (gx - min_sx_gx);
//            
//            		//        printf("oldsx: (%.2f) | sx: (%.2f) | gx: (%.2f) | min_sx_gx: (%.2f)\n",oldsx, sx,gx,min_sx_gx);
//            
//            		++tr_number;
//        	}while( (get_tr(tr_hdr,tr_data,nt,su_file)) && (sx==oldsx) );
//        
//        	//////////////////////////////////////////////////    
//        	/////////////////// MIGRATION ////////////////////
//        	//////////////////////////////////////////////////    
//
//        	printf("migrating shot [%d]\n",240-nxshot);
//        	--nxshot;
//    	}while(nxshot);
//

	//////////////////////////////////////////////////////////////////////////////////
        //				DEALOCATING VECTORS				//
	//////////////////////////////////////////////////////////////////////////////////

	free1int(shotidx);
	free1int(igz);
	free1int(igx);


	free1float(trace.tr_data);
	free1float(imagtrace.tr_data);

	free2float(psrc, nz, nx);
	free2float(prcv, nz, nx);
	free2float(imag, nz, nx);
	free2float(imag_filter, nz, nx);


	free1double(gama_x);
	free1double(gama_z);
	free1double(window_x);

    	free2double(vel, nzz, nxx );
//	free2double(sourcemask, nzz, nxx);

	free3double(p,nzz,nxx,2);


	free1su_trace(shotgather);

	printf("\n");

    	return 0;
}

