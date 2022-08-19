#include <stdio.h>
#include <stdlib.h>

#include "segy.h"
#include "string.h"


int get_set_hdr_data(int tr_index, int ns, FILE *su_file, FILE *su_file_out)
{
	fseek(su_file, tr_index * ( (sizeof(float) * ns) + HDRBYTES ), SEEK_SET);

    	int nread;    
    	su_header *tr_hdr;    
    	float *tr_data;
    	void *ptr;

    	tr_hdr = (su_header*) malloc(sizeof(su_header));
    	tr_data = (float*) calloc(sizeof(float), ns);
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
        		ptr = (float*) calloc(sizeof(float), ns);

        		fread((float*) ptr, sizeof(float) , ns, su_file);

       		    	memcpy(tr_data, (float*)ptr, sizeof(float) * ns);

                	free(ptr);

			fseek(su_file_out, tr_index * ( (sizeof(float) * ns) + HDRBYTES ), SEEK_SET);

                	// DADOS A EDITAR NO HEADER //
                	tr_hdr->gelev = (-25);
                	tr_hdr->sdepth = 10;
                	// DADOS A EDITAR NO HEADER //

			fwrite(tr_hdr, HDRBYTES, 1, su_file_out);
            		fwrite(tr_data, sizeof(float), ns, su_file_out);

	                return nread;
        	        break;
    	}
}


int main(int argc, char *argv[])
{
	FILE *su_file_in;		// complete su data - marmousi - 240 tiros - 23040 traços
	char *f_name1 = "MARMOUSI_xtLINUX_CORRECT.su";

	if( (su_file_in = fopen(f_name1, "rb")) == NULL )
   	{
       		printf("Error opening seismic data file\n");
        	return -1;
   	}

	FILE *su_file_out;		// complete su data - marmousi - 240 tiros - 23040 traços
	char *f_name2 = "MARMOUSI_xtLINUX_CORRECT2.su";

	if( (su_file_out = fopen(f_name2, "wb")) == NULL )
   	{
       		printf("Error opening seismic data file\n");
        	return -1;
   	}

    	int ntrc=0;
    	int ns=725;

	get_set_hdr_data(ntrc, ns, su_file_in, su_file_out);
    
    	do{
		ntrc++;
	}while( get_set_hdr_data(ntrc, ns, su_file_in, su_file_out) != 0 );

	fclose(su_file_in);
	fclose(su_file_out);

	printf("\n");
	//system("PAUSE");
   	return 0;
}
