/* 
	produces an "r8" file.  
	at the moment, just gives a slice through 3D data set.
*/

#include "decs.h"

#define DENSITY 	0
#define SDENSITY 	1
#define ENERGY		0
#define VORTICITY 	0 
#define CURRENT 	0

#define LOOPSLICE N1_LOOP N2_LOOP
//#define LOOPSLICE N2_LOOP N3_LOOP 

#if( MAKE_IMAGE ) 
void image(void)
{
        double iq,liq,a,lmax,lmin,Pmax,Pmin ;
	double j1,j2,j3,w1,w2,w3 ;
	int i,j,k ;
	FILE *im_file ;
	char dfnam[50];



	k = N3/2 ;

	sprintf(dfnam,"%s/im%04d.%04d",DIR_out[OUT_IMAGE],N_out[OUT_IMAGE],myid) ;
	im_file = fopen(dfnam,"w") ;
	if(im_file==NULL) {
	  fprintf(stderr,"image(): error opening image file = \n", dfnam) ;
	  fail( FAIL_BASIC,0 );
	}
 

	/* density mapping N1S logarithmic, in 255 steps
	   between e^lmax and e^lmin */
#if ENERGY
	Pmax = -1.e9 ;
	Pmin =  1.e9 ;
	LOOPSLICE {
		iq = p[i][j][k][UU] ;
		if(iq > Pmax) Pmax = iq ;
		if(iq < Pmin) Pmin = iq ;
	}
	lmax = Pmax ;
	lmin = Pmin ;
#endif

#if SDENSITY
	Pmax = -1.e9 ;
	Pmin =  1.e9 ;
	LOOPSLICE {
		iq = 0. ;
		for(k=N3S;k<=N3E;k++) {
			iq += p[i][j][k][RHO]*dV ;
		}
		if(iq > Pmax) Pmax = iq ;
		if(iq < Pmin) Pmin = iq ;
	}
	lmax = Pmax ;
	lmin = Pmin ;
#endif

#if CURRENT
	Pmax = -1.e9 ;
	Pmin =  1.e9 ;
	LOOPSLICE {
		j1 = (p[i][j+1][k][B3] - p[i][j-1][k][B3])/(2.*dx[2]) -
		     (p[i][j][k+1][B2] - p[i][j][k-1][B2])/(2.*dx[3]) ;
		j2 = (p[i][j][k+1][B1] - p[i][j][k-1][B1])/(2.*dx[3]) -
		     (p[i+1][j][k][B3] - p[i-1][j][k][B3])/(2.*dx[1]) ;
		j3 = (p[i+1][j][k][B2] - p[i-1][j][k][B2])/(2.*dx[1]) -
		     (p[i][j+1][k][B1] - p[i][j-1][k][B1])/(2.*dx[2]) ;
		iq = sqrt(j1*j1 + j2*j2 + j3*j3) ;
		if(iq > Pmax) Pmax = iq ;
		if(iq < Pmin) Pmin = iq ;
	}
	lmax = Pmax ;
	lmin = Pmin ;
#endif

#if VORTICITY
	Pmax = -1.e9 ;
	Pmin =  1.e9 ;
	LOOPSLICE {
		w1 = (p[i][j+1][k][U3] - p[i][j-1][k][U3])/(2.*dx[2]) -
		     (p[i][j][k+1][U2] - p[i][j][k-1][U2])/(2.*dx[3]) ;
		w2 = (p[i][j][k+1][U1] - p[i][j][k-1][U1])/(2.*dx[3]) -
		     (p[i+1][j][k][U3] - p[i-1][j][k][U3])/(2.*dx[1]) ;
		w3 = (p[i+1][j][k][U2] - p[i-1][j][k][U2])/(2.*dx[1]) -
		     (p[i][j+1][k][U1] - p[i][j-1][k][U1])/(2.*dx[2]) ;
		iq = sqrt(w1*w1 + w2*w2 + w3*w3) ;
		if(iq > Pmax) Pmax = iq ;
		if(iq < Pmin) Pmin = iq ;
	}
	lmax = Pmax ;
	lmin = Pmin ;
#endif

#if DENSITY
	Pmax = -1.e9 ;
	Pmin =  1.e9 ;
	LOOPSLICE {
		iq = p[i][j][k][RHO] ;
		if(iq > Pmax) Pmax = iq ;
		if(iq < Pmin) Pmin = iq ;
	}
	lmax = Pmax ;
	lmin = Pmin ;
#endif

	a = 255./(lmax - lmin) ;

	/* thN1S must be changed if loopslice N1S changed */
	//LOOPSLICE {
	for(j=N2E;j>=N2S;j--)
	for(i=N1S;i<=N1E;i++) {
#if DENSITY
		iq = p[i][j][k][RHO] ;
#endif
#if SDENSITY
		iq = 0. ;
		for(k=N3S;k<=N3E;k++)
			iq += p[i][j][k][RHO]*dV ;
#endif
#if ENERGY
		iq = p[i][j][k][UU] ;
#endif
#if VORTICITY
		w1 = (p[i][j+1][k][U3] - p[i][j-1][k][U3])/(2.*dx[2]) -
		     (p[i][j][k+1][U2] - p[i][j][k-1][U2])/(2.*dx[3]) ;
		w2 = (p[i][j][k+1][U1] - p[i][j][k-1][U1])/(2.*dx[3]) -
		     (p[i+1][j][k][U3] - p[i-1][j][k][U3])/(2.*dx[1]) ;
		w3 = (p[i+1][j][k][U2] - p[i-1][j][k][U2])/(2.*dx[1]) -
		     (p[i][j+1][k][U1] - p[i][j-1][k][U1])/(2.*dx[2]) ;
		iq = sqrt(w1*w1 + w2*w2 + w3*w3) ;
#endif
#if CURRENT
		j1 = (p[i][j+1][k][B3] - p[i][j-1][k][B3])/(2.*dx[2]) -
		     (p[i][j][k+1][B2] - p[i][j][k-1][B2])/(2.*dx[3]) ;
		j2 = (p[i][j][k+1][B1] - p[i][j][k-1][B1])/(2.*dx[3]) -
		     (p[i+1][j][k][B3] - p[i-1][j][k][B3])/(2.*dx[1]) ;
		j3 = (p[i+1][j][k][B2] - p[i-1][j][k][B2])/(2.*dx[1]) -
		     (p[i][j+1][k][B1] - p[i][j-1][k][B1])/(2.*dx[2]) ;
		iq = sqrt(j1*j1 + j2*j2 + j3*j3) ;
#endif

		liq = a*(iq - lmin) ;

		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}


	fclose(im_file) ;
	return;
}
#else
void image(void) { return; } 
#endif
