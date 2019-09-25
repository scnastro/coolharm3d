
#include "decs.h"



/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

 -- Linear wave test by Gammie etal 2003. 

***********************************************************************************/

void init(void)
{
  double ftmp;
  void init_data(void) ;

  /* Do some checks on the pre-compiler definitions : */ 
  ftmp = (double) N1TOT;  ftmp *= N2TOT  ;  ftmp *= N3TOT  ;  ftmp *= NP; 
  if( ftmp > ( (double) UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Max. no. of elements too large !!! N = %g \n", ftmp);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

      
  init_base();  /* Set global grid parameters and static work arrays */

  init_data();  /* Set MHD grid functions */

  return;
}

/***********************************************************************************/
/***********************************************************************************
  init_base(): 
  ---------
   -- initializes grid structure and constant grid functions; 

***********************************************************************************/
void init_base(void)
{
	int i,j,k ;
	double X[NDIM];

	void calc_all_geom(void) ;

	/* Set the parameters that define the grid and other constants : */ 
	gam = (4./3.) ;
	cour = 0.4 ;

	t = 0. ;                     /* Initial time */ 
	GridLength[0] = 0.9 ;        /* Length of X0 dimension (evolution period) */ 
	GridLength[1] = 1. ;        /* Length of X1 dimension */ 
	GridLength[2] = 1. ;        /* Length of X2 dimension */ 
	GridLength[3] = 1. ;        /* Length of X3 dimension */ 

	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	dx[0] = 1.e-5;               /* Discretization size in X0 direction */
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];

	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
	startx[1] =  0.;               /* Set the Physical Minimum boundary  */
	startx[2] =  0.;               /* Set the Physical Minimum boundary  */
	startx[3] =  0.;               /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	DT_out[OUT_ASCII]   = GridLength[0]/10. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/40. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = 0.05 ;	                /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	n_restart = 100;                /* number of time steps between restart dumps */
	

	calc_all_geom() ;  /* Calculate the static metric grid functions */

}

#define ENTROPY_WAVE		1
#define SLOW		2
#define ALFVEN		3
#define FAST		4


/***********************************************************************************/
/***********************************************************************************
  init_data(): 
  ---------
   -- calculates the initial distribution of the MHD fields;
   -- driver routine for various initial data prescriptions;
   -- this routine taken from Linear wave test of harm2d

***********************************************************************************/
void init_data(void) 
{
  int i, j, k,l ; 
  double X[NDIM], r;

  int branch ;
  double rmax,rmin,pmax,pmin ;
  double dp[NP] ;
  double eigs(double Bx0, int branch, double *dp) ;
  double AMP, tf ;
  double rho0,uu0,Bx0,kx,ky,w,delta,eph ;

  /* linear theory quantities */
  rho0 = 1. ;
  uu0 = 1./(gam - 1.) ;
  delta = 1. ;
  Bx0 = sqrt(delta) ;
  kx = 2.*M_PI ;
  ky = 2.*M_PI ;
  AMP = 1.e-4 ;
  branch = SLOW ;

  w = eigs(Bx0,branch,dp) ;

  tf = 2.*M_PI/w ;
  GridLength[0] = tf;
  DT_out[OUT_HDF5] = tf/10.;

  LOOP {
    coord(i,j,k,CENT,X) ; 
    eph = cos(kx*X[1] + ky*X[2]) ;

    p[i][j][k][RHO] = rho0 + dp[0]*AMP*eph ;
    p[i][j][k][UU]  = uu0  + dp[1]*AMP*eph ;
    p[i][j][k][U1]  = 0.   + dp[2]*AMP*eph ;
    p[i][j][k][U2]  = 0.   + dp[3]*AMP*eph ;
    p[i][j][k][U3]  = 0.   + dp[4]*AMP*eph ;
    p[i][j][k][B1]  = Bx0  + dp[5]*AMP*eph ;
    p[i][j][k][B2]  = 0.   + dp[6]*AMP*eph ;
    p[i][j][k][B3]  = 0.   + dp[7]*AMP*eph ;

  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

}


/* simplified eigensolvers */
double eigs(double Bx, int branch, double *dp) 
{
	int k,l ;
	double Bx2,aux1,r1,r2 ;

	PLOOP dp[l] = 0. ;

	if(branch == ENTROPY_WAVE) {
		dp[UU] = 1. ;
		return(1.) ;
	}
	else if(branch == SLOW) {
		Bx2 = Bx*Bx ;
		aux1 = 400. + 80.*Bx2 + 169.*Bx2*Bx2 ;
		r1 = 2.*M_PI*sqrt((20. + 17.*Bx2 - sqrt(aux1))/
			(75. + 15.*Bx2)) ;
		r2 = 2.*M_PI*sqrt((20. + 17.*Bx2 + sqrt(aux1))/
			(75. + 15.*Bx2)) ;
		if(r2 > r1) {
			dp[0] = (-20. + 13*Bx2 + sqrt(aux1))/20. ;
			dp[1] = (-20. + 13*Bx2 + sqrt(aux1))/5. ;
			dp[2] = sqrt((20. + 17.*Bx2 - sqrt(aux1))/
				(75.+15.*Bx2))*(13.*Bx2 +
				sqrt(aux1))/20. ;
			dp[3] = -sqrt((20. + 17.*Bx2 - sqrt(aux1))/
				(75.+15.*Bx2)) ;
			dp[5] = -Bx ;
			dp[6] = Bx ;
			return(r1) ;
		}
		else {
			dp[0] = (-20. + 13*Bx2 - sqrt(aux1))/20. ;
			dp[1] = (-20. + 13*Bx2 - sqrt(aux1))/5. ;
			dp[2] = sqrt((20. + 17.*Bx2 + sqrt(aux1))/
				(75.+15.*Bx2))*(13.*Bx2 -
				sqrt(aux1))/20. ;
			dp[3] = -sqrt((20. + 17.*Bx2 + sqrt(aux1))/
				(75.+15.*Bx2)) ;
			dp[5] = -Bx ;
			dp[6] = Bx ;
			return(r2) ;
		}
	}
	else if(branch == ALFVEN) {
		Bx2 = Bx*Bx ;
		dp[U3] = -1./sqrt(5. + Bx2) ;
		dp[B3] = 1. ;
		return(2.*M_PI*Bx/sqrt(5. + Bx2)) ;
	}
	else if(branch == FAST) {
		Bx2 = Bx*Bx ;
		aux1 = 400. + 80.*Bx2 + 169.*Bx2*Bx2 ;
		r1 = 2.*M_PI*sqrt((20. + 17.*Bx2 - sqrt(aux1))/
			(75. + 15.*Bx2)) ;
		r2 = 2.*M_PI*sqrt((20. + 17.*Bx2 + sqrt(aux1))/
			(75. + 15.*Bx2)) ;
		if(r1 > r2) {
			dp[0] = (-20. + 13*Bx2 + sqrt(aux1))/20. ;
			dp[1] = (-20. + 13*Bx2 + sqrt(aux1))/5. ;
			dp[2] = sqrt((20. + 17.*Bx2 - sqrt(aux1))/
				(75.+15.*Bx2))*(13.*Bx2 +
				sqrt(aux1))/20. ;
			dp[3] = -sqrt((20. + 17.*Bx2 - sqrt(aux1))/
				(75.+15.*Bx2)) ;
			dp[5] = -Bx ;
			dp[6] = Bx ;
			return(r1) ;
		}
		else {
			dp[0] = (-20. + 13*Bx2 - sqrt(aux1))/20. ;
			dp[1] = (-20. + 13*Bx2 - sqrt(aux1))/5. ;
			dp[2] = sqrt((20. + 17.*Bx2 + sqrt(aux1))/
				(75.+15.*Bx2))*(13.*Bx2 -
				sqrt(aux1))/20. ;
			dp[3] = -sqrt((20. + 17.*Bx2 + sqrt(aux1))/
				(75.+15.*Bx2)) ;
			dp[5] = -Bx ;
			dp[6] = Bx ;
			return(r2) ;
		}
	}
	else { 
	  fprintf(stderr,"Should not be here !! \n");
	  fflush(stderr);
	  return(-1.); 
	}
}

#undef ENTROPY_WAVE
#undef SLOW
#undef ALFVEN	
#undef FAST
