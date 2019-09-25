
#include "decs.h"



/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

Generic Alfven test of Giacomazzo and Rezzolla (2006)

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
	gam = (5./3.) ;
	cour = 0.8 ;

	t = 0. ;                     /* Initial time */ 
	GridLength[0] = 1.5 ;        /* Length of X0 dimension (evolution period) */ 
	GridLength[1] = 4. ;        /* Length of X1 dimension */ 
	GridLength[2] = 4. ;        /* Length of X2 dimension */ 
	GridLength[3] = 4. ;        /* Length of X3 dimension */ 

	dx[0] = 1.e-5;               /* Discretization size in X0 direction */
	dx[1] = GridLength[1]/N1 ;   /* Discretization size in X1 direction */
	dx[2] = GridLength[2]/N2 ;   /* Discretization size in X2 direction */
	dx[3] = GridLength[3]/N3 ;   /* Discretization size in X3 direction */

	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
	startx[1] = -2.;               /* Set the Physical Minimum boundary  */
	startx[2] = -2.;               /* Set the Physical Minimum boundary  */
	startx[3] = -2.;               /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	DT_out[OUT_ASCII]   = GridLength[0]/10. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/40. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = GridLength[0]/200. ;      /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	n_restart = 100;                /* number of time steps between restart dumps */
	

	calc_all_geom() ;  /* Calculate the static metric grid functions */

}


/***********************************************************************************/
/***********************************************************************************
  init_data(): 
  ---------
   -- calculates the initial distribution of the MHD fields;
   -- driver routine for various initial data prescriptions;

***********************************************************************************/
void init_data(void) 
{
  int i, j, k ; 
  double X[NDIM], r, vx1, vy1, vz1, vx2, vy2, vz2, gamma1, gamma2, v1, v2;
  double fr(double R,double min,double max);

  vx1 = 0.0;
  vy1 = 0.3;
  vz1 = 0.4;
  v1 = sqrt( vx1*vx1 + vy1*vy1 + vz1*vz1 ) ; 
  gamma1 = 1./sqrt((1.-v1)*(1.+v1));

  vx2 = 0.;
  vy2 = 0.;
  vz2 = 0.;
  v2 = sqrt( vx2*vx2 + vy2*vy2 + vz2*vz2 ) ; 
  gamma2 = 1./sqrt((1.-v2)*(1.+v2));


  /* Komissarov's explosion problem */
  LOOP {
    coord(i,j,k,CENT,X) ; 

    //r = sqrt(X[1]*X[1] + X[2]*X[2]) ;

    if( X[1] < 0. ) { 
      p[i][j][k][U1] = vx1*gamma1 ;
      p[i][j][k][U2] = vy1*gamma1 ;
      p[i][j][k][U3] = vz1*gamma1 ;
      p[i][j][k][B1] = 1. ;
      p[i][j][k][B2] = 6. ;
      p[i][j][k][B3] = 2. ;
      p[i][j][k][RHO] = 1.;
      p[i][j][k][UU] = 5./(gam - 1.) ;
    }
    else { 
      p[i][j][k][U1] = vx2*gamma2 ;
      p[i][j][k][U2] = vy2*gamma2 ;
      p[i][j][k][U3] = vz2*gamma2 ;
      p[i][j][k][B1] = 1. ;
      p[i][j][k][B2] = 5.;
      p[i][j][k][B3] = 2. ;
      p[i][j][k][RHO] = 0.9;
      p[i][j][k][UU] = 5.3/(gam - 1.);
    }
  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

}

double fr(double R,double min_r,double max_r)
{
        double r_outer,r_inner,dR ;

        r_outer = GridLength[1]/12. ;
        r_inner = 0.8*r_outer ;
        dR = (r_outer - r_inner) ;

        if(R > r_outer) return(min_r) ;
        if(R < r_inner) return(max_r) ;
        else {
                return( exp(
                        log(max_r)*(r_outer - R)/dR +
                        log(min_r)*(R - r_inner)/dR)) ;
        }
}