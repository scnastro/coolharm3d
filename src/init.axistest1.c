
#include "decs.h"



/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

 Unique test to check the boundary conditions along the axis in spherical coordinates.
  This version of the test is the "true solution" and is done in cartesian (x,y)  
  coordinates since the test should be "easy" in them.  The initial data is an
  azimuthally-symmetric Gaussian density and pressure profile initially at rest.  
  The pressure makes it expand outward and should reflect through the axis 
  symmetrically. 

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
    fail(FAIL_BASIC,0,0);
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
	cour = 0.8 ;

	t = 0. ;                     /* Initial time */ 
	GridLength[0] = 20. ;         /* Length of X0 dimension (evolution period) */ 
	GridLength[1] = 10. ;        /* Length of X1 dimension */ 
	GridLength[2] = 10. ;        /* Length of X2 dimension */ 
	GridLength[3] = 10. ;        /* Length of X3 dimension */ 

	dx[0] = 2.e-2;               /* Discretization size in X0 direction */
	dx[1] = GridLength[1]/N1 ;   /* Discretization size in X1 direction */
	dx[2] = GridLength[2]/N2 ;   /* Discretization size in X2 direction */
	dx[3] = GridLength[3]/N3 ;   /* Discretization size in X3 direction */

	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
	startx[1] = -0.5*GridLength[1];               /* Set the Physical Minimum boundary  */
	startx[2] = -0.5*GridLength[2];               /* Set the Physical Minimum boundary  */
	startx[3] = -0.5*GridLength[3];               /* Set the Physical Minimum boundary  */
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
  double X[NDIM], r;
  double fr(double R,double min,double max);


  /* Loop over the grid  */
  LOOP {
    coord(i,j,k,CENT,X) ; 
    r = sqrt(X[1]*X[1] + X[2]*X[2]) ;

    p[i][j][k][U1] = 0. ;
    p[i][j][k][U2] = 0. ;
    p[i][j][k][U3] = 0. ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
    p[i][j][k][RHO] = fr(r,RHOMIN,(RHOMIN*1.e4)) ;
    p[i][j][k][UU] = fr(r,UUMIN,(RHOMIN)*1.e5)/(gam - 1.) ;
  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

}

double fr(double R,double min_f,double max_f)
{
        double r_0, width,dR ;

	width = GridLength[1]/10. ;
        r_0 = GridLength[1]/4. ;
        dR = (r_0 - R)/width ;
	return( min_f + max_f * exp(-dR*dR) ); 
}
