
#include "decs.h"


/* Determine whether we specify the z-component of vector potential, else we specify B-field (and generate monopoles) */
/* The B-field option is to follow the test done in the Einstein Toolkit's MHD paper.  */
#define SET_AZ (1) 

/***********************************************************************************/
/* Parameter of the run : */ 
/***********************************************************************************/
static const double v0          = (1./2.) ;     /* Magnitude of the velocity component of the field loop: */
static const double loop_radius = 0.3  ;     /* Radius of loop */
static const double Az0         = 1.e-3;     /* Magnitude of the vector potential */
static double theta; /* Angle between the x & y velocity components */



/* static void set_mag_field( void ) ;  */
static void set_mag_field( double vcon[NDIM] ) ; 


/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

  Advected field loop test.  Assumes Cartesian coordinates and Minkowski. 

  See Gardiner & Stone (2005) and references therein for history. 

  Should setup it up so that N2 = N1/2

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

  if( n_cells_glob > ( UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Num. of global points too large !!! N = %ld \n", 
	    n_cells_glob);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

#if( TOP_TYPE_CHOICE != TOP_CARTESIAN ) 
  this-initial-data-routine-not-setup-for-noncartesian-grids
#endif 
      
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
	cour = 0.4 ;

	t = 0. ;                     /* Initial time */ 

	SET_GEOM_COORD_TIME_FUNCS(t,1);

	if( v0 < SMALL ) { 
	  GridLength[0] = 10. ;         /* Length of X0 dimension (evolution period) */ 
	}
	else { 
	  GridLength[0] = 25.;         /* Length of X0 dimension (evolution period) */ 
	} 
	dx[0] = 1.e-5;               /* Discretization size in X0 direction */

#if( COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN )
	GridLength[1] = 1. ;        /* Length of X1 dimension */ 
	GridLength[2] = 1. ;        /* Length of X2 dimension */ 
	startx[1] = 0.;               /* Set the Physical Minimum boundary  */
	startx[2] = 0.;              /* Set the Physical Minimum boundary  */
#else
	GridLength[1] = 6. ;        /* Length of X1 dimension */ 
	GridLength[2] = 6. ;        /* Length of X2 dimension */ 
	startx[1] = -3.;               /* Set the Physical Minimum boundary  */
	startx[2] = -3;              /* Set the Physical Minimum boundary  */
#endif

        dx[1] = GridLength[1] / totalsize[1];
        dx[2] = GridLength[2] / totalsize[2];
        dx[3] = dx[2] ;

	GridLength[3] = N3*dx[3] ;        /* Length of X3 dimension */ 

	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
	startx[3] = -0.5*dx[3];        /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	DT_out[OUT_ASCII]   = GridLength[0]/10. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/40. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = 0.5 ;      /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	n_restart = 1000;                /* number of time steps between restart dumps */
	

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
  double *ploc;
  double vcon[NDIM], vconp[NDIM], bcon[NDIM];

  double x,y,r;
  struct of_geom *geom;  
  struct of_coord *coords;  

  double press = 3.;

  double mink[NDIM][NDIM] = {{-1.,0.,0.,0.},
			     { 0.,1.,0.,0.},
			     { 0.,0.,1.,0.},
			     { 0.,0.,0.,1.} };

  /* flatspace components  of the 4-vel. :  */
  theta = atan2(1.,2.);
  double vx = v0;
  double vy = 0.5*vx;
  double vz = 0.;
  double gamma = 1./sqrt(1. - vx*vx - vy*vy - vz*vz); 

  vcon[TT] = gamma; 
  vcon[XX] = gamma*vx;
  vcon[YY] = gamma*vy;
  vcon[ZZ] = gamma*vz;

  /* Loop over the grid  */
  LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);
    get_geometry(i,j,k,CENT,ncurr,geom); 
    x = coords->x[1];
    y = coords->x[2];
    r = sqrt(x*x + y*y);
    ploc = p[i][j][k];

    ploc[RHO] = 1.;
    ploc[UU]  = press/(gam - 1.) ;

#if( SET_AZ ) 
      ploc[B1] = 0.;
      ploc[B2] = 0.;
      ploc[B3] = 0.;
#else
    if ( vx > SMALL || vy > SMALL || vz > SMALL) {
      fprintf(stderr,"\nWARNING! \n if the loop is moving, you should set SET_AZ to 1 ! \n\n"); 
      fflush(stderr);
    }
    if( r < loop_radius ) { 
      bcon[0] = 0.;
      bcon[1] = -Az0*y/r ;
      bcon[2] =  Az0*x/r ;
      bcon[3] = 0.;
      transform_rank1con2(coords->dxp_dx,bcon);
      ploc[B1] = bcon[1] ;
      ploc[B2] = bcon[2] ;
      ploc[B3] = bcon[3] ;
    } else {
      ploc[B1] = 0.;
      ploc[B2] = 0.;
      ploc[B3] = 0.;
    }
#endif

    vconp[0] = vcon[0];
    vconp[1] = vcon[1];
    vconp[2] = vcon[2];
    vconp[3] = vcon[3];

    transform_rank1con2(coords->dxp_dx,vconp);
    ucon2pr(ploc,vconp,geom->gcon);

  }

#if( SET_AZ ) 
  set_mag_field(vcon);
#endif

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

}

/***********************************************************************************/
/***********************************************************************************
  set_mag_field():
  ------------------
   -- Sets the magnetic field; 
   -- Currently only supports poloidal fields derived from the azimuthal component 
      of the vector potential; 
   -- Before the calculation, does a fixup() and sets the boundaries since 
      the derivatives of A_\phi require ghost cells (field follows density contours); 
   -- We assume below that the magnetic field is independent of X3 ; 
***********************************************************************************/
static void set_mag_field( double vcon[NDIM] )
{
  int i,j,k,l;

  double beta_act, norm, bsq_int, u_int ; 
  double x,y,r ; 
  double *ploc;
  struct of_geom *geom;  
  struct of_coord *coords;  
  extern void curl_of_A( double ****A , double ****prim ) ;

  double gamma = vcon[TT];
  double vx    = vcon[XX]/gamma ;
  double vy    = vcon[YY]/gamma ;
  double vz    = vcon[ZZ]/gamma ;
  double v2    = vx*vx + vy*vy + vz*vz ;

  double xx[NDIM] ;

  // coord transformation to the boosted frame
  double Lambda[NDIM][NDIM] ;

  Lambda[0 ][0 ]                  =  gamma ;
  Lambda[0 ][XX] = Lambda[XX][0 ] = -gamma*vx ;
  Lambda[0 ][YY] = Lambda[YY][0 ] = -gamma*vy ;
  Lambda[0 ][ZZ] = Lambda[ZZ][0 ] = -gamma*vz ;
  Lambda[XX][XX]                  = 1. + (gamma - 1.)*vx*vx/v2 ;
  Lambda[XX][YY] = Lambda[YY][XX] =      (gamma - 1.)*vx*vy/v2 ;
  Lambda[XX][ZZ] = Lambda[ZZ][XX] =      (gamma - 1.)*vx*vz/v2 ;
  Lambda[YY][YY]                  = 1. + (gamma - 1.)*vy*vy/v2 ;
  Lambda[YY][ZZ] = Lambda[YY][ZZ] =      (gamma - 1.)*vz*vy/v2 ;
  Lambda[ZZ][ZZ]                  = 1. + (gamma - 1.)*vz*vz/v2 ;


#if(0)
  // and its inverse (which we don't need)
  double invLambda[NDIM][NDIM] ;
  invLambda[0 ][0 ]                     = gamma ;
  invLambda[0 ][XX] = invLambda[XX][0 ] = gamma*vx ;
  invLambda[0 ][YY] = invLambda[YY][0 ] = gamma*vy ;
  invLambda[0 ][ZZ] = invLambda[ZZ][0 ] = gamma*vz ;
  invLambda[XX][XX]                     = 1. + (gamma - 1.)*vx*vx/v2 ;
  invLambda[XX][YY] = invLambda[YY][XX] =      (gamma - 1.)*vx*vy/v2 ;
  invLambda[XX][ZZ] = invLambda[ZZ][XX] =      (gamma - 1.)*vx*vz/v2 ;
  invLambda[YY][YY]                     = 1. + (gamma - 1.)*vy*vy/v2 ;
  invLambda[YY][ZZ] = invLambda[YY][ZZ] =      (gamma - 1.)*vz*vy/v2 ;
  invLambda[ZZ][ZZ]                     = 1. + (gamma - 1.)*vz*vz/v2 ;

  double  delta[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}} ;
  for (i=0; i<4 ; i++) {
    for (j=0; j<4 ; j++) {
      for (k=0; k<4 ; k++) {
        delta[i][j] += Lambda[i][k] * invLambda[k][j] ;
      }
      printf("delta^%d_%d = %.18lf \n",i,j, delta[i][j] );
    }
  }
#endif

  /*************************************************************************************
    Calculate the vector potential, A_mu: 
       --  Need to calc. A_\mu on the last ghost cells or all the corners surrounding the 
           physical cells' centers for finite diff. calculation for B^i 
       -- We first calculate A_\phi then we transform to xp;
       -- To this end we implicitly assume that the density distribution is 
           phi-independent (but not necessaryily xp3 independent);
  *************************************************************************************/
  double Acov[NDIM];

  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) for(k=N3S; k<=(N3E+1); k++) { 
	get_coord(i,j,k,CORN,ncurr,coords);

	xx[0] = coords->x[0];
	xx[1] = coords->x[1];
	xx[2] = coords->x[2];
	xx[3] = coords->x[3];

        transform_rank1con2(Lambda,xx);

	r = sqrt(xx[1]*xx[1] + xx[2]*xx[2]);

	Acov[0] = 0.;
	Acov[1] = 0.;
	Acov[2] = 0.;
	Acov[3] = 0.;

	if( r <= loop_radius ) { 
	  Acov[3] = Az0 * (loop_radius - r); 
	}
        transform_rank1cov2(Lambda,Acov);
        transform_rank1cov2(coords->dx_dxp,Acov);

	ploc = ph[i][j][k];
	ploc[B1] = Acov[1];
	ploc[B2] = Acov[2];
	ploc[B3] = Acov[3];
  }
    
  /*************************************************************************************
    Calculate the poloidal field components given A_\phi
       --  Note that we can take derivatives w.r.t. numerical coordinates here because 
           the connection is symmetric in its lower indices, we just have to use 
           the numerical metric determinant; 
  *************************************************************************************/
  curl_of_A(ph,p);


  /*************************************************************************************
    Calculate global internal energy and magnetic field energy:
  *************************************************************************************/
  bsq_int = u_int = 0.;

  N1_LOOP  N2_LOOP   N3_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); 
    ploc = p[i][j][k];
    bsq_int += geom->g * bsq_calc(ploc,geom) ;
    u_int   += geom->g * ploc[UU]; 
  }

  /************************************************************************************
     Normalize magnetic field to match specified value of beta (using umax);
   ************************************************************************************/
  mpi_global_sum( &bsq_int );
  fprintf(stdout,"initial bsq_int = %28.18e \n", bsq_int); fflush(stdout); 

  mpi_global_sum( &u_int );
  fprintf(stdout,"initial u_int = %28.18e \n", u_int); fflush(stdout); 

  beta_act = (gam - 1.)*u_int/(0.5*bsq_int) ;
  fprintf(stdout,"initial beta: %28.18e \n",beta_act); fflush(stdout);

  return;

}

