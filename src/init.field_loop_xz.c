
#include "decs.h"


/* Determine whether we specify the z-component of vector potential, else we specify B-field (and generate monopoles) */
/* The B-field option is to follow the test done in the Einstein Toolkit's MHD paper.  */
#define SET_AY (1) 

/***********************************************************************************/
/* Parameter of the run : */ 
/***********************************************************************************/
//static const double v0          = (1./12.) ;     /* Magnitude of the velocity component of the field loop: */
static const double v0          = 0. ;     /* Magnitude of the velocity component of the field loop: */
static const double loop_radius = 0.3  ;     /* Radius of loop */
static const double Ay0         = 1.e-3;     /* Magnitude of the vector potential */
static double theta; /* Angle between the x & y velocity components */



static void set_mag_field( void ) ; 


/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

  Advected field loop test.  Assumes Cartesian coordinates and Minkowski. 
  Like init.field_loop.c  but  rotated to lie in x-z plane instead of x-y plane. 

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
	  GridLength[0] = 24.;         /* Length of X0 dimension (evolution period) */ 
	} 

	dx[0] = 1.e-5;               /* Discretization size in X0 direction */

#if( COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN )
	GridLength[1] = 1. ;        /* Length of X1 dimension */ 
	GridLength[3] = 1. ;        /* Length of X2 dimension */ 
	startx[1] = 0.;               /* Set the Physical Minimum boundary  */
	startx[3] = 0.;              /* Set the Physical Minimum boundary  */
#else
	GridLength[1] = 2. ;        /* Length of X1 dimension */ 
	GridLength[3] = 1. ;        /* Length of X2 dimension */ 
	startx[1] = -1.;               /* Set the Physical Minimum boundary  */
	startx[3] = -0.5;              /* Set the Physical Minimum boundary  */
#endif

	dx[1] = GridLength[1]/N1 ;   /* Discretization size in X1 direction */
	dx[3] = GridLength[3]/N3 ;   /* Discretization size in X2 direction */
	dx[2] = dx[3] ;              /* Discretization size in X3 direction */

	GridLength[2] = N2*dx[2] ;        /* Length of X3 dimension */ 

	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
	startx[2] = -0.5*dx[2];        /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	DT_out[OUT_ASCII]   = GridLength[0]/10. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/40. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = GridLength[0]/201. ;      /* hdf frequency */
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
  double *ploc;
  double vcon[NDIM], vconp[NDIM], bcon[NDIM];

  double x,z,r;
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
  double vz = 0.5*vx;
  double vy = 0.;
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
    z = coords->x[3];
    r = sqrt(x*x + z*z);
    ploc = p[i][j][k];

    ploc[RHO] = 1.;
    ploc[UU]  = press/(gam - 1.) ;

#if( SET_AZ ) 
      ploc[B1] = 0.;
      ploc[B2] = 0.;
      ploc[B3] = 0.;
#else
    if( r < loop_radius ) { 
      bcon[0] = 0.;
      bcon[1] = -Ay0*z/r ;
      bcon[3] =  Ay0*x/r ;
      bcon[2] = 0.;
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
  set_mag_field();
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
static void set_mag_field( void )
{
  int i,j,k,l;

  double beta_act, norm, bsq_int, u_int ; 
  double x,z,r ; 
  double *ploc;
  struct of_geom *geom;  
  struct of_coord *coords;  
  extern void curl_of_A( double ****A , double ****prim ) ;

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
	x = coords->x[1];
	z = coords->x[3];
	r = sqrt(x*x + z*z);

	Acov[0] = 0.;
	Acov[1] = 0.;
	Acov[2] = 0.;
	Acov[3] = 0.;

	if( r <= loop_radius ) { 
	  Acov[2] = Ay0 * (loop_radius - r); 
	  transform_rank1cov2(coords->dx_dxp,Acov);
	}

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

