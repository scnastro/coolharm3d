
#include "decs.h"

static double rhomax=0., umax=0., bsq_max=0.;

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

   -- initial conditions to recreate the distributions of the primitive variables
      seen in the funnel at  t=2400M  of the run in 
             /data1/gauss/harm3d/devel/2d/harm3d-kd9-2d-hi-ee/

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
	double X[NDIM], x[NDIM];

	void calc_all_geom(void) ;
	void set_special_coord(void);

	/* Set the parameters that define the grid and other constants : */ 
	gam = (5./3.) ;
	cour = 0.4;


	/* Coordinate dependent quantities : */ 
	set_special_coord();

	t = 0. ;                     /* Initial time */ 

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 10000. ;         /* Length of X0 dimension (evolution period) */ 
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	GridLength[1] = Rout-Rin  ;   /* Length of X1 dimension */ 
	GridLength[2] = th_length ;   /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
	//	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ;              /* Length of X1 dimension */ 
	GridLength[2] = xi_diag2[n_diag2_lines]-xi_diag2[0]  ; /* Length of X2 dimension */ 
#endif
	GridLength[3] = M_PI/2. ;                 /* Length of X3 dimension */ 

	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];

	/**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
	**************************************************************************/
	startx[0] =  0.;                /* Set the Physical Minimum boundary  */
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	startx[1] = Rin;                /* Set the Physical Minimum boundary  */
	startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	//	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = 10. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	n_restart = 1000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	dx[0] = 1.e-5;
	calc_all_geom() ;  
	dt_global_min = find_min_dt();

	/**************************************************************************
	  Print out basic grid parameters: 
	**************************************************************************/
	fprintf(stdout,"Length: %10.4e %10.4e %10.4e %10.4e\n", 
		GridLength[0], GridLength[1], GridLength[2], GridLength[3]);
	fprintf(stdout,"dx    : %10.4e %10.4e %10.4e %10.4e\n", 
		dx[0], dx[1], dx[2], dx[3]);
	fprintf(stdout,"startx: %10.4e %10.4e %10.4e %10.4e\n", 
		startx[0], startx[1], startx[2], startx[3]);

	fprintf(stdout,"rhomin,minlimt = %28.18e %28.18e \n", RHOMIN, RHOMINLIMIT);
	fprintf(stdout,"uumin ,minlimt = %28.18e %28.18e \n", UUMIN,  UUMINLIMIT);
	
	fflush(stdout);


	/* in case you want to test the metric and coord. transf. routines */
	//	test_geom();   

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
  int i, j, k, l, g ; 
  double X[NDIM];
  
  void init_prim( int i, int j, int k, double *xp, double *pr);
  void init_prim2(int i, int j, int k, double *xp, double *pr);
  void set_mag_field( void ) ; 

  if( METRIC_TYPE_CHOICE != METRIC_KS_SPHERICAL ) { 
    fprintf( stderr, "This initial data type was written for KS spherical coordinates !\n");
    fflush(stderr); 
    fail( FAIL_BASIC,0 ) ; 
  }

  /*******************************************************************************
     Set the hydrodynamic quantities (rho,uu,v^i) : 
  *******************************************************************************/
  LOOP {
    coord(i,j,k,CENT,X) ; 
    
    init_prim( i, j, k, X, p[i][j][k]);

  }

  /*******************************************************************************
    Normalize the densities: 
  *******************************************************************************/
  mpi_global_max(&rhomax); 
  mpi_global_max(&umax); 
  fprintf(stdout, "init_data(): orig rhomax   = %28.18e \n", rhomax); 
  fprintf(stdout, "init_data(): orig   umax   = %28.18e \n",   umax);   fflush(stdout); 
  
//  LOOP {
//    p[i][j][k][RHO]  /=  rhomax; 
//    p[i][j][k][ UU]  /=  rhomax; 
//  }
//  umax /= rhomax; 
//  rhomax = 1. ; 
//
//  fprintf(stdout, "init_data(): new  rhomax = %28.18e \n", rhomax); 
//  fprintf(stdout, "init_data(): new    umax = %28.18e \n", umax);   fflush(stdout); 


  /*******************************************************************************
    Set the magnetic field :  
  *******************************************************************************/
  /* Correct bad points and setup boundary values since we will require them for B^i */
  fixup(p) ;
  bounds(p,0) ;

}


/***********************************************************************/
/***********************************************************************
  set_special_coord(): 
  --------------------
   -- set coordinate specific (global) quantities;
***********************************************************************/
void set_special_coord(void)
{
  int i;
  
  a = 0.9;          /* Spin of the black hole in units of M */ 
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);

  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  //  R0       = 0.9*r_horizon;  /* Offset in Radius from 1.  (GRMHD-like)       */
  Rout     = 120.;         /* Radial extent of the grid           */
  n_within_horizon = 5;   /* number of cells within horizon */
  Rin       = Rin_calc();


  h_slope   = 0.14;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.04 * M_PI;     /* Angular size to excise from the axis */
                                  /* 0.02655 gives a cutout of ~0.045Pi for h_slope = 0.35 */
  th_cutout = SMALL;
  
  /* New way : */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */


#if( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  /* Parameters are set so poles are extra coarse, then match to slow, linear reduction 
      of discretization size as you tend to the equator.  Symmetric about xp2=xp2_3  */

  /*********************************************************************************/
  /* Collection of straight lines: */
//  n_diag2_lines = 6;      /* Number of segments */
//  d0_diag2      = 0.085;  /* Sets the floating scale of the discretization, similar to the value of dx2 at xp2=xi_xp2[0] */
//  alloc_diag2();          
//
//  /* Locations of transitions : (assume xp2 goes from 0 to 1 , use th_length to scale thinsg) */
//  xi_diag2[0] = SMALL;
//  xi_diag2[1] = (2./192.) ;
//  xi_diag2[2] = 4. * xi_diag2[1];
//  xi_diag2[3] = 0.5 ;
//  xi_diag2[4] = 2*xi_diag2[3] - xi_diag2[2];
//  xi_diag2[5] = 2*xi_diag2[3] - xi_diag2[1];
//  xi_diag2[6] = 2*xi_diag2[3] - xi_diag2[0];
//
//  /* Constant values of d^2 th / dx2^2   */
//  ci_diag2[0] =  0. ;   ci_diag2[1] = -2.;   ci_diag2[2] = -0.02; 
//  ci_diag2[3] = -ci_diag2[2];  ci_diag2[4] = -ci_diag2[1];   ci_diag2[5] = -ci_diag2[0]; 
//  
//  for(i=0;i<n_diag2_lines;i++)  { si_diag2[i] = 1.e4; }  /* make all strong transitions */
//  setup_diag2(0);

  /*********************************************************************************/
  /* Collection of fewer straight lines with smoother transitions: 
     (si[1] was tuned so that the dtheta is approximately symmetric about equator) */
  n_diag2_lines = 4;      /* Number of segments */
  alloc_diag2();          

  /* Locations of transitions : (assume xp2 goes from 0 to 1 , use th_length to scale thinsg) */
  xi_diag2[0] = SMALL;
  xi_diag2[1] = (4./192.) ;
  xi_diag2[2] = 0.5 ;
  xi_diag2[3] = 2*xi_diag2[2] - xi_diag2[1];
  xi_diag2[4] = 2*xi_diag2[2] - xi_diag2[0];

  /* Constant values of d^2 th / dx2^2   */
  ci_diag2[0] =  0.;    ci_diag2[1] = -0.25;  ci_diag2[2] = -ci_diag2[1];  ci_diag2[3] = -ci_diag2[0];
  di_diag2[0] =  1.4;   
  di_diag2[1] =  0.25;  
  di_diag2[2] =  di_diag2[1]+2*xi_diag2[2]*ci_diag2[1];  
  di_diag2[3] =  di_diag2[0];
  si_diag2[0] =  1.e4;  si_diag2[1] =  4.893617;  si_diag2[2] = 1.e2 ;  si_diag2[3] = 2.3e2;
  
  setup_diag2(1);

#endif

  fprintf(stdout,"a                = %28.18e \n", a);
  fprintf(stdout,"R0               = %28.18e \n", R0               );
  fprintf(stdout,"Rin              = %28.18e \n", Rin              );
  fprintf(stdout,"Rout             = %28.18e \n", Rout             );
  fprintf(stdout,"n_within_horizon = %6d     \n", n_within_horizon );
  fprintf(stdout,"h_slope          = %28.18e \n", h_slope          );
  fprintf(stdout,"X1_slope         = %28.18e \n", X1_slope         );
  fprintf(stdout,"X1_0             = %28.18e \n", X1_0             );
  fprintf(stdout,"r_isco           = %28.18e \n", r_isco           );
  fprintf(stdout,"r_horizon        = %28.18e \n", r_horizon        );
  fflush(stdout);

  return;
}


/***********************************************************************************/
/***********************************************************************************
  init_prim(): 
  ------------------
   -- This is supposed to follow the routine in inits.f of the Hawley et al. code; 
       -- inconsistencies are due to not knowing what various grid functions represent 
          in inits.f  (there is very little documentation in that file); 
   -- when in doubt, I used   
             De Villiers, Hawley, Krolik  ApJ 599 (2003)

   -- The disk parameters mentioned in this paper are used below; 
   -- Uses BL coordinates until the very end and then transforms the velocities to KS; 
***********************************************************************************/
void init_prim(int i, int j, int k, double *xp, double *pr)
{
  int    l;
  double x[NDIM];
  double sth,cth,cth5;

  static int first_call = 1;


  x_of_xp(x,xp);
  sincos(x[TH],&sth,&cth);
  cth5 = cth*cth*cth*cth*cth;
  
  pr[RHO] = 2.e-9;
  pr[UU]  = 1.e-5*pow(x[RR],-2.5);
  pr[B1]  = 1.2e-1*pow(x[RR],-2.5)*cth5;
  pr[B2]  = -1.e-2*pow(x[RR],-2.5)*sth;
  pr[B3]  = -3.e-2*pow(x[RR],-1.5)*cth5;
  pr[U1]  = 3.e-1/sqrt(x[RR]) - 5.e-3/sth;
  pr[U2]  = -5.e-3*sth*cth;
  pr[U3]  = 1.e-1/(x[RR]*sth);

  if( pr[RHO] > rhomax )  { rhomax = pr[RHO] ; } 
  if( pr[ UU] >   umax )  {   umax = pr[ UU] ; } 

  return;

}


