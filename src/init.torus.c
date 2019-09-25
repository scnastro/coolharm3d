
#include "decs.h"


double rin, rdiskout, beta, rhomax=0, umax=0, bsq_max=0;

double  ranc(int iseed); 

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

 Fishbone Moncrief Torus (w/ or w/o Magnetic field, specify below) 

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
	gam = (4./3.) ;
	cour = 0.8;


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
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	//	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
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
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
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
	DT_out[OUT_HDF5]    = 1. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = 1.;         /* radiative flux dumps       */
	n_restart = 1000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	dx[0] = 1.e-5;
	calc_all_geom() ;  
	dt_global_min = find_min_dt();

	//	dx[0] = cour * dt_global_min; /* Discretization size in X0 direction, time step*/

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
  struct of_coord *coords;
  
  void init_prim(int i, int j, int k, struct of_coord *coords, double *pr);
  void set_mag_field( void ) ; 

//  if( METRIC_TYPE_CHOICE != METRIC_KS_SPHERICAL ) { 
//    fprintf( stderr, "This initial data type was written for KS spherical coordinates !\n");
//    fflush(stderr); 
//    fail( FAIL_BASIC,0 ) ; 
//  }

  /*******************************************************************************
     Set the hydrodynamic quantities (rho,uu,v^i) : 
  *******************************************************************************/
  ranc(0);
  LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);	
    
    init_prim( i, j, k, coords, p[i][j][k]);
  }
  mpi_global_max(&rdiskout); 
  fprintf(stdout,"rdiskout = %28.18e \n", rdiskout) ; fflush(stdout);

  /*******************************************************************************
    Normalize the densities: 
  *******************************************************************************/
  mpi_global_max(&rhomax); 
  mpi_global_max(&umax); 
  fprintf(stdout, "init_data(): orig rhomax = %28.18e \n", rhomax); 
  fprintf(stdout, "init_data(): orig   umax = %28.18e \n",   umax);   fflush(stdout); 

  LOOP {
    p[i][j][k][RHO]  /=  rhomax; 
    p[i][j][k][ UU]  /=  rhomax; 
  }
  umax /= rhomax; 
  rhomax = 1. ; 

  fprintf(stdout, "init_data(): new  rhomax = %28.18e \n", rhomax); 
  fprintf(stdout, "init_data(): new    umax = %28.18e \n", umax);   fflush(stdout); 


  /*******************************************************************************
    Set the magnetic field :  
  *******************************************************************************/
  /* Correct bad points and setup boundary values since we will require them for B^i */
  fixup(p) ;
  bounds(p,0) ;

#if( HYDRO_ONLY ) 
  LOOP { 
    p[i][j][k][B1] = p[i][j][k][B2] = p[i][j][k][B3] =  0.; 
  }
#else 
  set_mag_field(); 
#endif 

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
  
  //  a = 0.9375;          /* Spin of the black hole in units of M */ 
  a = 0.;          /* Spin of the black hole in units of M */ 
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);

  R0       = 0.;          /* Offset in Radius from 1.            */
  Rout     = 40.;         /* Radial extent of the grid           */
  n_within_horizon = 5;   /* number of cells within horizon */
  Rin       = Rin_calc();

  h_slope  = 0.35;      /* Severity of the focusing            */
  X1_slope = 1.;       /* Severity of transition              */
  X1_0     = log(1.e6);     /* Location of transition in X1 units  */ 
  
  th_cutout = SMALL * M_PI;       /* Angular size to excise from the axis */
  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */

  /* New way : */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */

  /* Old way: */
//  th_beg    = 0.;        /* beg/end set the limits of theta  and */
//  th_end    = M_PI;   /*  can override th_cutout              */
//  th_length = M_PI;    /* Angular size in theta of the grid    */

  /* COORD_DIAGONAL3 specific parameters (note h_slope is overloaded for dtheta_min/(pi/2) */
  diag3_exponent = 7 ;                            /* note that this should be an odd integer greater than 1  */
  diag3_factor   = 1. - 2*th_cutout/M_PI - h_slope;  /* Temporary variable */
  if( (diag3_exponent < 3) || (diag3_exponent % 2 == 0 )  ) { fprintf(stderr,"set_special_coord(): bad value for  diag3_exponent  = %d \n",diag3_exponent); fflush(stderr); exit(4112); }


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


/***********************************************************************/
/***********************************************************************
        SPECIAL TORUS FUNCTIONS 
***********************************************************************/



/***********************************************************************************/
/***********************************************************************************
  init_prim(): 
  ------------------
   This is pretty much the same routine as in HARM2D to set Fishbone & Moncrief 
    initial data. 
 ***********************************************************************************/
void init_prim(int i, int j, int k, struct of_coord *coords, double *pr);
{
  int  prograde;
  double r,th;
  double randfact;
  double sth, cth;
  double ur, uh, up, u, rho;
  double ucon[NDIM], ucon2[NDIM];
  struct of_geom geom_bl;
  struct of_geom *geom;

  /* for disk interior */
  double l, lnh, expm2chi, up1;
  double DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  double kappa, hm1;
  double rmax;

  const double slowdown_factor = 1.0;      /* reduce u_phi by this amount */

  double  lfish_calc2(double  r, int prograde );
  double  lfish_calc(double rmax);

  /* Set the parameters that define the torus : */
  prograde = 1;                      /* Set to non-zero for prograde disks, or set to 0 */
  rin = 10. ;                         /* Radius of inner edge of torus */
  rmax = 16. ;                       /* Radius of pressure maximum    */

  l = lfish_calc2(rmax,prograde) ;   /* Ang. mom. at pressure maximum used throughout disk*/
  //  l = lfish_calc(rmax) ;


  kappa = 1.e-3 ;       /* Constant entropy factor in polytropic EOS used at t=0 */
  beta = 1.e2 ;         /* Plasma beta parameter */
  //  randfact = 4.e-2;     /* Magnitude of random perturbation to internal energy density */
  randfact = 0.;     /* Magnitude of random perturbation to internal energy density */

  r = coords->x[RR];  th = coords->x[TH];
  sincos( th, &sth, &cth );

  /* calculate lnh */
  DD = r * r - 2. * r + a * a;  // 26.7
  AA = (r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth;  // 1500 
  SS = r * r + a * a * cth * cth;  // 38.2
  
  thin = M_PI / 2.;
  //  sincos( thin, &sthin, &cthin );
  sthin = 1.;    cthin = 0.; 

  DDin = rin * rin - 2. * rin + a * a;  // 24.9 
  AAin = (rin * rin + a * a) * (rin * rin + a * a) 
    - DDin * a * a * sthin * sthin;                  // 1360 
  SSin = rin * rin + a * a * cthin * cthin; // 37 
  
  if (r >= rin) {
    lnh = 0.5*log((1. + sqrt(1. + 4.*(l*l*SS*SS)*DD/
			     (AA*sth*AA*sth)))/(SS*DD/AA)) 
      - 0.5*sqrt(1. + 4.*(l*l*SS*SS)*DD/(AA*AA*sth*sth))
      - 2.*a*r*l/AA 
      - (0.5*log((1. + sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
			    (AAin*AAin*sthin*sthin)))/(SSin*DDin/AAin)) 
	 - 0.5*sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
		    (AAin*AAin*sthin*sthin)) 
	 - 2.*a*rin*l/AAin ) ;
  } 
  else { 
    lnh = 1.;
  }
  

  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {

    rho = -RHOMIN * 1.E-14;   /* Make sure that these values are replaced by floor model */
    u = UUMIN * 1.E-14;

    pr[RHO] = rho;
    pr[UU] = u;

    /* Set coordinate velocity in Kerr-Schild coordinates to zero : */
    // First find Kerr-Schild coordinate velocity: 
    pr[U1] = 0.;    pr[U2] = 0.;    pr[U3] = 0.;

//    // Now find prim. velocity in numerical coordinates: 
//    ucon_calc( pr, &geom, ucon ); 
//    transform_rank1con( x, xp, ucon ); 
//    ucon2pr( pr, ucon, gcon[i][j][CENT] );
  }
  /* region inside magnetized torus; u^i is calculated in
     Boyer-Lindquist coordinates, as per Fishbone & Moncrief, so it
     needs to be transformed at the end */
  else {

    if( r > rdiskout  )  rdiskout = r; 

    hm1 = exp(lnh) - 1.;
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));
    u = kappa * pow(rho, gam) / (gam - 1.);
    ur = 0.;
    uh = 0.;
    
    /* calculate u^phi */
    expm2chi = SS * SS * DD / (AA * AA * sth * sth);
    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0) - 0.5));
    ucon[1] = ur ;
    ucon[2] = uh ;
    ucon[3] = slowdown_factor * up;
    
    /* We start with 4-velocity in BL coordinates */
    get_special_geometry( coords, &geom_bl, METRIC_BL_SPHERICAL );  
    setutcon( ucon, geom_bl.gcov );             /* set t-component of our find BL 4-velocity */
    bl_to_ks_con( coords->x, ucon, ucon2 );          /* find KS 4-vel. from BL 4-vel.   */
    transform_rank1con2(coords->dxp_dx,ucon2); // transform from (r,th,ph) to (x1,x2,x3)
    get_geometry(i,j,k,CENT,ncurr,geom); 
    ucon2pr( pr, ucon2, geom->gcon );  /* finally set prim. velocities */

    if(  rho > rhomax )             { rhomax = rho ; } 
    if( ( u  > umax) && (r > rin) ) { umax   =  u  ; } 

//    fprintf(stdout,"blah : %g  %g  %g  %g  %g  %g  \n", r,rin,rmax,th,cth,sth); 
//    fprintf(stdout,"blah2: %g  %g  %g  %g  %g  %g  \n", l,lnh,DD,AA,SS,DDin); 
//    fprintf(stdout,"blah3: %g  %g  %g  %g  %g  %g  \n", AAin,SSin,ucon[0],ucon[1], 
//	    ucon[2],ucon[3]); 
//    fprintf(stdout,"blah4: %g  %g  %g  %g  %g  %g  %g  \n", rho,u,ur,uh,up,rhomax,umax);
//    fprintf(stdout,"blah5: %g  %g  %g  %g  %g  %g  %g  \n", rho,u,ur,uh,up,rhomax,umax);
//    fprintf(stdout,"blah6: %g  %g  %g  \n", thin,cthin,sthin); 
//    fflush(stdout);

  }

  return;

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
void set_mag_field( void )
{
  int i,j,k,l;

  double rho_avg, q, bsq_ij, beta_act, norm ; 
  struct of_geom *geom;  
  
  /*************************************************************************************
    Calculate the azimuthal component of the vector potential, A_\phi: 
       --  Use F[] as a work array to hold A_\phi . 
       --  We assume here that rho is initially independent of X3 
       --  Need to calc. A_\phi on first ghost cells for finite diff. calculation for B^i 
   *************************************************************************************/
  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) { 
    
    rho_avg =  0.25*( p[i  ][j  ][N3S][RHO] +
                      p[i-1][j  ][N3S][RHO] +
                      p[i  ][j-1][N3S][RHO] +
                      p[i-1][j-1][N3S][RHO]   ) ;

    q = rho_avg / rhomax - 0.2 ; 
    
    F[i][j][0][0][0] = ( q > 0. ) ?  q  :  0. ; 

  }
    
  /*************************************************************************************
    Calculate the poloidal field components given A_\phi
       --  Note that we can take derivatives w.r.t. numerical coordinates here because 
           the connection is symmetric in its lower indices, we just have to use 
           the numerical metric determinant; 
  *************************************************************************************/
  N1_LOOP  N2_LOOP { 
    get_geometry(i,j,0,CENT,ncurr,geom); 
    /* flux-ct */
    p[i][j][0][B1] = -( F[i  ][j  ][0][0][0] 
		      - F[i  ][j+1][0][0][0] 
		      + F[i+1][j  ][0][0][0] 
                      - F[i+1][j+1][0][0][0] 
		       )/(2.*dx[2]*geom->g) ;

    p[i][j][0][B2] = (  F[i  ][j  ][0][0][0] 
		      + F[i  ][j+1][0][0][0] 
		      - F[i+1][j  ][0][0][0] 
		      - F[i+1][j+1][0][0][0] 
			)/(2.*dx[1]*geom->g) ;
    
    p[i][j][0][B3] = 0. ;
    
    bsq_ij = bsq_calc(p[i][j][0],geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }

  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) {  F[i][j][0][0][0] = 0.; }


  mpi_global_max( &bsq_max );
  fprintf(stdout,"initial bsq_max = %28.18e \n", bsq_max); fflush(stdout); 


  /************************************************************************************
     Normalize magnetic field to match specified value of beta (using umax);
  ************************************************************************************/
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
  fprintf(stdout,"initial beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);

  norm = sqrt(beta_act/beta) ;
  bsq_max = 0. ;
  
  N1_LOOP  N2_LOOP { 
    get_geometry(i,j,0,CENT,ncurr,geom) ;

    p[i][j][0][B1] *= norm ;
    p[i][j][0][B2] *= norm ;

    bsq_ij = bsq_calc(p[i][j][0],geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }

  mpi_global_max( &bsq_max );
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
  fprintf(stdout,"final beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);

  /************************************************************************************
     Copy B^i  from X3=0 slice to other slices X3 slices : 
  ************************************************************************************/
  ALL_LOOP  for( l = B1 ; l <= B3 ; l++ ) { 
    p[i][j][k][l] = p[i][j][0][l] ; 
  }
    

  return;

}

/**********************************************************************/
/**********************************************************************
  lfish_calc2(): 
  ---------------
    -- simplified and more general version of lfish_calc() in that it 
       allows one to specify whether the disk is prograde (prograde != 0 )
       or retrograde (prograde == 0).
 
**********************************************************************/
double  lfish_calc2(double  r, int prograde )
{
  double f_sign, sqrt_r, rsq;

  f_sign = ( prograde ) ? 1. : -1. ;
  sqrt_r = sqrt(r);
  rsq = r*r;
  asq = a*a;

  return( f_sign * ( rsq*( rsq + asq ) - 2*r*asq  - f_sign*a*sqrt_r*(rsq-asq))
	       / ( sqrt_r*r*(rsq - 3.*r + 2.*f_sign*a*sqrt_r))   
	    );

}


/**********************************************************************/
/**********************************************************************
  lfish_calc(): 
  ---------------
    -- returns the angular momentum parameter for the Fishbone & Moncrief
       disk give the radius of the torus's pressure maximum (r). 
 
**********************************************************************/
double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
	   ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
	    sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
	    ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) +
					  pow(a,
					      2) * (2. + r))) / sqrt(1 +
								     (2.
								      *
								      a)
								     /
								     pow
								     (r,
								      1.5)
								     -
								     3.
								     /
								     r)))
	  / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
	     (pow(a, 2) + (-2. + r) * r))
	  );
}
