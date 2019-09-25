
#include "decs.h"


/* Set the configuration of the field used at t=0: 
   (see set_mag_field() for details)                */
#define DIPOLE_FIELD         (0) 
#define QUADRUPOLE_FIELD     (1)
#define FIELD_CONFIG_CHOICE  (DIPOLE_FIELD)


#define R_DISK_MAX  (1.e6)   /* Maximum extent of a disk */
#define L_DISK_MAX  (1.e4)   /* Maximum specific angular momentum in disk */

const double l_ang_precision = 1.e-6;  /* precision to determine l_ang or rpmax */

static int disktype; 
static double rin, beta, rho_pmax, rdiskout=0., rhomax=0., umax=0., bsq_max=0.;
static double eta, eta_2q, utcov_in, sign_l, lambda_in, f_l_in, k_eta;
static double adiab, rin, rpmax, rout, GRscale, l_ang_in, q;
static double exp_1, exp_2, exp_3, exp_4, alpha;
static struct of_geom geom_in;


static double find_l_ang_in( void );
static double get_rpmax( double l );
static void  rho_func( double x[NDIM], struct of_geom *geom, 
		       double l_ang_in, double *rho, double *uu, 
		       double *utcov, double *l_ang, double *utcon, double *upcon ) ;

double  ranc(int iseed); 

double height_denom[N1TOT] = {0.}; 
double height_numer[N1TOT] = {0.}; 
double height_avg;

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

"KD" or Keplerian disk initial data ala De Villiers and Hawley 2002-2004
  B-field  follows density contours. 

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
	DT_out[OUT_HDF5]    = 20. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = 1.;         /* radiative flux dumps       */
	n_restart = 1000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	calc_all_geom() ;  
	dt_global_min = find_min_dt();

	//	dx[0] = cour * dt_global_min; /* Discretization size in X0 direction, time step*/
	dx[0] = 1.e-5;

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
  double r;
  
  void init_prim(int i, int j, int k, struct of_coord *coords, double *pr);
  void set_mag_field( void ) ; 

  if( METRIC_TYPE_CHOICE != METRIC_KS_SPHERICAL ) { 
    fprintf( stderr, "This initial data type was written for KS spherical coordinates !\n");
    fflush(stderr); 
    fail( FAIL_BASIC,0 ) ; 
  }

  /*******************************************************************************
     Set the hydrodynamic quantities (rho,uu,v^i) : 
  *******************************************************************************/
  ranc(myid+1);
  LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);	
    
    init_prim( i, j, k, coords, p[i][j][k]);
  }
  mpi_global_max(&rdiskout); 
  fprintf(stdout,"rdiskout = %28.18e \n", rdiskout) ; fflush(stdout);

  /*******************************************************************************
    Print the disk height information:
  *******************************************************************************/
  j = NG + N2/2; 
  k = 0; 
  g = 0;
  height_denom[0] = height_numer[0] = height_avg = 0.;
  fflush(stdout);  fprintf(stdout,"\n"); fflush(stdout);
  N1_LOOP  if(height_denom[i] != 0.)  { 
    coord(i,j,k,CENT,X) ; 
    r = exp(X[1]);
    if( fabs( (r - rpmax)/(rpmax*dx[1]) ) < 0.5 ) { 
      fprintf(stdout,"disk-thickness   %26.16e   %26.16e \n",r,height_numer[i]/height_denom[i]/r); fflush(stdout);
    }
  }

  //  exit(100);

  /*******************************************************************************
    Normalize the densities: 
  *******************************************************************************/
  mpi_global_max(&rhomax); 
  mpi_global_max(&umax); 
  fprintf(stdout, "init_data(): orig rho_pmax = %28.18e \n", rho_pmax); 
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
  int i;
  
  a = 0.9;          /* Spin of the black hole in units of M */ 
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);

  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  //  R0       = 0.9*r_horizon;  /* Offset in Radius from 1.  (GRMHD-like)       */
  Rout     = 80.;         /* Radial extent of the grid           */
  n_within_horizon = 25;   /* number of cells within horizon */
  Rin       = Rin_calc();


  h_slope   = 0.0729419909586;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.0* M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
  th_cutout = M_PI*SMALL;     /* Angular size to excise from the axis */
                                  /* 0.02655 gives a cutout of ~0.045Pi for h_slope = 0.35 */

  /* New way : */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */

  /* Old way: */
//  th_beg    = 0.;               /* beg/end set the limits of theta  and */
//  th_end    = M_PI;             /*  can override th_cutout              */
//  th_length = M_PI;             /* Angular size in theta of the grid    */
  
  /* COORD_DIAGONAL3 specific parameters (note h_slope is overloaded for dtheta_min/(pi/2) */
  diag3_exponent = 9;                            /* note that this should be an odd integer greater than 1  */
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

#define KEPLERIAN (0) 
#define FAT       (1)


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
void init_prim(int i, int j, int k, struct of_coord *coords, double *pr)
{
  int    l;
  double r,th, rho, uu;
  struct of_geom geom_bl;
  struct of_geom *geom;
  struct of_coord coords_loc;

  double lambda, Omega, h, f_l, utcov, l_ang, epsilon;
  double ucon[NDIM], ucov[NDIM], ucov_bl[NDIM], ucon_bl[NDIM];
  double utcon, upcon;

  static int first_call = 1;


  /***********************************************************************************************
    Set parameters used for all points : 
  **********************************************************************************************/  
  
  if( first_call ) { 
    /* for disk interior */
    disktype = KEPLERIAN;
    adiab    = 0.01;    /* constant "K" in P = K rho**gam , assume this EOS in disk */
//    rin      = 15.;     /* Radius of inner edge of the disk ?? */
//    rpmax    = 25.;     /* Radius of the pressure maximum ?? */
    rin      = 20.;     /* Radius of inner edge of the disk ?? */
    rpmax    = 35.;     /* Radius of the pressure maximum ?? */
    rout     = 0.0;     /* Not used */
    GRscale  = 0.01;    /* Relative magnitude of random perturbations to u */
    beta     = 1.e2 ;   /* Plasma beta parameter */
    //    q        = 1.68;    /* exponent in definition of Omega  */
    q        = 1.61871;    /* exponent in definition of Omega  */
  
    alpha = q/(q-2.);
    exp_1=(2.*(q-1.))/q;
    exp_4 = alpha + 1.;   // 2*(q-1)/(q-2)
    exp_2 = 1./exp_4;      // (q-2)/(2*(q-1))
    exp_3 = 1./alpha;     // (q-2)/q

    /* Set the parameters defined at rin, the radial coordinate of innermost part of the disk */
    coord_of_r(rin, &coords_loc);
    get_special_geometry(&coords_loc, &geom_in, METRIC_BL_SPHERICAL );

    l_ang_in = find_l_ang_in();    /* ang. mom. at rin, this sets other rin values too  */

    /* Find the values at rpmax so that we know that the real pressures max. is at rpmax : */
    coord_of_r(rpmax, &coords_loc);
    get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );

    rho_func( coords_loc.x, &geom_bl, l_ang_in, &(pr[RHO]), &(pr[UU]), &utcov, &l_ang, &utcon, &upcon );

    rho_pmax = pr[RHO];

    fprintf(stdout,"\n##################################################\n");
    fprintf(stdout,"  KD TORUS PARAMETERS \n------------------------------------\n");
    fprintf(stdout,"\t  disktype  =  %d    (Keplerian = %d,  Fat = %d) \n",disktype,KEPLERIAN,FAT); 
    fprintf(stdout,"\t  adiab     =  %28.18e \n",adiab   ); 
    fprintf(stdout,"\t  rin       =  %28.18e \n",rin     ); 
    fprintf(stdout,"\t  rpmax     =  %28.18e \n",rpmax   ); 
    fprintf(stdout,"\t  rout      =  %28.18e \n",rout    ); 
    fprintf(stdout,"\t  GRscale   =  %28.18e \n",GRscale ); 
    fprintf(stdout,"\t  l_ang_in  =  %28.18e \n",l_ang_in); 
    fprintf(stdout,"\t  beta      =  %28.18e \n",beta    ); 
    fprintf(stdout,"\t  q         =  %28.18e \n",q       ); 
    fprintf(stdout,"  Values at rin:\n------------------------------------\n");
    fprintf(stdout,"\t  eta       =  %28.18e \n",eta      ); 
    fprintf(stdout,"\t  eta_2q    =  %28.18e \n",eta_2q   ); 
    fprintf(stdout,"\t  k_eta     =  %28.18e \n",k_eta    ); 
    fprintf(stdout,"\t  lambda_in =  %28.18e \n",lambda_in); 
    fprintf(stdout,"\t  sign_l    =  %28.18e \n",sign_l   ); 
    fprintf(stdout,"\t  f_l_in    =  %28.18e \n",f_l_in   ); 
    fprintf(stdout,"\t  utcov_in  =  %28.18e \n",utcov_in ); 
    fprintf(stdout,"  Values at rpmax:\n------------------------------------\n");
    fprintf(stdout,"\t  l_ang     =  %28.18e \n",l_ang   ); 
    fprintf(stdout,"\t  utcov     =  %28.18e \n",utcov   ); 
    fprintf(stdout,"\t  rho       =  %28.18e \n",pr[RHO] ); 
    fprintf(stdout,"\t  uu        =  %28.18e \n",pr[UU ] ); 

    first_call = 0;
  }
  
  /***********************************************************************************************
    Set the values in the disk :
  **********************************************************************************************/  
  get_special_geometry( coords, &geom_bl, METRIC_BL_SPHERICAL );

  /* get rho and uu */ 
  rho_func( coords->x, &geom_bl, l_ang_in, &(pr[RHO]), &(pr[UU]), &utcov, &l_ang, &utcon, &upcon );

  /* set the coordinate velocity */
  if( pr[RHO] > 0. ) { 
//    ucov_bl[TT] = utcov;  
//    ucov_bl[RR] = 0.;
//    ucov_bl[TH] = 0.;
//    ucov_bl[PH] = -l_ang*utcov; 
//    bl_to_ks_cov(x,ucov_bl,ucov);  // transform from BL to KS coordinates 
//    transform_rank1cov(x,xp,ucov); // transform from (r,th,ph) to (x1,x2,x3)
    get_geometry(i,j,k,CENT,ncurr,geom); // get (x1,x2,x3) metric, from here on it's only (x1-3)
//    raise(ucov,&geom,ucon);        
//
//    fprintf(stdout,"ucon1 = %26.20e %26.20e %26.20e %26.20e \n",ucon[0],ucon[1],ucon[2],ucon[3]); 

    ucon_bl[TT] = utcon;
    ucon_bl[RR] = 0.;
    ucon_bl[TH] = 0.;
    ucon_bl[PH] = upcon;
    bl_to_ks_con(coords->x,ucon_bl,ucon);  // transform from BL to KS coordinates 
    transform_rank1con2(coords->dxp_dx,ucon); // transform from (r,th,ph) to (x1,x2,x3)

    //    fprintf(stdout,"ucon2 = %26.20e %26.20e %26.20e %26.20e \n",ucon[0],ucon[1],ucon[2],ucon[3]); 

    ucon2pr( pr, ucon, geom->gcon );  // get prim. velocities;
    if( coords->x[RR] > rdiskout ) { rdiskout = x[RR]; }
    rho = geom->g * pr[RHO] ;
    height_numer[i] += sqrt(geom_bl.gcov[2][2]) * fabs(coords->x[TH]-0.5*M_PI) * rho;
    height_denom[i] += rho;
  }
  else { 
    PLOOP pr[l] = 0.;
  }
    
  if( pr[RHO] > rhomax )  { rhomax = pr[RHO] ; } 
  if( pr[ UU] >   umax )  {   umax = pr[ UU] ; } 

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

  double rho_avg, rho_cut, Aphi, beta_act, norm, bsq_int, u_int ; 
  double x[NDIM], xp[NDIM];
  struct of_geom *geom;  

#if(   FIELD_CONFIG_CHOICE == QUADRUPOLE_FIELD )
  fprintf(stdout,"set_mag_field():  Using QUADRUPOLE_FIELD %d \n",FIELD_CONFIG_CHOICE); 
#elif( FIELD_CONFIG_CHOICE ==     DIPOLE_FIELD )
  fprintf(stdout,"set_mag_field():  Using DIPOLE_FIELD %d \n",FIELD_CONFIG_CHOICE); 
#else 
  fprintf(stderr,"set_mag_field():  Invalid value of FIELD_CONFIG_CHOICE =  %d \n",
	  FIELD_CONFIG_CHOICE);  fflush(stderr); 
  fail(FAIL_BASIC,0); 
#endif
  fflush(stdout);

  /*************************************************************************************
    Calculate the azimuthal component of the vector potential, A_\phi: 
       --  Use F[] as a work array to hold A_\phi . 
       --  We assume here that rho is initially independent of X3 
       --  Need to calc. A_\phi on first ghost cells for finite diff. calculation for B^i 
  *************************************************************************************/
  k = N3S; 

  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) { 

    rho_avg =  0.25*( p[i  ][j  ][N3S][RHO] +
                      p[i-1][j  ][N3S][RHO] +
                      p[i  ][j-1][N3S][RHO] +
                      p[i-1][j-1][N3S][RHO]   ) ;

#if(   FIELD_CONFIG_CHOICE ==     DIPOLE_FIELD )
    rho_cut = 0.25 * rhomax; 
    Aphi = rho_avg - rho_cut ;    /* dipole field */

#elif( FIELD_CONFIG_CHOICE == QUADRUPOLE_FIELD )
    rho_cut = 0.25 * rhomax; 
    get_coord(i,j,k,CORN,ncurr,coords);
    Aphi = (rho_avg - rho_cut) * cos(coords->x[TH]) / rhomax ;
#endif


    /* Truncate A_phi at zero so that field lines lie within disk */
    F[i][j][0][0][0] = ( rho_avg > rho_cut ) ?  Aphi  :  0. ; 

  }
    
  /*************************************************************************************
    Calculate the poloidal field components given A_\phi
       --  Note that we can take derivatives w.r.t. numerical coordinates here because 
           the connection is symmetric in its lower indices, we just have to use 
           the numerical metric determinant; 
  *************************************************************************************/
  bsq_int = u_int = 0.;

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

    bsq_int += geom->g * bsq_calc(p[i][j][0],geom) ;
    u_int   += geom->g * p[i][j][0][UU]; 
  }

  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) {  F[i][j][0][0][0] = 0.; }

  /************************************************************************************
     Normalize magnetic field to match specified value of beta (using umax);
   ************************************************************************************/
  mpi_global_sum( &bsq_int );
  fprintf(stdout,"initial bsq_int = %28.18e \n", bsq_int); fflush(stdout); 

  mpi_global_sum( &u_int );
  fprintf(stdout,"initial u_int = %28.18e \n", u_int); fflush(stdout); 

  beta_act = (gam - 1.)*u_int/(0.5*bsq_int) ;
  fprintf(stdout,"initial beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);

  norm = sqrt(beta_act/beta) ;
  bsq_max = 0. ;
  
  N1_LOOP  N2_LOOP { 
    p[i][j][0][B1] *= norm ;
    p[i][j][0][B2] *= norm ;
  }

  bsq_int = 0.;
  N1_LOOP  N2_LOOP { 
    get_geometry(i,j,0,CENT,ncurr,geom); 
    bsq_int += geom->g * bsq_calc(p[i][j][0],geom) ;
  }

  mpi_global_sum( &bsq_int );
  fprintf(stdout,"final bsq_int = %28.18e \n", bsq_int); fflush(stdout); 

  beta_act = (gam - 1.)*u_int/(0.5*bsq_int) ;
  fprintf(stdout,"final beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);
  

  /************************************************************************************
     Copy B^i  from X3=0 slice to other slices X3 slices : 
  ************************************************************************************/
  ALL_LOOP  for( l = B1 ; l <= B3 ; l++ ) { 
    p[i][j][k][l] = p[i][j][0][l] ; 
  }
    

  return;

}

/***********************************************************************************/
/***********************************************************************************
  find_l_bracket():
  ------------------

***********************************************************************************/
int find_l_bracket( int prograde,  double *l_lo, double *l_hi ) 
{
  int bad_bracket = 1;
  double lmid; 
  double r_lo, r_hi;
  double x,at,xsq;

  /* Set initial bounds for l_ang_in  */
  if( prograde ) { 
     at = a/M;
     x = sqrt(rin/M);  xsq = x*x;
     *l_lo = M * ( xsq * ( xsq - 2.*at ) + at*at ) / ( x * ( xsq - 2. ) + at );  /* Circular equatorial timelike geodesic */
     
     x = sqrt(rpmax/M);  xsq = x*x;
     *l_hi = M * ( xsq * ( xsq - 2.*at ) + at*at ) / ( x * ( xsq - 2. ) + at );  
  }
  else { 
    at = -a/M;
     x = sqrt(rin/M);  xsq = x*x;
     *l_lo = -M * ( xsq * ( xsq - 2.*at ) + at*at ) / ( x * ( xsq - 2. ) + at );  /* Circular equatorial timelike geodesic */

     x = sqrt(rpmax/M);  xsq = x*x;
     *l_hi = -M * ( xsq * ( xsq - 2.*at ) + at*at ) / ( x * ( xsq - 2. ) + at );  
  }

  r_lo = get_rpmax(*l_lo); 
  r_hi = get_rpmax(*l_hi); 

  if( (r_lo < 0.) && (r_hi < 0.) ) { 
    fprintf(stdout,"find_l_bracket():  Estimates failed!!  "); fflush(stdout); 
    exit(231);
  }

  if( (r_lo > 0.) && (r_lo < rpmax) && (r_hi > rpmax) ) { 
    return(0); 
  }

  /* One of the solutions were impossible, we need to determine which one and then find the other : */
  if( r_lo < 0. ) {   /* Bad l_lo   but  Good  l_hi */
    r_lo = r_hi;
    *l_lo = *l_hi;
  }
  else {           /* Bad l_hi   but  Good  l_lo */
    r_hi = r_lo;
    *l_hi = *l_lo;
  }

  if( r_hi > rpmax ) { 
    while( bad_bracket ) { 
      *l_lo -= 0.01*(*l_lo); 
      r_lo = get_rpmax(*l_lo); 
      if( r_lo  < rpmax )                {  bad_bracket = 0; } 
      if( fabs(*l_lo) > (L_DISK_MAX*M) ) {  bad_bracket = 0; } 
      if( ((*l_lo) * (*l_hi)) < 0. )     {  bad_bracket = 0; } 
    }
    if( r_lo > rpmax ) { 
      fprintf(stdout,"find_l_bracket():  Estimates failed  2 !!  l_lo = %26.16e ", *l_lo); fflush(stdout); 
      return(2);
    }
  }
  else { 
    while( bad_bracket ) { 
      *l_hi += 0.01*(*l_hi); 
      r_hi = get_rpmax(*l_hi); 
      if( r_hi  > rpmax )                {  bad_bracket = 0; } 
      if( fabs(*l_hi) > (L_DISK_MAX*M) ) {  bad_bracket = 0; } 
      if( ((*l_lo) * (*l_hi)) < 0. )     {  bad_bracket = 0; } 
    }
    if( r_hi < rpmax ) { 
      fprintf(stdout,"find_l_bracket():  Estimates failed  3 !!  l_hi = %26.16e ", *l_hi); fflush(stdout); 
      return(3);
    }
  }

  return(0);
}


/***********************************************************************************/
/***********************************************************************************
  find_l_ang_in():
  ------------------
   -- Does a bisection search in order to  find the value of l_ang_in that puts 
      the location of the pressure maximum at "rpmax"; 
 
   -- given the adiabatic EOS, the pressure max. is also the density max.; 

***********************************************************************************/
double find_l_ang_in( void )
{
  int do_search; 
  double error, rpmax_hi, rpmax_lo, r, l, l_lo, l_hi, dl, dr;

//  if( (REL_DIFF_FUNC(rin,15.) > 1.e-14)  || (REL_DIFF_FUNC(rpmax,25.) > 1.e-14) ) { 
//    fprintf(stdout,"find_l_ang_in(): Will most likely have to find new values for l_lo,hi !!\n");
//    fflush(stdout);
//    myexit(1);
//  }

  do_search = 1; 

  if( find_l_bracket(1,&l_lo,&l_hi) ) { 
    fprintf(stdout,"find_l_ang_in():  Cannot find bracket !!  ");  fflush(stdout);
    exit(21);
  }
  
//  l_lo = 4.;
//  l_hi = 4.5;

  dl = l_hi - l_lo;

  rpmax_lo =   get_rpmax(l_lo);
  rpmax_hi =   get_rpmax(l_hi);

  if( (rpmax_lo-rpmax)*(rpmax_hi-rpmax) > 0. ) { 
    fprintf(stderr,"find_l_ang_in(): Need a good bracket before starting!! \n");
    fprintf(stderr,"find_l_ang_in():  l_lo,l_hi = %g %g\n", l_lo, l_hi);
    fprintf(stderr,"find_l_ang_in():  r_lo,r_hi = %g %g\n", rpmax_lo,rpmax_hi);
  }

  fprintf(stdout,"#######################################################\n");
  fprintf(stdout," SEARCHING FOR l_ang_in....    \n");
  fprintf(stdout," l_ang_in bracket =  [ %26.16e  ,   %26.16e  ]     \n",l_lo,l_hi);
  fflush(stdout);

  while( do_search ) { 
    l = 0.5*(l_lo + l_hi); 
    r = get_rpmax(l);

    if( r > rpmax ) {   rpmax_hi = r;   l_hi = l;  }
    else            {   rpmax_lo = r;   l_lo = l;  }

    dr = fabs((rpmax_hi - rpmax_lo)/rpmax); 
    dl = fabs(2.*(l_hi - l_lo)/(l_hi+l_lo)); 
    
    if( MIN(dl, dr) < l_ang_precision )   {  do_search = 0; } 
    fprintf(stdout,"l,r = %28.18e  %28.18e \n", l,r); fflush(stdout);
  }

  l = 0.5*(l_lo + l_hi); 
  r = get_rpmax(l);

  fprintf(stdout,"l_ang_in      = %28.18e \n", l ); 
  fprintf(stdout,"Real rpmax    = %28.18e \n", r ); 
  fprintf(stdout,"Desired rpmax = %28.18e \n", rpmax ); 
  fflush(stdout);

  if( REL_DIFF_FUNC(r,rpmax) > 1.e-3 ) {
    fprintf(stderr,"Real and Desired rpmax  are too far apart!! \n  Try adjusting  l_lo and l_hi  \n "); 
    fflush(stderr); exit(200); 
  }

  return( l );
}

/***********************************************************************************/
/***********************************************************************************
  get_rpmax():
  ------------------
   -- return the location of rho's maximum given a value of l_ang_in; 
   -- uses bisection search on drho/dr;
***********************************************************************************/
double get_rpmax( double l )
{
  int finding_r_hi, finding_rpmax; 
  double r, r_lo, r_hi, rho, rho1, rho2, drho_lo, drho_hi, drho, dr, dr_big, dr_frac;
  double uu, utcov, l_ang, utcon, upcon;
  struct of_geom geom_bl;
  struct of_geom coords_loc;


  r_hi = r_lo = rin; 
  dr_frac = 1.e-3;
  dr_big = 0.5;
  dr = dr_big * dr_frac;
  
  r = r_lo ;

  coord_of_r( r, &coords_loc );  get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );
  rho_func( coords_loc.x, &geom_bl, l, &rho1, &uu, &utcov, &l_ang, &utcon, &upcon );

  coord_of_r( r+dr, &coords_loc );  get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );
  rho_func( coords_loc.x, &geom_bl, l, &rho2, &uu, &utcov, &l_ang, &utcon, &upcon);
  drho_lo = (rho2 - rho1)/dr;

  finding_r_hi = 1;

  fprintf(stdout,"r_lo, drho_lo, dr = %28.18e  %28.18e  %28.18e \n", r_lo, drho_lo,dr);
  fflush(stdout);

  while( finding_r_hi ) { 
    dr_big *= 2; 
    r_hi = r_lo + dr_big;
    r = r_hi ;
    dr = dr_big * dr_frac;

    coord_of_r( r-dr, &coords_loc );  get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );
    rho_func( coords_loc.x, &geom_bl, l, &rho1, &uu, &utcov, &l_ang, &utcon, &upcon );

    coord_of_r( r+dr, &coords_loc );  get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );
    rho_func( coords_loc.x, &geom_bl, l, &rho2, &uu, &utcov, &l_ang, &utcon, &upcon);

    drho_hi = 0.5*(rho2 - rho1)/dr;

    if( drho_lo*drho_hi < 0. ) {  finding_r_hi  = 0 ; } 

    fprintf(stdout,"r_hi, drho_hi, dr = %28.18e  %28.18e  %28.18e \n", r, drho_hi,dr);
    fflush(stdout);

    if( r_hi > (M*R_DISK_MAX) ) {  return( -1. );  }
  }

  finding_rpmax = 1; 

  while( finding_rpmax ) { 
    dr_big = r_hi - r_lo; 
    r = 0.5*(r_hi+r_lo) ;
    dr = dr_big * dr_frac;

    coord_of_r( r-dr, &coords_loc );  get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );
    rho_func( coords_loc.x, &geom_bl, l, &rho1, &uu, &utcov, &l_ang, &utcon, &upcon );

    coord_of_r( r+dr, &coords_loc );  get_special_geometry( &coords_loc, &geom_bl, METRIC_BL_SPHERICAL );
    rho_func( coords_loc.x, &geom_bl, l, &rho2, &uu, &utcov, &l_ang, &utcon, &upcon);

    drho = 0.5*(rho2 - rho1)/dr;
    
    if( drho >= 0. ) {  r_lo = r; }
    else             {  r_hi = r; }

    /* Need to find this a more precisely than l_ang_in so that rpmax 
      changes smoothly with l_ang_in : */
    if( drho == 0. ) { finding_rpmax = 0 ; } 
    if( (dr_big/r) <  0.1*l_ang_precision ) {  finding_rpmax = 0 ; } 

//    fprintf(stdout,"r_lo, r_hi, drho_max, dr = %28.18e  %28.18e  %28.18e  %28.18e \n", 
//	    r_lo,r_hi, drho,dr);
//    fflush(stdout);
  }    

  r = r_lo; 

  return( r ); 
}

/***********************************************************************************/
/***********************************************************************************
  rho_func():
  ------------------
   -- calculate the density at a specific "r" and "l_ang_in"
***********************************************************************************/
void  rho_func( double x[NDIM], struct of_geom *geom, 
		double l_ang_in, double *rho, double *uu, double *utcov, double *l_ang, 
		double *utcon, double *upcon) 
{
  double lambda, Omega, Omega_in, f_l, h, epsilon;

  /* Quantities defined at rin: */
  sign_l = (l_ang_in >= 0.)   ?  1.  :  -1. ;
  lambda_in = sqrt( -geom_in.gcon[TT][TT] / geom_in.gcon[PH][PH] );  
  eta = fabs(l_ang_in) * pow(lambda_in, (q-2.));
  eta_2q = sign_l * pow( eta, 2./q );
  k_eta = pow( eta , (-2./(q-2.)) );
  Omega_in = eta * pow( lambda_in, -q );

  //  f_l_in = pow( fabs(1. - k_eta*pow(l_ang_in,exp_4)) , exp_2 );
  f_l_in = pow( fabs(1. - k_eta*pow(l_ang_in,exp_4)) , exp_2 );
  utcov_in = -1./sqrt(-geom_in.gcon[TT][TT] 
		      + l_ang_in*(2*geom_in.gcon[TT][PH] - l_ang_in*geom_in.gcon[PH][PH]));

  /* Quantities defined at x[] : */
  lambda = sqrt( -geom->gcon[TT][TT] / geom->gcon[PH][PH] ); 
  Omega = eta * pow( lambda, -q );
  
  if( disktype == KEPLERIAN ) { 
    *l_ang = eta_2q * pow( Omega, exp_3 );
    //    f_l   = pow( fabs(1. - k_eta*pow((*l_ang),exp_4)) , exp_2 );
    f_l   = pow( fabs(1. - k_eta*pow((*l_ang),exp_4)) , exp_2 );
  }
  else { 
    *l_ang = l_ang_in;
    f_l   = f_l_in; 
  }

  *utcov = -1./sqrt(-geom->gcon[TT][TT] 
		    + (*l_ang)*(2*geom->gcon[TT][PH] - (*l_ang)*geom->gcon[PH][PH]));

  *utcon = 1./sqrt( -(geom->gcov[TT][TT] + Omega * ( 2.*geom->gcov[TT][PH] + Omega*geom->gcov[PH][PH] ) ) );
  *upcon = Omega * (*utcon);

  h = utcov_in * f_l_in / ( (*utcov) * f_l ); 
  if( (h > 1.) && (x[RR] > rin) ) {
    epsilon = (h - 1.) / gam; 
    *rho =  pow( (epsilon * (gam - 1.) / adiab), (1./(gam-1.)) );
    *uu  = (*rho) * epsilon * (1. + GRscale * (ranc(0) - 0.5))   ;
  }
  else { 
    *uu = *rho = -1.;
  }

  return;
}


#undef KEPLERIAN
#undef FAT      

