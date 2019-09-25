
#include "decs.h"
#include "metric.h"

/* Set the configuration of the field used at t=0: 
   (see set_mag_field() for details)                */
#define DIPOLE_FIELD         (0) 
#define QUADRUPOLE_FIELD     (1)
#define FIELD_CONFIG_CHOICE  (DIPOLE_FIELD)


/* Parameters for Newton-Raphson procedure used here */ 
#define MAX_NEWT_ITER 30     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-10    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER 2
#define NEWT_DIM 1 


static double beta, rdiskout=0., rhomax=0., umax=0.;
static double eta, utcov_in, lambda_in, lambda_p, f_l_in, k_eta;
static double adiab, rin, rpmax, GRscale, Omega_in, Omega_p, q;
static double l_ang_in, l_ang_p, l_ang_in_min, l_ang_in_max;
static double exp_1, exp_2, exp_3, exp_4, exp_5, alpha;
static double gtp, gtt, gpp;
static int metric_type;

static double f_l_func(double l_ang) ;

static int general_newton_raphson( double x[], int n, 
				   void (*funcd) (double [], double [], double [], 
						  double [][NEWT_DIM], double *, 
						  double *, int) );

static void resid_func(double x[], double dx[], double resid[], 
		       double jac[][NEWT_DIM], double *f, double *df, int n);
  
static void general_keplerian_functions(int prograde, double r, double *Omega_K, double *l_K, double *ucov_K, double *ucon_K); 

static void set_mag_field( void ) ; 
static void init_all_prims(int doing_height_survey, double l_factor);
static double point_interp(int n, double x[], double f[], double xpt);

double  ranc(int iseed); 


static double *height_denom, *height_numer, *r_loc;

int output_now = 0; 
double r_glob;

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

   -- Initial data for a cicumbinary disk in an evolving spacetime determined by 
       post-Newtonian approximations.  The initial data consists of a magnetized 
       disk surrounding a binary.  The disk starts at approximately twice the initial 
       separation of the binary.  The numerical grid starts outside the binary's orbit, 
       so the black holes live outside the grid.  The inner boundary condition is imposed
       to being inflow conditions.  The coordinates are spherical coordinates. 

    -- The disk model is a generalization of that in  init.kd.c  or "KD" or Keplerian disk initial data 
          ala De Villiers and Hawley 2002-2004 B-field  follows density contours, but follows more
          the solution procedure described in Chakrabarti (1985). 

***********************************************************************************/

void init(void)
{
  double ftmp;
  void init_data(void) ;

  TRACE_BEG;

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

#if( TOP_TYPE_CHOICE != TOP_SPHERICAL ) 
  this-initial-data-routine-not-setup-for-nonspherical-grids
#endif 


  TRACE_END;
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

  TRACE_BEG;

	/* Set the parameters that define the grid and other constants : */ 
	gam = (5./3.) ;
	cour = 0.8;

	t = 0. ;                     /* Initial time */ 

	/* Coordinate dependent quantities : */ 
	set_special_coord();


	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 20000. ;         /* Length of X0 dimension (evolution period) */ 
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	GridLength[1] = Rout-Rin  ;   /* Length of X1 dimension */ 
	GridLength[2] = th_length ;   /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
	//GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ;              /* Length of X1 dimension */ 
	GridLength[2] = xi_diag2[n_diag2_lines]-xi_diag2[0]  ; /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	GridLength[1] = log((coord_params->upsilon_r-coord_params->f0_r)/(1.-coord_params->f0_r)) ;              /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	//	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 

//Miguel
#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[1] = 1.  ;              /* Length of X1 dimension */
        GridLength[2] = 1.  ;                    /* Length of X2 dimension */

#endif
	GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
//Miguel
#if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[3] = 1. ;                 /* Length of X3 dimension */ 
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
	//startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	startx[1] = 0.;        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */

//Miguel
#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	startx[1] = 1.e-10;        /* Set the Physical Minimum boundary  */
        //	startx[3] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* history file frequency */
	DT_out[OUT_SURFACE] = 30. ;                	/* surface file frequency */
	DT_out[OUT_HDF5]    = 20. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = 1.;         /* radiative flux dumps       */
	DT_out[OUT_PHOTOSPHERE] = 1.;         /* radiative flux dumps       */
	DT_out[OUT_MIN_DT]  = DT_out[OUT_HISTORY];
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


#if( MAKE_PHOTOSPHERE )
	extern void  create_photosphere_grid(void);
	create_photosphere_grid();
#endif
	

	/* in case you want to test the metric and coord. transf. routines */
	//	test_geom();   

  TRACE_END;

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
  int i, j, k, l; 
  double xp[NDIM], x[NDIM];

  const int doing_height_survey = 0; 

  int n_l_factor;
  const double lf_max = 50./200.;
  //  const double lf_max = 1.;
  double l_factor, lf_min, dlf;
  struct of_coord *coords;

  TRACE_BEG;

  if( doing_height_survey )  { 
    n_l_factor = 100;
    lf_min = 40./200.;
    //    lf_min = 0.;
    dlf = ((lf_max-lf_min)/(n_l_factor-1));
  }
  else { 
    n_l_factor = 1;
    lf_min = 0.5;
    dlf = 0.;
  }


  /*******************************************************************************
     Set up utilities:
  *******************************************************************************/
  ranc(myid+1);

  ALLOC_ARRAY(height_denom,N1TOT);
  ALLOC_ARRAY(height_numer,N1TOT);
  ALLOC_ARRAY(r_loc,N1TOT);
  N1ALL_LOOP {  height_denom[i] = 0.;  }
  N1ALL_LOOP {  height_numer[i] = 0.;  }

  j = N2TOT/2;   k = N3S;
  N1ALL_LOOP {   get_coord(i,j,k,CENT,ncurr,coords) ; r_loc[i] = coords->x[RR];   }

  /*******************************************************************************
     Set up prim var arrays:
  *******************************************************************************/

  for(l=0; l<n_l_factor; l++) { 
    l_factor = lf_min + dlf*l;
      
    init_all_prims(doing_height_survey,l_factor); 

    fprintf(stdout,"l_ang_in  h_o_r(rpmax) q  eta = %26.16e %26.16e  %26.16e %26.16e \n",
	    l_ang_in, point_interp(N1TOT,r_loc,height_numer,rpmax),q,eta); 
    fflush(stdout);
    N1_LOOP { fprintf(stdout,"r height/r = %26.16e  %26.16e \n", r_loc[i], height_numer[i]); }
    fflush(stdout);
  }

  mpi_global_max(&rdiskout); 
  fprintf(stdout,"rdiskout = %28.18e \n", rdiskout) ; fflush(stdout);

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

#if( HYDRO_ONLY ) 
  LOOP { 
    p[i][j][k][B1] = p[i][j][k][B2] = p[i][j][k][B3] =  0.; 
  }
#else 
  set_mag_field(); 
#endif 

  fixup(p) ;    
  bounds(p,0) ;

  FREE(r_loc);
  FREE(height_denom);
  FREE(height_numer);

  TRACE_END;
  return;
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
  
  /* Spacetime parameters : */
  M = m_bh_tot = 1.;

  /* BBH parameters: */
#if( BBH_SPACETIME ) 
  m_bh1 = 0.5*M;                            /* Mass of left-most BH */
  m_bh2 = 0.5*M;                            /* Mass of right-most BH */
  initial_bbh_separation = 20.*m_bh_tot;  /* Initial separtion of BHs in units of total mass */

  r_horizon1 = m_bh1;             /* Radius of the horizon about BH1 in Harmonic coordinates */
  r_horizon2 = m_bh2;             /* Radius of the horizon about BH2 in Harmonic coordinates */

  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;

  a = 0.0;          /* Spin of the black hole in units of M */ 
  asq = a*a;
  r_isco    = 0.;
  r_horizon = 0.;

  Rout     = 13.*initial_bbh_separation;        /* Radial extent of the grid           */
  n_within_horizon = 0;            /* number of cells within horizon */
  //  Rin       = Rin_calc();

  Rin = 0.75 * (initial_bbh_separation) ; /*  put inner boundary at 2. times the binary orbit's radius */
  //  R0       = 0.9*Rin;
  R0 = 0.;

  /* Time at which to start shrinking the binary,  set to 0. if you want to shrink from the beginning and 
     set to an impossibly large number if you never want to shrink : */
  t_shrink_bbh = 3.e4; 
  phi_0_bbh    = 0.;

#else 
  /* Single BH Parameters: */
  a = 0.; 
  asq = a*a;
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  //  R0       = 0.9*r_horizon;  /* Offset in Radius from 1.  (GRMHD-like)       */
  Rout     = 70.;         /* Radial extent of the grid           */
  n_within_horizon = 25;   /* number of cells within horizon */
  Rin       = Rin_calc();

  r_horizon1 = r_horizon2 = r_horizon;
  m_bh1 = m_bh2 = M;
  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;
  initial_bbh_separation = 0.;

#endif

  h_slope   = 0.0729419909586;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.0* M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
  th_cutout = M_PI*SMALL;     /* Angular size to excise from the axis */
                                  /* 0.02655 gives a cutout of ~0.045Pi for h_slope = 0.35 */

  /* New way : */
//  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
//  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
//  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */

  /* Old way: */
  th_beg    = 0.;               /* beg/end set the limits of theta  and */
  th_end    = M_PI;             /*  can override th_cutout              */
  th_length = M_PI;             /* Angular size in theta of the grid    */
  
  /* COORD_DIAGONAL3 specific parameters (note h_slope is overloaded for dtheta_min/(pi/2) */
  diag3_exponent = 9;                            /* note that this should be an odd integer greater than 1  */
  diag3_factor   = 1. - 2*th_cutout/M_PI - h_slope;  /* Temporary variable */
  if( (diag3_exponent < 3) || (diag3_exponent % 2 == 0 )  ) { fprintf(stderr,"set_special_coord(): bad value for  diag3_exponent  = %d \n",diag3_exponent); fflush(stderr); exit(4112); }

  n_phi_avg = 100; 
  inv_n_phi_avg = 1./n_phi_avg; 
  d_phi_avg = 2.*M_PI * inv_n_phi_avg; 

#if( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  this-initial-data-is-not-setup-for-COORD_DIAGONAL2
#endif 

  SET_GEOM_COORD_TIME_FUNCS(t,1);

  if( myid == printer_pid ) { 
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
  }

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
   -- This is supposed to follow the routine in inits.f of the Hawley et al. code; 
       -- inconsistencies are due to not knowing what various grid functions represent 
          in inits.f  (there is very little documentation in that file); 
   -- when in doubt, I used   
             De Villiers, Hawley, Krolik  ApJ 599 (2003)

   -- The disk parameters mentioned in this paper are used below; 
   -- Uses BL coordinates until the very end and then transforms the velocities to KS; 
***********************************************************************************/
void init_all_prims(int doing_height_survey, double l_factor)
{
  int   i,j,k,l, prograde;
  int   j_glob_eq, j_glob, jbeg, jend;
  int   retval;

  double rho, uu;
  double *pr;

  double *lambda_eq, *l_eq, *h_eq;
  double xsol[1];

  struct of_geom geom_eq;
  struct of_geom *geom;
  struct of_coord *coords;
  struct of_coord coords_loc;

  double lambda,  h, f_l, l_ang, epsilon;
  double ucon[NDIM], ucov[NDIM], ucon_bl[NDIM];
  double utcov, upcov, utcon, upcon, Omega; 


  /***********************************************************************************************
    BEGIN:  Code block to integrate the equations from the equator over theta/X2: 
  **********************************************************************************************/  
#define X2_INTEGRATION_CODE_BLOCK                                                                            \
  {                                                                                                          \
    j = N2S + j_glob - jbeg ;                                                                                \
    if( (N2S <= j) && (j <= N2E) )  {                                                                        \
      get_coord(i,j,k,CENT,ncurr,coords);                                                                    \
      get_special_geometry( coords, &geom_eq, metric_type);                                                  \
      gtt = geom_eq.gcon[TT][TT];   gtp = geom_eq.gcon[TT][PH];   gpp = geom_eq.gcon[PH][PH];                \
      xsol[0] = lambda;                                                                                      \
      retval = general_newton_raphson( xsol, 1, resid_func);                                                 \
      lambda = xsol[0];                                                                                      \
      if( retval ) {                                                                                         \
	fprintf(stderr,"init_gen_axi_disk.c: init_all_prims():  ");                                          \
	fprintf(stderr,"Failed to find a solution at [i=%6d, j=%6d, j_glob=%6d] [r=%26.16e] [th=%26.16e]\n", \
                i, j, j_glob, coords->x[RR], coords->x[TH]);                                                 \
        fflush(stdout);                                                                                      \
      }                                                                                                      \
      else{                                                                                                  \
        l_ang  = eta * pow( lambda , (2.-q) );                                                               \
        f_l = f_l_func(l_ang);                                                                               \
        utcov = -1./sqrt(-gtt + l_ang*(2*gtp - l_ang*gpp));                                                  \
        upcov = -l_ang * utcov;                                                                              \
        utcon = gtt*utcov + gtp*upcov;                                                                       \
        upcon = gtp*utcov + gpp*upcov;                                                                       \
        h = utcov_in * f_l_in / ( utcov * f_l );                                                             \
        if( h > 1. ) {                                                                                       \
          epsilon = (h - 1.) / gam;                                                                          \
          rho =  pow( (epsilon * (gam - 1.) / adiab), (1./(gam-1.)) );                                       \
          uu  = rho * epsilon;                                                                               \
          pr = p[i][j][k];                                                                                   \
          pr[RHO] = rho;                                                                                     \
          pr[UU ] = uu;                                                                                      \
          pr[U1 ] = utcon;   /* store the phi-independent components here for now  */                        \
          pr[U3 ] = upcon;   /* to be transformed in a phi-dependent way later     */                        \
          get_geometry(i,j,k,CENT,ncurr,geom); /* get (x1,x2,x3) metric,   */                                \
          rho = geom->g * pr[RHO] ;                                                                          \
          height_numer[i] += sqrt(geom_eq.gcov[2][2]) * fabs(coords->x[2]-0.5*M_PI) * rho;                   \
          height_denom[i] += rho;                                                                            \
        }                                                                                                    \
      }                                                                                                      \
    }                                                                                                        \
  }                                                                                                          

  /***********************************************************************************************
    END:  Code block
  **********************************************************************************************/  


  /***********************************************************************************************
    Set parameters used for all points : 
  **********************************************************************************************/  
  adiab    = 1.e-2;    /* constant "K" in P = K rho**gam , assume this EOS in disk */
#if( BBH_SPACETIME ) 
  rin      = 1.5*(initial_bbh_separation+m_bh_tot); /* Radius of inner edge of the disk in BL coordinates */
  rpmax    = 2.5*(initial_bbh_separation+m_bh_tot); /* Radius of the pressure maximum in BL coordinates */
#else
  rin      = 20.;     /* Radius of inner edge of the disk ?? */
  rpmax    = 35.;     /* Radius of the pressure maximum ?? */
#endif
  prograde = 1;       /* non-zero for prograde disks */ 
  //  GRscale  = 1.e-2;    /* Relative magnitude of random perturbations to u */
  GRscale  = 0.;    /* Relative magnitude of random perturbations to u */
  beta     = 1.e2 ;   /* Plasma beta parameter */

#if( (METRIC_TYPE_CHOICE==METRIC_GENERAL_STATIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG2)&&(METRIC_DYNAMIC_TYPE_CHOICE!=METRIC_DYNAMIC_KS_SPHERICAL) )
  metric_type = METRIC_GENERAL_PHI_AVG;
#else
  metric_type = METRIC_BL_SPHERICAL;
#endif

#if(USE_SIMPLE_EOS) 
  N1ALL_LOOP { Katm[i] = adiab; }
#endif

  general_keplerian_functions(prograde, rin  , &Omega_in, &l_ang_in_min, ucov, ucon);
  general_keplerian_functions(prograde, rpmax, &Omega_p , &l_ang_p , ucov, ucon);
  lambda_p  = sqrt( l_ang_p  / Omega_p  );

  l_ang_in_max = l_ang_p;

  //  if( (l_factor < 0.)  ) { l_factor = 0.5; } 

  /* for when you want to search for desired h/r : */
  //  l_ang_in = l_factor*l_ang_in_max + (1.-l_factor)*l_ang_in_min; 

  /* for a single BH  h/r=0.1 disk with rin=20, rpmax=35  :  */
  //  l_ang_in = 5.2196077516713117;

  /* for a BBH  h/r=0.1 disk with rin=3*a, rpmax=5*a, a=20M  :  */
  l_ang_in = 8.7431737221816430e0;

  coord_of_r(rin, &coords_loc);
  get_special_geometry( &coords_loc, &geom_eq, metric_type);
  gtt = geom_eq.gcon[TT][TT];
  gtp = geom_eq.gcon[TT][PH];
  gpp = geom_eq.gcon[PH][PH];
  utcov_in = -1./sqrt(-gtt + l_ang_in*(2*gtp - l_ang_in*gpp));
  Omega_in = (gtp - gpp*l_ang_in)/(gtt - gtp*l_ang_in); 
  lambda_in = sqrt( l_ang_in / Omega_in );
  q = 2. - log(l_ang_in/l_ang_p) / log(lambda_in/lambda_p)  ; 
  eta = l_ang_p * pow( lambda_p , q-2.); 
  k_eta = pow( eta , (-2./(q-2.)) );

  alpha = q/(q-2.);
  exp_1=(2.*(q-1.))/q;
  exp_4 = alpha + 1.;   // 2*(q-1)/(q-2)
  exp_2 = 1./exp_4;      // (q-2)/(2*(q-1))
  exp_3 = 1./alpha;     // (q-2)/q
  exp_5 = 1./(2.-q);
  f_l_in = f_l_func(l_ang_in);

  if( fabs(alpha-1.) < 1.e-3 ) { 
    fprintf(stdout,"init_prim(): alpha may be getting too close to one :   %26.16e \n", alpha); 
    fail(FAIL_BASIC,0); 
  }

  if( !doing_height_survey ) { 
    fprintf(stdout,"\n##################################################\n");
    fprintf(stdout,"  AXI-SYM TORUS PARAMETERS \n------------------------------------\n");
    fprintf(stdout,"\t  adiab       =  %28.18e \n",adiab    ); 
    fprintf(stdout,"\t  rin         =  %28.18e \n",rin      ); 
    fprintf(stdout,"\t  rpmax       =  %28.18e \n",rpmax    ); 
    fprintf(stdout,"\t  GRscale     =  %28.18e \n",GRscale  ); 
    fprintf(stdout,"\t  beta        =  %28.18e \n",beta     ); 
    fprintf(stdout,"\t  q           =  %28.18e \n",q        ); 
    fprintf(stdout,"\t  alpha (exp) =  %28.18e \n",alpha    ); 
    fprintf(stdout,"\t  eta         =  %28.18e \n",eta      ); 
    fprintf(stdout,"\t  k_eta       =  %28.18e \n",k_eta    ); 
    fprintf(stdout,"\t  metric_type =  %d      \n",metric_type); 
    fprintf(stdout,"  Values at rin:\n------------------------------------\n");
    fprintf(stdout,"\t  l_ang_in    =  %28.18e \n",l_ang_in ); 
    fprintf(stdout,"\t  l_ang_in_min=  %28.18e \n",l_ang_in_min ); 
    fprintf(stdout,"\t  lambda_in   =  %28.18e \n",lambda_in); 
    fprintf(stdout,"\t  f_l_in      =  %28.18e \n",f_l_in   ); 
    fprintf(stdout,"\t  utcov_in    =  %28.18e \n",utcov_in ); 
    fprintf(stdout,"\t  Omega_in    =  %28.18e \n",Omega_in ); 
    fprintf(stdout,"  Values at rpmax:\n------------------------------------\n");
    fprintf(stdout,"\t  l_ang_p     =  %28.18e \n",l_ang_p  ); 
    fprintf(stdout,"\t  lambda_p    =  %28.18e \n",lambda_p ); 
    fprintf(stdout,"\t  utcov       =  %28.18e \n",ucov[TT] ); 
    fprintf(stdout,"\t  Omega_p     =  %28.18e \n",Omega_p  ); 
    fprintf(stdout,"##################################################\n");
  }

  /**********************************************************************************************
    Find the values in the equator:
  **********************************************************************************************/  
  k = N3S;   j = totalsize[2]/2;

  ALLOC_ARRAY( lambda_eq, N1TOT); 
  ALLOC_ARRAY(      l_eq, N1TOT); 
  ALLOC_ARRAY(      h_eq, N1TOT); 

  N1ALL_LOOP { lambda_eq[i] = -1. ;  } 
  N1ALL_LOOP {      l_eq[i] = -1. ;  } 
  N1ALL_LOOP {      h_eq[i] = -1. ;  } 
      
  N1_LOOP { 
    coord_of_r(r_loc[i], &coords_loc);

    general_keplerian_functions(prograde, coords_loc.x[RR]  , &Omega, &l_ang, ucov, ucon);

    fprintf(stdout,"Keplerian functions = %26.16e %26.16e %26.16e %26.16e %26.16e \n",coords_loc.x[RR], Omega, l_ang, ucov[TT], ucon[TT]); fflush(stdout); 

    if( coords_loc.x[RR] >= rin ) { 
      if( (i <= N1S)  || (lambda_eq[i-1] < 0.) ) { 
	l_eq[i] = l_ang_in; 
	lambda_eq[i] = pow(  (l_eq[i]/eta), exp_5 );
      } 
      else { 
	lambda_eq[i] = lambda_eq[i-1]; 
      }

      get_special_geometry( &coords_loc, &geom_eq, metric_type);
      gtt = geom_eq.gcon[TT][TT];
      gtp = geom_eq.gcon[TT][PH];
      gpp = geom_eq.gcon[PH][PH];

      xsol[0] = lambda_eq[i];
      retval = general_newton_raphson( xsol, 1, resid_func);

      if( retval ) { 
	fprintf(stderr,"init_gen_axi_disk.c: init_all_prims():  Failed to find a solution at [i=%6d]   [r=%26.16e]  \n", i, coords_loc.x[RR]); 
	fflush(stdout); 
      }
      else{ 
	lambda_eq[i] = xsol[0];
	l_ang  = eta * pow( lambda_eq[i] , (2.-q) ); 
	l_eq[i] = l_ang;
	f_l = f_l_func(l_ang);
	utcov = -1./sqrt(-gtt + l_ang*(2*gtp - l_ang*gpp));
	h_eq[i] = utcov_in * f_l_in / ( utcov * f_l ); 
//	fprintf(stdout,"ut_in f_in l_in lambda_in = %20.10e %20.10e %20.10e %20.10e\n", 
//		utcov_in, f_l_in, l_ang_in, lambda_in); 
	Omega = (gtp - gpp*l_ang)/(gtt - gtp*l_ang); 
//	fprintf(stdout,"r ut f l lambda h Omega = %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", 
//		coords_loc.x[RR],utcov, f_l, l_ang, lambda_eq[i], h_eq[i], Omega); 
//	fflush(stdout);

      }

      if( (h_eq[i] > 1.) && (coords_loc.x[RR] > rdiskout) ) { rdiskout = coords_loc.x[RR] ; }
      
      fprintf(stdout,"equator-values: myid i r lambda l h = %6d %6d %26.16e %26.16e %26.16e %26.16e \n",
	      myid, i,coords_loc.x[RR],lambda_eq[i],l_eq[i],h_eq[i]); 
      fflush(stdout);


    }   /*   if( coords_loc.x[RR] >= rin )   */
  }

  
  /**********************************************************************************************
    Now we set the phi-independent data; 
    Set the primitives to zeros as default values to be filled in by floor
  **********************************************************************************************/
  k = N3S; 
  N1_LOOP N2_LOOP PLOOP  {  p[i][j][k][l] = 0.; }
  N1_LOOP { height_denom[i] = 0. ; } 
  N1_LOOP { height_numer[i] = 0. ; } 

  /**********************************************************************************************
    Set the values in the disk :
  **********************************************************************************************/  
  j_glob_eq = totalsize[2]/2;
  get_global_ijk(0,0,0,&i,&jbeg,&k);
  get_global_ijk(N1E-N1S,N2E-N2S,N3E-N3S,&i,&jend,&k);
  k = N3S;

  N1_LOOP { 
    h  = h_eq[i] ;     
    if( h > 1. ) { 
      lambda = lambda_eq[i] ; 
      
      if( jend >= j_glob_eq ) { 
	for( j_glob = j_glob_eq; j_glob < totalsize[2]; j_glob++ ) { 
	  X2_INTEGRATION_CODE_BLOCK;
	}
      }

      lambda = lambda_eq[i] ; 

      if( jbeg <= j_glob_eq ) { 
	for( j_glob = j_glob_eq-1; j_glob >= 0; j_glob-- ) { 
	  X2_INTEGRATION_CODE_BLOCK;
	}
      }
    }
  }

  N1_LOOP { height_numer[i] /= height_denom[i]*r_loc[i]; } 
  N1_LOOP { if( !finite(height_numer[i]) ) { height_numer[i] = 0.; } }
  

  /**********************************************************************************************
    Now copy the data over all phi, transforming to numerical coordinatea and adding random 
     perturbations along the way:
  **********************************************************************************************/  
  if( doing_height_survey ) { return; } 
  N1_LOOP N2_LOOP { 
    k = N3S; 
    pr = p[i][j][k];
    rho         = pr[RHO]; 
    uu          = pr[UU ]; 
    ucon_bl[TT] = pr[U1 ];
    ucon_bl[RR] = 0.;
    ucon_bl[TH] = 0.;
    ucon_bl[PH] = pr[U3 ];
   

    if( rho > 0. ) { 
      fprintf(stdout,"KATM = %26.16e  %26.16e  %26.16e \n", rho,uu,((gam - 1.) * uu / pow( rho, gam ))) ; fflush(stdout);

      if( rho > rhomax )  { rhomax = rho ; } 
      if(  uu >   umax )  {   umax =  uu ; } 

      N3_LOOP { 

	pr = p[i][j][k];
	get_coord(i,j,k,CENT,ncurr,coords);	
	get_geometry(i,j,k,CENT,ncurr,geom); // get (x1,x2,x3) metric, from here on it's only (x1-3)

	/* Add perturbation to internal energy density : */
	pr[RHO] = rho;
	pr[UU]  = uu * (1. + GRscale * (ranc(0) - 0.5));

	/* Transform 4-velocity to numerical coordinates and into primitive variable velocity */
#if( (METRIC_TYPE_CHOICE==METRIC_KS_SPHERICAL) || \
     ((METRIC_TYPE_CHOICE==METRIC_GENERAL_STATIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC)&&\
      (METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_KS_SPHERICAL) ) )
	bl_to_ks_con(coords->x,ucon_bl,ucon);  // transform from BL to KS coordinates 
#else
	for(l=0;l<NDIM;l++) { ucon[l] = ucon_bl[l]; }
#endif
	transform_rank1con2(coords->dxp_dx,ucon); // transform from (r,th,ph) to (x1,x2,x3)
	ucon2pr( pr, ucon, geom->gcon );  // get prim. velocities;
      }      
    }
  }

  FREE( lambda_eq); 
  FREE(      l_eq); 
  FREE(      h_eq); 
    
#undef X2_INTEGRATION_CODE_BLOCK

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
static void set_mag_field( void )
{
  int i,j,k,l;

  double rho_avg, rho_cut, Aphi, beta_act, norm, bsq_int, u_int ; 
  struct of_geom *geom;  
  struct of_coord *coords;  

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

  N1_LOOP  N2_LOOP   N3_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); 

    /* flux-ct */
    p[i][j][k][B1] = -( F[i  ][j  ][0][0][0] 
		      - F[i  ][j+1][0][0][0] 
		      + F[i+1][j  ][0][0][0] 
                      - F[i+1][j+1][0][0][0] 
		       )/(2.*dx[2]*geom->g) ;

    p[i][j][k][B2] = (  F[i  ][j  ][0][0][0] 
		      + F[i  ][j+1][0][0][0] 
		      - F[i+1][j  ][0][0][0] 
		      - F[i+1][j+1][0][0][0] 
			)/(2.*dx[1]*geom->g) ;
    
    p[i][j][k][B3] = 0. ;

    bsq_int += geom->g * bsq_calc(p[i][j][k],geom) ;
    u_int   += geom->g * p[i][j][k][UU]; 
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
  
  N1_LOOP  N2_LOOP N3_LOOP { 
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
  }

  bsq_int = 0.;
  N1_LOOP  N2_LOOP N3_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); 
    bsq_int += geom->g * bsq_calc(p[i][j][k],geom) ;
  }

  mpi_global_sum( &bsq_int );
  fprintf(stdout,"final bsq_int = %28.18e \n", bsq_int); fflush(stdout); 

  beta_act = (gam - 1.)*u_int/(0.5*bsq_int) ;
  fprintf(stdout,"final beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);
  

  return;

}

/***********************************************************************************/
/***********************************************************************************
  general_keplerian_functions():
  ------------------
   -- calculates the velocity and associated functions for a circular equatorial orbit 
      in a general axi-symmetric (of BL form) spacetime; 
***********************************************************************************/
static void general_keplerian_functions(int prograde, double r, double *Omega_K, double *l_K, double *ucov_K, double *ucon_K)
{
  
#define DELR  (1.e-5)   /* Resolution in relative to r to perform the derivatives */  

  int i,j,k,l;
  double dg_tt, dg_tp, dg_pp, dr, inv_dr, Omega;
  struct of_geom  geom_m, geom_p;
  struct of_coord  coords_loc, coords_m, coords_p;

  /**********************************************************************************  
     First, calculate the metric derivatives for the orbital frequency:
  **********************************************************************************/
  dr = DELR * r; 

  coord_of_r(r   , &coords_loc);
  coord_of_r(r+dr, &coords_p  );
  coord_of_r(r-dr, &coords_m  );

  inv_dr = 0.5/dr;

  get_special_geometry( &coords_m, &geom_m, metric_type);
  get_special_geometry( &coords_p, &geom_p, metric_type);

  dg_tt = (geom_p.gcov[TT][TT] - geom_m.gcov[TT][TT])  * inv_dr;
  dg_tp = (geom_p.gcov[TT][PH] - geom_m.gcov[TT][PH])  * inv_dr;
  dg_pp = (geom_p.gcov[PH][PH] - geom_m.gcov[PH][PH])  * inv_dr;

  if( prograde ) { 
    Omega =  - ( dg_tp - sqrt( dg_tp * dg_tp  - dg_pp * dg_tt ) ) / dg_pp ; 
  }
  else { 
    Omega =  - ( dg_tp + sqrt( dg_tp * dg_tp  - dg_pp * dg_tt ) ) / dg_pp ; 
  }

  //  fprintf(stdout,"r dg_tt dg_tp dg_pp = %26.16e  %26.16e  %26.16e  %26.16e \n", r, dg_tt, dg_tp, dg_pp);

  if( !finite( Omega ) ) { 
    fprintf(stderr,"general_keplerian_functions(): Omega_K is not finite :   %26.16e   at   r= %26.16e \n", (Omega), r); 
    fflush(stderr); 
    DLOOP1 {  ucov_K[i] = 0.; }
    DLOOP1 {  ucon_K[i] = 0.; }
    *l_K = *Omega_K = 0.;
    return;
  }

  *Omega_K = Omega;

  get_special_geometry( &coords_loc, &geom_p, metric_type);

  ucon_K[TT] = sqrt(-1. / (geom_p.gcov[TT][TT] + (Omega)*(2.*geom_p.gcov[TT][PH] + (Omega)*geom_p.gcov[PH][PH]) ) ) ;
  ucon_K[RR] = 0.;
  ucon_K[TH] = 0.;
  ucon_K[PH] = (Omega)*ucon_K[TT];

  lower(ucon_K,&geom_p, ucov_K); 

  *l_K = - ucov_K[PH] / ucov_K[TT];

  return;

#undef DELR

}

/***********************************************************************************/
/***********************************************************************************
  f_l_func():
  ------------------
   -- function used to find enthalpy solution;
***********************************************************************************/
static double f_l_func(double l_ang) 
{
  return( 
	 pow( fabs(  1. - k_eta * pow(l_ang, exp_4) ) , exp_2 ) 
  ) ; 
}


/**********************************************************************/
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

n       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( double x[], int n, 
				   void (*funcd) (double [], double [], double [], 
						  double [][NEWT_DIM], double *, 
						  double *, int) )
{
  double f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  double errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;

  int   keep_iterating, i_increase;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  //-fast  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;

    x_old[0] = x[0] ;

    /* don't use line search : */
    x[0] += dx[0]  ;

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = fabs(x[0]);

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (!finite(f)) || (!finite(df)) || (!finite(x[0]))  ) {
#if(LTRACE)
    fprintf(stderr,"\ngnr not finite, f,df,x,x_o = %26.20e %26.20e %26.20e %26.20e \n", f,df,x[0],x_old[0]); fflush(stderr); 
#endif
    return(2);
  }

  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  } 

  return(0);

}

/**********************************************************************/
/*********************************************************************************
   resid_func():

        -- residual/jacobian routine to calculate lambda

        resid =  g^{t phi} l^2 - g^{tt} l  + g^{t phi} lambda^2  - g^{phi phi} l lambda^2 = 0 

        resid =  g^{t phi} eta^2 lambda^{4 - 2*q} - g^{tt} eta lambda^{2-q} 
                    + g^{t phi} lambda^2  - g^{phi phi} eta lambda^{4-q} = 0 


     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void resid_func(double x[], double dx[], double resid[], 
			 double jac[][NEWT_DIM], double *f, double *df, int n)
{
  double  t1 ,  t100,  t9 ,  t3 ,  t4 ,  t7 ,  t23,  t24,  t32,  t33,  t15,  t16, lambda;

  
  lambda = x[0];

  t1 = eta*eta;
  t100 = pow(lambda,-q);
  t9 = lambda*lambda;
  //      t3 = pow(lambda,2.0-q);
  t3 = t9*t100;  
  t4 = t3*t3;    // pow(lambda,4-2*q);
  t7 = t3*eta;
  resid[0] = gtp*t4*t1-gtt*t7+(gtp-gpp*t7)*t9;
  //      t23 = pow(lambda,1.0-q);
  t23 = t100*lambda;
  t24 = eta*t23;
  //      t32 = pow(lambda,3.0-q);
  t32 = t23 * t9;
  t33 = t32*eta;
  //      t15 = pow(lambda,3.0-2.0*q);
  t15 = t32 * t100;
  t16 = t1*t15;

  jac[0][0] = 4.0*t16*gtp-2.0*t16*gtp*q-2.0*t24*gtt+t24*gtt*q+2.0*
    lambda*gtp-4.0*t33*gpp+t33*gpp*q;

  dx[0] = -resid[0]/jac[0][0];
  *df = - resid[0]*resid[0];
  *f = -0.5*(*df);


//  if( output_now ) { 
//    fprintf(stdout,"r_eq lambda l resid = %26.16e  %26.16e  %26.16e  %26.16e  \n", r_glob, lambda, (eta * pow( lambda, (2.-q))), resid[0] ); 
//    fflush(stdout); 
//  }

  return;

}


#undef MAX_NEWT_ITER 
#undef NEWT_TOL   
#undef MIN_NEWT_TOL
#undef EXTRA_NEWT_ITER 
#undef NEWT_DIM

/*********************************************************************************/
/*********************************************************************************
   point_interp():
  ------------------

    -- returns with the linear interpolated value of function f[] over 
        at point xpt;
    -- assumes x[] is monotonic; 
*********************************************************************************/
static double point_interp(int n, double x[], double f[], double xpt)
{
  int i, imin ;
  double fmin;
  
  /* find imin : */
  /* If point lies beyond x[]'s domain, then we must extrapolate: */
  if( xpt > x[n-1] ) { 
    imin = n-2;
  }
  else { 
    for(i=0;i<n;i++) { 
      if( xpt < x[i] ) { 
	imin = i-1;  
	break; 
      }
    }
    if( imin < 0   ) { imin = 0;   }
    if( imin > n-2 ) { imin = n-2; }
  }

  fmin = (xpt-x[imin])/(x[imin+1]-x[imin]); 
  return(  f[imin+1]*fmin + f[imin]*(1.-fmin) );

}
