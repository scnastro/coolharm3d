
#include "decs.h"
#include "metric.h"
#if( USEMPI )
#include "mpi.h"
#endif

/* Define whether to set Mini-Disks ID */
#define INITIALIZE_MINI_DISKS (1)
#if( INITIALIZE_MINI_DISKS )
#define BHS_ON_GRID (1)
#define NO_CIRC (0)
#define SUPERIMPOSE_INTERP_DATA (0)
#if( (SUPERIMPOSE_INTERP_DATA && !NO_CIRC) || RUN_FROM_INTERP_DATA )
  --Must-choose-between-circ.-initial-data-or-interpolated-data;
#endif
double ****pd1, ****pd2; // primitive-like arrays for disk1/2
static double rdisk_max, rdisk_min;
static double r1min,r2min,r1max,r2max;
static double beta_disk1, beta_disk2;
/************************************************************************************************
    Special size of arrays for disk solutions in memory space and loop aliases
**********************************************************************************************/
#define N1TOT_MINI (200)
#define N2TOT_MINI (200)
#define N3TOT_MINI (1)
#define NCELLS_MINI ( N1TOT_MINI * N2TOT_MINI * N3TOT_MINI )
#define NFUNCS_MINI (12) // Number of functions to store 
/* Mnemonics for pd1 pd2 */
#define RHOM (0) //mini-disk array density
#define UUM  (1) //mini-disk array internal energy
#define UTTM (2) //mini-disk array u^t
#define URRM (3) //mini-disk array u^r
#define UTHM (4) //mini-disk array u^th
#define UPHM (5) //mini-disk array u^ph
#define ATTM (6) //mini-disk array A_t
#define ARRM (7) //mini-disk array A_r
#define ATHM (8) //mini-disk array A_th
#define APHM (9) //mini-disk array A_ph
#define XRRM (10)//mini-disk array r-coord
#define XTHM (11)//mini-disk array th-coord

#define LOOP_MINI   for(i=0;i<N1TOT_MINI;i++) for(j=0;j<N2TOT_MINI;j++) for(k=0;k<N3TOT_MINI;k++)
#define N1LOOP_MINI for(i=0;i<N1TOT_MINI;i++)
#define N2LOOP_MINI for(j=0;j<N2TOT_MINI;j++)
#define N3LOOP_MINI for(k=0;k<N3TOT_MINI;k++)
#define PLOOP_MINI  for(l=0;l<NFUNCS_MINI;l++)

# if( !(BHS_ON_GRID) )
  --dont-compile-need-bhog-with-mini-disks
# endif
#endif

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

static int first_call_disk1=1, first_call_disk2=1;
static int first_lookup1=1, first_lookup2=1;

static double beta, rdiskout=0., rhomax=0., umax=0.;
static double eta, utcov_in, lambda_in, lambda_p, f_l_in, k_eta;
static double adiab, rin, rpmax, GRscale, Omega_in, Omega_p, q;
static double l_ang_in, l_ang_p, l_ang_in_min, l_ang_in_max;
static double exp_1, exp_2, exp_3, exp_4, exp_5, alpha;
static double gtp, gtt, gpp;
static int metric_type, prograde;
static double h_o_r;
static double r_glob;
static int n_in; 

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
static void calc_eq_values(double *lambda_eq, double *l_eq, double *h_eq, struct of_coord *coords_r);

double  ranc(int iseed); 


double *height_denom, *height_numer, *r_beg, *r_end, *height_numer_bins, *height_denom_bins;


int output_now = 0; 


/* Special Mini-Disk functions */
void init_mini_disks(int i, int j, int k, struct of_coord *coords, double *pr);
void create_disk_solution(int disk, int doing_height_survey, double l_factor);
static void calc_eq_values_mini(double *lambda_eq, double *l_eq, double *h_eq, struct of_coord *coords_r);
void lookup_mini_disk( int disk, double xBL[NDIM],double *pr, double ucon_bl[NDIM], double Acov_bl[NDIM] );
void renormalize_ucon( double ucon[NDIM], double gcov[NDIM][NDIM] );
static void set_vect_potential_mini( double ****prl );
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
       to being inflow conditions.  

    -- The coordinate are NOT assumed to be 1-to-1 related to numerical coordinate systems, i.e. the coordinates
        can be warped spherical or completely arbitrary; 
        --this is the key difference between this routine and that found in init.gen_axi_disk.c 

    -- The disk model is a generalization of that in  init.kd.c  or "KD" or Keplerian disk initial data 
          ala De Villiers and Hawley 2002-2004 B-field  follows density contours, but follows more
          the solution procedure described in Chakrabarti (1985). 

    -- this version does not assume that "physical" coordinates (r,th,phi) are independent; 
        in other words, this routine treats  r=r(i,j,k), th=th(i,j,k), phi=phi(i,j,k);

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

#if( COORD_RADIUS_DIM == 1 )
//--HERE  fprintf(stdout,"init.gen_disk.c:  we should not be here.... use init.gen_axi_disk.c instead \n"); fflush(stdout); fail(FAIL_BASIC,0);
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
	cour = 0.45;

	t = 0. ;                     /* Initial time */ 

	/* Coordinate dependent quantities : */ 
	set_special_coord();

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 7.e3 ;         /* Length of X0 dimension (evolution period) */ 
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
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	GridLength[1] = log((coord_params->upsilon_r-coord_params->f0_r)/(1.-coord_params->f0_r)) ;              /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	//	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[1] = 1.  ;              /* Length of X1 dimension */
        GridLength[2] = 1.  ;                    /* Length of X2 dimension */

#endif
	GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
#if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[3] = 1. ;                 /* Length of X3 dimension */ 
#endif

	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];
#if(EQUATORIAL_RUN)
 	dx[2] = 1./160.;
#endif

	/**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
	**************************************************************************/
	startx[0] =  0.;                /* Set the Physical Minimum boundary  */
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	startx[1] = Rin;                /* Set the Physical Minimum boundary  */
	startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5*GridLength[2]-0.5*dx[2] ; 
# else 
 	startx[2] = th_beg      ;
# endif 

#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
        //	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5 - 0.5*dx[2] ;
# else 
 	startx[2] = 0.            ; 
# endif
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	startx[1] = 0.;        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */

#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	startx[1] = 1.e-10;        /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5-0.5*dx[2] ; 
# else 
 	startx[2] = 0.            ; 
# endif 
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 10. ;                	/* history file frequency */
	DT_out[OUT_SURFACE] = 30. ;                	/* surface file frequency */
	DT_out[OUT_HDF5]    = 10. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 4.*DT_out[OUT_HDF5];      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_STAT];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	DT_out[OUT_MIN_DT]  = DT_out[OUT_HISTORY];
	n_restart = 10000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	dx[0] = 1.e-2;
	calc_all_geom() ;  
	dt_global_min = find_min_dt();

	//	dx[0] = cour * dt_global_min; /* Discretization size in X0 direction, time step*/

	/**************************************************************************
	  Print out basic grid parameters: 
	**************************************************************************/
	if( myid == printer_pid ) {  
	  fprintf(stdout,"Length: %10.4e %10.4e %10.4e %10.4e\n", 
		  GridLength[0], GridLength[1], GridLength[2], GridLength[3]);
	  fprintf(stdout,"dx    : %10.4e %10.4e %10.4e %10.4e\n", 
		  dx[0], dx[1], dx[2], dx[3]);
	  fprintf(stdout,"startx: %10.4e %10.4e %10.4e %10.4e\n", 
		  startx[0], startx[1], startx[2], startx[3]);

	  fprintf(stdout,"rhomin,minlimt = %28.18e %28.18e \n", RHOMIN, RHOMINLIMIT);
	  fprintf(stdout,"uumin ,minlimt = %28.18e %28.18e \n", UUMIN,  UUMINLIMIT);
	
	  fflush(stdout);
	}

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

  double l_factor,reldiff;

  struct of_coord *coords;

  double lf_min = 0.;   /* use maximum range to start bisection search */
  double lf_max = 1.;

  const int doing_height_survey = 0; 
  const double tol              = 1e-13;                    /* tolerance for the height survey bisection search */
  const double h_o_r_target     = 1.0044771758781224e-01;   /* set this to be your desired scaleheight when "doing_height_survey" */
  const int max_iterations = 200;    /* maximum number of iterations allowed in bisection before giving up */


  TRACE_BEG;

  /*******************************************************************************
     Set up utilities:
  *******************************************************************************/
  ranc(myid+1);

  /*******************************************************************************
     Set up prim var arrays:
  *******************************************************************************/
  if( doing_height_survey )  { 
    i=0;
    do { 
      l_factor = 0.5*(lf_min+lf_max);

      init_all_prims(doing_height_survey,l_factor); 

      if( h_o_r <= h_o_r_target ) { lf_min = l_factor;  }  
      else                        { lf_max = l_factor;  }

      reldiff = fabs(lf_max-lf_min)/(l_factor);
      if( myid == printer_pid ) { 
	fprintf(stdout,"bisection  %10d %26.16e  %26.16e %26.16e %26.16e  %26.16e %26.16e %26.16e \n", i,l_ang_in, lf_min,lf_max,l_factor,reldiff,h_o_r,h_o_r_target); fflush(stdout);
      }
      i++;
    } while( (reldiff > tol) || (i > max_iterations) ); 
  }
  else { 
    init_all_prims(doing_height_survey,0.5); 
  }

  mpi_global_max(&rdiskout); 
  if( myid == printer_pid ) {  fprintf(stdout,"rdiskout = %28.18e \n", rdiskout) ; fflush(stdout); }


  /*******************************************************************************
     Set up mini-disks:
  *******************************************************************************/
  //DEBUG_OUT("herehere beg");
#if( INITIALIZE_MINI_DISKS )
  LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    init_mini_disks(i,j,k,coords,p[i][j][k]);
  }
#endif
#if( USEMPI )
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  //DEBUG_OUT("herehere end");

  /*******************************************************************************
    Normalize the densities: 
  *******************************************************************************/
  mpi_global_max(&rhomax); 
  mpi_global_max(&umax); 
  if( myid == printer_pid ) { 
    fprintf(stdout, "init_data(): orig rhomax   = %28.18e \n", rhomax); 
    fprintf(stdout, "init_data(): orig   umax   = %28.18e \n",   umax);   fflush(stdout); 
  }
  
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


#if( SUPERIMPOSE_INTERP_DATA )
/*******************************************************************************
    if reading from interpolated data we have to transform vector components to
    the numerical coordinates
 *******************************************************************************/

  double *rho_tmp, *uu_tmp ;
  double *v1_tmp, *v2_tmp, *v3_tmp;
  double *B1_tmp, *B2_tmp, *B3_tmp;

  struct of_coord *coords;
  double Bcon[NDIM], vcon[NDIM];
  char hdf_name[200];
  FILE *fp;
  
  /* Set the name of the restart file : */
  sprintf(hdf_name,"%s/rdump_start_interp_data.h5", DIR_out[OUT_HDF5]);

  /* Open existing file : wrapper handles error checking */ 
  file_id = myH5_Fopen(hdf_name);


  ALLOC_ARRAY(rho_tmp, NCELLS);
  ALLOC_ARRAY(uu_tmp , NCELLS);
  ALLOC_ARRAY(v1_tmp , NCELLS);
  ALLOC_ARRAY(v2_tmp , NCELLS);
  ALLOC_ARRAY(v3_tmp , NCELLS);
  ALLOC_ARRAY(B1_tmp , NCELLS);
  ALLOC_ARRAY(B2_tmp , NCELLS);
  ALLOC_ARRAY(B3_tmp , NCELLS);

  myH5_read_gfunc(file_id, "/rho_dest", rho_tmp ) ;
  myH5_read_gfunc(file_id, "/uu_dest" , uu_tmp ) ;
  myH5_read_gfunc(file_id, "/v1_dest" , v1_tmp ) ;
  myH5_read_gfunc(file_id, "/v2_dest" , v2_tmp ) ;
  myH5_read_gfunc(file_id, "/v3_dest" , v3_tmp ) ;
  myH5_read_gfunc(file_id, "/B1_dest" , B1_tmp ) ;
  myH5_read_gfunc(file_id, "/B2_dest" , B2_tmp ) ;
  myH5_read_gfunc(file_id, "/B3_dest" , B3_tmp ) ;

  ind = 0 ;
  LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);
    if( coords->r > initial_bbh_separation ) {
      
      p[i][j][k][RHO] = rho_tmp[ind]  ;
      p[i][j][k][UU]  =  uu_tmp[ind]  ;
      
      vcon[0] = 0 ;
      vcon[1] = v1_tmp[ind] ;
      vcon[2] = v2_tmp[ind] ;
      vcon[3] = v3_tmp[ind] ;
      Bcon[0] = 0 ;
      Bcon[1] = B1_tmp[ind] ;
      Bcon[2] = B2_tmp[ind] ;
      Bcon[3] = B3_tmp[ind] ;
      
      transform_rank1con2(coords->dxp_dx,vcon);
      transform_rank1con2(coords->dxp_dx,Bcon);
      
      p[i][j][k][U1]  =  vcon[1]  ;
      p[i][j][k][U2]  =  vcon[2]  ;
      p[i][j][k][U3]  =  vcon[3]  ;
      p[i][j][k][B1]  =  Bcon[1]  ;
      p[i][j][k][B2]  =  Bcon[2]  ;
      p[i][j][k][B3]  =  Bcon[3]  ;
    }
    ind++;
  }

  DEALLOC_ARRAY(rho_tmp, NCELLS);
  DEALLOC_ARRAY(uu_tmp , NCELLS);
  DEALLOC_ARRAY(v1_tmp , NCELLS);
  DEALLOC_ARRAY(v2_tmp , NCELLS);
  DEALLOC_ARRAY(v3_tmp , NCELLS);
  DEALLOC_ARRAY(B1_tmp , NCELLS);
  DEALLOC_ARRAY(B2_tmp , NCELLS);
  DEALLOC_ARRAY(B3_tmp , NCELLS);
#endif


  TRACE_END;
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

  /* BB parameters: */
#if( BBH_SPACETIME ) 
  m_bh1 = 0.5*M;                            /* Mass of left-most BH */
  m_bh2 = 0.5*M;                            /* Mass of right-most BH */
  a_bh1 = 0.;                               /* Spin of BH1 */
  a_bh2 = 0.;                               /* Spin of BH2 */
  initial_bbh_separation = 20.*m_bh_tot;  /* Initial separtion of BHs in units of total mass */

  r_horizon1 = m_bh1;             /* Radius of the horizon about BH1 in Harmonic coordinates */
  r_horizon2 = m_bh2;             /* Radius of the horizon about BH2 in Harmonic coordinates */

  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;

  a = 0.0;          /* Spin of the black hole in units of M */ 
  asq = a*a;
  r_isco    = 0.;
  r_horizon = 0.;

//  Rout     = 13.*initial_bbh_separation;        /* Radial extent of the grid           */
  Rout = 300. ;
  n_within_horizon = 0;            /* number of cells within horizon */
  //  Rin       = Rin_calc();

  Rin = 0.75 * (initial_bbh_separation) ; /*  put inner boundary at 2. times the binary orbit's radius */
  //  R0       = 0.9*Rin;
  R0 = 0.;

#if( SET_TSHRINK )
  /* Time at which to start shrinking the binary,  set to 0. if you want to shrink from the beginning and 
     set to an impossibly large number if you never want to shrink : */
  t_shrink_bbh = 0.;//1.e12; 
  phi_0_bbh    = 0.;
#endif

#else  /*   #if( BBH_SPACETIME )  */
  /* Single BH Parameters: */
  a = 0.; 
  asq = a*a;
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  Rout     = 300.;         /* Radial extent of the grid           */
  n_within_horizon = 5;   /* number of cells within horizon */
  Rin       = Rin_calc();
  
  r_horizon1 = r_horizon2 = r_horizon;
  m_bh1 = m_bh2 = M;
  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;
  initial_bbh_separation = 0.;

#endif

  h_slope   = 0.13;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

#if( EQUATORIAL_RUN ) 
  th_cutout = 0.;
#else 
  th_cutout = 2.*SMALL;     
  /* 0.02655 gives a cutout of ~0.045Pi for h_slope = 0.35 */
  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */ 
  //  th_cutout = 0.0* M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
#endif                 

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
  int   i,j,k,l,n;
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
  struct of_coord *coords_r;

  double lambda,  h, f_l, l_ang, epsilon;
  double ucon[NDIM], ucov[NDIM], ucon_bl[NDIM];
  double utcov, upcov, utcon, upcon, Omega; 

  extern void  coord_of_xp2( double *xp, struct of_coord *coords ) ;
  extern void bin_data_in_r(int bin_state, int n_in, double *data_in, double *r_beg, double *r_end, double *data_out);

  /***********************************************************************************************
    BEGIN:  Code block to integrate the equations from the equator over theta/X2: 
  **********************************************************************************************/  
#define X2_INTEGRATION_CODE_BLOCK                                                                            \
  {                                                                                                          \
    j = N2S + j_glob - jbeg ;                                                                                \
    if( (N2S <= j) && (j <= N2E) )  {                                                                        \
      r_beg[n] = r_end[n] = 1.e-20;                                                                          \
      get_coord(i,j,k,CENT,ncurr,coords);                                                                    \
      get_special_geometry( coords, &geom_eq, metric_type);                                                  \
      gtt = geom_eq.gcon[TT][TT];   gtp = geom_eq.gcon[TT][PH];   gpp = geom_eq.gcon[PH][PH];                \
      xsol[0] = lambda;                                                                                      \
      retval = general_newton_raphson( xsol, 1, resid_func);                                                 \
      lambda = xsol[0];                                                                                      \
      if( retval ) {                                                                                         \
	fprintf(stderr,"init_gen_disk.c: init_all_prims():  ");                                              \
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
        if( (h > 1.) && finite(h) && finite(utcon) && finite(upcon) ) {	                                     \
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
          height_numer[n] = sqrt(geom_eq.gcov[2][2]) * fabs(coords->x[2]-0.5*M_PI) * rho;                    \
          height_denom[n] = rho;                                                                             \
          get_coord(i,j,k,FACE1,ncurr,coords);                                                               \
	  r_beg[n] = coords->r;                                                                              \
	  get_coord((i+1),j,k,FACE1,ncurr,coords);                                                           \
          r_end[n] = coords->r;	                                                                             \
        }                                                                                                    \
      }                                                                                                      \
      n++;                                                                                                   \
    }                                                                                                        \
  }                                                                                                          

  /***********************************************************************************************
    END:  Code block
  **********************************************************************************************/  

  TRACE_BEG;

  /***********************************************************************************************
    Set parameters used for all points : 
  **********************************************************************************************/  
  adiab    = 1.e-2;    /* constant "K" in P = K rho**gam , assume this EOS in disk */
#if( BBH_SPACETIME ) 
  rin      = 3.*(initial_bbh_separation+m_bh_tot); /* Radius of inner edge of the disk in BL coordinates */
  rpmax    = 5.*(initial_bbh_separation+m_bh_tot); /* Radius of the pressure maximum in BL coordinates */
#else
  rin      = 20.;     /* Radius of inner edge of the disk ?? */
  rpmax    = 35.;     /* Radius of the pressure maximum ?? */
#endif
  prograde = 1;       /* non-zero for prograde disks */ 
  //  GRscale  = 1.e-2;    /* Relative magnitude of random perturbations to u */
  GRscale  = 0.;    /* Relative magnitude of random perturbations to u */
  beta     = 1.e2 ;   /* Plasma beta parameter */

#if( (METRIC_TYPE_CHOICE==METRIC_GENERAL_STATIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG2)&&(METRIC_DYNAMIC_TYPE_CHOICE!=METRIC_DYNAMIC_KS_SPHERICAL) )
  metric_type = METRIC_GENERAL_PHI_AVG2;
  int metric_type2 = metric_type;
#else
  metric_type = METRIC_BL_SPHERICAL;
  int metric_type2 = METRIC_KS_SPHERICAL;
#endif

#if(USE_SIMPLE_EOS) 
  N1ALL_LOOP { Katm[i] = adiab; }
#endif

  general_keplerian_functions(prograde, rin  , &Omega_in, &l_ang_in_min, ucov, ucon);
  general_keplerian_functions(prograde, rpmax, &Omega_p , &l_ang_p     , ucov, ucon);
  lambda_p  = sqrt( l_ang_p  / Omega_p  );

  l_ang_in_max = l_ang_p;

  if( doing_height_survey ) { 
    /* for when you want to search for desired h/r : */
    l_ang_in = l_factor*l_ang_in_max + (1.-l_factor)*l_ang_in_min; 
  }
  else { 
    /* for a single BH  h/r=0.1 disk with rin=30, rpmax=50  :  */
    //  l_ang_in = 6.2162995320638608 ;

    /* for a single BH  h/r=0.1 disk with rin=30, rpmax=45  :  */
    //  l_ang_in = 6.3454464305836780;

    /* for a single BH  h/r=0.1 disk with rin=20, rpmax=35  :  */
    //  l_ang_in = 5.2196077516713117;

    /* for a BBH  h/r=0.1 disk with rin=3*a, rpmax=5*a, a=20M  :  */
    l_ang_in = 8.7431737221816430e0;

    /* for a BBH h/r=0.1 disk with rin=3*a, rpmax=5*a, a=30M : */
    //l_ang_in = 1.0536939733523131e+01;

    /* for a BBH h/r=0.1 disk with rin=3*a rpmax=5*a, a=40M : */
    //l_ang_in = 1.2064466188780221e+01;

    /* for a BBH  h/r=0.1 disk with rin=3*a, rpmax=5*a, a=50M  :  */
    // l_ang_in = 1.3391110387962730e+01 ;

    /* for a BBH  h/r=0.1 disk with rin=2.5*a, rpmax=4.2*a, a=50M  :  */
    // l_ang_in = 1.2213162202203518e+01 ;

    /* for a BBH  h/r=0.1 disk with rin=2.5*a, rpmax=4.2*a, a=100M  :  */
    //l_ang_in = 1.7078426447122045e+01 ; 

    /* for a BBH  h/r=0.1 disk with rin=2.5*a, rpmax=4.2*a, a=1000M  :  */
    //  l_ang_in =   
  }

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

  if( myid == printer_pid ) { 
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
  }

  /**********************************************************************************************
    Set the primitives to zeros as default values to be filled in by floor
  **********************************************************************************************/
  N1_LOOP N2_LOOP N3_LOOP PLOOP  {  p[i][j][k][l] = 0.; }

  /**********************************************************************************************
    Setup the arrays for calculating the scaleheight :
  **********************************************************************************************/
  set_bin_arrays(); 

  n_in = NCELLS;
  ALLOC_ARRAY(height_denom_bins, n_r_bins); 
  ALLOC_ARRAY(height_numer_bins, n_r_bins); 
  ALLOC_ARRAY(height_denom     , n_in    ); 
  ALLOC_ARRAY(height_numer     , n_in    ); 
  ALLOC_ARRAY(r_beg            , n_in    ); 
  ALLOC_ARRAY(r_end            , n_in    ); 

  /**********************************************************************************************
    Set the values in the disk :
  **********************************************************************************************/  
  ALLOC_ARRAY( lambda_eq, N1TOT); 
  ALLOC_ARRAY(      l_eq, N1TOT); 
  ALLOC_ARRAY(      h_eq, N1TOT); 
  ALLOC_ARRAY(  coords_r, N1TOT);

  N1ALL_LOOP { lambda_eq[i] = -1. ;  } 
  N1ALL_LOOP {      l_eq[i] = -1. ;  } 
  N1ALL_LOOP {      h_eq[i] = -1. ;  } 

  j_glob_eq = totalsize[2]/2;
  get_global_ijk(0,0,0,&i,&jbeg,&k);
  get_global_ijk(N1E-N1S,N2E-N2S,N3E-N3S,&i,&jend,&k);

  double x[NDIM],xp[NDIM];

  n = 0; 

  N3_LOOP { 
    j = j_glob_eq - globalpos[2];   /* this is so we can use coord() to correctly find global xp[] values with local function  */
    N1_LOOP { 
      coord(i,j,k,CENT,xp);
      coords_r[i].i   = i; 
      coords_r[i].j   = j; 
      coords_r[i].k   = k; 
      coords_r[i].pos = CENT; 
      coord_of_xp2(xp,&(coords_r[i]));
    }
    calc_eq_values(lambda_eq, l_eq, h_eq, coords_r);      

    if( k == N3S ) { 
      N1_LOOP{ 
	fprintf(stdout,"eq-values [%d,%d,%d,%d]  %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e\n",i,j,k,myid,coords_r[i].r,coords_r[i].x[TH],coords_r[i].x[PH],lambda_eq[i],l_eq[i],h_eq[i]); 
      }
      fflush(stdout);
    }

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
  }

  if( n > NCELLS ) { 
    fprintf(stdout,"init_all_prims(): mismatch with the counting:  %d %d \n", n, NCELLS); fflush(stdout); 
    fail(FAIL_BASIC,0) ; 
  }

  bin_data_in_r(1, n, height_numer, r_beg, r_end, height_numer_bins); 
  bin_data_in_r(0, n, height_denom, r_beg, r_end, height_denom_bins); 

  FREE(height_numer); 
  FREE(height_denom); 
  ALLOC_ARRAY(height_denom, n_r_bins); 
  ALLOC_ARRAY(height_numer, n_r_bins); 

  mpi_sum_to_head(height_numer_bins,height_numer,n_r_bins);
  mpi_sum_to_head(height_denom_bins,height_denom,n_r_bins);

  if( myid==printer_pid ) { 
    for(i=0;i<n_r_bins;i++) { 
      if( height_denom[i] != 0. ) { 
	height_numer[i] /= height_denom[i]*(r_bins[i]+0.5*dr_bins[i]); 
      }
      if( !finite(height_numer[i]) ) { height_numer[i] = 0.; }
    }
    h_o_r = point_interp(n_r_bins,r_bins,height_numer,rpmax);
    fprintf(stdout,"l_ang_in  h_o_r(rpmax) q  eta = %26.16e %26.16e  %26.16e %26.16e \n",
	    l_ang_in, h_o_r,q,eta); 
    fflush(stdout);
    for(i=0;i<n_r_bins;i++) { fprintf(stdout,"r height/r = %26.16e  %26.16e \n", r_bins[i], height_numer[i]); }
    fflush(stdout);
  }
  sync_val(&h_o_r);

  FREE(height_numer); 
  FREE(height_denom); 
  FREE(height_numer_bins); 
  FREE(height_denom_bins); 
  FREE(r_beg); 
  FREE(r_end); 

  if( doing_height_survey ) {   TRACE_END;  return; } 

  /**********************************************************************************************
    Now copy the data over all phi, transforming to numerical coordinatea and adding random 
     perturbations along the way:
  **********************************************************************************************/  
  N1_LOOP N2_LOOP N3_LOOP { 
    pr  = p[i][j][k];
    rho = pr[RHO]; 

    get_coord(i,j,k,CENT,ncurr,coords);
    get_geometry(i,j,k,CENT,ncurr,geom); // get (x1,x2,x3) metric, from here on it's only (x1-3)
   
    if( rho > 0. ) { 
      uu          = pr[UU ]; 
      ucon_bl[TT] = pr[U1 ];
      ucon_bl[RR] = 0.;
      ucon_bl[TH] = 0.;
      ucon_bl[PH] = pr[U3 ];

      //      fprintf(stdout,"KATM = %26.16e  %26.16e  %26.16e \n", rho,uu,((gam - 1.) * uu / pow( rho, gam ))) ; fflush(stdout); 

      if( rho > rhomax )  { rhomax = rho ; } 
      if(  uu >   umax )  {   umax =  uu ; } 

      /* Add perturbation to internal energy density : */
      pr[UU]  *= (1. + GRscale * (ranc(0) - 0.5));

#if( (METRIC_TYPE_CHOICE==METRIC_KS_SPHERICAL) ||			\
     ((METRIC_TYPE_CHOICE==METRIC_GENERAL_STATIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC)&&	\
      (METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_KS_SPHERICAL) ) )
    bl_to_ks_con(coords->x,ucon_bl,ucon);  // transform from BL to KS coordinates 
#else
    for(l=0;l<NDIM;l++) { ucon[l] = ucon_bl[l]; }
#endif
    /* Transform 4-velocity to numerical coordinates and into primitive variable velocity */
    transform_rank1con2(coords->dxp_dx,ucon); // transform from (r,th,ph) to (x1,x2,x3)
    ucon2pr( pr, ucon, geom->gcon );  // get prim. velocities;

    }
    else { 
      /* Set the "floor" state: -- we can keep the density and pressure zero, but we need to transfer the velocity 
       in order to handle arbitrarily warped coordinates: */
      pr[U1] = 0. ; 
      pr[U2] = 0. ; 
      pr[U3] = 0. ; 

#if( !BHS_ON_GRID )       
      get_special_geometry( coords, &geom_eq, metric_type2);
      ucon_calc(pr,&geom_eq,ucon); 
    /* Transform 4-velocity to numerical coordinates and into primitive variable velocity */
    transform_rank1con2(coords->dxp_dx,ucon); // transform from (r,th,ph) to (x1,x2,x3)
    ucon2pr( pr, ucon, geom->gcon );  // get prim. velocities;
#endif

    }
  }
 
  FREE( lambda_eq); 
  FREE(      l_eq); 
  FREE(      h_eq); 
  FREE(  coords_r); 
#undef X2_INTEGRATION_CODE_BLOCK
   
  TRACE_END;

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
   -- Dennis extend this for INITIALIZE_MINI_DISKS to use A_\mu from init_mini-disks
      for cells within mini-disk
***********************************************************************************/
static void set_mag_field( void )
{
  int i,j,k,l;

  double rho_avg, rho_cut, Aphi, beta_act, norm, bsq_int, u_int ; 
  struct of_geom *geom;  
  struct of_coord *coords;  
  extern void curl_of_A( double ****A , double ****prim ) ;

#if( INITIALIZE_MINI_DISKS )
  double bsq_int_disk1, bsq_int_disk2, u_int_disk1, u_int_disk2;
  double norm_disk1, norm_disk2;
  double beta_act_disk1, beta_act_disk2;
#endif

  TRACE_BEG;

  if( myid == printer_pid ) { 
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
  }

#if( INITIALIZE_MINI_DISKS )
  int in_mini_disk = 0;
  double x,y,z,xbh1,xbh2,ybh1,ybh2,zbh1,zbh2;
  double r1,r2;
  
  /*************************************************************************************
    Calculate the vector potential, A_mu for all cells not within mini-disk: 
       --  Need to calc. A_\mu on the last ghost cells or all the corners surrounding the 
           physical cells' centers for finite diff. calculation for B^i 
       -- We first calculate A_\phi then we transform to xp;
       -- To this end we implicitly assume that the density distribution is 
           phi-independent (but not necessaryily xp3 independent);
  *************************************************************************************/
  double Acov[NDIM];

  /* Get BH Traj Data */ 
  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  struct of_bbh_traj bbh_traj;
  get_bbh_traj_data(&bbh_traj);
  xbh1 = bbh_traj.xi1x;  ybh1 = bbh_traj.xi1y;  zbh1 = bbh_traj.xi1z;
  xbh2 = bbh_traj.xi2x;  ybh2 = bbh_traj.xi2y;  zbh2 = bbh_traj.xi2z;	
	
  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) for(k=N3S; k<=(N3E+1); k++) { 
	
	/* Are we within a mini-disk? */
	get_coord(i,j,k,CENT,ncurr,coords);
	x = coords->xcart[XX]; y  = coords->xcart[YY]; z  = coords->xcart[ZZ];
	
	/* Set the distance to each BH from center of cell in PN harmonic coordinates*/
	r1   = sqrt( (x - xbh1)*(x - xbh1) + (y - ybh1)*(y - ybh1) + (z - zbh1)*(z - zbh1) );
	r2   = sqrt( (x - xbh2)*(x - xbh2) + (y - ybh2)*(y - ybh2) + (z - zbh2)*(z - zbh2) );
	if( r1min <= r1 && r1 <= r1max )      { in_mini_disk = 1; }
	else if( r2min <= r2 && r2 <= r2max ) { in_mini_disk = 1; }
	else                                  { in_mini_disk = 0; }
	
	if( in_mini_disk ) {
	  /* Move p[i][j][k][B1-3] to ph for curl_of_A calculation later */
	  ph[i][j][k][B1] = p[i][j][k][B1];
	  ph[i][j][k][B2] = p[i][j][k][B2];
	  ph[i][j][k][B3] = p[i][j][k][B3];
	  /* Reset those p compts to 0 */
	  p[i][j][k][B1] = p[i][j][k][B2] = p[i][j][k][B3] = 0.;
	}	
	else {
	  /* If not within a mini-disk calculate A_{\mu} */
	  rho_avg =   
	    p[i  ][j  ][k  ][RHO] +
	    p[i-1][j  ][k  ][RHO] +
	    p[i  ][j-1][k  ][RHO] +
	    p[i-1][j-1][k  ][RHO] +
	    p[i  ][j  ][k-1][RHO] +
	    p[i-1][j  ][k-1][RHO] +
	    p[i  ][j-1][k-1][RHO] +
	    p[i-1][j-1][k-1][RHO] ;
	  
	  rho_avg *= 0.125;
	  get_coord(i,j,k,CORN,ncurr,coords);
	  
#if(   FIELD_CONFIG_CHOICE ==     DIPOLE_FIELD )
	  rho_cut = 0.25 * rhomax; 
	  Aphi = rho_avg - rho_cut ;    /* dipole field */
	  
#elif( FIELD_CONFIG_CHOICE == QUADRUPOLE_FIELD )
	  rho_cut = 0.25 * rhomax; 
	  Aphi = (rho_avg - rho_cut) * cos(coords->x[TH]) / rhomax ;
#endif
	  
	  /* Truncate A_\mu at zero so that field lines lie within disk */
	  Acov[0] = 0.;
	  Acov[1] = 0.;
	  Acov[2] = 0.;
	  Acov[3] = 0.;
	  
	  if( rho_avg > rho_cut ) { 
	    Acov[3] = Aphi;
	    transform_rank1cov2(coords->dx_dxp,Acov);
	  }
	  
	  ph[i][j][k][B1] = Acov[1];
	  ph[i][j][k][B2] = Acov[2];
	  ph[i][j][k][B3] = Acov[3];
	}
      }

  /*************************************************************************************
    Calculate the poloidal field components given A_\phi
       --  Note that we can take derivatives w.r.t. numerical coordinates here because 
           the connection is symmetric in its lower indices, we just have to use 
           the numerical metric determinant; 
  *************************************************************************************/
  curl_of_A(ph,p);


  /*************************************************************************************
    Calculate global internal energy and magnetic field energy: --Dennis need a separate sum for each disk to normalize
  *************************************************************************************/
  bsq_int = u_int = 0.;
  bsq_int_disk1 = u_int_disk1 = 0.;
  bsq_int_disk2 = u_int_disk2 = 0.;

  N1_LOOP  N2_LOOP   N3_LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    get_geometry(i,j,k,CENT,ncurr,geom); 
    x = coords->xcart[XX]; y  = coords->xcart[YY]; z  = coords->xcart[ZZ];
    
    /* Set the distance to each BH from center of cell in PN harmonic coordinates*/
    r1   = sqrt( (x - xbh1)*(x - xbh1) + (y - ybh1)*(y - ybh1) + (z - zbh1)*(z - zbh1) );
    r2   = sqrt( (x - xbh2)*(x - xbh2) + (y - ybh2)*(y - ybh2) + (z - zbh2)*(z - zbh2) );
    if( r1min <= r1 && r1 <= r1max ) { 
      bsq_int_disk1 += geom->g * bsq_calc(p[i][j][k],geom) ;
      u_int_disk1   += geom->g * p[i][j][k][UU]; 
    }
    else if( r2min <= r2 && r2 <= r2max ) {
      bsq_int_disk2 += geom->g * bsq_calc(p[i][j][k],geom) ;
      u_int_disk2   += geom->g * p[i][j][k][UU]; 
    }
    else {
      bsq_int += geom->g * bsq_calc(p[i][j][k],geom) ;
      u_int   += geom->g * p[i][j][k][UU]; 
    }
  }

  /************************************************************************************
     Normalize magnetic field to match specified value of beta (using umax);
   ************************************************************************************/
  mpi_global_sum( &bsq_int       );
  mpi_global_sum( &bsq_int_disk1 );
  mpi_global_sum( &bsq_int_disk2 );
  if( myid == printer_pid ) {  fprintf(stdout,"initial bsq_int       = %28.18e \n", bsq_int      ); fflush(stdout);  }
  if( myid == printer_pid ) {  fprintf(stdout,"initial bsq_int_disk1 = %28.18e \n", bsq_int_disk1); fflush(stdout);  }
  if( myid == printer_pid ) {  fprintf(stdout,"initial bsq_int_disk2 = %28.18e \n", bsq_int_disk2); fflush(stdout);  }

  mpi_global_sum( &u_int       );
  mpi_global_sum( &u_int_disk1 );
  mpi_global_sum( &u_int_disk2 );
  if( myid == printer_pid ) {  fprintf(stdout,"initial u_int       = %28.18e \n", u_int      ); fflush(stdout); }
  if( myid == printer_pid ) {  fprintf(stdout,"initial u_int_disk1 = %28.18e \n", u_int_disk1); fflush(stdout); }
  if( myid == printer_pid ) {  fprintf(stdout,"initial u_int_disk2 = %28.18e \n", u_int_disk2); fflush(stdout); }

  beta_act       = (gam - 1.)*u_int      /(0.5*bsq_int      ) ;
  beta_act_disk1 = (gam - 1.)*u_int_disk1/(0.5*bsq_int_disk1) ;
  beta_act_disk2 = (gam - 1.)*u_int_disk2/(0.5*bsq_int_disk2) ;
  if( myid == printer_pid ) {  fprintf(stdout,"initial beta      : %28.18e (should be %28.18e)\n",beta_act      ,beta      ); fflush(stdout); }
  if( myid == printer_pid ) {  fprintf(stdout,"initial beta_disk1: %28.18e (should be %28.18e)\n",beta_act_disk1,beta_disk1); fflush(stdout); }
  if( myid == printer_pid ) {  fprintf(stdout,"initial beta_disk2: %28.18e (should be %28.18e)\n",beta_act_disk2,beta_disk2); fflush(stdout); }

  norm       = sqrt(beta_act      /beta      ) ;
  norm_disk1 = sqrt(beta_act_disk1/beta_disk1) ;
  norm_disk2 = sqrt(beta_act_disk2/beta_disk2) ;
  
  N1_LOOP  N2_LOOP N3_LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);
    x = coords->xcart[XX]; y  = coords->xcart[YY]; z  = coords->xcart[ZZ];
    
    /* Set the distance to each BH from center of cell in PN harmonic coordinates*/
    r1   = sqrt( (x - xbh1)*(x - xbh1) + (y - ybh1)*(y - ybh1) + (z - zbh1)*(z - zbh1) );
    r2   = sqrt( (x - xbh2)*(x - xbh2) + (y - ybh2)*(y - ybh2) + (z - zbh2)*(z - zbh2) );
    if( r1min <= r1 && r1 <= r1max ) { 
      p[i][j][k][B1] *= norm_disk1 ;
      p[i][j][k][B2] *= norm_disk1 ;
      p[i][j][k][B3] *= norm_disk1 ;
    }
    else if( r2min <= r2 && r2 <= r2max ) {
      p[i][j][k][B1] *= norm_disk2 ;
      p[i][j][k][B2] *= norm_disk2 ;
      p[i][j][k][B3] *= norm_disk2 ;
    }
    /* /\* If not in a mini-disk or circumbinary zero out the B field *\/ */
    /* else if( coords->x[RR] < 1.5 * (initial_bbh_separation + m_bh_tot) ) { */
    /*   p[i][j][k][B1] *= 0.; */
    /*   p[i][j][k][B2] *= 0.; */
    /*   p[i][j][k][B3] *= 0.; */
    /* } */
    else {
      p[i][j][k][B1] *= norm ;
      p[i][j][k][B2] *= norm ;
      p[i][j][k][B3] *= norm ;
    }
  }

  bsq_int       = 0.;
  bsq_int_disk1 = 0.;
  bsq_int_disk2 = 0.;
  N1_LOOP  N2_LOOP N3_LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    get_geometry(i,j,k,CENT,ncurr,geom); 
    x = coords->xcart[XX]; y  = coords->xcart[YY]; z  = coords->xcart[ZZ];
    
    /* Set the distance to each BH from center of cell in PN harmonic coordinates*/
    r1   = sqrt( (x - xbh1)*(x - xbh1) + (y - ybh1)*(y - ybh1) + (z - zbh1)*(z - zbh1) );
    r2   = sqrt( (x - xbh2)*(x - xbh2) + (y - ybh2)*(y - ybh2) + (z - zbh2)*(z - zbh2) );
    
    if( r1min <= r1 && r1 <= r1max )      { bsq_int_disk1 += geom->g * bsq_calc(p[i][j][k],geom) ; }
    else if( r2min <= r2 && r2 <= r2max ) { bsq_int_disk2 += geom->g * bsq_calc(p[i][j][k],geom) ; }
    else                                  { bsq_int       += geom->g * bsq_calc(p[i][j][k],geom) ; }
  }

  mpi_global_sum( &bsq_int );
  mpi_global_sum( &bsq_int_disk1 );
  mpi_global_sum( &bsq_int_disk2 );
  if( myid == printer_pid ) {   fprintf(stdout,"final bsq_int       = %28.18e \n", bsq_int      ); fflush(stdout); }
  if( myid == printer_pid ) {   fprintf(stdout,"final bsq_int_disk1 = %28.18e \n", bsq_int_disk1); fflush(stdout); }
  if( myid == printer_pid ) {   fprintf(stdout,"final bsq_int_disk2 = %28.18e \n", bsq_int_disk2); fflush(stdout); }

  beta_act       = (gam - 1.)*u_int      /(0.5*bsq_int      ) ;
  beta_act_disk1 = (gam - 1.)*u_int_disk1/(0.5*bsq_int_disk1) ;
  beta_act_disk2 = (gam - 1.)*u_int_disk2/(0.5*bsq_int_disk2) ;
  if( myid == printer_pid ) {  fprintf(stdout,"final beta      : %28.18e (should be %28.18e)\n",beta_act      ,beta      ); fflush(stdout); }
  if( myid == printer_pid ) {  fprintf(stdout,"final beta_disk1: %28.18e (should be %28.18e)\n",beta_act_disk1,beta_disk1); fflush(stdout); }
  if( myid == printer_pid ) {  fprintf(stdout,"final beta_disk2: %28.18e (should be %28.18e)\n",beta_act_disk2,beta_disk2); fflush(stdout); }

#else //#if(INITIALIZE_MINI_DISKS)
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

	rho_avg =   
	  p[i  ][j  ][k  ][RHO] +
	  p[i-1][j  ][k  ][RHO] +
	  p[i  ][j-1][k  ][RHO] +
	  p[i-1][j-1][k  ][RHO] +
	  p[i  ][j  ][k-1][RHO] +
	  p[i-1][j  ][k-1][RHO] +
	  p[i  ][j-1][k-1][RHO] +
	  p[i-1][j-1][k-1][RHO] ;
	
	rho_avg *= 0.125;

	get_coord(i,j,k,CORN,ncurr,coords);

#if(   FIELD_CONFIG_CHOICE ==     DIPOLE_FIELD )
	rho_cut = 0.25 * rhomax; 
	Aphi = rho_avg - rho_cut ;    /* dipole field */

#elif( FIELD_CONFIG_CHOICE == QUADRUPOLE_FIELD )
	rho_cut = 0.25 * rhomax; 
	Aphi = (rho_avg - rho_cut) * cos(coords->x[TH]) / rhomax ;
#endif

	/* Truncate A_\mu at zero so that field lines lie within disk */
	Acov[0] = 0.;
	Acov[1] = 0.;
	Acov[2] = 0.;
	Acov[3] = 0.;

	if( rho_avg > rho_cut ) { 
	  Acov[3] = Aphi;
	  transform_rank1cov2(coords->dx_dxp,Acov);
	}

	ph[i][j][k][B1] = Acov[1];
	ph[i][j][k][B2] = Acov[2];
	ph[i][j][k][B3] = Acov[3];
  }

  /*************************************************************************************
    Calculate the poloidal field components given A_\phi
       --  Note that we can take derivatives w.r.t. numerical coordinates here because 
           the connection is symmetric in its lower indices, we just have to use 
           the numerical metric determinant; 
  *************************************************************************************/
  curl_of_A(ph,p);


  /*************************************************************************************
    Calculate global internal energy and magnetic field energy: --Dennis need a separate sum for each disk to normalize
  *************************************************************************************/
  bsq_int = u_int = 0.;

  N1_LOOP  N2_LOOP   N3_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); 

    bsq_int += geom->g * bsq_calc(p[i][j][k],geom) ;
    u_int   += geom->g * p[i][j][k][UU]; 
  }

  /************************************************************************************
     Normalize magnetic field to match specified value of beta (using umax);
   ************************************************************************************/
  mpi_global_sum( &bsq_int );
  if( myid == printer_pid ) {  fprintf(stdout,"initial bsq_int = %28.18e \n", bsq_int); fflush(stdout);  }

  mpi_global_sum( &u_int );
  if( myid == printer_pid ) {  fprintf(stdout,"initial u_int = %28.18e \n", u_int); fflush(stdout); }

  beta_act = (gam - 1.)*u_int/(0.5*bsq_int) ;
  if( myid == printer_pid ) {  fprintf(stdout,"initial beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout); }

  norm = sqrt(beta_act/beta) ;
  
  N1_LOOP  N2_LOOP N3_LOOP { 
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
    p[i][j][k][B3] *= norm ;
  }

  bsq_int = 0.;
  N1_LOOP  N2_LOOP N3_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); 
    bsq_int += geom->g * bsq_calc(p[i][j][k],geom) ;
  }

  mpi_global_sum( &bsq_int );
  if( myid == printer_pid ) {   fprintf(stdout,"final bsq_int = %28.18e \n", bsq_int); fflush(stdout); }

  beta_act = (gam - 1.)*u_int/(0.5*bsq_int) ;
  if( myid == printer_pid ) {  fprintf(stdout,"final beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout); }

#endif //#if( INITIALIZE_MINI_DISKS )  

  TRACE_END;
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

       -- x[]   = set of independent variables to solve for;
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


// if( output_now ) { 
//   fprintf(stdout,"r_eq lambda l resid = %26.16e  %26.16e  %26.16e  %26.16e  \n", r_glob, lambda, (eta * pow( lambda, (2.-q))), resid[0] ); 
//   fprintf(stdout,"resid_func  %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e\n",r_glob,eta,q,gtt,gtp,gpp);
//   fflush(stdout); 
// }

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

/*********************************************************************************/
/*********************************************************************************
   calc_eq_values():
  ------------------

    -- calculates the relevant disk functions along the equator;
*********************************************************************************/
static void calc_eq_values(double *lambda_eq, double *l_eq, double *h_eq, struct of_coord *coords_r)
{
  int i,j,k,retval;
  double Omega, l_ang, ucon[NDIM], ucov[NDIM], xsol[1],f_l,utcov;
  struct of_geom geom_eq;


  N1ALL_LOOP { lambda_eq[i] = l_eq[i] = h_eq[i] = -1. ;  } 
      
  N1_LOOP { 
    general_keplerian_functions(prograde, coords_r[i].r , &Omega, &l_ang, ucov, ucon);

    //    fprintf(stdout,"Keplerian functions = %26.16e %26.16e %26.16e %26.16e %26.16e \n",coords_r[i].r, Omega, l_ang, ucov[TT], ucon[TT]); fflush(stdout); 

    if( coords_r[i].r >= rin ) { 
      if( (i <= N1S)  || (lambda_eq[i-1] < 0.) ) { 
	l_eq[i] = l_ang_in; 
	lambda_eq[i] = pow(  (l_eq[i]/eta), exp_5 );
      } 
      else { 
	lambda_eq[i] = lambda_eq[i-1]; 
      }

      get_special_geometry( &(coords_r[i]), &geom_eq, metric_type);
      gtt = geom_eq.gcon[TT][TT];
      gtp = geom_eq.gcon[TT][PH];
      gpp = geom_eq.gcon[PH][PH];

      //      r_glob = coords_r[i].r;
      //      fprintf(stdout,"gtttp %26.16e %26.16e %26.16e %26.16e \n",r_glob,gtt,gtp,gpp); fflush(stdout); 

      xsol[0] = lambda_eq[i];
      retval = general_newton_raphson( xsol, 1, resid_func);

      if( retval ) { 
	fprintf(stderr,"init_gen_disk.c: calc_eq_values():  Failed to find a solution at [i=%6d]   [r=%26.16e]  \n", i, coords_r[i].r); 
	fflush(stdout); 
      }
      else { 
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
//		coords_r[i].r,utcov, f_l, l_ang, lambda_eq[i], h_eq[i], Omega); 
//	fflush(stdout);

      }

      if( (h_eq[i] > 1.) && (coords_r[i].r > rdiskout) ) { rdiskout = coords_r[i].r ; }
      
//      fprintf(stdout,"equator-values: myid i r lambda l h = %6d %6d %26.16e %26.16e %26.16e %26.16e \n",
//	      myid, i,coords_r[i].r,lambda_eq[i],l_eq[i],h_eq[i]); 
//      fflush(stdout);


    }   /*   if( coords_r[i].r >= rin )   */
  }
  
  return;
}

#if( INITIALIZE_MINI_DISKS )
/****************************************/
/* M I N I - D I S K  F U N C T I O N S */
/****************************************/
/***********************************************************************************/
/***********************************************************************************
  init_mini_disks(): 
  ------------------
   -- Initialize mini-disks around each BH ( all hydro quantities )
   -- Calculates A_{\mu} wrt numerical coordinates and stores A_{1-3} in pr[B1-B3]
***********************************************************************************/
void init_mini_disks(int i, int j, int k, struct of_coord *coords, double *pr)
{
  int mu, nu,l;
  double r,th, ph;
  double x,y,z;
  double xbh1, xbh2, r1, r2, Mass, ybh1, ybh2, zbh1, zbh2;
  double vbx1, vbx2, vby1, vby2, vbz1, vbz2, vsq;
  double rdisk, phdisk;
  double norm, center, stddev;
  double rho, uu, OmegaK, v_r;
  double WW;
  double Spin,a_bh1,a_bh2;
  struct of_geom geom_bl;
  struct of_geom *geom;

  double ucon[NDIM], ucon_cart[NDIM]; // Center of mass reference frame
  double ucon_prime[NDIM], ucon_cart_prime[NDIM]; //BH reference frame
  double xNZ[NDIM], xIZ[NDIM], xBL[NDIM];
  double dxIZ_dxBL[NDIM][NDIM], dxNZ_dxIZ[NDIM][NDIM], dxIZ_dxNZ[NDIM][NDIM],dxs_dxc[NDIM][NDIM],dxBL_dxIZ[NDIM][NDIM],dxc_dxs[NDIM][NDIM];
  double gcov_special[NDIM][NDIM];

  double coord_cart_CM[NDIM], coord_cart_BH[NDIM];
  double ucon_BL[NDIM];
  
  double vboost[NDIM];
  double Lambda2BH[NDIM][NDIM], Lambda2CM[NDIM][NDIM];

  /* Definitions for stolen bits from init_all_prims */
  int retval;
  double *lambda_eq, *l_eq, *h_eq;
  double xsol[1];

  struct of_geom geom_eq;
  //struct of_geom *geom;
  struct of_coord *coords_disk;
  struct of_coord coords_loc;
  struct of_coord *coords_r;

  double lambda,  h, f_l, l_ang, epsilon;
  double ucov[NDIM], ucon_bl[NDIM], Acov_bl[NDIM];
  double utcov, upcov, utcon, upcon, Omega; 

  int in_mini_disk = 0;

  int count;
  double l_factor, reldiff;
  double lf_min = 0.;   /* use maximum range to start bisection search */
  double lf_max = 1.;

  const int doing_height_survey = 1; 
  const double tol              = 1e-13;                    /* tolerance for the height survey bisection search */
  const double h_o_r_target     = 1.0044771758781224e-01;   /* set this to be your desired scaleheight when "doing_height_survey" */
  const int max_iterations = 200;    /* maximum number of iterations allowed in bisection before giving up */


  TRACE_BEG;

  /* Set the radial extent of each mini-disk */
  r1max = 0.4  * initial_bbh_separation;
  r2max = 0.4  * initial_bbh_separation;
  r1min = 2.5;//0.05 * initial_bbh_separation;
  r2min = 2.5;//0.05 * initial_bbh_separation;

  /* Set the coordinates of cell and BHs */
  r = coords->x[RR];     th = coords->x[TH];     ph = coords->x[PH];
  x = coords->xcart[XX]; y  = coords->xcart[YY]; z  = coords->xcart[ZZ];
  /* Get BH Traj Data */ 
  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  struct of_bbh_traj bbh_traj;
  get_bbh_traj_data(&bbh_traj);
  xbh1 = bbh_traj.xi1x;  ybh1 = bbh_traj.xi1y;  zbh1 = bbh_traj.xi1z;
  xbh2 = bbh_traj.xi2x;  ybh2 = bbh_traj.xi2y;  zbh2 = bbh_traj.xi2z;

  /* Set the distance to each BH from center of cell in PN harmonic coordinates*/
  r1   = sqrt( (x - xbh1)*(x - xbh1) + (y - ybh1)*(y - ybh1) + (z - zbh1)*(z - zbh1) );
  r2   = sqrt( (x - xbh2)*(x - xbh2) + (y - ybh2)*(y - ybh2) + (z - zbh2)*(z - zbh2) );

  if( r1 <= r1max && r1 >= r1min ) {
    in_mini_disk = 1;
    /* Set desired disk parameters for disk1 */    
    MASS_TYPE = MBH1;
    Mass  = m_bh1;
    Spin  = 0.;
    /* Set extent of disk in BL coordinates --dennis here we really need to be careful
       that our disk contains rpmax and rdisk_max etc. calculate rdisk_max and rdisk_min
       by transforming the PN harmonic r1max/r1min to BL coordinates!!!!!*/
    

    if( first_call_disk1 ) {
      lf_min = 0.;   /* use maximum range to start bisection search */
      lf_max = 1.;
  
      first_call_disk1 = 0;
      ALLOC_4D_ARRAY(pd1,N1TOT_MINI,N2TOT_MINI,N3TOT_MINI,NFUNCS_MINI); 
      /* Set inner and outter extent in BL coordinates to create solution */
      rdisk_min = 2.;
      rdisk_max = 60.;

      adiab    = 1.e-2;    /* constant "K" in P = K rho**gam , assume this EOS in disk */
      rin      = 2.5;     /* Radius of inner edge of the disk ??  units of Mtot*/
      rpmax    = 4.2;     /* Radius of the pressure maximum ??  units of Mtot*/
      prograde = 1;       /* non-zero for prograde disks */ 
      //  GRscale  = 1.e-2;    /* Relative magnitude of random perturbations to u */
      GRscale  = 0.;    /* Relative magnitude of random perturbations to u */
      beta_disk1  = 1.e2 ;   /* Plasma beta parameter */
      
      metric_type = METRIC_BL_SPHERICAL; /* Metric type to use for mini-disk solution */

      if( doing_height_survey ) {
	count=0;
	do {
	  l_factor = 0.5*(lf_min+lf_max);
	  
	  create_disk_solution(MBH1,doing_height_survey,l_factor);
	  if( h_o_r <= h_o_r_target ) { lf_min = l_factor; }
	  else                        { lf_max = l_factor; }

	  reldiff = fabs(lf_max - lf_min)/(l_factor);
	  if( myid == printer_pid ) { 
	    fprintf(stdout,"bisection  %10d %26.16e  %26.16e %26.16e %26.16e  %26.16e %26.16e %26.16e \n", count,l_ang_in, lf_min,lf_max,l_factor,reldiff,h_o_r,h_o_r_target); fflush(stdout);
	  }
	  count++;
	} while( (reldiff > tol) && ( count < max_iterations) );
	//if( myid == printer_pid ) {
	fprintf(stdout,"Final Mini-Disk 1 solution: rin %20.16e rpmax %20.16e l_ang_in %20.16e h_o_r %20.16e\n",rin,rpmax,l_ang_in,h_o_r); fflush(stdout);
	  //}
      }
      else {      
      /* For now we assume solution already solved for
	 -- Improve this in future to automate using height_survey */
      /* for m=0.5,a=0,rin=10,rpmax=15,rout=60 */
      l_ang_in = 2.6787442569269468;
      /* for a single BH h/r=0.1 disk with rin=20, rpmax=35, mass=0.5,spin=0 */
      //l_ang_in =3.4969144610144682 ; //wrong value
      
      /* Allocate and store mini-disk solution */
      create_disk_solution(MBH1, 0, 0.5);     
      }
    }    
  }
  else if( r2 <= r2max && r2 >= r2min ) {
    in_mini_disk = 1;
    /* Set desired disk parameters for disk2 */
    MASS_TYPE = MBH2;
    Mass = m_bh2;    
    Spin = 0.;
   
    if( first_call_disk2 ) {
      lf_min = 0.;   /* use maximum range to start bisection search */
      lf_max = 1.;
  
      first_call_disk2 = 0;
      ALLOC_4D_ARRAY(pd2,N1TOT_MINI,N2TOT_MINI,N3TOT_MINI,NFUNCS_MINI); 
      /* Set inner and outter extent in BL coordinates to create solution */
      rdisk_min = 2.;
      rdisk_max = 60.;

      adiab    = 1.e-2;    /* constant "K" in P = K rho**gam , assume this EOS in disk */
      rin      = 2.5;     /* Radius of inner edge of the disk ??  units of Mtot*/
      rpmax    = 4.2;     /* Radius of the pressure maximum ??  units of Mtot*/
      prograde = 1;       /* non-zero for prograde disks */ 
      //  GRscale  = 1.e-2;    /* Relative magnitude of random perturbations to u */
      GRscale  = 0.;    /* Relative magnitude of random perturbations to u */
      beta_disk2 = 1.e2 ;   /* Plasma beta parameter */
      
      metric_type = METRIC_BL_SPHERICAL; /* Metric type to use for mini-disk solution */
      
      if( doing_height_survey ) {
      	count=0;
      	do {
      	  l_factor = 0.5*(lf_min+lf_max);
	  
      	  create_disk_solution(MBH2,doing_height_survey,l_factor);
      	  if( h_o_r <= h_o_r_target ) { lf_min = l_factor; }
      	  else                        { lf_max = l_factor; }

      	  reldiff = fabs(lf_max - lf_min)/(l_factor);
      	  if( myid == printer_pid ) {
      	    fprintf(stdout,"bisection  %10d %26.16e  %26.16e %26.16e %26.16e  %26.16e %26.16e %26.16e \n", count,l_ang_in, lf_min,lf_max,l_factor,reldiff,h_o_r,h_o_r_target); fflush(stdout);
      	  }
      	  count++;
      	} while( (reldiff > tol) && ( count < max_iterations) );
	//if( myid == printer_pid ) {
	  fprintf(stdout,"Final Mini-Disk 2 solution: rin %20.16e rpmax %20.16e l_ang_in %20.16e h_o_r %20.16e\n",rin,rpmax,l_ang_in,h_o_r); fflush(stdout);
	  //}
      }
      else {
	/* For now we assume solution already solved for
	   -- Improve this in future to automate using height_survey */
	/* for m=0.5,a=0,rin=10,rpmax=15,rout=60 */
	l_ang_in = 2.6787442569269468;
	/* for a single BH h/r=0.1 disk with rin=20, rpmax=35, mass=0.5,spin=0 */
	//l_ang_in =3.4969144610144682 ; //wrong value
	
	/* for a=30 h/r=0.0395, rin=4, rpmax=5.5, rout=15 */
	//l_ang_in = 2.0246198830860416;
	create_disk_solution(MBH2,0,0.5);
      }
    }
  }
  else{ 
#if( NO_CIRC )
    pr[RHO] = coords->rhoflr;
    pr[UU ] = coords->uuflr;
    pr[U1 ] = pr[U2] = pr[U3] = 0.;
#endif
    MASS_TYPE = MTOT;   
    TRACE_END;   
    return; 
  } // Leave init_mini_disks if cell not in mini-disk


  if( in_mini_disk ) {    
    /* Convert NZ coordinates to IZ coordinates */
    xNZ[TT] = t; xNZ[XX] = x;
    xNZ[YY] = y; xNZ[ZZ] = z;
    coordsIZ_of_coordsNZ( MASS_TYPE, xNZ, xIZ );
    
    /* Convert IZ coordinates to BL coordinates */
    coordsBL_of_coordsIZ(Mass,Spin,xIZ,xBL);
    
    /* Perform look up and interpolation */
    if     ( MASS_TYPE == MBH1 ) { lookup_mini_disk( MBH1, xBL, pr, ucon_bl, Acov_bl)     ;                 }
    else if( MASS_TYPE == MBH2 ) { lookup_mini_disk( MBH2, xBL, pr, ucon_bl, Acov_bl )    ;                 }
    else                         { fprintf(stderr,"error in init_mini_disks MASS_TYPE \n"); fflush(stderr); }

    /* Solution now specified in BL coordinates
       -- Must convert back to NZ coordinates */
    double *det_tmp;
    ALLOC_ARRAY(det_tmp,1);
    
    /* Convert 4-velocity to IZ coords */
    jac_xBL_to_xIZ(Mass,Spin,xBL,dxIZ_dxBL);  //Get Jacobian
    transform_rank1con2( dxIZ_dxBL, ucon_bl ); //Transform
    /* Convert 4-velocity to NZ coords */
#if( SET_TSHRINK )
    if( t < t_shrink_bbh ) { jac_xNZ_to_xIZ_r12const(MASS_TYPE, xNZ, dxIZ_dxNZ); } // Get Jacobian 
    else                   { jac_xNZ_to_xIZ(MASS_TYPE,xNZ,dxIZ_dxNZ)           ; } // Get Jacobian    
#else
    jac_xNZ_to_xIZ(MASS_TYPE,xNZ,dxIZ_dxNZ);
#endif
    invert_matrix2( dxIZ_dxNZ, dxNZ_dxIZ, det_tmp ); //Invert for Jacobian we need
    transform_rank1con2( dxNZ_dxIZ, ucon_bl ); //Transform
    if( MASS_TYPE == MBH2 ) { ucon_bl[XX] *= -1.; ucon_bl[YY] *= -1.; } //correct Jac error
    
    /* Convert to Spherical Coordinates */
    xspher_of_xcart(coords->x, coords->xcart, dxs_dxc ); // convert to spherical coordinates    
    transform_rank1con2( dxs_dxc, ucon_bl );
    
    for(l=0;l<NDIM;l++) { ucon[l] = ucon_bl[l]; }    
    
    /* Convert 4-velocity to prims renormalize etc.*/
    transform_rank1con2(coords->dxp_dx,ucon); // transform from (r,th,ph) to (x1,x2,x3)
    get_geometry(i,j,k,CENT,ncurr,geom); /* get (x1,x2,x3) metric,   */ 
    renormalize_ucon(ucon,geom->gcov); // Renormalize the 4-velocity
    ucon2pr( pr, ucon, geom->gcon ); // get prim. velocities
    
#if( HYDRO_ONLY )
    pr[B1] = pr[B2] = pr[B3] = 0.;
#else
    /* Convert A_{\mu} to IZ coords */
    invert_matrix2( dxIZ_dxBL, dxBL_dxIZ, det_tmp ); //Get Jacobian from dxIZ_dxBL
    transform_rank1cov2( dxBL_dxIZ, Acov_bl ); //Transform

    /* Convert A_{\mu} to NZ coords */
    transform_rank1cov2( dxIZ_dxNZ, Acov_bl ); //Transform
    if( MASS_TYPE == MBH2 ) { Acov_bl[XX] *= -1.; Acov_bl[YY] *= -1.; } //correct Jac error

    /* Convert to Spherical Coordinates */
    invert_matrix2( dxs_dxc, dxc_dxs, det_tmp ); //Get Jacobian from dxs_dxc
    transform_rank1cov2( dxc_dxs, Acov_bl );

    /* Convert to numerical Grid */
    transform_rank1cov2(coords->dx_dxp, Acov_bl); //Vector potential to use for curl calculation
    
    /* Store Spacial Parts of A in B compts of pr */
    pr[B1] = Acov_bl[RR];
    pr[B2] = Acov_bl[TH];
    pr[B3] = Acov_bl[PH];
#endif

    /* Reset MASS_TYPE to be safe */
    MASS_TYPE = MTOT;
    FREE(det_tmp);
  }
 
  TRACE_END;
  return;
}

/***********************************************************************************/
/***********************************************************************************
  create_disk_solution(): 
  ------------------
   -- This is supposed to follow init_all_prim() except that we will use preset
      array sizes rather than those specified at runtime
   -- We will set the solution in BL coordinates so that we can later access 
      rho, u , u^{alpha} do not need primitive velocities here
   -- Assumes that the disk parameters have been set in init_mini_disks()

   Arguments:
      -- int disk specifies where to store the solution in memory space ( disk1 or disk2)
      -- others are same as init_all_prims() currently not implementing heigh_survey
***********************************************************************************/
void create_disk_solution(int disk, int doing_height_survey, double l_factor)
{
  int   i,j,k,l,n;
  int   j_glob_eq, j_glob, jbeg, jend;
  int   retval;

  double r,th,ph;
  double rho, uu;
  double *pr;

  double *lambda_eq, *l_eq, *h_eq;
  double xsol[1];
  double *r_loc;

  int j_glob_eq_mini,jbeg_mini,jend_mini;
  double dr_mini,dth_mini,dph_mini;
  double beta_mini;

  struct of_geom geom_eq;
  struct of_geom *geom;
  struct of_coord *coords;
  struct of_coord coords_loc;
  struct of_coord *coords_r;

  double lambda,  h, f_l, l_ang, epsilon;
  double ucon[NDIM], ucov[NDIM], ucon_bl[NDIM],ucon_floor[NDIM];
  double utcov, upcov, utcon, upcon, Omega; 

  extern void  coord_of_xp2( double *xp, struct of_coord *coords ) ;
  extern void bin_data_in_r(int bin_state, int n_in, double *data_in, double *r_beg, double *r_end, double *data_out);


  /***********************************************************************************************
    BEGIN:  Code block to integrate the equations from the equator over theta/X2: --Need to generalize for mini-disks
  **********************************************************************************************/  
#define X2_INTEGRATION_CODE_BLOCK_MINI                                                                       \
  {                                                                                                          \
    if( (0 <= j) && (j < N2TOT_MINI) )  {				                                     \
      /*r_beg[n] = r_end[n] = 1.e-20;*/					\
      th = SINGSMALL + j * dth_mini;                                                                         \
      ph = 0.;                                                                                               \
      coords_r[i].x[TH] = th;                                                                                    \
      coords_r[i].x[PH] = ph;                                                                                    \
      get_special_geometry( &(coords_r[i]), &geom_eq, metric_type);	\
      gtt = geom_eq.gcon[TT][TT];   gtp = geom_eq.gcon[TT][PH];   gpp = geom_eq.gcon[PH][PH];                \
      xsol[0] = lambda;                                                                                      \
      retval = general_newton_raphson( xsol, 1, resid_func);                                                 \
      lambda = xsol[0];                                                                                      \
      if( retval ) {                                                                                         \
	fprintf(stderr,"init.hydro_station_mini_disks.c: create_disk_solutions():  ");                       \
	fprintf(stderr,"Failed to find a solution at [i=%6d, j=%6d] [r=%26.16e] [th=%26.16e]\n", \
                i, j, coords_r[i].x[RR], coords_r[i].x[TH]);                                                 \
        fflush(stderr);                                                                                      \
      }                                                                                                      \
      else{                                                                                                  \
        l_ang  = eta * pow( lambda , (2.-q) );                                                               \
        f_l = f_l_func(l_ang);                                                                               \
        utcov = -1./sqrt(-gtt + l_ang*(2*gtp - l_ang*gpp));                                                  \
        upcov = -l_ang * utcov;                                                                              \
        utcon = gtt*utcov + gtp*upcov;                                                                       \
        upcon = gtp*utcov + gpp*upcov;                                                                       \
        h = utcov_in * f_l_in / ( utcov * f_l );                                                             \
        if( (h > 1.) && finite(h) && finite(utcon) && finite(upcon) ) {	                                     \
          epsilon = (h - 1.) / gam;                                                                          \
          rho =  pow( (epsilon * (gam - 1.) / adiab), (1./(gam-1.)) );                                       \
          uu  = rho * epsilon;                                                                               \
          if( disk == MBH1 )      { pr = pd1[i][j][k]; }                                                     \
	  else if( disk == MBH2 ) { pr = pd2[i][j][k]; }                                                     \
	  else { fprintf(stderr,"integration block create_disk_solution why are we here \n?"); }             \
          pr[RHOM] = rho;                                                                                    \
          pr[UUM ] = uu;                                                                                     \
          pr[UTTM] = utcon;   /* store the phi-independent components here for later  */                     \
	  pr[URRM] = 0.;                                                                                     \
	  pr[UTHM] = 0.;                                                                                     \
          pr[UPHM] = upcon;   /* to be transformed in a phi-dependent way later     */                       \
	  pr[XRRM] = coords_r[i].x[RR];                                                                                      \
	  pr[XTHM] = th;                                                                                     \
          rho = geom_eq.g * pr[RHOM] ;                                                                       \
          height_numer[i] += sqrt(geom_eq.gcov[2][2]) * fabs(th-0.5*M_PI) * rho; \
	  height_denom[i] += rho;					\
	  /*Solutions here are nan-free */\
	  /*fprintf(stderr,"nan in disk i,j pr? %d %d %d %e %e %e %e %e %e %e %e\n",disk,i,j,pr[RHOM],pr[UUM],pr[UTTM],pr[URRM],pr[UTHM],pr[UPHM],pr[XRRM],pr[XTHM]);fflush(stderr);*/ \
	}								                                     \
      }                                                                                                      \
      n++;                                                                                                   \
      }                                         							     \
  }                                                                                                          

  /***********************************************************************************************
    END:  Code block
  **********************************************************************************************/  

  TRACE_BEG;

  /***********************************************************************************************
    Don't Set parameters used for all points (assumes parameters set in init_mini_disks() ) : 
  **********************************************************************************************/  
  metric_type  = METRIC_BL_SPHERICAL;
  int metric_type2 = METRIC_BL_SPHERICAL;

  general_keplerian_functions(prograde, rin  , &Omega_in, &l_ang_in_min, ucov, ucon);
  general_keplerian_functions(prograde, rpmax, &Omega_p , &l_ang_p     , ucov, ucon);
  lambda_p  = sqrt( l_ang_p  / Omega_p  );

  l_ang_in_max = l_ang_p;

  if( doing_height_survey ) {
    /* for when you want to search for desired h/r : */
    l_ang_in = l_factor*l_ang_in_max + (1.-l_factor)*l_ang_in_min; 
  }


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
    fprintf(stdout,"create_disk_solution(): alpha may be getting too close to one :   %26.16e \n", alpha); 
    fail(FAIL_BASIC,0); 
  }


  /**********************************************************************************************
    Allocate memory for mini-disk solutions and initialize the primitive-like array to 0.
  **********************************************************************************************/
  if( disk == MBH1 ) { 
    LOOP_MINI PLOOP_MINI {  pd1[i][j][k][l] = 0.; }
    beta_mini = beta_disk1;
  }
  else if( disk == MBH2 ) {
    LOOP_MINI PLOOP_MINI {  pd2[i][j][k][l] = 0.; }
    beta_mini = beta_disk2;
  }
  else { fprintf(stderr,"Error in create_disk_solution() wrong disk val\n"); }

  if( myid == printer_pid ) { 
    if( !doing_height_survey ) { 
      fprintf(stdout,"\n##################################################\n");
      fprintf(stdout,"  AXI-SYM TORUS PARAMETERS \n------------------------------------\n");
      fprintf(stdout,"\t  adiab       =  %28.18e \n",adiab    ); 
      fprintf(stdout,"\t  rin         =  %28.18e \n",rin      ); 
      fprintf(stdout,"\t  rpmax       =  %28.18e \n",rpmax    ); 
      fprintf(stdout,"\t  GRscale     =  %28.18e \n",GRscale  ); 
      fprintf(stdout,"\t  beta        =  %28.18e \n",beta_mini); 
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
  }


  /**********************************************************************************************
    Set the values in the disk :
  **********************************************************************************************/  
  ALLOC_ARRAY(height_denom, N1TOT_MINI);
  ALLOC_ARRAY(height_numer, N1TOT_MINI);
  ALLOC_ARRAY(       r_loc, N1TOT_MINI);
  ALLOC_ARRAY(   lambda_eq, N1TOT_MINI); 
  ALLOC_ARRAY(        l_eq, N1TOT_MINI); 
  ALLOC_ARRAY(        h_eq, N1TOT_MINI); 
  ALLOC_ARRAY(    coords_r, N1TOT_MINI);

  N1LOOP_MINI {    lambda_eq[i] = -1. ;  } 
  N1LOOP_MINI {         l_eq[i] = -1. ;  } 
  N1LOOP_MINI {         h_eq[i] = -1. ;  } 
  N1LOOP_MINI { height_denom[i] = 0.  ;  }
  N1LOOP_MINI { height_numer[i] = 0.  ;  }

  j_glob_eq_mini = N2TOT_MINI/2;
  jbeg_mini      = 0;
  jend_mini      = N2TOT_MINI - 1;
  dr_mini        = ( rdisk_max     - rdisk_min    ) / (double)( N1TOT_MINI - 1 ) ;
  dth_mini       = ( M_PI          - 2.*SINGSMALL ) / (double)( N2TOT_MINI - 1 ) ;
  if( N3TOT_MINI == 1 ) { dph_mini = 1.e-10; }
  else                  { dph_mini = ( GridLength[3] - 0. ) / (double)( N3TOT_MINI - 1 ) ; }


  n = 0; 

  N3LOOP_MINI { 
    j = j_glob_eq_mini;   /* this is so we can use coord() to correctly find global xp[] values with local function  */
    N1LOOP_MINI { 
      /* Set up local mini-disk coordinates at equator phi coord. doesn't matter */
      coords_r[i].i   = i;
      coords_r[i].j   = j;
      coords_r[i].k   = k;
      coords_r[i].pos = CENT;
      coords_r[i].x[RR] = rdisk_min + i * dr_mini;
      coords_r[i].x[TH] = 0.5 * M_PI;
      coords_r[i].x[PH] = 0. + k * dph_mini;
      coords_r[i].r     = coords_r[i].x[RR];
    }
    calc_eq_values_mini(lambda_eq, l_eq, h_eq, coords_r);      

    if( k == 0 ) { 
      N1LOOP_MINI{ 
	fprintf(stdout,"eq-values [%d,%d,%d,%d]  %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e\n",i,j,k,myid,coords_r[i].r,coords_r[i].x[TH],coords_r[i].x[PH],lambda_eq[i],l_eq[i],h_eq[i]); 
      }
      fflush(stdout);
    }

    N1LOOP_MINI { 
      h  = h_eq[i] ;     
      if( h > 1. ) { 
	lambda = lambda_eq[i] ; 
	/* If above equator index integrate */
	for( j = j_glob_eq_mini; j < N2TOT_MINI; j++ ) { 
	    X2_INTEGRATION_CODE_BLOCK_MINI;	                                                 
	}
	
	lambda = lambda_eq[i] ; 
	/* If below equator index integrate */
	for( j = j_glob_eq_mini-1; j >= 0; j-- ) { 
	  X2_INTEGRATION_CODE_BLOCK_MINI;
	}
      }
    }
  }


  if( n > NCELLS_MINI ) { 
    fprintf(stdout,"create_disk_solutions(): mismatch with the counting:  %d %d \n", n, NCELLS); fflush(stdout); 
    fail(FAIL_BASIC,0) ; 
  }

  N1LOOP_MINI { height_numer[i] /= height_denom[i]*coords_r[i].r; } 
  N1LOOP_MINI { if( !finite(height_numer[i]) ) { height_numer[i] = 0.; } }
  N1LOOP_MINI { r_loc[i] = coords_r[i].r ; }
  h_o_r = point_interp( N1TOT_MINI, r_loc, height_numer, rpmax );

  fprintf(stdout,"mini-disk %d l_ang_in h_o_r(rpmax) q eta = %26.16e %26.16e %26.16e %26.16e \n", disk, l_ang_in, h_o_r,q,eta);
  fflush(stdout);
  N1LOOP_MINI { fprintf(stdout, "mini-disk %d r height/r = %26.16e %26.16e \n",disk,r_loc[i], height_numer[i]); }
  fflush(stdout);

  /**********************************************************************************************
    Perform a copy in upper and lower most th cells -- issues finding solutions here ( should be floor )
  **********************************************************************************************/
  N1LOOP_MINI N3LOOP_MINI PLOOP_MINI {
    if( disk == MBH1 ) { 
      pd1[i][0             ][k][l] = pd1[i][1             ][k][l];
      pd1[i][N2TOT_MINI - 1][k][l] = pd1[i][N2TOT_MINI - 2][k][l];
    }
    else if( disk == MBH2 ) {
      pd2[i][0             ][k][l] = pd2[i][1             ][k][l];
      pd2[i][N2TOT_MINI - 1][k][l] = pd2[i][N2TOT_MINI - 2][k][l];
    }
  }

  //if( doing_height_survey ) {   TRACE_END;  return; } 

  /**********************************************************************************************
    Check if solution in floor state... if not then add perturbation to UU
  **********************************************************************************************/  
  N1LOOP_MINI N2LOOP_MINI N3LOOP_MINI { 
    if( disk == MBH1 )      { pr = pd1[i][j][k]; }
    else if( disk == MBH2 ) { pr = pd2[i][j][k]; }

    rho = pr[RHOM]; 
   
    if( rho > 0. ) { 
      /* if( rho > rhomax )  { rhomax = rho ; }  */
      /* if(  uu >   umax )  {   umax =  uu ; }  */

      /* Add perturbation to internal energy density : */
      pr[UUM]  *= (1. + GRscale * (ranc(0) - 0.5));

    }
    else { 
      /* Set the "floor" state: -- we can keep the density and pressure zero, but we need to 
	 set the coordinates and 4-velocity */
      pr[XRRM] = rdisk_min + i * dr_mini;
      pr[XTHM] = SINGSMALL + j * dth_mini;
      ucon_floor[TT] = -1.;
      ucon_floor[RR] =  0.;
      ucon_floor[TH] =  0.;
      ucon_floor[PH] =  0.;
      renormalize_ucon( ucon_floor, geom_eq.gcov );
      pr[UTTM] = ucon_floor[TT];
      pr[URRM] = ucon_floor[RR];
      pr[UTHM] = ucon_floor[TH];
      pr[UPHM] = ucon_floor[PH];
      
    }
  }
  
  /**********************************************************************************************
    Set the Vector Poential for mini-disk solution 
  **********************************************************************************************/  
#if( !HYDRO_ONLY )
  if( disk == MBH1 )      { set_vect_potential_mini(pd1); }
  else if( disk == MBH2 ) { set_vect_potential_mini(pd2); }
#endif



  FREE(height_numer);
  FREE(height_denom);
  FREE(       r_loc);
  FREE(   lambda_eq); 
  FREE(        l_eq); 
  FREE(        h_eq); 
  FREE(    coords_r); 
#undef X2_INTEGRATION_CODE_BLOCK_MINI

  TRACE_END;
  return;

}


/*********************************************************************************/
/*********************************************************************************
   calc_eq_values_mini():
  ------------------

    -- calculates the relevant disk functions along the equator;
*********************************************************************************/
static void calc_eq_values_mini(double *lambda_eq, double *l_eq, double *h_eq, struct of_coord *coords_r)
{
  int i,j,k,retval;
  double Omega, l_ang, ucon[NDIM], ucov[NDIM], xsol[1],f_l,utcov;
  struct of_geom geom_eq;


  //N1ALL_LOOP { lambda_eq[i] = l_eq[i] = h_eq[i] = -1. ;  } 
  N1LOOP_MINI { lambda_eq[i] = l_eq[i] = h_eq[i] = -1. ; }
      
  //N1_LOOP { 
  N1LOOP_MINI {
    general_keplerian_functions(prograde, coords_r[i].r , &Omega, &l_ang, ucov, ucon);

    //    fprintf(stdout,"Keplerian functions = %26.16e %26.16e %26.16e %26.16e %26.16e \n",coords_r[i].r, Omega, l_ang, ucov[TT], ucon[TT]); fflush(stdout); 

    /* Need to fix this condition for mini-disks */
    if( coords_r[i].r >= rdisk_min ) { 
      if( (i==0) || (lambda_eq[i-1] < 0.) ) { 
	l_eq[i] = l_ang_in; 
	lambda_eq[i] = pow(  (l_eq[i]/eta), exp_5 );
      } 
      else { 
	lambda_eq[i] = lambda_eq[i-1]; 
      }

      get_special_geometry( &(coords_r[i]), &geom_eq, metric_type);
      gtt = geom_eq.gcon[TT][TT];
      gtp = geom_eq.gcon[TT][PH];
      gpp = geom_eq.gcon[PH][PH];

      //      r_glob = coords_r[i].r;
      //      fprintf(stdout,"gtttp %26.16e %26.16e %26.16e %26.16e \n",r_glob,gtt,gtp,gpp); fflush(stdout); 

      xsol[0] = lambda_eq[i];
      retval = general_newton_raphson( xsol, 1, resid_func);

      if( retval ) { 
	fprintf(stderr,"init_gen_disk.c: calc_eq_values():  Failed to find a solution at [i=%6d]   [r=%26.16e]  \n", i, coords_r[i].r); 
	fflush(stdout); 
      }
      else { 
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
//		coords_r[i].r,utcov, f_l, l_ang, lambda_eq[i], h_eq[i], Omega); 
//	fflush(stdout);

      }

      if( (h_eq[i] > 1.) && (coords_r[i].r > rdiskout) ) { rdiskout = coords_r[i].r ; }
      
//      fprintf(stdout,"equator-values: myid i r lambda l h = %6d %6d %26.16e %26.16e %26.16e %26.16e \n",
//	      myid, i,coords_r[i].r,lambda_eq[i],l_eq[i],h_eq[i]); 
//      fflush(stdout);


    }   /*   if( coords_r[i].r >= rin )   */
  }
  
  return;
}



/**********************************************************************/ 
/*********************************************************************************
  lookup_mini_disk():
         -- Looks up the solutions stores in pd1 and pd2
         -- Interpolates these solutions onto the numerical grid
         -- Specifies prim array and 4-vel. in bl coordinates
         -- Interpolation based off pages 132-133 of Numerical Recipes
             y_{ij} = y(x_{1i}, x_{2j}) 
                  t = (x_1 - x_{1i})/(x_{1(i+1)} -x_{1i})
                  u = (x_2 - x_{2j})/(x_{2(j+1)} -x_{2j})
           y(x1,x2) = (1-t)(1-u)y_{ij} + t(1-u)y_{(i+1)j} + tuy_{(i+1)(j+1)} + (1-t)uy_{i(j+1)}

         -- In this example x1->r,x2->th,y->(rho,uu,u^{\alpha},A_{\alpha})
     Arguments:
         disk: Determines which mini-disk the cell belongs to ( MBH1 or MBH2 )
          xBL: Coordinates of cell for mini-disk in BL coordinates
           pr: Primitive array of cell
      ucon_bl: returns 4-velocity for cell in BL coordinates for the cell
 *********************************************************************************/
void lookup_mini_disk( int disk, double xBL[NDIM],double *pr, double ucon_bl[NDIM], double Acov_bl[NDIM] )
{
  double ****prl;
  int i,j,k,icrit,jcrit;
  double rdiff, rdiff_min, thdiff, thdiff_min;
  double r,th,tt,uu;
  double ww,rho;

  TRACE_BEG;

  /* First we need to specify which solution we want to match to */
  if     ( disk == MBH1 ) { prl = pd1; }
  else if( disk == MBH2 ) { prl = pd2; }

  rdiff_min = thdiff_min = 1.e6;
  /* Let us find the index of prl most closely related to our radial coordinate --note r > r_{icrit} */
  N1LOOP_MINI {
    if( xBL[RR] > prl[i][N2TOT_MINI/2][0][XRRM] ) {
      rdiff = fabs(xBL[RR] - prl[i][N2TOT_MINI/2][0][XRRM]);
      if( rdiff < rdiff_min ) {
	rdiff_min = rdiff;
	icrit     = i;
      }
    }
  }
  /* Now search along that radial coordinate for closest theta --note th > th_{jcrit}*/
  N2LOOP_MINI {
    if( xBL[TH] > prl[icrit][j][0][XTHM] ) {
      thdiff = fabs(xBL[TH] - prl[icrit][j][0][XTHM]);
      if( thdiff < thdiff_min ) {
	thdiff_min = thdiff;
	jcrit      = j;
      }
    }
  }

  /* Need to be careful not to index outside of bounds in our stored solution */
  if( icrit == N1TOT_MINI-1 ) { icrit -= 1 ; }
  if( jcrit == N2TOT_MINI-1 ) { jcrit -= 1 ; }
  

  /* Now need to bi-linearly interpolate onto our evolution grid */
  r  = xBL[RR];
  th = xBL[TH];
  tt = ( r  - prl[icrit][jcrit][0][XRRM] ) / ( prl[icrit+1][jcrit  ][0][XRRM] - prl[icrit][jcrit][0][XRRM] );
  ww = ( th - prl[icrit][jcrit][0][XTHM] ) / ( prl[icrit  ][jcrit+1][0][XTHM] - prl[icrit][jcrit][0][XTHM] );

  
  rho = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][RHOM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][RHOM] + tt*ww*prl[icrit+1][jcrit+1][0][RHOM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][RHOM];
  if( rho > 0 ) { 
    if( ( tt < 0 ) || ( tt > 1 ) ) { fprintf(stderr,"ERROR tt < 0 OR tt > 1 icrit %d jcrit %d tt %e \n",icrit,jcrit,tt); fflush(stderr); }
    if( ( ww < 0 ) || ( ww > 1 ) ) { fprintf(stderr,"ERROR ww < 0 OR ww > 1 icrit %d jcrit %d ww %e \n",icrit,jcrit,ww); fflush(stderr); }

    uu          = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][UUM ] + tt*(1.-ww)*prl[icrit+1][jcrit][0][UUM ] + tt*ww*prl[icrit+1][jcrit+1][0][UUM ] + (1.-tt)*ww*prl[icrit][jcrit+1][0][UUM ];
    ucon_bl[TT] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][UTTM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][UTTM] + tt*ww*prl[icrit+1][jcrit+1][0][UTTM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][UTTM];
    ucon_bl[RR] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][URRM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][URRM] + tt*ww*prl[icrit+1][jcrit+1][0][URRM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][URRM];
    ucon_bl[TH] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][UTHM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][UTHM] + tt*ww*prl[icrit+1][jcrit+1][0][UTHM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][UTHM];
    ucon_bl[PH] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][UPHM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][UPHM] + tt*ww*prl[icrit+1][jcrit+1][0][UPHM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][UPHM];
  }
  else { 
    rho = 0;
    uu  = 0;
    ucon_bl[TT] = prl[icrit][jcrit][0][UTTM];
    ucon_bl[RR] = prl[icrit][jcrit][0][URRM];
    ucon_bl[TH] = prl[icrit][jcrit][0][UTHM];
    ucon_bl[PH] = prl[icrit][jcrit][0][UPHM];
  }


#if( HYDRO_ONLY )
  pr[B1] = 0.;
  pr[B2] = 0.;
  pr[B3] = 0.;
#else
  /* /\*Need to figure out vector potential stuff*\/ */
  /* fprintf(stderr,"INITIALIZE_MINI_DISKS NOT COMPATIBLE WITH MHD\n"); */
  /* fflush(stderr); */
  /* fail(FAIL_BASIC,0); */

  if( rho > 0 ) {
    if( ( tt < 0 ) || ( tt > 1 ) ) { fprintf(stderr,"ERROR tt < 0 OR tt > 1 icrit %d jcrit %d tt %e \n",icrit,jcrit,tt); fflush(stderr); }
    if( ( ww < 0 ) || ( ww > 1 ) ) { fprintf(stderr,"ERROR ww < 0 OR ww > 1 icrit %d jcrit %d ww %e \n",icrit,jcrit,ww); fflush(stderr); }

    Acov_bl[TT] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][ATTM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][ATTM] + tt*ww*prl[icrit+1][jcrit+1][0][ATTM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][ATTM];
    Acov_bl[RR] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][ARRM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][ARRM] + tt*ww*prl[icrit+1][jcrit+1][0][ARRM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][ARRM];
    Acov_bl[TH] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][ATHM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][ATHM] + tt*ww*prl[icrit+1][jcrit+1][0][ATHM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][ATHM];
    Acov_bl[PH] = (1.-tt)*(1.-ww)*prl[icrit][jcrit][0][APHM] + tt*(1.-ww)*prl[icrit+1][jcrit][0][APHM] + tt*ww*prl[icrit+1][jcrit+1][0][APHM] + (1.-tt)*ww*prl[icrit][jcrit+1][0][APHM];
  }
  else {
    Acov_bl[TT] = prl[icrit][jcrit][0][ATTM];
    Acov_bl[RR] = prl[icrit][jcrit][0][ARRM];
    Acov_bl[TH] = prl[icrit][jcrit][0][ATHM];
    Acov_bl[PH] = prl[icrit][jcrit][0][APHM];
  }
#endif

  /* Set primitive array */
  pr[RHO] = rho;
  pr[UU ] = uu;

  TRACE_END;
  return;
}


/**********************************************************************/ 
/*********************************************************************************
  renormalize_ucon(): 
         --Takes a contravariant 4-velocity and 4-metric then renormalizes
         --Solves quadratic equation for ucon[TT] which has normalized 4-velocity
         --Rejects solutions with ucon[TT] < 0.
         --Takes solution which is closes to current ucon[TT]
         -- g_{tt}u^t u^t + 2 g_{t i} u^t u^i + g_{ij} u^i u^j = -1
            -- A = g_{tt}, B = 2g_{t i} u^i, C = g_{ij}u^i u^j + 1
     Arguments:
         ucon: 4-velocity to renormalize
         gcov: Metric to normalize wrt
 *********************************************************************************/
void renormalize_ucon( double ucon[NDIM], double gcov[NDIM][NDIM] )
{
  int mu,nu,i,j;
  double old_ucontt, sol1, sol2;
  double AA,BB,CC;

  /* Store orignal value before renormalizing */
  old_ucontt = ucon[TT];
  
  /* Coefficients for quadratic equation */
  AA = BB = 0.;
  CC = 1.;
  AA = gcov[TT][TT];
  for( i = 1; i < NDIM; i++ ) { BB += 2. * gcov[TT][i] * ucon[i] ; }
  for( i = 1; i < NDIM; i++ ) for( j = 1; j < NDIM; j++ ) { CC += gcov[i][j] * ucon[i] * ucon[j] ; }

  sol1 = ( -BB + sqrt(BB*BB - 4.*AA*CC) ) / ( 2. * AA );
  sol2 = ( -BB - sqrt(BB*BB - 4.*AA*CC) ) / ( 2. * AA );

  /* Reject negative ucon[TT] solutions */
  if( sol1 <= 0 )      { ucon[TT] = sol2;                                                      }
  else if( sol2 <= 0 ) { ucon[TT] = sol1;                                                      }
  else                 { ucon[TT] = fmin( fabs( ucon[TT] - sol1 ) , fabs( ucon[TT] - sol2 ) ); }

  if( ucon[TT] <= 0 ) {
    fprintf(stderr," WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING renormalize_ucon found only ucon[TT] < 0 \n");
    ucon[TT] = old_ucontt; // Set back to old value just to be safe... really should kill run if this happens though 
  }
    
  return;
}

/***********************************************************************************/
/***********************************************************************************
  set_vect_potential_mini():
  ------------------
   -- Sets the vector potential A_{\mu} in stored disk solution ( pd1 || pd2 ); 
   -- Implicitly assume that N3TOT_MINI == 0!!!
   -- Currently only supports poloidal fields derived from the azimuthal component 
      of the vector potential; 
   -- Before the calculation, does a fixup() and sets the boundaries since 
      the derivatives of A_\phi require ghost cells (field follows density contours); 
   -- We assume below that the magnetic field is independent of X3 ; 
***********************************************************************************/
static void set_vect_potential_mini( double ****prl )
{
  int i,j,k,l;

  double rho_avg, rho_cut, Aphi, beta_act, norm, bsq_int, u_int ; 
  double dr_mini, dth_mini, dph_mini;
  double rcorn, thcorn;
  double rhomax_loc;
  struct of_geom *geom;  
  struct of_coord *coords;  
  extern void curl_of_A( double ****A , double ****prim ) ;

  TRACE_BEG;

  if( myid == printer_pid ) { 
#if(   FIELD_CONFIG_CHOICE == QUADRUPOLE_FIELD )
    fprintf(stdout,"set_vect_potential_mini():  Using QUADRUPOLE_FIELD %d \n",FIELD_CONFIG_CHOICE); 
#elif( FIELD_CONFIG_CHOICE ==     DIPOLE_FIELD )
    fprintf(stdout,"set_vect_potential_mini():  Using DIPOLE_FIELD %d \n",FIELD_CONFIG_CHOICE); 
#else 
    fprintf(stderr,"set_vect_potential_mini():  Invalid value of FIELD_CONFIG_CHOICE =  %d \n",
	    FIELD_CONFIG_CHOICE);  fflush(stderr); 
    fail(FAIL_BASIC,0); 
#endif
    fflush(stdout);
  }

  /*************************************************************************************
    Calculate the vector potential, A_mu: 
       --  Need to calc. A_\mu on the last ghost cells or all the corners surrounding the 
           physical cells' centers for finite diff. calculation for B^i 
       -- We first calculate A_\phi then we transform to xp;
       -- To this end we implicitly assume that the density distribution is 
           phi-independent (but not necessaryily xp3 independent);
  *************************************************************************************/
  double Acov[NDIM];

  //for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) for(k=N3S; k<=(N3E+1); k++) { 
  /* Set k=0 assuming N3TOT_MINI = 1 and reset dr_mini,dth_mini to find coords at corner */
  dr_mini        = ( rdisk_max     - rdisk_min    ) / (double)( N1TOT_MINI - 1 ) ;
  dth_mini       = ( M_PI          - 2.*SINGSMALL ) / (double)( N2TOT_MINI - 1 ) ;
  if( N3TOT_MINI == 1 ) { dph_mini = 1.e-10; }
  else                  { dph_mini = ( GridLength[3] - 0. ) / (double)( N3TOT_MINI - 1 ) ; }
  
  /* Find rhomax in the stored solution --may not be necessary*/
  rhomax_loc = 0.;
  LOOP_MINI { if( prl[i][j][k][RHOM] > rhomax_loc ) { rhomax_loc = prl[i][j][k][RHOM] ; } }

  k = 0;
  for( i = 1; i < N1TOT_MINI; i++) for( j = 1; j < N2TOT_MINI; j++ ) {

	rho_avg =   
	  prl[i  ][j  ][k][RHOM] +
	  prl[i-1][j  ][k][RHOM] +
	  prl[i  ][j-1][k][RHOM] +
	  prl[i-1][j-1][k][RHOM] +
	  prl[i  ][j  ][k][RHOM] +
	  prl[i-1][j  ][k][RHOM] +
	  prl[i  ][j-1][k][RHOM] +
	  prl[i-1][j-1][k][RHOM] ;
	
	rho_avg *= 0.125;

	//get_coord(i,j,k,CORN,ncurr,coords); // coords on corner?
	/* Coords at bottom left corner */
	rcorn  = prl[i][j][k][XRRM] - 0.5*dr_mini;
	thcorn = prl[i][j][k][XTHM] - 0.5*dth_mini;

#if(   FIELD_CONFIG_CHOICE ==     DIPOLE_FIELD )
	rho_cut = 0.25 * rhomax_loc; 
	Aphi = rho_avg - rho_cut ;    /* dipole field */

#elif( FIELD_CONFIG_CHOICE == QUADRUPOLE_FIELD )
	rho_cut = 0.25 * rhomax_loc; 
	Aphi = (rho_avg - rho_cut) * cos(thcorn) / rhomax_loc ;
#endif

	/* Truncate A_\mu at zero so that field lines lie within disk */
	Acov[0] = 0.;
	Acov[1] = 0.;
	Acov[2] = 0.;
	Acov[3] = 0.;

	if( rho_avg > rho_cut ) { 
	  Acov[3] = Aphi;
	  //transform_rank1cov2(coords->dx_dxp,Acov);
	}
	prl[i][j][k][ATTM] = Acov[0];
	prl[i][j][k][ARRM] = Acov[1];
	prl[i][j][k][ATHM] = Acov[2];
	prl[i][j][k][APHM] = Acov[3];
    }
  /* We couldn't set rho_avg for i or j = 0 
     So copy next result for A_{\mu} */
  N2LOOP_MINI {
    prl[0][j][k][ATTM] = prl[1][j][k][ATTM];
    prl[0][j][k][ARRM] = prl[1][j][k][ARRM];    
    prl[0][j][k][ATHM] = prl[1][j][k][ATHM];
    prl[0][j][k][APHM] = prl[1][j][k][APHM];
  }
  N1LOOP_MINI {
    prl[i][0][k][ATTM] = prl[i][1][k][ATTM];
    prl[i][0][k][ARRM] = prl[i][1][k][ARRM];    
    prl[i][0][k][ATHM] = prl[i][1][k][ATHM];
    prl[i][0][k][APHM] = prl[i][1][k][APHM];
  }


  return;
}

#endif
