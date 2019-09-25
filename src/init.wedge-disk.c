
#include "decs.h"
#include "metric.h"

static void init_all_prims(void);

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
	GridLength[0] = 5.e4 ;         /* Length of X0 dimension (evolution period) */ 
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
	//GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
	//GridLength[3] = 0.5*M_PI ;                 /* Length of X3 dimension */ 
	GridLength[3] = 2.*M_PI/1000. ;                 /* Length of X3 dimension */ 
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
	DT_out[OUT_HDF5]    = 1.e-6 ;                     /* hdf frequency */
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


  TRACE_BEG;

  /*******************************************************************************
     Set up utilities:
  *******************************************************************************/
  init_all_prims(); 
  
  /*******************************************************************************
    Set the magnetic field :  
  *******************************************************************************/
  /* Correct bad points and setup boundary values since we will require them for B^i */
  fixup(p) ;
  bounds(p,0) ;

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
  t_shrink_bbh = 0.; 
  phi_0_bbh    = 0.;
#endif

#else 
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

  //  h_slope   = 0.13;           /* Severity of the focusing             */
  h_slope   = (1.- 4.*th_cutout/M_PI) ;           /* Severity of the focusing             */
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
  diag3_exponent = 3;                            /* note that this should be an odd integer greater than 1  */
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
void init_all_prims(void)
{
  int   i,j,k,l,n;
  double dthloc;
  double *prim; 
  struct of_coord *coords;

  LOOP {
    prim = p[i][j][k]; 
    get_coord(i,j,k,CENT,ncurr,coords);
    dthloc = fabs( coords->x[TH] - 0.5 * M_PI ) ;
    if( (dthloc < 0.1)  && (coords->x[RR] > r_isco) ) {
      prim[RHO] = 1.e-2;
      prim[ UU] = 1.e-2;
    }
    else {
      prim[RHO] = 0.;
      prim[ UU] = 0.;
    }
    prim[U1] = 0. ; 
    prim[U2] = 0. ; 
    prim[U3] = 0. ; 
    prim[B1] = 0. ; 
    prim[B2] = 0. ; 
    prim[B3] = 0. ; 
  }
  
  return;

}
