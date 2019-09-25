
#include "decs.h"
#include "defs.h"

#include <time.h>

static void set_initial_globals(void) ;
static void set_arrays( void ) ;
static void global_init( void );
static void macro_warnings(void); 

extern void check_for_crash( void );
#if(MAKE_TIMERS)
extern void dump_timers(void) ;
#endif

int main(int argc,char *argv[])
{
  int i,j,k,l,d;

  time_t  time1, time2; 
  double nsecs;

  set_initial_globals();

#if( ALLOW_NEGATIVE_DENSITIES )
	treat_floor_as_failure = 0;
#else
	//	treat_floor_as_failure = 1;
	treat_floor_as_failure = 0;
#endif

	/* Initialize MPI structures, domain decomp., set other global var's */
	setup_mpi(argc,argv);

	macro_warnings(); 

	time1 = time(NULL);  
	if( myid == printer_pid ) { fprintf(stdout,"TIME: start time %d = %s \n",myid, asctime(localtime(&time1))); fflush(stdout); }

	set_arrays(); 

#if(USE_NUMREL_DATA)
   read_numrel_data();
   test_numrel_data();
#endif

	/* perform initializations */
	if( !restart_init() ) { 
	  init() ;
	  t = t_old = startx[0];

	  diag(OUT_INIT) ;	  /* do initial diagnostics */
#if( 1 || (COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL) )
     OUT_LOOP  diag(i); 
#endif
	}

	global_init();  /* set various global var's based on init()/restart_init() */


	//	myexit(0);  /* used to quit out sometimes to tset initial data */

	while(t < GridLength[0]) {
		BEG_TIMING(TIMER_TYPE_TOTAL);

		time2 = time(NULL);  nsecs = difftime(time2,time1);

  	        if( myid == printer_pid ) { 
		  fprintf(stdout,"timestep %10d %20g %16.6e  %16.6e  :  %26.16e %26.16e %26.16e \n",
			  nstep,nsecs,t,dx[0],M_tot,E_tot,L_tot) ;
		  fflush(stdout); 
		}

#if( RUNTIME_SECONDS > 0 )
		sync_val(&nsecs);
		if( nsecs > RUNTIME_SECONDS ) { break; }
#endif

		/* step variables forward in time */
		step_ch() ;

		nstep++ ;

		/* Output  events  periodic in time : */
		BEG_TIMING(TIMER_TYPE_DIAG);
		diag(OUT_EVERY);  /* Needs to be first since some dump*() routines reset grid functions */
		OUT_LOOP  if(t >= T_out[i])  diag(i); 
		END_TIMING(TIMER_TYPE_DIAG);


		/* Write restart file every n_restart time steps : */
		if( (nstep % n_restart) == 0 ) 	restart_write() ;

//		time2 = time(NULL);  fprintf(stdout,"TIME: curr  time   = %s \n", asctime(localtime(&time2))); 
		fflush(stdout);


		/*  See if run has reached a bad state  */ 
		check_for_crash(); 

		exit_status();   /* Make sure that all nodes are still running */

		//--testing
		//		myexit(0);
		//		if( nstep == 20 ) break;


		END_TIMING(TIMER_TYPE_TOTAL);

#if(MAKE_TIMERS)
		if( (nstep % n_timer_interval) == 0 ) { dump_timers(); } 
#endif

	}

	if( myid == printer_pid ) { fprintf(stdout,"ns,ts: %d %ld\n",nstep,((long int) nstep*NCELLS)) ; fflush(stdout); }

	/* do final diagnostics */
	restart_write() ;
	diag(OUT_FINAL) ;

#if(USE_NUMREL_DATA)
   free_numrel_data();
#endif

	myexit(0) ;
}

/******************************************************************************************/
/******************************************************************************************
  set_arrays(): 
 --------------
      -- Initializes all arrays to zero , and initializes recon_type[] to default 
         reconstruction method specified by RECON_TYPE_CHOICE ;
      -- should not require anything to be set (e.g. via init()) 
           but the pre-compiler macros ; 
******************************************************************************************/
void set_arrays( void ) 
{ 
  int i,j,k,d,l,g; 
  unsigned long int n ;
  extern void alloc_grid_arrays(void) ;

  TRACE_BEG;

  alloc_grid_arrays(); 

  set_levi_civita();

  t_old = t = 0.;

  x0_bh = y0_bh = z0_bh = 0.;


  for(l=0; l<N_NFAIL; l++) { nfailtot[l] = 0; }

  ALL_LOOP PLOOP {        p[i][j][k][l] = 0.; }
  ALL_LOOP PLOOP {       ph[i][j][k][l] = 0.; }
  ALL_LOOP PLOOP {  U_gf[0][i][j][k][l] = 0.; }
  ALL_LOOP PLOOP {  U_gf[1][i][j][k][l] = 0.; }

    //   ALL_LOOP PLOOP { p_old[i][j][k][l] = 0.; }

  ALL_LOOP FACE_LOOP PLOOP  {  p_L[i][j][k][d][l] = 0.; } 
  ALL_LOOP FACE_LOOP PLOOP  {  p_R[i][j][k][d][l] = 0.; } 
  ALL_LOOP FACE_LOOP PLOOP  {    F[i][j][k][d][l] = 0.; } 
  ALL_LOOP FACE_LOOP        { c_gf[i][j][k][d][0] = c_gf[i][j][k][d][1] = 0.; }

#if( CALC_CURRENT ) 
  ALL_LOOP {
    for(l=0;l<NDIM;l++)  for(d=0;d<NDIM;d++) { 
      faraday[0][i][j][k][l][d] = faraday[1][i][j][k][l][d] = 0.;
    }
    for(l=0;l<NDIM;l++) jcon[i][j][k][l] = 0.;
  }
#endif 

  ALL_LOOP { pflag[i][j][k] = 0; }

  for(l=0; l<N_NFAIL; l++) ALL_LOOP { nfail[l][i][j][k] = 0; }

#if( USE_LOCAL_RECON_TYPE )
  ALL_LOOP FACE_LOOP { recon_type[i][j][k][d] = RECON_TYPE_CHOICE ; }
#endif 


  /* Set names of directories in which we will store data : */
//  sprintf(DIR_out[OUT_HISTORY], "history");
//  sprintf(DIR_out[OUT_ASCII  ], "dumps"  );
//  sprintf(DIR_out[OUT_SDF    ], "dumps"  );
//  sprintf(DIR_out[OUT_HDF5   ], "dumps"  );
//  sprintf(DIR_out[OUT_IMAGE  ], "images" );
//  sprintf(DIR_out[OUT_STAT   ], "dumps" );
//  sprintf(DIR_out[OUT_STAT2  ], "dumps" );
//  sprintf(DIR_out[OUT_RADFLUX], "dumps" );

  sprintf(DIR_out[OUT_HISTORY    ], ".");
  sprintf(DIR_out[OUT_ASCII      ], ".");
  sprintf(DIR_out[OUT_SDF        ], ".");
  sprintf(DIR_out[OUT_HDF5       ], ".");
  sprintf(DIR_out[OUT_IMAGE      ], ".");
  sprintf(DIR_out[OUT_STAT       ], ".");
  sprintf(DIR_out[OUT_STAT2      ], ".");
  sprintf(DIR_out[OUT_RADFLUX    ], ".");
  sprintf(DIR_out[OUT_PHOTOSPHERE], ".");
  sprintf(DIR_out[OUT_MIN_DT     ], ".");
  sprintf(DIR_out[OUT_SURFACE    ], ".");

  if( (totalsize[2] != 1) && EQUATORIAL_RUN ) { 
    fprintf(stdout,"set_arrays(): in order to use EQUATORIAL_RUN consistently, you must have totalsize[2] = 1 :  %d %d \n",
	    totalsize[2], EQUATORIAL_RUN); fflush(stdout);  
    fail(FAIL_BASIC,0);
  }


#if( MAKE_TIMERS ) 
  n_timer_interval = 20;
  for(i=0;i<N_TIMER_TYPES;i++) { elapsed_times[i][0] = elapsed_times[i][1] = 0.; }
#endif


  TRACE_END;

  return;
}


/**********************************************************************************************/
/**********************************************************************************************
  free_global_arrays():
 ------------------
    -- frees all the global, dynamically allocated arrays ;
**********************************************************************************************/
void free_global_arrays(void)
{
  int i,j,k,d,l,g; 
  time_t  time1;

  TRACE_BEG;

#if( MAKE_STAT && DUMP_ALL_STAT ) 
  PLOOP { FREE( U_out[l] ); }
  PLOOP { FREE( U_pre[l] ); }
  NPH_LOOP for(g=0; g<N_STAT; g++)  { FREE( dU_stat[l][g] );  }
#endif

#if( MAKE_STAT2 ) 
  FACE_LOOP  PLOOP { FREE( FL_stat2[d][l] ); } 
  FACE_LOOP  PLOOP { FREE( FR_stat2[d][l] ); } 
  FACE_LOOP  PLOOP { FREE( UL_stat2[d][l] ); } 
  FACE_LOOP  PLOOP { FREE( UR_stat2[d][l] ); }
  FACE_LOOP        { FREE( ctop_stat2[d]  ); }
  for(l=0;l<2;l++) { FREE( S_stat2[l]     ); } 
#endif 

#if( MAKE_RADFLUX ) 
  for(l=0;l<NDIM;l++) { FREE( coolflux[l]     ); } 
#endif

#if( COORD_TYPE_CHOICE == COORD_DIAGONAL2 )
  dealloc_diag2();
#endif

  time1 = time(NULL);  
  if( myid == printer_pid ) {  fprintf(stdout,"TIME: final time %d = %s \n",myid, asctime(localtime(&time1)));  fflush(stdout); }


  TRACE_END;

  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  global_init():
 ------------------
    -- sets various global variables and arrays that should be after init() or 
          restart_init() has been called; 
    -- this should hold all the miscellaneous operations that need to be done 
         from main() before integration starts;
**********************************************************************************************/
void global_init( void )
{
  int i,j,k,d,l,g; 
  struct of_geom *geom;
  struct of_state q;

  TRACE_BEG;

#if( KEEP_CONSERVED_VARS ) 
  /* Need to set the conserved variable array after init(): */
  LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); get_state(  p[i][j][k], geom, &q );
    primtoflux( p[i][j][k], &q, 0, geom, U_gf[0][i][j][k]) ;
  }
#endif

  gam_m1_o_gam = (gam-1.)/gam; 


  /* Warn if M!=1. as some routines may assume this */
  if( myid == printer_pid ) { 
    if( REL_DIFF_FUNC(M,1.) > SMALL ) { 
      fprintf(stderr,"\nWARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n");
      fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n");
      fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
      fprintf(stderr,"   M  is assumed to be unity several places, use non-unity value at your own risk   M=%26.16e  !!! \n\n",M);
      fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n");
      fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n");
      fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
      fflush(stderr);
    }
  }



/* Clean up divergence constraint violations (aka magnetic monopoles) before we
   start integrating in time: */
#if( CLEAN_MONOPOLES )

  diag(OUT_HDF5); 

  double restot           = 1.e10 ;

  const double tol_CORN   = 1.e-8 ;
  const int itmax_CORN    = 50000 ;
  /* const int itmax_CORN    = 500 ; */
  const int global_it_max = 2000000 ;

  char h5filename[50];

  int it = 0 ;
  while( restot > tol_CORN && it < global_it_max ) {

    restot =  clean_monopoles(tol_CORN, itmax_CORN);

    sprintf(h5filename, "%s/monopole_cleaner_it%09d.h5", DIR_out[OUT_HDF5], it);
    dump_monopole_cleaner(h5filename) ;

    it += itmax_CORN ;
  }


#endif


#if( INTERPOLATE_DATA )
  //  extern void read_and_interp_all_funcs(void);
  //  read_and_interp_all_funcs();

  extern void read_and_interp_all_funcs_parallel(void);
  read_and_interp_all_funcs_parallel();

# if( !RUN_FROM_INTERP_DATA )
  DEBUG_OUT("  ending run from here because we do not want to run with interpolated data");
  myexit(0);  /* exit if we do not wish to run with interpolated data, we just wanted to perform the interpolation */
# endif

#endif


  TRACE_END;

  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  macro_warnings():
 ------------------
   --  report possible problems with compile-time macro definitions;
**********************************************************************************************/
void macro_warnings(void)
{

  if( myid == printer_pid ) { 

#if( (METRIC_TYPE_CHOICE == METRIC_GENERAL_STATIC) || (METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC) )
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    fprintf(stderr,"FLOOR IS GENERALLY GAUGE DEPENDENT ---  BE CAREFUL !!! \n\n");
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    fflush(stderr);
#endif

#if( USE_COOLING_FUNCTION ) 
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    fprintf(stderr,"COOLING FUNCTION IS GAUGE DEPENDENT --- BE CAREFUL !!! \n\n");
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    fflush(stderr);
#endif

#if( TOP_TYPE_CHOICE == TOP_CARTESIAN ) 
#if( (BC_TYPE_CHOICE != BC_CARTESIAN_OUTFLOW) &&	\
     (BC_TYPE_CHOICE != BC_CARTESIAN_PERIODIC) &&	\
     (BC_TYPE_CHOICE != BC_STATIC_ALL        ) &&	\
     (BC_TYPE_CHOICE != BC_CARTESIAN_REFLECT) )
    fprintf(stderr,"ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n\n");
    fprintf(stderr," CARTESIAN COORDINATES MUST HAVE CARTESIAN BOUNDARY CONDITIONS!!!!!  \n"); 
    fprintf(stderr," TOP_TYPE_CHOICE = %d        BC_TYPE_CHOICE = %d  \n", TOP_TYPE_CHOICE, BC_TYPE_CHOICE); 
    fprintf(stderr,"ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n\n");
#endif  
#endif

#if( TOP_TYPE_CHOICE == TOP_SPHERICAL ) 
#if( (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW      ) &&	\
     (BC_TYPE_CHOICE != BC_CYLINDRICAL_OUTFLOW    ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_CUTOUT       ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_CONSINTERP   ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW2     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_CUTOUT2      ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_CUTOUT3      ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW3     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW4     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW5     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW6     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW7     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW8     ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_OUTFLOW9     ) &&	\
     (BC_TYPE_CHOICE != BC_STATIC_ALL             ) &&	\
     (BC_TYPE_CHOICE != BC_SPHERICAL_STATIC       ) )
    fprintf(stderr,"ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n\n");
    fprintf(stderr," SPHERICAL COORDINATES MUST HAVE SPHERICAL BOUNDARY CONDITIONS!!!!!  \n"); 
    fprintf(stderr," TOP_TYPE_CHOICE = %d        BC_TYPE_CHOICE = %d  \n", TOP_TYPE_CHOICE, BC_TYPE_CHOICE); 
    fprintf(stderr,"ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n\n");
#endif  
#endif

#if( TOP_TYPE_CHOICE == TOP_CYLINDRICAL ) 
#if( (BC_TYPE_CHOICE != BC_CYLINDRICAL_OUTFLOW) )
    fprintf(stderr,"ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n\n");
    fprintf(stderr," CYLINDRICAL COORDINATES MUST HAVE CYLINDRICAL BOUNDARY CONDITIONS!!!!!  \n"); 
    fprintf(stderr," TOP_TYPE_CHOICE = %d        BC_TYPE_CHOICE = %d  \n", TOP_TYPE_CHOICE, BC_TYPE_CHOICE); 
    fprintf(stderr,"ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n\n");
#endif  
#endif

#if( USE_LIGHT_SPEED && !(MAKE_MIN_DT) ) 
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    fprintf(stderr,"IF YOU REALLY WANT TO USE THE LIGHT SPEED TO SET TIME STEP YOU SHOULD TURN ON  MAKE_MIN_DT  !!! \n\n");
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    dtmin_glob = 1.e200;   /* set it to be big so that it is always larger than characteristic MHD speed, i.e. revert to MHD method */
#endif

#if( FAST_AND_FURIOUS ) 
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
    fprintf(stderr," WE ARE RUNNING  WITH  --->  FAST_AND_FURIOUS  <----  TURNED ON    !!! \n\n");
    fprintf(stderr,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING \n\n");
#endif

    fflush(stderr);
  }

/******************************************************************************************  
  Macro logic verification. 
********************************************************************************************/

#if( USE_MPI_IO && !(USEMPI) )
  USE_MPI_IO should only be set when USEMPI is set --dont-compile
#endif  

#if( USEMPI && USE_MPI_IO_SURF && MAKE_SURFACE )
  MPI-IO not yet setup for MAKE_SURFACE  --dont-compile
#endif  

#if( USE_MASK ) 
# if( BUFFER_WIDTH < 1 )
  BUFFER_WIDTH-is-too-small
# endif

# if( (BUFFER_WIDTH != 3) && USE_FLUX_CT_PARA )
  BUFFER_WIDTH-is-wrong
# endif

# if( (BUFFER_WIDTH != 2) && (!USE_FLUX_CT_PARA) )
  BUFFER_WIDTH-is-wrong
# endif
#endif 

#if((METRIC_TYPE_CHOICE==METRIC_GENERAL_STATIC || METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC)&&(TOP_TYPE_CHOICE==TOP_CYLINDRICAL)) 
   general-metrics-not-compatible-with-cylindrical-topology-yet
#endif

#if((METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG)&&(TOP_TYPE_CHOICE==TOP_CARTESIAN)) 
   phi-avg-metric-cannot-be-cartesian 
#endif


#if(USE_SIMPLE_EOS && USE_ENTROPY_EQ )
  cannot-have-both-simple-eos-and-entropy-eq
#endif

#if(USE_NUMREL_DATA && !USE_GSL )
  need-gsl-with-USE_NUMREL_DATA-----wont-compile
#endif

#if( (METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG2) && (COORD_TYPE_CHOICE != COORD_IDENTITY))
  this-metric-does-not-transform-to-xp-coordinates-so-it-will-not-evolve-correctly------wont-compile
#endif

#if(USE_SIMPLE_EOS )
  -need-to-fix-Katm-assignment-for-warped-coordinates
#endif

#if(METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC)
#if( !(SET_TSHRINK) && !((METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ) ||(METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_DROP)|| (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND)|| (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND) || (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP) )) 
  set-tshrink-switch-off-currently-only-supported-for-METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ-and-METRIC_DYNAMIC_FULLPN_NZ_DROP-and-METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND-and-METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND-and-METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP
#endif
#endif


#if( (USE_NUMREL_DATA) && !((METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC) && (METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NR_FZ)))
  wrong-settings-for-USE_NUMREL_DATA
#endif

#if( DYNAMIC_COORDINATES && !DYNAMIC_SPACETIME )
  dont-know-how-we can have dynamic coordinates without a dynamic metric --dont-compile	
#endif

#if( USE_SIMPLE_EOS ) 
# if(						\
     (COORD_TYPE_CHOICE != COORD_IDENTITY ) &&	\
     (COORD_TYPE_CHOICE != COORD_DIAGONAL ) &&	\
     (COORD_TYPE_CHOICE != COORD_DIAGONAL2) &&	\
     (COORD_TYPE_CHOICE != COORD_DIAGONAL3) )
   --cannot-use-nondiagonal-coord-system-with---simple_eos----wont-compile
# endif
#endif

#if( USE_MASK && (METRIC_TYPE_CHOICE != METRIC_GENERAL_STATIC) && (METRIC_TYPE_CHOICE != METRIC_GENERAL_DYNAMIC) && (METRIC_TYPE_CHOICE != METRIC_GENERAL_PHI_AVG) && (METRIC_TYPE_CHOICE != METRIC_GENERAL_PHI_AVG2) )
  the-mask-is-set-in-set_general_geometry-so-it-will-not-be-set
#endif

#if( USE_COOLING_FUNCTION && (TOP_TYPE_CHOICE != TOP_SPHERICAL) )
 non-spherical cooling function not implemented yet. 
#endif


  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  set_initial_globals():
 ------------------
   --  sets initial values of global variables before anything else is called; 
**********************************************************************************************/
void set_initial_globals(void) 
{

  nstep =  ncurr = 0 ;
  t_old = t = 0.;

  M_tot = E_tot = L_tot = 0.;
  M = 1.;

  x0_bh = y0_bh = z0_bh = 0.;
  m_bh1 = m_bh2 = initial_bbh_separation = 0.;
  t_shrink_bbh = 1.e200; 
  phi_0_bbh    = 0.;

  sibling_failed = failed = 0; 
  boundary_mpi_pflag = boundary_phys_pflag = failure_exists = 0;

  N_hist_dump = 0; 
  rdump_cnt = 0; 
  using_restart = 0; 

  dt_conn = 9.e-6;   half_inv_dt_conn = 0.5/dt_conn;
	
  n_r_bins = n_phi_bins = 0; 

  return;
}

