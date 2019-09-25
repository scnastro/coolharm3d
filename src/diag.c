
#include "decs.h"
#include "metric.h"

/*
#include <sys/stat.h>
#include <sys/types.h>
*/


/* List of output routines, should have the same call structure */
extern void image(void);
extern void dump_ascii(void);
extern void dump_hdf5(void);
extern void dump_history(void);
extern void dump_surface(void);
extern void gdump(void);
extern void dump_sdf(void);
extern void dump_stat(void);
extern void dump_stat2(void);
extern void dump_radflux(void);
extern void dump_photosphere(void);
extern void calc_dU_stat(int fail_type, int state);
extern void dump_min_dt(void);

void dump_ener(int call_code);
void dump_bbh_trajectories(int call_code);

/* Array of pointers to these functions, so that we can make the code simpler, 
   this needs to be in the order of the OUT_* macros in decs.h */
const int mask[]       = {MAKE_HISTORY, MAKE_ASCII, MAKE_SDF, MAKE_HDF5, MAKE_IMAGE, MAKE_STAT, MAKE_STAT2, MAKE_RADFLUX, MAKE_PHOTOSPHERE, MAKE_MIN_DT, MAKE_SURFACE };
void (*dump[])( void ) = {dump_history, dump_ascii, dump_sdf, dump_hdf5, image, dump_stat, dump_stat2, dump_radflux, dump_photosphere, dump_min_dt, dump_surface } ;

/******************************************************************************************/
/******************************************************************************************
  diag(): 
  -------

    -- routine that manages the calls to the various output routines; 
    -- It is responsible for setting the output directories and updating the 
        time of the next output for each type; 

    -- Note: assigns the output directory for restart files, but IS NOT responsible
       for the calls to the restart dump routines;  that is done in main() since 
       it is best to dump restart files between a certain number of time steps. 
     
******************************************************************************************/

void diag(int call_code)
{
  int i;
  char sys_command[40] ;

  TRACE_BEG;


  stat_call_code = call_code; 

  
  switch( call_code ) { 

    /******************************************************************************
      Perform initial configuration of output : 
    *******************************************************************************/
  case OUT_INIT : 

    /* initialize counters and open directories : */
    OUT_LOOP if( mask[i] ) { 
      N_out[i] = 0; 
      T_out[i] = t + SMALL ; 
      //      sprintf(sys_command,"mkdir -p %s ",); 
      //      mkdir(DIR_out[i], S_IRWXU); /* Create directory with only user RWX access */
      //      system(sys_command);
    }

    gdump();   /* dump metric functions once only at the beginning */

    OUT_LOOP if( mask[i] ) { 
      fprintf(stdout,"DT_out[%d] = %28.18e \n", i, DT_out[i]); fflush(stdout);
    }
    dump_ener(call_code);
#if( MAKE_STAT )
    dump_stat();  	/* Tally this timestep's statistical information */
#endif
#if(BBH_SPACETIME)
    dump_bbh_trajectories(call_code);
#endif

#if( MAKE_HISTORY || MAKE_SURFACE )
    set_bin_arrays();
#endif

    break;


    /******************************************************************************
      Make sure that everything is written at beginning and end of the run :
    *******************************************************************************/
  case OUT_FINAL : 
//    OUT_LOOP  if( mask[i] ) { 
//      dump[i]();
//      N_out[i]++;
//      T_out[i] +=  DT_out[i];
//    }

//    gdump();   /* dump metric functions once only at the beginning */

    dump_ener(call_code);
    break;


    /******************************************************************************
      Output that is done every timestep : 
    *******************************************************************************/
  case OUT_EVERY : 
    if( (nstep % 10) == 0 ) { 
      dump_ener(call_code);

#if( MAKE_STAT )
    dump_stat();  	/* Tally this timestep's statistical information */
#endif

#if(DYNAMIC_SPACETIME)
    dump_bbh_trajectories(call_code);
#endif
    }

    break;

    /******************************************************************************
      Otherwise, write what the call intended 
    *******************************************************************************/
  default : 
    i = call_code ; 
    if( mask[i] ) { 
      dump[i]();
      N_out[i]++;
      T_out[i] +=  DT_out[i];
    }

  } /* End  switch(call_code) */

  TRACE_END;
  return;
}

/***************************************************************************************************/
/***************************************************************************************************
  fail(): 
 --------
    -- Routine that handles the failure state of the code, e.g. handles file closures, 
       memory deallocations, error reports, etc. if any;
    -- Responsible for any error monitoring, recording (for statistical analysis), etc.
    -- "state" controls how to record or flag errors differently depending on the 
        caller's situation;

***************************************************************************************************/

void fail(int fail_type, int state)
{
  int ptmp;
  const  int max_fails = 10000;
  static int fail_count[N_FAIL]={0}; 
  static const int pflag_response[3] = {0,PFLAG_INTERP_PRIM,PFLAG_INTERP_V};

  /* Compile the statistical information only at the last temporal substep and only 
     when we do not intend to interpolate over the cell's value : */
#if( MAKE_STAT && DUMP_ALL_STAT ) 
  if( (n_substep == (N0-1)) && (state == 0) ) {   calc_dU_stat(fail_type,state); }
#endif 

  switch( fail_type ) { 

    /*****  Utoprim() returned with densities below the floor : ***********/
  case  FAIL_FLOOR   :
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_FLOOR][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_FLOOR]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: utoprim hit floor (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//	fprintf(stderr,"fail: FLOOR  no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break;

    /*****  evolved to a state with excessive temperature,  u > rho * Tmax *********/
  case  FAIL_TMAX   :
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_TMAX][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_TMAX]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: too hot, u > rho Tmax (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//	fprintf(stderr,"fail: TMAX  no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break;

    /*****  Reached the ceiling for gamma  : ******************/
  case  FAIL_GAMMA_MAX     :
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_GAMMA_MAX][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_GAMMA_MAX]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: hit gamma ceiling (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//	fprintf(stderr,"fail: GAMMA_MAX  no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break;

    /*****  Problem with the calculation of gamma, probably superluminal velocity *****/
  case  FAIL_GAMMA_CALC    :
    if( state )   { 
      failure_exists = 1 ; failure_type[pflag_response[state]] = 1;
      switch( pcurr ) { 
      case CENT  :  check_boundary_pflag(icurr  ,jcurr  ,kcurr  ); pflag[icurr  ][jcurr  ][kcurr  ] = pflag_response[state]; break;
      case FACE1 :  check_boundary_pflag(icurr-1,jcurr  ,kcurr  ); pflag[icurr-1][jcurr  ][kcurr  ] = pflag_response[state]; break ; 
      case FACE2 :  check_boundary_pflag(icurr  ,jcurr-1,kcurr  ); pflag[icurr  ][jcurr-1][kcurr  ] = pflag_response[state]; break ; 
      case FACE3 :  check_boundary_pflag(icurr  ,jcurr  ,kcurr-1); pflag[icurr  ][jcurr  ][kcurr-1] = pflag_response[state]; break ; 
      }
    }
    nfail[NFAIL_GAMMA_CALC][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_GAMMA_CALC]++;
    if( fail_count[fail_type-1] < max_fails ) { 
      struct of_coord *coords;
      get_coord(icurr,jcurr,kcurr,pcurr,ncurr,coords);
      fprintf(stderr,"\n\n fail: gamma_calc failure(pid=%d): %d : %d %d %d %d %28.18e %28.18e %28.18e %28.18e \n",
	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t,coords->x[1],coords->x[2],coords->x[3]) ;    
      fail_count[fail_type-1]++;
      if( fail_count[fail_type-1] >= max_fails ) { 
	fprintf(stderr,"fail: GAMMA_CALC no longer printing error messages!!  \n");
      }
      fflush(stderr);
    }
    break;

    /*****  using the full inversion solution  *************/
  case  FAIL_USE_FULL_INV  :
    nfail[NFAIL_USE_FULL_INV][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_USE_FULL_INV]++;
    break;

    /*****  using the entropy equation fix  *************/
  case  FAIL_USE_EE  :
    nfail[NFAIL_USE_EE][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_USE_EE]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: using entropy equation  (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//	fprintf(stderr,"fail: USE_EE   no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break;

    /*****  Utoprim() failed to converge to a reasonable solution : *************/
  case  FAIL_UTOPRIM  :
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_UTOPRIM][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_UTOPRIM]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: utoprim did not converge (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//	fprintf(stderr,"fail: UTOPRIM  no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break;

    /*****  Utoprim_ee() for entropy equation failed to converge to a reasonable solution : *************/
  case  FAIL_UTOPRIM_EE  :
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_UTOPRIM_EE][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_UTOPRIM_EE]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: utoprim did not converge (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//	fprintf(stderr,"fail: UTOPRIM_EE  no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break;

    /*****  fixup_interp_prim()'s interpolation procedure failed to find adequate stencil : */
  case  FAIL_USE_INTERP_PRIM         : 
    nfail[NFAIL_USE_INTERP_PRIM][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_USE_INTERP_PRIM]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: using interpolated values from fixup_interp_prim() (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: USE_INTERP_PRIM no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  fixup_interp_v()'s interpolation procedure failed to find adequate stencil : */
  case  FAIL_USE_INTERP_V        : 
    nfail[NFAIL_USE_INTERP_V][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_USE_INTERP_V]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: using interpolated values from fixup_interp_prim() (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: USE_INTERP_V no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  fixup_interp_prim()'s interpolation procedure failed to find adequate stencil : */
  case  FAIL_INTERP_PRIM         : 
    nfail[NFAIL_INTERP_PRIM][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_INTERP_PRIM]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: cannot find good stencil at (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: INTERP_PRIM no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  fixup_interp_v()'s interpolation procedure failed to find adequate stencil : */
  case  FAIL_INTERP_V        : 
    nfail[NFAIL_INTERP_V][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_INTERP_V]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: interp_v cannot find good stencil at (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: INTERP_V no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  whether we smoothed the pressure near the cutout boundary : */
  case  FAIL_CUTOUT_PRESSURE       : 
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_CUTOUT_PRESSURE][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_CUTOUT_PRESSURE]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: pressure change near cutout at (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: CUTOUT_PRESSURE no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  whether we smoothed the Lorentz factor  near the cutout boundary : */
  case  FAIL_CUTOUT_GAMMA       : 
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_CUTOUT_GAMMA][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_CUTOUT_GAMMA]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: gamma change near cutout at (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: CUTOUT_GAMMA no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  whether set_cutout_boundary*() changed anything */
  case  FAIL_CUTOUT_BC          :
    if( state ) { failure_exists = 1 ; check_boundary_pflag(icurr,jcurr,kcurr); pflag[icurr][jcurr][kcurr] = pflag_response[state]; failure_type[pflag_response[state]] = 1;}
    nfail[NFAIL_CUTOUT_BC][icurr][jcurr][kcurr]++;
    nfailtot[NFAIL_CUTOUT_BC]++;
//    if( fail_count[fail_type-1] < max_fails ) { 
//      fprintf(stderr,"\n\n fail: cutout BC change at (pid=%d): %d : %d %d %d %d %28.18e \n",
//	      myid,fail_type,icurr,jcurr,kcurr,pcurr,t) ;  
//      fail_count[fail_type-1]++;
//      if( fail_count[fail_type-1] >= max_fails ) { 
//    	fprintf(stderr,"fail: CUTOUT_BC no longer printing error messages!!  \n");
//      }
//      fflush(stderr);
//    }
    break; 

    /*****  Negative magnetosonic speed ****/
  case  FAIL_VCHAR_NEG     : 
    if( fail_count[fail_type-1] < max_fails ) { 
      fprintf(stderr,"\n\n fail: vchar neg (pid=%d): %d : [%d,%d,%d] [%d,%d,%d]g  %d %28.18e \n",
	      myid,fail_type,icurr,jcurr,kcurr,
	      globalpos[1]-NG+icurr,globalpos[2]-NG+jcurr,globalpos[3]-NG+kcurr,pcurr,t) ;
      fail_count[fail_type-1]++;
      if( fail_count[fail_type-1] >= max_fails ) { 
	fprintf(stderr,"fail: VCHAR_NEG no longer printing error messages!!  \n");
      }
      fflush(stderr);
    }
    break; 


    /*****  Superluminal magnetosonic speed ****/
  case  FAIL_VCHAR_SUPER   : 
    if( fail_count[fail_type-1] < max_fails ) { 
      fprintf(stderr,"\n\n fail: vchar super (pid=%d): %d : [%d,%d,%d] [%d,%d,%d]g  %d %28.18e \n",
	      myid,fail_type,icurr,jcurr,kcurr,
	      globalpos[1]-NG+icurr,globalpos[2]-NG+jcurr,globalpos[3]-NG+kcurr,pcurr,t) ;
      fail_count[fail_type-1]++;
      if( fail_count[fail_type-1] >= max_fails ) { 
	fprintf(stderr,"fail: VCHAR_SUPER no longer printing error messages!!  \n");
      }
      fflush(stderr);
    }
    break; 


    /*****  Complex characteristic speed  ****/
  case  FAIL_VCHAR_DISCR   : 
    if( fail_count[fail_type-1] < max_fails ) { 
      fprintf(stderr,"\n\n fail: vchar complex (pid=%d): %d : [%d,%d,%d] [%d,%d,%d]g  %d %28.18e \n",
	      myid,fail_type,icurr,jcurr,kcurr,
	      globalpos[1]-NG+icurr,globalpos[2]-NG+jcurr,globalpos[3]-NG+kcurr,pcurr,t) ;
      fail_count[fail_type-1]++;
      if( fail_count[fail_type-1] >= max_fails ) { 
	fprintf(stderr,"fail: VCHAR_DISCR no longer printing error messages!!  \n");
      }
      fflush(stderr);
    }
    break; 

    /*********************************************************************************
      These are fatal failures: 
    *********************************************************************************/
  case  FAIL_RESTART       :
  case  FAIL_MPI_BASIC     : 
  case  FAIL_BASIC         : 
  case  FAIL_METRIC        : 
  case  FAIL_HDF           :
    /* Need to first check if we have been here before since diag() below may cause failures */
    if( failed ) {  myexit(fail_type) ;  }

    /* Only output this message if this process is at fault : */
    if( !sibling_failed ) { 
      failed = fail_type ;
      fprintf(stderr,"\n\n fail: fundamental(pid=%d): %d : %d %d %d %d  %28.18e \n",
	      myid,fail_type, icurr,jcurr,kcurr,pcurr,t) ;    fflush(stderr);

      // area_map(icurr,jcurr,kcurr,p) ;

      /* Make sure that all non-failed processes get to this point too, so put exit_status() before each MPI call : */
      exit_status() ;  

    }

    /* Dump final output:  */
    diag(OUT_FINAL) ;

    myexit(fail_type) ;  /* this must be here */

    break;

  default :
    fprintf(stderr,"\n\n fail(): unknown type = %d \n",fail_type) ;   
    fflush(stderr);
    myexit(fail_type) ;
    break;

  }

}

/***********************************************************************************************/
/***********************************************************************************************
  divb_calc(): 
 --------------
    -- divergence of B^i  residual  based on stencil used in Flux CT scheme, so it should be 
        very small;
    -- the differencing is centered at the corner (i.e. CORN position) of cell (i,j,k)  ;
    -- note that we impicitly difference via negative signs in the Btmp[1-3] assigmnments;
***********************************************************************************************/
double divb_calc( int i, int j, int k ) 
{
  unsigned int id;
  double divb, det_g[8], Btmp1[8], Btmp2[8], Btmp3[8], ftmp;
  
  struct of_geom *geom;
  
  id = 0;
  get_geometry(i-1,j-1,k-1,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i-1,j-1,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i-1,j  ,k-1,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i-1,j  ,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j-1,k-1,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j-1,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j  ,k-1,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j  ,k  ,CENT,ncurr,geom);  det_g[id  ] = geom->g;

  id = 0;
  Btmp1[id++] = -p[i-1][j-1][k-1][B1];
  Btmp1[id++] = -p[i-1][j-1][k  ][B1];
  Btmp1[id++] = -p[i-1][j  ][k-1][B1];
  Btmp1[id++] = -p[i-1][j  ][k  ][B1];
  Btmp1[id++] =  p[i  ][j-1][k-1][B1];
  Btmp1[id++] =  p[i  ][j-1][k  ][B1];
  Btmp1[id++] =  p[i  ][j  ][k-1][B1];
  Btmp1[id  ] =  p[i  ][j  ][k  ][B1];

  id = 0;
  Btmp2[id++] = -p[i-1][j-1][k-1][B2];
  Btmp2[id++] = -p[i-1][j-1][k  ][B2];
  Btmp2[id++] =  p[i-1][j  ][k-1][B2];
  Btmp2[id++] =  p[i-1][j  ][k  ][B2];
  Btmp2[id++] = -p[i  ][j-1][k-1][B2];
  Btmp2[id++] = -p[i  ][j-1][k  ][B2];
  Btmp2[id++] =  p[i  ][j  ][k-1][B2];
  Btmp2[id  ] =  p[i  ][j  ][k  ][B2];

  id = 0;
  Btmp3[id++] = -p[i-1][j-1][k-1][B3];
  Btmp3[id++] =  p[i-1][j-1][k  ][B3];
  Btmp3[id++] = -p[i-1][j  ][k-1][B3];
  Btmp3[id++] =  p[i-1][j  ][k  ][B3];
  Btmp3[id++] = -p[i  ][j-1][k-1][B3];
  Btmp3[id++] =  p[i  ][j-1][k  ][B3];
  Btmp3[id++] = -p[i  ][j  ][k-1][B3];
  Btmp3[id  ] =  p[i  ][j  ][k  ][B3];

  divb = 0.;

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += Btmp1[id]*det_g[id]; } 
  divb += ftmp*invdx[1];

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += Btmp2[id]*det_g[id]; } 
  divb += ftmp*invdx[2];

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += Btmp3[id]*det_g[id]; } 
  divb += ftmp*invdx[3];

  divb *= 0.25;
  divb = fabs(divb);


//--old   divb = fabs( 
//--old 	      0.25*(
//--old 		      p[i  ][j  ][k  ][B1] * gdet[i  ][j  ][CENT] 
//--old 		    + p[i  ][j-1][k  ][B1] * gdet[i  ][j-1][CENT] 
//--old 		    + p[i  ][j  ][k-1][B1] * gdet[i  ][j  ][CENT] 
//--old 		    + p[i  ][j-1][k-1][B1] * gdet[i  ][j-1][CENT]
//--old 		    - p[i-1][j  ][k  ][B1] * gdet[i-1][j  ][CENT] 
//--old 		    - p[i-1][j-1][k  ][B1] * gdet[i-1][j-1][CENT] 
//--old 		    - p[i-1][j  ][k-1][B1] * gdet[i-1][j  ][CENT] 
//--old 		    - p[i-1][j-1][k-1][B1] * gdet[i-1][j-1][CENT]
//--old 		    )/dx[1] +
//--old 	      0.25*(
//--old 		      p[i  ][j  ][k  ][B2] * gdet[i  ][j  ][CENT] 
//--old 		    + p[i-1][j  ][k  ][B2] * gdet[i-1][j  ][CENT] 
//--old 		    + p[i  ][j  ][k-1][B2] * gdet[i  ][j  ][CENT] 
//--old 		    + p[i-1][j  ][k-1][B2] * gdet[i-1][j  ][CENT]
//--old 		    - p[i  ][j-1][k  ][B2] * gdet[i  ][j-1][CENT] 
//--old 		    - p[i-1][j-1][k  ][B2] * gdet[i-1][j-1][CENT] 
//--old 		    - p[i  ][j-1][k-1][B2] * gdet[i  ][j-1][CENT] 
//--old 		    - p[i-1][j-1][k-1][B2] * gdet[i-1][j-1][CENT]
//--old 		    )/dx[2] +
//--old 	      0.25*(
//--old 		      p[i  ][j  ][k  ][B3] * gdet[i  ][j  ][CENT] 
//--old 		    + p[i-1][j  ][k  ][B3] * gdet[i-1][j  ][CENT] 
//--old 		    + p[i  ][j-1][k  ][B3] * gdet[i  ][j-1][CENT] 
//--old 		    + p[i-1][j-1][k  ][B3] * gdet[i-1][j-1][CENT]
//--old 		    - p[i  ][j  ][k-1][B3] * gdet[i  ][j  ][CENT] 
//--old 		    - p[i-1][j  ][k-1][B3] * gdet[i-1][j  ][CENT] 
//--old 		    - p[i  ][j-1][k-1][B3] * gdet[i  ][j-1][CENT] 
//--old 		    - p[i-1][j-1][k-1][B3] * gdet[i-1][j-1][CENT]
//--old 		    )/dx[3] 
//--old 	      ) ;

  return( divb  );
}

/***********************************************************************************************/
/***********************************************************************************************
  divb_cen_calc(): 
 --------------
    -- divergence of B^i  residual  based on the stencil centered at the center (i.e. the CENT
       position)  of cell (i,j,k);
    -- this residual need not be small as it is not based on the stencil used in our 
          constraint transport scheme;
***********************************************************************************************/
double divb_cen_calc( int i, int j, int k ) 
{
  double divb_cen;
  unsigned int id;
  double det_g[6], Btmp[6];
  struct of_geom *geom;
  
  id = 0;
  get_geometry(i+1,j  ,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i-1,j  ,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j+1,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j-1,k  ,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j  ,k+1,CENT,ncurr,geom);  det_g[id++] = geom->g;
  get_geometry(i  ,j  ,k-1,CENT,ncurr,geom);  det_g[id  ] = geom->g;

  id = 0;
  Btmp[id] = p[i+1][j  ][k  ][B1] * det_g[id];    id++ ; 
  Btmp[id] = p[i-1][j  ][k  ][B1] * det_g[id];    id++ ; 
  Btmp[id] = p[i  ][j+1][k  ][B2] * det_g[id];    id++ ; 
  Btmp[id] = p[i  ][j-1][k  ][B2] * det_g[id];    id++ ; 
  Btmp[id] = p[i  ][j  ][k+1][B3] * det_g[id];    id++ ; 
  Btmp[id] = p[i  ][j  ][k-1][B3] * det_g[id];

  id = 0;
  divb_cen = 0.5*fabs( 
		  (Btmp[0]-Btmp[1])*invdx[1]  + 
		  (Btmp[2]-Btmp[3])*invdx[2]  + 
		  (Btmp[4]-Btmp[5])*invdx[3]  
		  );

//--old   divb_cen = fabs(
//--old 		  (
//--old 		     p[i+1][j][k][B1]*gdet[i+1][j][CENT] 
//--old 		   - p[i-1][j][k][B1]*gdet[i-1][j][CENT]
//--old 		   )/(2.*dx[1]) +
//--old 		  (
//--old 		     p[i][j+1][k][B2]*gdet[i][j+1][CENT] 
//--old 		   - p[i][j-1][k][B2]*gdet[i][j-1][CENT]
//--old 		   )/(2.*dx[2]) +
//--old 		  (  
//--old 		     p[i][j][k+1][B3]*gdet[i][j][CENT] 
//--old 		   - p[i][j][k-1][B3]*gdet[i][j][CENT] 
//--old 		   )/(2.*dx[3])
//--old 		  ) ;

  return( divb_cen  );
}

/****************************************************************************************/
/*****************************************************************************************
  current_time():
  --------------
     -- since the calculation of the current is time-consuming, this routine 
        is designed to tell harm3d when it should start saving the faraday tensor 
        and when it should calculate the current; 
     -- the current calculation for t^n is done during the calculation for t^{n+1}; 
        the current requires the faraday tensor at the previous and present time steps, 
        so we need to return with (1) when we are at a timestep before the next output 
        time; 
*****************************************************************************************/
int current_time(void) 
{
  int i;
  
  double tnext;

  tnext = GridLength[0];

  OUT_LOOP { 
    if( mask[i] && (tnext > T_out[i]) )  { tnext = T_out[i]; }  
  }

  if( t >= (tnext - 2.1*SAFE*dx[0]) ) {
    return( 1 ) ; 
  }
  else { 
    return( 0 ) ; 
  }

}


/****************************************************************************************/
/*****************************************************************************************
  dump_ener():
  --------------
     -- write out per timestep information  that is in  <RUN_TAG>_ener.out 
     -- output format : 
         nstep  t  dt   U_tot[0-(NP-1)]  M_tot  F_tot[0-2][0-1][0-(NP-1)]

*****************************************************************************************/
void dump_ener(int call_code)
{
  int i,j,k,l,d,g,nvars;
  double v_send[100], v_recv[100];
  double U_tot[NP]={0.}, F_tot[NDIM-1][2][NP]={0.};
  
  char ener_filename[200];
  
  static int local_first_time = 1 ; 
  static FILE *ener_file; 

  TRACE_BEG;
  //  fprintf(stdout,"dump_ener begin pid=%d   nstep=%d\n",myid,nstep); fflush(stdout);

  /* Open the file if this is the first time, never clobber the file :  */  
  if( local_first_time ) { 
    if( myid == ener_out_pid ) { 
      sprintf(ener_filename,"%s_ener.out",RUN_TAG);
      ener_file = fopen(ener_filename,"at");
      if(ener_file==NULL) {
	fprintf(stderr,"dump_ener():  error opening %s file\n",ener_filename) ;
	if( call_code != OUT_FINAL ) { fail( FAIL_BASIC,0 ); }
      }
    }
    local_first_time = 0; 
  }

  if( call_code == OUT_FINAL ) { 
    if( myid == ener_out_pid ) { fclose(ener_file); }   
    TRACE_END;
    return;
  } 

  /* The fluxes and conserved variables have not been set yet so do not dump : */
  if( call_code == OUT_INIT ) {  
    TRACE_END;
    return; 
  }   

  /**************************************************************************************
    Calculated integrated quantities: 
  ***************************************************************************************/
#if( USE_MASK )  
  LOOP if( evol_mask[ncurr][i][j][k] ) {  
#else
  LOOP  {  
#endif
      PLOOP {  U_tot[l] +=  U_gf[n_U][i][j][k][l] ; }
  }

  for(l=0;l<=NP;l++) { U_tot[l] *= dV; } 
  
  if(bc_pid[1][BCDN] == BC_PHYS) { N2_LOOP N3_LOOP PLOOP { F_tot[0][BCDN][l] += F[N1S  ][j    ][k    ][0][l] ; } }
  if(bc_pid[1][BCUP] == BC_PHYS) { N2_LOOP N3_LOOP PLOOP { F_tot[0][BCUP][l] += F[N1E+1][j    ][k    ][0][l] ; } }
  if(bc_pid[2][BCDN] == BC_PHYS) { N1_LOOP N3_LOOP PLOOP { F_tot[1][BCDN][l] += F[i    ][N2S  ][k    ][1][l] ; } }
  if(bc_pid[2][BCUP] == BC_PHYS) { N1_LOOP N3_LOOP PLOOP { F_tot[1][BCUP][l] += F[i    ][N2E+1][k    ][1][l] ; } }
  if(bc_pid[3][BCDN] == BC_PHYS) { N1_LOOP N2_LOOP PLOOP { F_tot[2][BCDN][l] += F[i    ][j    ][N3S  ][2][l] ; } }
  if(bc_pid[3][BCUP] == BC_PHYS) { N1_LOOP N2_LOOP PLOOP { F_tot[2][BCUP][l] += F[i    ][j    ][N3E+1][2][l] ; } }

  PLOOP F_tot[0][BCDN][l] *= dx[2]*dx[3]; 
  PLOOP F_tot[0][BCUP][l] *= dx[2]*dx[3]; 
  PLOOP F_tot[1][BCDN][l] *= dx[1]*dx[3]; 
  PLOOP F_tot[1][BCUP][l] *= dx[1]*dx[3]; 
  PLOOP F_tot[2][BCDN][l] *= dx[1]*dx[2]; 
  PLOOP F_tot[2][BCUP][l] *= dx[1]*dx[2]; 

  /**************************************************************************************
    Integrate quantities over subdomains: 
  ***************************************************************************************/
  i=0;
  PLOOP { v_send[i++] = U_tot[l];  } 
  v_send[i++] = M_tot; 
  for(d=0;d<NDIM-1;d++) for(g=0;g<2;g++) PLOOP {  v_send[i++] = F_tot[d][g][l]; } 

  nvars = i;
  if( nvars > 100 ) { 
    fprintf(stderr,"dump_ener(): Need to increase dimensions of v_send and v_recv !!   nvars = %d \n",nvars);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  mpi_sum_to_head(v_send,v_recv,nvars);

  /**************************************************************************************
    Print out the integrated quantities: 
  ***************************************************************************************/
  if( myid == ener_out_pid ) { 
    fprintf(ener_file," %d %28.18e %28.18e ",nstep,t,dx[0]); 
    for(i=0;i<nvars;i++) { fprintf(ener_file," %28.18e ",v_recv[i]); }
    fprintf(ener_file,"\n"); 
    fflush(ener_file);

    E_tot = -(v_recv[UU]-v_recv[RHO]);
    L_tot = v_recv[U3];
  }


  //  fprintf(stdout,"dump_ener end pid=%d   nstep=%d\n",myid,nstep); fflush(stdout);

  TRACE_END;
  return;
}

/****************************************************************************************/
/*****************************************************************************************
  trace_message():
  --------------
     -- prints a message to stdout useful in debugging the code ; 
     -- should be called at the beginning and end of each routine in order to ensure 
         adequate tracing of bugs;
*****************************************************************************************/
void trace_message(char *routine_name, int instance)
{
  if( instance==0 ) { fprintf(stdout,"%s beg : %d %d %d \n",routine_name,myid,nstep,n_substep); }
  else              { fprintf(stdout,"%s end : %d %d %d \n",routine_name,myid,nstep,n_substep); }
  fflush(stdout);
  return;
}

/****************************************************************************************/
/*****************************************************************************************
  dump_bbh_trajectory():
  --------------
     -- prints the trajectories of the pair of black holes over time; 
*****************************************************************************************/
void dump_bbh_trajectories(int call_code)
{
  int i,j,k,l,d,g,nvars;
  char traj_filename[200];
  static int local_first_time = 1 ; 
  static FILE *traj_file; 

  TRACE_BEG;

  SET_GEOM_ONLY_TIME_FUNCS(t,1);

#if( BBH_SPACETIME )

  //  fprintf(stdout,"dump_ener begin pid=%d   nstep=%d\n",myid,nstep); fflush(stdout);

  /* Open the file if this is the first time, never clobber the file :  */  
  if( local_first_time ) { 
    if( myid == traj_out_pid ) { 
      sprintf(traj_filename,"%s_bbh_trajectory.out",RUN_TAG);
      traj_file = fopen(traj_filename,"at");
      if(traj_file==NULL) {
	fprintf(stderr,"dump_bbh_trajectory():  error opening %s file\n",traj_filename) ;
	fail( FAIL_BASIC,0 );
      }

      
#if( (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_FAST ) ||	\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ   ) ||	\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ) ||	\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NR_FZ)    ||  \
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_DROP)  ||  \
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_SIMPLE_NEWTONIAN)||  \
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW) )
      fprintf(traj_file,"# LEGEND: \n"); 
      fprintf(traj_file,"# ------------ \n"); 
      fprintf(traj_file,"#\t    0): t     ==  current time                                                            \n");
      fprintf(traj_file,"#\t    1): m1    ==  mass of BH1								\n");
      fprintf(traj_file,"#\t    2): m2    ==  mass of BH2								\n");
      fprintf(traj_file,"#\t    3): b	    ==  initial separation							\n");
      fprintf(traj_file,"#\t    4): xi1x  ==  x coordinate of current position of BH1					\n");
      fprintf(traj_file,"#\t    5): xi1y  ==  y coordinate of current position of BH1					\n");
      fprintf(traj_file,"#\t    6): xi2x  ==  x coordinate of current position of BH2					\n");
      fprintf(traj_file,"#\t    7): xi2y  ==  y coordinate of current position of BH2					\n");
      fprintf(traj_file,"#\t    8): v1x   ==  x component of velocity of BH1						\n");
      fprintf(traj_file,"#\t    9): v1y   ==  y component of velocity of BH1						\n");
      fprintf(traj_file,"#\t   10): v2x   ==  x component of velocity of BH2						\n");
      fprintf(traj_file,"#\t   11): v2y   ==  y component of velocity of BH2						\n");
      fprintf(traj_file,"#\t   12): v1    ==  magnitude of velocity of BH1						\n");
      fprintf(traj_file,"#\t   13): v2    ==  magnitude of velocity of BH2						\n");
      fprintf(traj_file,"#\t   14): v12x  ==  x component of v1-v2							\n");
      fprintf(traj_file,"#\t   15): v12y  ==  y component of v1-v2							\n");
      fprintf(traj_file,"#\t   16): v21x  ==  x component of v2-v1							\n");
      fprintf(traj_file,"#\t   17): v21y  ==  y component of v2-v1							\n");
      fprintf(traj_file,"#\t   18): v12   ==  magnitude of v1-v2							\n");
      fprintf(traj_file,"#\t   19): v21   ==  magnitude of v2-v1							\n");
      fprintf(traj_file,"#\t   20): v1v2  ==  dot product of v1 and v2  						\n");
      fprintf(traj_file,"#\t   21): v2v1  ==  dot product of v2 and v1 (different from v1v2 due to PN approximation)	\n");
      fprintf(traj_file,"#\t   22): t_c   ==  time to merger from t=0							\n");
      fprintf(traj_file,"#\t   23): phi   ==  current orbital phase (w.r.t. phi(t=0) = 0 )				\n");
      fprintf(traj_file,"#\t   24): omega ==  current rate of change of orbital phase					\n");
      fprintf(traj_file,"#\t   25): r12   ==  current separation							\n");
      fprintf(traj_file,"#\t   26): r21   ==  current separation (from BH2's perspective)				\n");
      fprintf(traj_file,"#\t   27): xi1z  ==  z coordinate of current position of BH1					\n");
      fprintf(traj_file,"#\t   28): xi2z  ==  z coordinate of current position of BH2					\n");
      fprintf(traj_file,"#\t   29): v1z   ==  z component of velocity of BH1						\n");
      fprintf(traj_file,"#\t   30): v2z   ==  z component of velocity of BH2						\n");
      fprintf(traj_file,"#\t   31): v12z  ==  z component of v1-v2							\n");
      fprintf(traj_file,"#\t   32): v21z  ==  z component of v2-v1                                                    \n"); 
      fprintf(traj_file,"#\t   33): r1T0  ==  Inner1-Near radius trans. func. parameter                               \n"); 
      fprintf(traj_file,"#\t   34): w1T0  ==  Inner1-Near width trans. func. parameter                               \n"); 
      fprintf(traj_file,"#\t   35): r2T0  ==  Inner2-Near radius trans. func. parameter                               \n"); 
      fprintf(traj_file,"#\t   36): w2T0  ==  Inner2-Near width trans. func. parameter                               \n"); 
      fprintf(traj_file,"#\t   37): xNT0  ==  Near-Inner radius trans. func. parameter                               \n"); 
      fprintf(traj_file,"#\t   38): wNT0  ==  Near-Inner width trans. func. parameter                               \n"); 
      fprintf(traj_file,"#\t   39): lambda  ==  Near-Far trans. func. parameter                               \n"); 
#elif( METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND ||\
       METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND )
      fprintf(traj_file,"# LEGEND: \n");
      fprintf(traj_file,"# ------------ \n");
      fprintf(traj_file,"#\t    0): t     ==  current time                                                            \n");
      fprintf(traj_file,"#\t    1): m1    ==  mass of BH1                        \n");
      fprintf(traj_file,"#\t    2): m2    ==  mass of BH2                        \n");
      fprintf(traj_file,"#\t    3): b      ==  initial separation                   \n");
      fprintf(traj_file,"#\t    4): xi1x  ==  x coordinate of current position of BH1              \n");
      fprintf(traj_file,"#\t    5): xi1y  ==  y coordinate of current position of BH1              \n");
      fprintf(traj_file,"#\t    6): xi2x  ==  x coordinate of current position of BH2              \n");
      fprintf(traj_file,"#\t    7): xi2y  ==  y coordinate of current position of BH2              \n");
      fprintf(traj_file,"#\t    8): v1x   ==  x component of velocity of BH1                 \n");
      fprintf(traj_file,"#\t    9): v1y   ==  y component of velocity of BH1                 \n");
      fprintf(traj_file,"#\t   10): v2x   ==  x component of velocity of BH2                 \n");
      fprintf(traj_file,"#\t   11): v2y   ==  y component of velocity of BH2                 \n");
      fprintf(traj_file,"#\t   12): v1    ==  magnitude of velocity of BH1                \n");
      fprintf(traj_file,"#\t   13): v2    ==  magnitude of velocity of BH2                \n");
      fprintf(traj_file,"#\t   14): v12x  ==  x component of v1-v2                     \n");
      fprintf(traj_file,"#\t   15): v12y  ==  y component of v1-v2                     \n");
      fprintf(traj_file,"#\t   16): v12   ==  magnitude of v1-v2                    \n");
      fprintf(traj_file,"#\t   17): v1v2  ==  dot product of v1 and v2                    \n");
      fprintf(traj_file,"#\t   18): n12x  ==  x component of n1-n2 \n");
      fprintf(traj_file,"#\t   19): n12y  ==  y component of n1-n2 \n");
      fprintf(traj_file,"#\t   20): n12v12  ==  dot product: (n1-n2) . (v1-v2) \n");
      fprintf(traj_file,"#\t   21): n12v1  ==  dot product: (n1-n2) . v1 \n");
      fprintf(traj_file,"#\t   22): n12v2  ==  dot product: (n1-n2) . v2 \n");
      fprintf(traj_file,"#\t   23): t_c   ==  time to merger from t=0                     \n");
      fprintf(traj_file,"#\t   24): phi   ==  current orbital phase (w.r.t. phi(t=0) = 0 )            \n");
      fprintf(traj_file,"#\t   25): omega ==  current rate of change of orbital phase              \n");
      fprintf(traj_file,"#\t   26): r12   ==  current separation                    \n");
      fprintf(traj_file,"#\t   27): r12dot ==  current separation time derivative                     \n");
      fprintf(traj_file,"#\t   28): xi1z  ==  z coordinate of current position of BH1              \n");
      fprintf(traj_file,"#\t   29): xi2z  ==  z coordinate of current position of BH2              \n");
      fprintf(traj_file,"#\t   30): v1z   ==  z component of velocity of BH1                 \n");
      fprintf(traj_file,"#\t   31): v2z   ==  z component of velocity of BH2                 \n");
      fprintf(traj_file,"#\t   32): v12z  ==  z component of v1-v2                     \n");
      fprintf(traj_file,"#\t   33): n12z  ==  z component of n1-n2                     \n");
      fprintf(traj_file,"#\t   34): r1T0  ==  Inner1-Near radius trans. func. parameter                               \n");
      fprintf(traj_file,"#\t   35): w1T0  ==  Inner1-Near width trans. func. parameter                               \n");
      fprintf(traj_file,"#\t   36): r2T0  ==  Inner2-Near radius trans. func. parameter                               \n");
      fprintf(traj_file,"#\t   37): w2T0  ==  Inner2-Near width trans. func. parameter                               \n");
      fprintf(traj_file,"#\t   38): xNT0  ==  Near-Inner radius trans. func. parameter                               \n");
      fprintf(traj_file,"#\t   39): wNT0  ==  Near-Inner width trans. func. parameter                               \n");
      fprintf(traj_file,"#\t   40): lambda  ==  Near-Far trans. func. parameter                               \n");
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )
      fprintf(traj_file,"# LEGEND: \n"); 
      fprintf(traj_file,"# ------------ \n"); 
      fprintf(traj_file,"#\t    SYSTEM:\n");
      fprintf(traj_file,"#\t    ------\n"); 
      fprintf(traj_file,"#\t    Title: Black-Hole binary with spins aligned to orbital angular momentum \n");
      fprintf(traj_file,"#\t    m1    ==  mass of BH1 = %26.16e \n", m_bh1);
      fprintf(traj_file,"#\t    m2    ==  mass of BH2 = %26.16e \n", m_bh2);
      fprintf(traj_file,"#\t    chi1  ==  z spin of BH1 = %26.16e \n", chi_bh1);
      fprintf(traj_file,"#\t    chi2  ==  z spin of BH2 = %26.16e \n", chi_bh2);
      fprintf(traj_file,"#\t    b     ==  initial separation = %26.16e \n", initial_bbh_separation);
      fprintf(traj_file,"#\t    t_c   ==  time to merger from t=0 = %26.6e \n",bbh_params->t_c);
      fprintf(traj_file,"#\t    ------\n"); 
      fprintf(traj_file,"#\t    COLUMNS:\n");
      fprintf(traj_file,"#\t    ------\n"); 
      fprintf(traj_file,"#\t    0): t     ==  current time  \n");
      fprintf(traj_file,"#\t    1): xi1x  ==  x coordinate of current position of BH1     \n");
      fprintf(traj_file,"#\t    2): xi1y  ==  y coordinate of current position of BH1     \n");
      fprintf(traj_file,"#\t    3): xi1z  ==  z coordinate of current position of BH1     \n");
      fprintf(traj_file,"#\t    4): xi2x  ==  x coordinate of current position of BH2     \n");
      fprintf(traj_file,"#\t    5): xi2y  ==  y coordinate of current position of BH2     \n");
      fprintf(traj_file,"#\t    6): xi2z  ==  z coordinate of current position of BH2     \n");
      fprintf(traj_file,"#\t    7): v1x   ==  x component of velocity of BH1      \n");
      fprintf(traj_file,"#\t    8): v1y   ==  y component of velocity of BH1      \n");
      fprintf(traj_file,"#\t    9): v1z   ==  z component of velocity of BH1      \n");
      fprintf(traj_file,"#\t   10): v2x   ==  x component of velocity of BH2      \n");
      fprintf(traj_file,"#\t   11): v2y   ==  y component of velocity of BH2      \n");
      fprintf(traj_file,"#\t   12): v2z   ==  z component of velocity of BH2      \n");
      fprintf(traj_file,"#\t   13): v1    ==  magnitude of velocity of BH1      \n");
      fprintf(traj_file,"#\t   14): v2    ==  magnitude of velocity of BH2      \n");
      fprintf(traj_file,"#\t   15): phi   ==  current orbital phase (w.r.t. phi(t=0) = 0 ) \n");
      fprintf(traj_file,"#\t   16): omega ==  current rate of change of orbital phase      \n");
      fprintf(traj_file,"#\t   17): r12   ==  current separation                        \n");
      fprintf(traj_file,"#\t   18): r12dot ==  current separation time derivative       \n");
      fprintf(traj_file,"#\t   19): r1T0  ==  Inner1-Near radius trans. func. parameter  \n"); 
      fprintf(traj_file,"#\t   20): w1T0  ==  Inner1-Near width trans. func. parameter   \n"); 
      fprintf(traj_file,"#\t   21): r2T0  ==  Inner2-Near radius trans. func. parameter  \n"); 
      fprintf(traj_file,"#\t   22): w2T0  ==  Inner2-Near width trans. func. parameter   \n"); 
      fprintf(traj_file,"#\t   23): xNT0  ==  Near-Inner radius trans. func. parameter   \n"); 
      fprintf(traj_file,"#\t   24): wNT0  ==  Near-Inner width trans. func. parameter    \n"); 
      fprintf(traj_file,"#\t   25): lambda  ==  Near-Far trans. func. parameter          \n"); 
#else
      fprintf(traj_file,"# LEGEND: \n"); 
      fprintf(traj_file,"# ------------ \n"); 
      fprintf(traj_file,"#\t      0): t     ==  current time                                                            \n");
      fprintf(traj_file,"#\t      1): m1    ==  mass of BH1								\n");
      fprintf(traj_file,"#\t      2): m2    ==  mass of BH2								\n");
      fprintf(traj_file,"#\t      3): b	    ==  initial separation							\n");
      fprintf(traj_file,"#\t      4): xi1x  ==  x coordinate of current position of BH1					\n");
      fprintf(traj_file,"#\t      5): xi1y  ==  y coordinate of current position of BH1					\n");
      fprintf(traj_file,"#\t      6): xi2x  ==  x coordinate of current position of BH2					\n");
      fprintf(traj_file,"#\t      7): xi2y  ==  y coordinate of current position of BH2					\n");
      fprintf(traj_file,"#\t      8): xi1z  ==  z coordinate of current position of BH1					\n");
      fprintf(traj_file,"#\t      9): xi2z  ==  z coordinate of current position of BH2					\n");
#endif
      fprintf(traj_file,"# \n");
      fflush(traj_file);
    }
    local_first_time = 0; 
  }


#if( (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_FAST )     ||\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ   )     ||\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ)     ||\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NR_FZ   )     ||\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_DROP)      ||\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_SIMPLE_NEWTONIAN)    ||\
     (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW ) )

    struct of_bbh_traj bbh_traj; 
    extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
    get_bbh_traj_data(&bbh_traj) ;

  if( myid == traj_out_pid ) { 
    fprintf(traj_file,"%26.16e %26.16e %26.16e %26.16e ", t, m_bh1, m_bh2, initial_bbh_separation); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi1x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi1y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi2x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi2y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1x   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1y   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2x   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2y   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1    ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2    ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v21x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v21y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v21   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1v2  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2v1  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.t_c   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.phi   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.omega ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r12   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r21   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi1z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi2z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1z   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2z   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v21z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r1T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.w1T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r2T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.w2T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xNT0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.wNT0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.lambda); 
    fprintf(traj_file,"\n"); 
    fflush(traj_file);
    if( call_code == OUT_FINAL ) { fclose(traj_file); }
  }

#elif( (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND)    ||\
       (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND) )

      struct of_bbh_traj bbh_traj; 
    extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
    get_bbh_traj_data(&bbh_traj) ;

  if( myid == traj_out_pid ) { 
    fprintf(traj_file,"%26.16e %26.16e %26.16e %26.16e ", t, m_bh1, m_bh2, initial_bbh_separation); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi1x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi1y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi2x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi2y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1x   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1y   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2x   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2y   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1    ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2    ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1v2  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.n12x  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.n12y  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.n12v12); 
    fprintf(traj_file,"%26.16e ",bbh_traj.n12v1 ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.n12v2 ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.t_c   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.phi   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.omega ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r12   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r12dot); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi1z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xi2z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v1z   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v2z   ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.v12z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.n12z  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r1T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.w1T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.r2T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.w2T0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.xNT0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.wNT0  ); 
    fprintf(traj_file,"%26.16e ",bbh_traj.lambda); 
    fprintf(traj_file,"\n"); 
    fflush(traj_file);
    if( call_code == OUT_FINAL ) { fclose(traj_file); }
  }
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )
  if( myid == traj_out_pid ) { 
    fprintf(traj_file,"%26.16e ", bbh_params->t); // Current simulation time
    // Binary BH trajectory coordinates:
    fprintf(traj_file,"%26.16e ", bbh_params->xi1x); 
    fprintf(traj_file,"%26.16e ", bbh_params->xi1y); 
    fprintf(traj_file,"%26.16e ", bbh_params->xi1z); 
    fprintf(traj_file,"%26.16e ", bbh_params->xi2x); 
    fprintf(traj_file,"%26.16e ", bbh_params->xi2y); 
    fprintf(traj_file,"%26.16e ", bbh_params->xi2z); 
    // Binary BH velocity components:
    fprintf(traj_file,"%26.16e ", bbh_params->v1x); 
    fprintf(traj_file,"%26.16e ", bbh_params->v1y); 
    fprintf(traj_file,"%26.16e ", bbh_params->v1z); 
    fprintf(traj_file,"%26.16e ", bbh_params->v2x); 
    fprintf(traj_file,"%26.16e ", bbh_params->v2y); 
    fprintf(traj_file,"%26.16e ", bbh_params->v2z); 
    // Binary BH velocity magnitudes:
    fprintf(traj_file,"%26.16e ", bbh_params->v1); 
    fprintf(traj_file,"%26.16e ", bbh_params->v2); 
    // Binary BH orbital phase, separation and time derivatives:
    fprintf(traj_file,"%26.16e ", bbh_params->phi);    // Current binary orbital phase (w.r.t. phi(t=0) = 0) 
    fprintf(traj_file,"%26.16e ", bbh_params->omega);  // Current binary orbital phase rate of change 
    fprintf(traj_file,"%26.16e ", bbh_params->r12);    // Current binary separation
    fprintf(traj_file,"%26.16e ", bbh_params->r12dot); // Current binary radial velocity (separation time derivative)
    // Buffer zone parameters:
    fprintf(traj_file,"%26.16e ", bbh_params->r1T0);   // Inner1-Near radius trans. func. parameter
    fprintf(traj_file,"%26.16e ", bbh_params->w1T0);   // Inner1-Near width trans. func. parameter
    fprintf(traj_file,"%26.16e ", bbh_params->r2T0);   // Inner2-Near radius trans. func. parameter
    fprintf(traj_file,"%26.16e ", bbh_params->w2T0);   // Inner2-Near width trans. func. parameter
    fprintf(traj_file,"%26.16e ", bbh_params->xNT0);   // Near-Inner radius trans. func. parameter
    fprintf(traj_file,"%26.16e ", bbh_params->wNT0);   // Near-Inner width trans. func. parameter
    fprintf(traj_file,"%26.16e ", bbh_params->lambda); // Near-Far trans. func. parameter
    fprintf(traj_file,"\n");
    fflush(traj_file);
    if( call_code == OUT_FINAL ) { fclose(traj_file); }
  }
#else
  if( myid == traj_out_pid ) { 
    fprintf(traj_file,"%26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e  \n",
	    bh1_traj[ncurr][TT], 
	    m_bh1, m_bh2, initial_bbh_separation,
	    bh1_traj[ncurr][XX], bh1_traj[ncurr][YY],
	    bh2_traj[ncurr][XX], bh2_traj[ncurr][YY], 
	    bh1_traj[ncurr][ZZ], bh2_traj[ncurr][ZZ]);
    fflush(traj_file);
    if( call_code == OUT_FINAL ) { fclose(traj_file); }
  }
#endif

#endif   /*  #if( BBH_SPACETIME )  */

  TRACE_END;
  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  set_bin_arrays();
 ------------------
   -- sets coordinate arrays used for various things including history and surface outputs;
   -- the reason why it is here is because the arrays are used by other routines (e.g. init data), 
      in addition to the history and surface routines;
   -- r_bins, phi_bins  are the locations of the left-most edge of the bins within which we bin up data 

   -- phi_bins = [0,2pi]/(2*totalsize[3])
   -- r = r0*exp(xp1)    -->  xp1 = log(r/r0) 
**********************************************************************************************/
void set_bin_arrays(void)
{
  int i,j,k; 
  struct of_coord *coords;

  TRACE_BEG;

#if( COORD_RADIUS_DIM == 1 )
  //--HERE  fprintf(stdout,"set_bin_arrays():  we should not be here.... \n"); fflush(stdout); fail(FAIL_BASIC,0);
#endif

  /* Set the parameters if we are starting from scratch and not a checkpoint file: */
  if( n_r_bins < 1 ) { 
    n_r_bins   = 2*totalsize[1];
    n_phi_bins = 2*totalsize[3];
    r0_bins = 1.;
    phi_min_bins = 0.;
    phi_max_bins = 2.*M_PI; 
    dphi_bins = (phi_max_bins-phi_min_bins) / n_phi_bins; 
    inv_dphi_bins = 1./dphi_bins;

    if( (n_r_bins < 1) || (n_phi_bins < 1) ) { 
      fprintf(stdout,"make_bin_coords(): bad value for n_r_bins or n_phi_bins :  %d  %d   \n", n_r_bins,n_phi_bins); fflush(stdout); fail(FAIL_BASIC,0); 
    }

    r_max_bins = -1.e200;
    r_min_bins =  1.e200;

#if( TOP_TYPE_CHOICE == TOP_SPHERICAL ) 
    /* use the global min and max radii as the bounds for the bins, we assume r is always monotonic with "i" : */
    i = N1S;
    N2_LOOP N3_LOOP { 
      get_coord(i,j,k,FACE1,ncurr,coords);
      if( coords->r < r_min_bins ) { r_min_bins = coords->r ; }
    }
    i = N1E+1;
    N2_LOOP N3_LOOP { 
      get_coord(i+1,j,k,FACE1,ncurr,coords);
      if( coords->r > r_max_bins ) { r_max_bins = coords->r ; }
    }
#else
    /* use the global min and max radii as the bounds for the bins, no assumption on r is made : */
    LOOP { 
      get_coord(i,j,k,CENT,ncurr,coords);
      if( coords->r < r_min_bins ) { r_min_bins = coords->r ; }
      if( coords->r > r_max_bins ) { r_max_bins = coords->r ; }
    }
#endif

    mpi_global_min(&r_min_bins); 
    mpi_global_max(&r_max_bins); 

    xp1_min_bins = log(r_min_bins / r0_bins); 
    xp1_max_bins = log(r_max_bins / r0_bins); 
    dxp1_bins = (xp1_max_bins - xp1_min_bins) / n_r_bins; 
    inv_dxp1_bins = 1./dxp1_bins;

    if( myid == printer_pid ) { 
      fprintf(stdout,"\n##################################################\n");
      fprintf(stdout,"  Bin arrays PARAMETERS \n------------------------------------\n");
      fprintf(stdout,"\t  n_r_bins      =  %10d \n"   ,n_r_bins   ); 
      fprintf(stdout,"\t  n_phi_bins    =  %10d \n"   ,n_phi_bins ); 
      fprintf(stdout,"\t  r0_bins       =  %28.18e \n",r0_bins       ); 
      fprintf(stdout,"\t  r_min_bins    =  %28.18e \n",r_min_bins    ); 
      fprintf(stdout,"\t  r_max_bins    =  %28.18e \n",r_max_bins    ); 
      fprintf(stdout,"\t  xp1_min_bins  =  %28.18e \n",xp1_min_bins  ); 
      fprintf(stdout,"\t  xp1_max_bins  =  %28.18e \n",xp1_max_bins  ); 
      fprintf(stdout,"\t  dxp1_bins     =  %28.18e \n",dxp1_bins     ); 
      fprintf(stdout,"\t  inv_dxp1_bins =  %28.18e \n",inv_dxp1_bins ); 
      fprintf(stdout,"\t  phi_min_bins  =  %28.18e \n",phi_min_bins  ); 
      fprintf(stdout,"\t  phi_max_bins  =  %28.18e \n",phi_max_bins  ); 
      fprintf(stdout,"\t  dphi_bins     =  %28.18e \n",dphi_bins     ); 
      fprintf(stdout,"\t  inv_dphi_bins =  %28.18e \n",inv_dphi_bins ); 
      fprintf(stdout,"\n##################################################\n");
      fflush(stdout);
    }
  }

  /* Now setup the arrays: */
  FREE(      r_bins );
  FREE(     dr_bins );
  FREE( inv_dr_bins );
  FREE(    xp1_bins );
  FREE(    phi_bins );

  ALLOC_ARRAY(      r_bins,  n_r_bins); 
  ALLOC_ARRAY(  r_bins_mid,  n_r_bins); 
  ALLOC_ARRAY(     dr_bins,  n_r_bins); 
  ALLOC_ARRAY( inv_dr_bins,  n_r_bins); 
  ALLOC_ARRAY(    xp1_bins,  n_r_bins); 
  ALLOC_ARRAY(    phi_bins,n_phi_bins); 
  ALLOC_ARRAY(phi_bins_mid,n_phi_bins); 

  for(i=0;i<n_r_bins;i++) { 
    xp1_bins[i] = xp1_min_bins + i*dxp1_bins;
    r_bins[i] = r0_bins * exp(xp1_bins[i]); 
  }
  for(i=0;i<(n_r_bins-1);i++) { 
    dr_bins[i] = r_bins[i+1] - r_bins[i];
  }
  i = n_r_bins-1;
  dr_bins[i] = r0_bins*exp(xp1_min_bins + (i+1)*dxp1_bins) - r_bins[i] ;

  for(i=0;i<n_r_bins;i++) { 
    r_bins_mid[i] = r_bins[i] + 0.5*dr_bins[i];
  }
  for(i=0;i<n_r_bins;i++) { 
    inv_dr_bins[i] = 1./dr_bins[i];
  }

  for(i=0;i<n_phi_bins;i++) { 
    phi_bins[i] = phi_min_bins + i * dphi_bins; 
  }
  for(i=0;i<n_phi_bins;i++) { 
    phi_bins_mid[i] = phi_min_bins + (i+0.5) * dphi_bins; 
  }


  TRACE_END;

  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  r_to_xp1_bins;
 ------------------
   -- calculates the xp1 coordainates from the radial coordinates, for binning purposes;
   -- responsible for the allocation of the xp1 array;
   -- returns with the pointer to the head of the array;  
      the new array has the same length "n_in" as the radial coordinates array;

**********************************************************************************************/
double *r_to_xp1_bins(int n_in, double *r_in) 
{
  int i; 
  double *xp1;

  ALLOC_ARRAY(xp1,n_in);

  for(i=0;i<n_in;i++) {  xp1[i] = log(r_in[i]/r0_bins);   }

  return(xp1);
}

/**********************************************************************************************/
/**********************************************************************************************
  set_bin_weights_r():
 --------------------
   -- calculates how the data overlaps with the bins and the weights needed to bin the data 
      correctly; 

**********************************************************************************************/
void set_bin_weights_r(int n_in, double *r_beg, double *r_end, double ***weights, int **ibins_beg, int **ibins_end) 
{
  int i,j,n,ibin1,ibin2; 
  double x1, x2, r1, r2, dr, rb, re;

  TRACE_BEG;


  ALLOC_ARRAY(  (*weights),  n_in);
  ALLOC_ARRAY( (*ibins_beg), n_in);
  ALLOC_ARRAY( (*ibins_end), n_in);

  const int last_index = n_r_bins-1;

  for(i=0;i<n_in;i++) { 
    r1 = r_beg[i];
    r2 = r_end[i];

    x1 = log(r1/r0_bins);
    x2 = log(r2/r0_bins);

    /* Skip the bin if the data falls outside the binning range: */
    if( (x2 < xp1_min_bins) || (x1 > xp1_max_bins) ) { 
      (*ibins_beg)[i] = 1;
      (*ibins_end)[i] = 0;
      (*weights)[i] = NULL; 
      continue; 
    }
    
    ibin1 = (int) ((x1 - xp1_min_bins)*inv_dxp1_bins);
    ibin2 = (int) ((x2 - xp1_min_bins)*inv_dxp1_bins);

    ibin1 = MAX(  0       ,ibin1); 
    ibin1 = MIN(last_index,ibin1); 
    ibin2 = MAX(  0       ,ibin2); 
    ibin2 = MIN(last_index,ibin2); 

    n = ibin2 - ibin1 + 1;
    
    ALLOC_ARRAY((*weights)[i], n);

    (*ibins_beg)[i] = ibin1; 
    (*ibins_end)[i] = ibin2; 

    //    fprintf(stdout,"set_bin1: %d %d %d %d %d %26.16e %26.16e %26.16e %26.16e\n",n_in,i,ibin1,ibin2,n,r1,r2,x1,x2); fflush(stdout);

    n = 0; 
    for(j=ibin1;j<=ibin2;j++ ) { 
      dr = dr_bins[j];
      rb = r_bins[j];
      re = rb + dr;

      rb = MAX(rb,r1); 
      re = MIN(re,r2);

      (*weights)[i][n]  = (re-rb)*inv_dr_bins[j];
      //      fprintf(stdout,"set_bin2: %d %d %d %26.16e %26.16e %26.16e %26.16e %26.16e\n",i,j,n,dr,r_bins[j],rb,re,(*weights)[i][n]); fflush(stdout);
      
      n++;
    }
  }

  TRACE_END;

  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  set_bin_weights_phi():
 --------------------
   -- calculates how the data overlaps with the bins and the weights needed to bin the data 
      correctly; 

**********************************************************************************************/
void set_bin_weights_phi(int n_in, double *phi_beg, double *phi_end, double ***weights, int **ibins_beg, int **ibins_end)
{
  int i,j,n,ibin1,ibin2; 
  double r1, r2, dr, rb, re;

  TRACE_BEG;


  ALLOC_ARRAY(   (*weights), n_in);
  ALLOC_ARRAY( (*ibins_beg), n_in);
  ALLOC_ARRAY( (*ibins_end), n_in);

  const int last_index = n_phi_bins-1;

  for(i=0;i<n_in;i++) { 
    r1 = phi_beg[i];
    r2 = phi_end[i];

    if( (r2 < phi_min_bins) || (r1 > phi_max_bins) ) { 
      fprintf(stdout,"set_bin_weights_phi(): data lies outside phi bin range, you may want to check how coordinates are set, %26.16e %26.16e %26.16e %26.16e ",
	      r1,r2,phi_min_bins,phi_max_bins); 
      fflush(stdout);
      fail(FAIL_BASIC,0);
    }

    ibin1 = (int) ((r1 - phi_min_bins)*inv_dphi_bins);
    ibin2 = (int) ((r2 - phi_min_bins)*inv_dphi_bins);

    ibin1 = MAX(  0       ,ibin1); 
    ibin1 = MIN(last_index,ibin1); 
    ibin2 = MAX(  0       ,ibin2); 
    ibin2 = MIN(last_index,ibin2); 

    n = ibin2 - ibin1 + 1;
    
    ALLOC_ARRAY((*weights)[i], n);

    (*ibins_beg)[i] = ibin1; 
    (*ibins_end)[i] = ibin2; 

    n = 0; 
    for(j=ibin1;j<=ibin2;j++ ) { 
      rb = phi_bins[j];
      re = rb + dphi_bins;

      rb = MAX(rb,r1); 
      re = MIN(re,r2);

      (*weights)[i][n]  = (re-rb)*inv_dphi_bins;
      n++;
    }
  }

  TRACE_END;

  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  bin_data_in_r();
 ------------------
   -- using the existing weights array, bin the data; 
   -- calculates ibins_beg[], ibins_end[], weights[] locally  and store them for future use;
   -- bin_state dictates what this routine does:

          bin_state =  0  :  normal binning operation, assumes local arrays are set; 
                    =  1  :  sets the local arrays first; 

**********************************************************************************************/
void bin_data_in_r(int bin_state, int n_in, double *data_in, double *r_beg, double *r_end, double *data_out)
{
  int i,j,n,ibin1,ibin2;
  static int n_in_old=0; 

  TRACE_BEG;

  /*************************************************************************************
    Handle different ways of calling this routine: 
   *************************************************************************************/
  switch( bin_state ) { 

    /* do nothing for special for this case :*/ 
  case 0 : break;


    /* setup the local arrays for present and future use: */
  case 1 : 
    if( bin_weights_r != NULL ) { 
      for(i=0;i<n_in_old;i++) { 
	if( bin_weights_r[i] != NULL ) { FREE(bin_weights_r[i]); }
      }
      FREE(bin_weights_r); 
    }
    if( ibins_beg != NULL ) { FREE( ibins_beg ); }
    if( ibins_end != NULL ) { FREE( ibins_end ); }

    n_in_old = n_in; 

    set_bin_weights_r(n_in, r_beg, r_end, &bin_weights_r, &ibins_beg, &ibins_end);

    break;
  }

  if( n_in_old != n_in ) {
    fprintf(stdout,"bin_data_in_r(): problem with sequence or mode of calls to this routine, n_in_old != n_in:  %d %d \n",
	    n_in_old,n_in);  fflush(stdout); fail(FAIL_BASIC,0); 
  }

  /*************************************************************************************
    Perform the binning procedure: 
   *************************************************************************************/
  for(i=0;i<n_in;i++) { 
    ibin1 = ibins_beg[i]; 
    ibin2 = ibins_end[i]; 

    n = 0;
    for(j=ibin1;j<=ibin2;j++) { 
      data_out[j] += bin_weights_r[i][n] * data_in[i];
      n++;
    }
  }

  TRACE_END;

  return;
}

/**********************************************************************************************/
/**********************************************************************************************
  bin_data_in_r_phi();
 ------------------
   -- using the existing weights array, bin the data; 
   -- calculates ibins_beg[], ibins_end[], weights[] locally  and store them for future use;
   -- bin_state dictates what this routine does:

          bin_state =  0  :  normal binning operation, assumes local arrays are set; 
                    =  1  :  sets the local arrays first; 

**********************************************************************************************/
void bin_data_in_r_phi(int bin_state, int n_in, double *data_in, 
		       double   *r_beg, double   *r_end, 
		       double *phi_beg, double *phi_end, 
		       double **data_out)
{
  int i,j,k,n1,n2,ibin1,ibin2,kbin1,kbin2;
  static int n_in_old=0; 

  TRACE_BEG;

  /*************************************************************************************
    Handle different ways of calling this routine: 
   *************************************************************************************/
  switch( bin_state ) { 

    /* do nothing for special for this case :*/ 
  case 0 : break;


    /* setup the local arrays for present and future use: */
  case 1 : 
    if( bin_weights_r != NULL ) { 
      for(i=0;i<n_in_old;i++) { 
	if( bin_weights_r[i] != NULL ) { FREE(bin_weights_r[i]); }
      }
      FREE(bin_weights_r); 
    }
    if( ibins_beg != NULL ) { FREE( ibins_beg ); }
    if( ibins_end != NULL ) { FREE( ibins_end ); }

    if( bin_weights_phi != NULL ) { 
      for(i=0;i<n_in_old;i++) { 
	if( bin_weights_phi[i] != NULL ) { FREE(bin_weights_phi[i]); }
      }
      FREE(bin_weights_phi); 
    }
    if( kbins_beg != NULL ) { FREE( kbins_beg ); }
    if( kbins_end != NULL ) { FREE( kbins_end ); }

    n_in_old = n_in; 

    set_bin_weights_r(  n_in,   r_beg,   r_end, &bin_weights_r,   &ibins_beg, &ibins_end);
    set_bin_weights_phi(n_in, phi_beg, phi_end, &bin_weights_phi, &kbins_beg, &kbins_end);

    break;
  }

  if( n_in_old != n_in ) {
    fprintf(stdout,"bin_data_in_r_phi(): problem with sequence or mode of calls to this routine, n_in_old != n_in:  %d %d \n",
	    n_in_old,n_in);  fflush(stdout); fail(FAIL_BASIC,0); 
  }

  /*************************************************************************************
    Perform the binning procedure: 
   *************************************************************************************/
  for(i=0;i<n_in;i++) { 
    ibin1 = ibins_beg[i]; 
    ibin2 = ibins_end[i]; 
    kbin1 = kbins_beg[i]; 
    kbin2 = kbins_end[i]; 

    n1 = 0;
    for(j=ibin1;j<=ibin2;j++) { 
      n2 = 0;
      for(k=kbin1;k<=kbin2;k++) { 
	data_out[j][k] += bin_weights_r[i][n1] * bin_weights_phi[i][n2] * data_in[i];
	n2++;
      }
      n1++;
    }
  }

  TRACE_END;

  return;
}

