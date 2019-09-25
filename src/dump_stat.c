
#include "decs.h"

#if( MAKE_STAT ) 
#include <hdf5.h>

#if( USEMPI )
#include "mpi.h"
#endif 

#define IND_LOOP for(ind=0;ind<NCELLS;ind++) 

static hid_t file_id;  /* ID variable for hdf5 statistics file */ 
  
extern void     myH5_write_gfunc( hid_t loc_id, char *name, double *value ) ;
extern void     myH5_write_int_gfunc( hid_t loc_id, char *name, int *value );
extern void myH5_write_scalar2( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
extern hid_t myH5_Fcreate(char *name,int overwrite);
extern void     setup_hdf5(void); 

/**********************************************************************************************/
/**********************************************************************************************
  dump_stat(): 
 -----------
  -- routine that dumps various gridfunctions that store statistical information 
        about the time evolution; 

  --  data dumped to an hdf5 file called  <RUN_TAB>.STAT.<timeid>.h5

  -- responsible for doing the necessary calculations required for our in situ 
        statistical analyses of the simulation;

   -- we clobber existing files;

   -- full dump file names take the form of : 

        RUN_TAG.STAT.<dump_id>.h5

          where  
 
            <dump_id> = time index number 

     -- use the USE_CHUNKS macro to control whether or not to use chunks in parallel IO;
         -- we assume that the chunk size is the memory of one gridfunction in the 
             local domain (i.e. N1*N2*N3);

    -- integrated quantities are dumped per timestep in file "<RUN_TAG>_stat.out" 
           -- these include the number of failures per type; 
           -- the amout of dU per failure 

    -- there are several types of error states that a cell can be in: 
         0) inversion successful
         1) inversion unsuccesful (bad state or non-convergent soln, use old P[])
         2) gamma > GAMMAMAX
         3) rho,uu < rhoflr,uuflr  or  pv[UU] > pv[RHO]*MAX_UU_RHO_RATIO
         4) cannot calculate physical gamma; 

         -- obviously  0 & 1  cannot both happen
         -- State #1 is exclusive to other states since it leaves the P[] unchanged; 
         -- 0 & 4 are exclusive since the inversion calculation is only successful if 
             a physical velocity is obtained (e.g. a good value of gamma can be calculated);  
         -- 2 & 4 are exclusive since gamma > GAMMAMAX is a valid value of gamma (unless gamma=inf which never happens)
         -- 4 has never before occured in recent simulations, so let's ignore it
         -- for a cell to have 2,3,4 to happen, it must also have 0
         -- so those State 0 cells can have 2,3,4 or neither. 
         -- currently, a fixup interpolation (ala fixup_utoprim()) is done for 1,2,4
         -- hence, we need only calculate U(P) after  0, 3 (is fixed) and after the 
             fixup interpolation; 
         -- the changes are only calculated where there have been changes, and the 
             cells that have changed will have non-zero values of U(P)[0]; 
         -- also calculates the integrated change in the conserved quantities 

    -- statistical gridfunctions (S = space is allocated, T = temporary, no space allocated ): 
              U_pre[0-7]   : conserved variable before timestep (to measure change from update);
            U_fin[0-4]     : conserved variables calculated after all fixups (B^i never change)
            U_avg[0-7]     : average of conserved variables (U_fin[]) since last dump;
           dU_stat[0-4]_0  : FAIL_FLOOR change from floor
           dU_stat[0-4]_1  : FAIL_TMAX change from Tmax
           dU_stat[0-4]_2  : FAIL_GAMMA_MAX change from gamma ceiling 
           dU_stat[0-4]_3  : FAIL_USE_FULL_INV change from full inversion error, i.e.  U - U(P(U))
           dU_stat[0-4]_4  : FAIL_USE_EE change from using fixup_entropy_eq()
           dU_stat[0-4]_5  : FAIL_USE_INTERP_PRIM change from using fixup_interp_prim()
           dU_stat[0-4]_6  : FAIL_USE_INTERP_V change from using fixup_interp_v()
           dU_stat[0-4]_7  : FAIL_INTERP_PRIM change from using old solution since fixup_interp_prim() failed to interpolate
           dU_stat[0-4]_8  : FAIL_CUTOUT_PRESSURE change from using cutout fix for the pressure 
           dU_stat[0-4]_9  : FAIL_CUTOUT_GAMMA change from using cutout fix for gamma
           dU_stat[0-4]_10 : FAIL_CUTOUT_BC change from using cutout BC condition

    -- integrations/time-averages (shown in order they appear per line of the "stat".out file):
	       nstep       : timestep iteration number 
	           t       : time 
	          dt       : timestep size 
     nfailtot[0:N_NFAIL-1] : volume-tallied total number of occurrences of each type of failure
              sint_U[0-7]  : spatially integrated U_fin[]; 
       sint_dU[0-4][0-10]  : spatially integrated dU_stat[][]

    -- full hdf5 dump is done whenever stat_call_code = OUT_STAT or OUT_FINAL
          -- stat_call_code set in diag()
             
**********************************************************************************************/
void dump_stat(void) 
{
  int i,j,k,l,d,g,ind;
  int itype;
  char var_name[200];

  int ncells_loc_old;

  int *itmp;

  const int nvars = NP+NPH*N_STAT;
  double sint_U[NP], sint_dU[NPH][N_STAT], sint_send[nvars], sint_recv[nvars];

  unsigned long int nfailtot_sum[N_NFAIL];

  const int dump_nfail[] = { 
    1,                                  /* NFAIL_FLOOR	         ( 0) */
    USE_TEMPERATURE_MAX,		/* NFAIL_TMAX            ( 1) */
    1,					/* NFAIL_GAMMA_MAX	 ( 2) */
    1,					/* NFAIL_GAMMA_CALC      ( 3) */
    1,					/* NFAIL_USE_FULL_INV    ( 4) */
    USE_ENTROPY_EQ,			/* NFAIL_USE_EE          ( 5) */
    1,					/* NFAIL_UTOPRIM         ( 6) */
    USE_ENTROPY_EQ,			/* NFAIL_UTOPRIM_EE      ( 7) */
    1,					/* NFAIL_USE_INTERP_PRIM ( 8) */
    1, 					/* NFAIL_USE_INTERP_V    ( 9) */
    1, 					/* NFAIL_INTERP_PRIM     (10) */
    1,					/* NFAIL_INTERP_V        (11) */
    FIX_CUTOUT_PRESSURE,		/* NFAIL_CUTOUT_PRESSURE (12) */
    FIX_CUTOUT_GAMMA,			/* NFAIL_CUTOUT_GAMMA    (13) */
    BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT /* NFAIL_CUTOUT_BC       (14) */
  };


  char stat_filename[200]; 
  static FILE *stat_file;

  static int nstep_last; 
  static double t_last;
  double inv_steps_passed, dt_stat, inv_dt_stat;

  static int local_first_time = 1; 

  void setup_stat_hdf5(void);

  TRACE_BEG;
  
  //  exit_status();  /* Make sure that all processes are alive */

  /******************************************************************************
     INITIALIZATIONS : 
  ******************************************************************************/ 

  /* Set up persistent data used for all time by all hdf5 routines: */
  setup_hdf5();    
    
  /* Initialize file for writing : */
  if( local_first_time ) { 
    /* Make sure that dump_nfail[] array is set correctly : */
    if( NFAIL_FLOOR	      !=  0 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_FLOOR	          %d\n",NFAIL_FLOOR	     ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_TMAX            !=  1 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_TMAX             %d\n",NFAIL_TMAX           ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_GAMMA_MAX	      !=  2 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_GAMMA_MAX	  %d\n",NFAIL_GAMMA_MAX	     ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_GAMMA_CALC      !=  3 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_GAMMA_CALC       %d\n",NFAIL_GAMMA_CALC     ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_USE_FULL_INV    !=  4 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_USE_FULL_INV     %d\n",NFAIL_USE_FULL_INV   ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_USE_EE          !=  5 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_USE_EE           %d\n",NFAIL_USE_EE         ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_UTOPRIM         !=  6 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_UTOPRIM          %d\n",NFAIL_UTOPRIM        ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_UTOPRIM_EE      !=  7 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_UTOPRIM_EE       %d\n",NFAIL_UTOPRIM_EE     ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_USE_INTERP_PRIM !=  8 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_USE_INTERP_PRIM  %d\n",NFAIL_USE_INTERP_PRIM); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_USE_INTERP_V    !=  9 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_USE_INTERP_V     %d\n",NFAIL_USE_INTERP_V   ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_INTERP_PRIM     != 10 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_INTERP_PRIM      %d\n",NFAIL_INTERP_PRIM    ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_INTERP_V        != 11 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_INTERP_V         %d\n",NFAIL_INTERP_V       ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_CUTOUT_PRESSURE != 12 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_CUTOUT_PRESSURE  %d\n",NFAIL_CUTOUT_PRESSURE); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_CUTOUT_GAMMA    != 13 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_CUTOUT_GAMMA     %d\n",NFAIL_CUTOUT_GAMMA   ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( NFAIL_CUTOUT_BC       != 14 ) { fprintf(stderr,"dump_stat(): bad value for NFAIL_CUTOUT_BC        %d\n",NFAIL_CUTOUT_BC      ); fflush(stderr); fail(FAIL_BASIC,0); }

    if( myid == out_pid[OUT_STAT] ) { 
      sprintf(stat_filename,"%s_stat.out",RUN_TAG);
      stat_file = fopen(stat_filename,"at");
      if(stat_file==NULL) {
	fprintf(stderr,"dump_stat():  error opening %s file\n",stat_filename) ;
	fail( FAIL_BASIC,0 );
      }
    }

#if( DUMP_ALL_STAT )
    PLOOP IND_LOOP {  U_avg[l][ind] = 0.; }
#endif

    nstep_last = nstep;
    t_last = t;
    local_first_time = 0;
  }

  if( stat_call_code == OUT_INIT ) {  
    TRACE_END;
    return; 
  }

  /*****************************************************************************************
    Compile this timestep's statistical output for the *_stat.out  file. 
  ******************************************************************************************/ 
#if( DUMP_ALL_STAT )
  PLOOP IND_LOOP { U_avg[l][ind] += U_out[l][ind]; }
  
  /* -- sint_U[] = spatially integrated U_fin[] or U_out; */
  PLOOP { 
    sint_U[l] = 0.; 
    IND_LOOP { sint_U[l] += U_out[l][ind]; }
    sint_U[l] *= dV; 
  }
  
  /* -- sint_dU[] = spatially integrated dU_stat[]; */
  NPH_LOOP for(g=0; g<N_STAT; g++) { 
    sint_dU[l][g] = 0.; 
    IND_LOOP { sint_dU[l][g] += dU_stat[l][g][ind]; }
    sint_dU[l][g] *= dV; 
  }

  /* Integrate these values over all subdomains and send only to the master node : */
  ind = 0; 
  PLOOP                            { sint_send[ind++] = sint_U[l];     } 
  NPH_LOOP for(g=0; g<N_STAT; g++) { sint_send[ind++] = sint_dU[l][g]; } 

  if( ind != nvars ) { 
    fprintf(stderr,"dump_stat(): problem with counting elements for sint_send/recv : %d %d \n",
	    ind, nvars);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  mpi_sum_to_head(sint_send,sint_recv,nvars);
#endif
  
  /******************************************************************************
      OUTPUT THE TIMESTEP INFORMATION: 
  ******************************************************************************/ 
  for(i=0;i<N_NFAIL;i++) {  nfailtot_sum[i] = 0; } 
  mpi_ulong_sum_to_head(nfailtot,nfailtot_sum,N_NFAIL);

  if( myid == out_pid[OUT_STAT] ) { 
    fprintf(stat_file,"%d %28.18e %28.18e ", nstep, t,dx[0]);   //  1 -  3 
    for(l=0;l<N_NFAIL;l++) { fprintf(stat_file," %ld ",nfailtot_sum[l]); }   // 4 - 18
#if( DUMP_ALL_STAT )
    for(i=0;i<nvars;i++) { fprintf(stat_file," %28.18e ",sint_recv[i]); }   // 19 - (19+nvars-1)
#endif
    fprintf(stat_file,"\n");
    fflush(stat_file);
  }

  /*****************************************************************************************
   Write the statistical gridfunctions if it is time : 
  *****************************************************************************************/
  if( stat_call_code == OUT_STAT || stat_call_code == OUT_FINAL ) { 
    setup_stat_hdf5();
    inv_dt_stat = (t == t_last) ?   1./dx[0]  :  1./(t - t_last); 
    inv_steps_passed  = (nstep == nstep_last) ?  1. : 1./((double) (nstep - nstep_last));

#if( DUMP_ALL_STAT )
    /* -- U_pre (U before the update)  */
    PLOOP { 
      sprintf(var_name,"U_pre%d",l);  myH5_write_gfunc(file_id, var_name, U_pre[l] ); 
    }
    
    /* -- U_out (U after all fixups)  */
    PLOOP { 
      sprintf(var_name,"U_fin%d",l);  myH5_write_gfunc(file_id, var_name, U_out[l] ); 
    }
    
    /* -- U_avg (average of U_fin since last full stat dump)  */
    PLOOP IND_LOOP { U_avg[l][ind] *= inv_dt_stat; } 
    PLOOP { 
      sprintf(var_name,"U_avg%d",l);  myH5_write_gfunc(file_id, var_name, U_avg[l] ); 
    }
    PLOOP IND_LOOP  U_avg[l][ind] = 0.;

    
    /* -- dU_stat[][]  */
    NPH_LOOP for(g=0; g<N_STAT; g++) {
      sprintf(var_name,"dU_stat%d_%d",l,g); myH5_write_gfunc(file_id, var_name, dU_stat[l][g]);
      IND_LOOP  dU_stat[l][g][ind] = 0.;
    }
#endif

    /* -- nfail[]  */
    ALLOC_ARRAY(itmp,NCELLS);
    for(l=0;l<N_NFAIL;l++) if(dump_nfail[l]) { 
      ind = 0;
      LOOP  { itmp[ind++] = nfail[l][i][j][k]; } 
      sprintf(var_name,"int_nfail%02d",l); myH5_write_int_gfunc(file_id, var_name, itmp);
    }
    for(l=0;l<N_NFAIL;l++)  ALL_LOOP { nfail[l][i][j][k] = 0; } 
    for(l=0;l<N_NFAIL;l++) { nfailtot[l] = 0; } 
    FREE(itmp);
    

    nstep_last = nstep; 
    t_last = t;

#if( GATHER_IO ) 
    if( myid == out_pid[OUT_STAT] ) { 
#else
    {
#endif
      H5Fclose(file_id);      /* Terminate access to the file. */
    }

  }

  /******************************************************************************
      FINALIZE: CLOSE FILES, ETC. 
  ******************************************************************************/ 

  /* Make sure we close the file if this is the last call */ 
  if( stat_call_code == OUT_FINAL ) { if( myid == out_pid[OUT_STAT] )  {  fclose( stat_file ); } }

#if( USEMPI ) 
  //  exit_status();
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 

  TRACE_END;

  return;

}



/**********************************************************************************************/
/**********************************************************************************************
  setup_stat_hdf5():
 ------------------
    -- creates the hdf5 file, writes the header information and closes the header;

**********************************************************************************************/
void setup_stat_hdf5(void)
{
  int n1,n2,n3;
  char hdf_name[200];  /* Name of current hdf5 statistics file */
  hid_t header_id; 

  /* Set the name of this time's hdf file : */
#if( USE_MPI_IO  || (!USEMPI) || GATHER_IO )
  sprintf(hdf_name,"%s/%s.%s.%06d.h5",DIR_out[OUT_STAT],RUN_TAG,"STAT",N_out[OUT_STAT]); 
#else 
  sprintf(hdf_name,"%s/%s.%s.%06d.p%05d.h5",DIR_out[OUT_STAT],RUN_TAG,"STAT",N_out[OUT_STAT],myid); 
#endif
  
  if( myid == printer_pid ) {   fprintf(stdout,"Dumping STAT %s .... \n",hdf_name); fflush(stdout); }
    
  /******************************************************************************
      Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Create a new, clobbering file using serial or parallel properties */
  file_id = myH5_Fcreate(hdf_name,0);
  
  /* Create a subgroup to which we attach attributes or header information */
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_STAT] ) { 
#else
   {
#endif
     header_id = H5Gcreate(file_id, "Header", 0);
     if( header_id  < 0 ) { 
       fprintf(stderr,"setup_stat_hdf5(): Cannot create Header group in file %s \n", hdf_name);
       fflush(stderr);     fail(FAIL_HDF,0);
     }
  
     /******************************************************************************
      Dump the grid parameters into the HDF file : 
     ******************************************************************************/ 
     n1 = N1;   n2 = N2;  n3 = N3; 
     myH5_write_scalar2(header_id, "N1"         , H5T_NATIVE_INT   ,  &n1             );
     myH5_write_scalar2(header_id, "N2"         , H5T_NATIVE_INT   ,  &n2             );
     myH5_write_scalar2(header_id, "N3"         , H5T_NATIVE_INT   ,  &n3             );
     myH5_write_scalar2(header_id, "totalsize0" , H5T_NATIVE_INT   ,  &(totalsize[0]) );
     myH5_write_scalar2(header_id, "totalsize1" , H5T_NATIVE_INT   ,  &(totalsize[1]) );
     myH5_write_scalar2(header_id, "totalsize2" , H5T_NATIVE_INT   ,  &(totalsize[2]) );
     myH5_write_scalar2(header_id, "totalsize3" , H5T_NATIVE_INT   ,  &(totalsize[3]) );
     myH5_write_scalar2(header_id, "cpupos0"    , H5T_NATIVE_INT   ,  &(cpupos[0])    );
     myH5_write_scalar2(header_id, "cpupos1"    , H5T_NATIVE_INT   ,  &(cpupos[1])    );
     myH5_write_scalar2(header_id, "cpupos2"    , H5T_NATIVE_INT   ,  &(cpupos[2])    );
     myH5_write_scalar2(header_id, "cpupos3"    , H5T_NATIVE_INT   ,  &(cpupos[3])    );
     myH5_write_scalar2(header_id, "startx0"    , H5T_NATIVE_DOUBLE,  &(startx[0])    );
     myH5_write_scalar2(header_id, "startx1"    , H5T_NATIVE_DOUBLE,  &(startx[1])    );
     myH5_write_scalar2(header_id, "startx2"    , H5T_NATIVE_DOUBLE,  &(startx[2])    );
     myH5_write_scalar2(header_id, "startx3"    , H5T_NATIVE_DOUBLE,  &(startx[3])    );
     myH5_write_scalar2(header_id, "dx0"        , H5T_NATIVE_DOUBLE,  &(dx[0])        );
     myH5_write_scalar2(header_id, "dx1"        , H5T_NATIVE_DOUBLE,  &(dx[1])        );
     myH5_write_scalar2(header_id, "dx2"        , H5T_NATIVE_DOUBLE,  &(dx[2])        );
     myH5_write_scalar2(header_id, "dx3"        , H5T_NATIVE_DOUBLE,  &(dx[3])        );
     myH5_write_scalar2(header_id, "gridlength0", H5T_NATIVE_DOUBLE,  &(GridLength[0]));
     myH5_write_scalar2(header_id, "gridlength1", H5T_NATIVE_DOUBLE,  &(GridLength[1]));
     myH5_write_scalar2(header_id, "gridlength2", H5T_NATIVE_DOUBLE,  &(GridLength[2]));
     myH5_write_scalar2(header_id, "gridlength3", H5T_NATIVE_DOUBLE,  &(GridLength[3]));
     myH5_write_scalar2(header_id, "t"          , H5T_NATIVE_DOUBLE,  &t              );
     myH5_write_scalar2(header_id,  "nstep"     , H5T_NATIVE_INT   ,  &nstep          );
  
     H5Gclose(header_id);    /* Terminate access to the group /Header. */
   }

  return;

}

/****************************************************************************************/
/*****************************************************************************************
  calc_dU_stat():
  --------------
     -- calculates the change in the conserved variables from the fixup that just
         took place; 
     -- only used if MAKE_STAT is set;
     -- should only be called at the last temporal substep (aka corrector step)
*****************************************************************************************/
void calc_dU_stat(int fail_type, int state)
{
  int i,j,k,l,g;
  int ind_U;  
  double U_new[NP];  
  struct of_state q;  
  struct of_geom *geom;

  static int index_fail_to_dU[N_FAIL] = {0};
  static int local_first_time=1;


  //  fprintf(stdout,"calc_dU_stat begin pid=%d   nstep=%d\n",myid,nstep); fflush(stdout);

  if( local_first_time ) { 
    /* First check to see that the FAIL_* macros have not changed since we explicitly use their values below :*/
    if( FAIL_FLOOR           !=  1 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_FLOOR            %d\n",FAIL_FLOOR          ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_TMAX            !=  2 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_TMAX             %d\n",FAIL_TMAX           ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_GAMMA_MAX       !=  3 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_GAMMA_MAX        %d\n",FAIL_GAMMA_MAX      ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_GAMMA_CALC      !=  4 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_GAMMA_CALC       %d\n",FAIL_GAMMA_CALC     ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_USE_FULL_INV    !=  5 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_USE_FULL_INV     %d\n",FAIL_USE_FULL_INV   ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_USE_EE          !=  6 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_USE_EE           %d\n",FAIL_USE_EE         ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_UTOPRIM         !=  7 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_UTOPRIM          %d\n",FAIL_UTOPRIM        ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_UTOPRIM_EE      !=  8 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_UTOPRIM_EE       %d\n",FAIL_UTOPRIM_EE     ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_USE_INTERP_PRIM !=  9 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_USE_INTERP_PRIM  %d\n",FAIL_USE_INTERP_PRIM); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_USE_INTERP_V    != 10 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_USE_INTERP_V     %d\n",FAIL_USE_INTERP_V   ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_INTERP_PRIM     != 11 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_INTERP_PRIM      %d\n",FAIL_INTERP_PRIM    ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_INTERP_V        != 12 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_INTERP_V         %d\n",FAIL_INTERP_V       ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_CUTOUT_PRESSURE != 13 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_CUTOUT_PRESSURE  %d\n",FAIL_CUTOUT_PRESSURE); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_CUTOUT_GAMMA    != 14 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_CUTOUT_GAMMA     %d\n",FAIL_CUTOUT_GAMMA   ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_CUTOUT_BC       != 15 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_CUTOUT_BC        %d\n",FAIL_CUTOUT_BC      ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_VCHAR_NEG       != 16 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_VCHAR_NEG        %d\n",FAIL_VCHAR_NEG      ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_VCHAR_SUPER     != 17 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_VCHAR_SUPER      %d\n",FAIL_VCHAR_SUPER    ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_VCHAR_DISCR     != 18 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_VCHAR_DISCR      %d\n",FAIL_VCHAR_DISCR    ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_BASIC           != 19 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_BASIC            %d\n",FAIL_BASIC          ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_METRIC          != 20 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_METRIC           %d\n",FAIL_METRIC         ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_RESTART         != 21 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_RESTART          %d\n",FAIL_RESTART        ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_MPI_BASIC       != 22 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_MPI_BASIC        %d\n",FAIL_MPI_BASIC      ); fflush(stderr); fail(FAIL_BASIC,0); }
    if( FAIL_HDF             != 23 ) { fprintf(stderr,"calc_dU_stat(): bad value for FAIL_HDF              %d\n",FAIL_HDF            ); fflush(stderr); fail(FAIL_BASIC,0); }

    index_fail_to_dU[0                    ] = -1;
    index_fail_to_dU[FAIL_FLOOR           ] =  0;
    index_fail_to_dU[FAIL_TMAX            ] =  1;
    index_fail_to_dU[FAIL_GAMMA_MAX       ] =  2;
    index_fail_to_dU[FAIL_GAMMA_CALC      ] = -1;
    index_fail_to_dU[FAIL_USE_FULL_INV    ] =  3;
    index_fail_to_dU[FAIL_USE_EE          ] =  4;
    index_fail_to_dU[FAIL_UTOPRIM         ] = -1;
    index_fail_to_dU[FAIL_UTOPRIM_EE      ] = -1;
    index_fail_to_dU[FAIL_USE_INTERP_PRIM ] =  5;
    index_fail_to_dU[FAIL_USE_INTERP_V    ] =  6;
    index_fail_to_dU[FAIL_INTERP_PRIM     ] =  7;
    index_fail_to_dU[FAIL_INTERP_V        ] = -1;
    index_fail_to_dU[FAIL_CUTOUT_PRESSURE ] =  8;
    index_fail_to_dU[FAIL_CUTOUT_GAMMA    ] =  9;
    index_fail_to_dU[FAIL_CUTOUT_BC       ] = 10;
    index_fail_to_dU[FAIL_VCHAR_NEG       ] = -1;
    index_fail_to_dU[FAIL_VCHAR_SUPER     ] = -1;
    index_fail_to_dU[FAIL_VCHAR_DISCR     ] = -1;
    index_fail_to_dU[FAIL_BASIC           ] = -1;
    index_fail_to_dU[FAIL_METRIC          ] = -1;
    index_fail_to_dU[FAIL_RESTART         ] = -1;
    index_fail_to_dU[FAIL_MPI_BASIC       ] = -1;
    index_fail_to_dU[FAIL_HDF             ] = -1;

    local_first_time = 0; 
  }

  if( (fail_type < 0)  || (fail_type >= N_FAIL) ) { 
    fprintf(stderr,"calc_dU_stat(): invalid value of fail_type :  %d \n",fail_type); 
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  g = index_fail_to_dU[fail_type];

  /* Leave now if we do not recalculate U[] after this type of "failure" */
  if( g < 0 )   return;  

  /* We are only concerned with the conserved variables at the cell centers : */
  if( pcurr != CENT ) return;

  /* Only count the change when it occurs on the physical grid : */
  if((icurr > N1E)||(jcurr > N2E)||(kcurr > N3E)||(icurr < N1S)||(jcurr < N2S)||(kcurr < N3S))  return;


  /* Now get down to calculating the change in the conserved variables : */ 
  ind_U = (kcurr-N3S) + N3*((jcurr-N2S) + N2*(icurr-N1S));
  get_geometry(icurr,jcurr,kcurr,pcurr,ncurr,geom);  /* Assuming axisymmetric spacetime */
  get_state(  p[icurr][jcurr][kcurr], geom, &q );
  primtoflux( p[icurr][jcurr][kcurr], &q, 0, geom, U_new) ;

  NPH_LOOP { dU_stat[l][g][ind_U] += U_new[l] - U_out[l][ind_U] ; }

  NPH_LOOP { U_out[l][ind_U] = U_new[l]  ; }

  //  fprintf(stdout,"calc_dU_stat end pid=%d   nstep=%d\n",myid,nstep); fflush(stdout);

  return;

}

///**********************************************************************************************/
///**********************************************************************************************
//  calc_U_of_P():
// ------------------
//    -- given a primitive variable array pv[], calculates 
//        U(P) for a given U_of_P[] element; 
//
//**********************************************************************************************/
//void calc_U_of_P(int i_elem, double (* pv)[N2TOT][N3TOT][NP])
//{
//  int i,j,k,l,p,g,ind;
//  double U[NP];
//  struct of_state q;
//  struct of_geom *geom; 
//
//  ind = 0; 
//  N1_LOOP N2_LOOP N3_LOOP { 
//    get_geometry(i,j,k,CENT,ncurr,geom) ; /* assume axisym. metric */
//      get_state(  pv[i][j][k], geom, &q );
//      primtoflux( pv[i][j][k], &q, 0, geom, U) ;
//      PLOOP  U_of_P[l][i_elem][ind] = U[l]; 
//      ind++;
//  }
//
//  return;
//}
//  

#undef IND_LOOP 

#else 
void dump_stat(void) { return; } 
#endif


