
/* restart functions; restart_init and restart_dump */

#include "decs.h"


#if( MAKE_HDF5 ) 
extern void dump_hdf5_gen(int do_restart);
extern int restart_check_hdf5( char *name ) ;
extern void restart_read_hdf5( char *name ) ;
#endif 

/***********************************************************************/
/***********************************************************************
  restart_write():
     -- writes current state of primitive variables to the 
        checkpointing/restart file. 
     -- uses ASCII text format ;
     -- when changing this routine, be sure to make analogous changes 
        in restart_read();
************************************************************************/
void restart_write(void)
{
  FILE *fp ;
  int idum,i,j,k,l ;
  char dfnam[50]; 

  TRACE_BEG;
  
#if( MAKE_HDF5 ) 
  dump_hdf5_gen(1); /* generate restart file in hdf5 format instead */

#else  

  if(rdump_cnt%2 == 0) {
    sprintf(dfnam,"rdump0.%04d",myid); 
  }
  else {
    sprintf(dfnam,"rdump1.%04d",myid); 
  }

  fp = fopen(dfnam,"wt") ;
  fprintf(stdout,"RESTART  file=%s\n",dfnam) ;
  if(fp == NULL) {
    fprintf(stderr,"Cannot open restart file on pid=%d \n",myid ) ;
    fail( FAIL_BASIC,0 ) ; 
  }

  /*************************************************************
	  Write the header of the restart file: 
  *************************************************************/
  fprintf(fp, FMT_DBL_OUT, t        );
  fprintf(fp, FMT_INT_OUT, N1       );
  fprintf(fp, FMT_INT_OUT, N2       );
  fprintf(fp, FMT_INT_OUT, N3       );
  fprintf(fp, FMT_DBL_OUT, GridLength[0]);
  fprintf(fp, FMT_INT_OUT, nstep    );
  fprintf(fp, FMT_DBL_OUT, gam      );
  fprintf(fp, FMT_DBL_OUT, cour     );
  OUT_LOOP  fprintf(fp, FMT_INT_OUT, N_out[i] ); 
  fprintf(fp, FMT_INT_OUT, rdump_cnt);
  fprintf(fp, FMT_DBL_OUT, dx[0]    );

  fprintf(fp,"\n");

  /*************************************************************
	  Write the body of the restart file: 
  *************************************************************/
  RDUMP_LOOP {
    PLOOP fprintf(fp, FMT_DBL_OUT, p[i][j][k][l]); 
    fprintf(fp, "\n"); 
  }

  fclose(fp) ;
#endif


  rdump_cnt++ ;
	
  TRACE_END;
  return;
}


/***********************************************************************/
/***********************************************************************
  restart_init():
     -- main driver for setting initial conditions from a checkpoint 
        or restart file. 
     -- determines if there are any restart files to use; if there are 
         two, then picks "rdump0";  can not prompt the user 
         for input when doing MPI, so we use one as our default rdump; 
     -- then initializes run with restart data;
     -- relies on init() to set coordinate system and metric functions
************************************************************************/
int restart_init(void)
{
  FILE *fp, *fp1, *fp0 ;
  char dfnam1[50],dfnam0[50],hdf_name[200]; 
  //int strncmp(char *s1, char *s2, int n) ;
  int i,j,k,l ;
  unsigned int ui;

  void restart_read(FILE *fp);

  TRACE_BEG;

#if( MAKE_HDF5 ) 
  if( !restart_check_hdf5(hdf_name) ) { 
    if( myid == printer_pid ) {
      fprintf(stdout,"restart(): Not using restart file!\n");
      fflush(stdout);
    }
    TRACE_END;
    return(0);
  }

  using_restart = 1 ;   /* Make sure all routines know that we are restarting */

  /* set up global arrays and parameters
     NOTE:  parameters for coordinates and geometry are all defined
            at compile time to avoid portability issues (e.g. some chips 
	    may calculate Rin or the metric functions differently)
   */
  init_base();

  restart_read_hdf5(hdf_name);

#else 
  /********************************************************************
   Check to see which restart files exist. 
   Use the only one that exists, else prompt user to decide 
     which one to use if we have a choice : 
  ********************************************************************/
  sprintf(dfnam0,"rdump0.%04d",myid); 
  sprintf(dfnam1,"rdump1.%04d",myid); 
  fp0 = fopen(dfnam0,"r") ;
  fp1 = fopen(dfnam1,"r") ;

  if( (fp0 == NULL) && (fp1 == NULL) ) {
    if( myid == printer_pid ) {  fprintf(stdout,"No restart file\n") ; fflush(stdout);  }
    TRACE_END;
    return(0) ;
  }
  else { 
    if( myid == printer_pid ) {  fprintf(stdout,"\nRestart file exists! \n") ; fflush(stdout); }
    if( fp0 == NULL ) { 
      if( myid == printer_pid ) { fprintf(stdout,"Using %s ... \n",dfnam1); fflush(stdout); }
      fp = fp1;
    }
    else if( fp1 == NULL ) { 
      if( myid == printer_pid ) { fprintf(stdout,"Using %s ... \n",dfnam0); fflush(stdout); }
      fp = fp0;
    }
    else { 
      fprintf(stderr,"restart_init(): Should not be here  ... \n");
      fail(FAIL_BASIC,0);
    }
  }

  using_restart = 1 ;   /* Make sure all routines know that we are restarting */

  /********************************************************************
   Now that we know we are restarting from a checkpoint file, then 
   we need to read in data, assign grid functions and define the grid: 
  ********************************************************************/
  init_base();
  restart_read(fp) ;
  fclose(fp) ;

#endif 

  if( myid == printer_pid ) {   fprintf(stdout,"restart_init():  using_restart =  %d\n", using_restart); fflush(stdout);  }


#if( DYNAMIC_SPACETIME )
  /* Need to recalculate metric as dx[0] has most likely changed after reading in checkpoint */
  init_general_metric();
#endif

#if( USE_SIMPLE_EOS )
  set_Katm();
#endif 

  /* bound */
  fixup(p);
  bounds(p,0) ;

  /* set half-step primitives everywhere */
  ALL_LOOP PLOOP { ph[i][j][k][l] = p[i][j][k][l] ;  }

  if( myid == printer_pid ) { fprintf(stdout,"done with restart init.\n") ; fflush(stdout);  }

  TRACE_END;
  return(1) ;

}

/***********************************************************************/
/***********************************************************************
  restart_read():
     -- reads in data from the restart file, which is specified in 
         restart_init() but is usually named "rdump[0,1]" 
************************************************************************/
void restart_read(FILE *fp)
{
  int idum,i,j,k,l ;

  TRACE_BEG;

  /*************************************************************
	  READ the header of the restart file: 
  *************************************************************/
  fscanf(fp, "%lf", &t        );

  fscanf(fp, "%d", &idum  );
  if(idum != N1) {
    fprintf(stderr,"error reading restart file; N1 differs\n") ;
    fail( FAIL_MPI_BASIC ,0) ; 
  }
  fscanf(fp, "%d", &idum  );
  if(idum != N2) {
    fprintf(stderr,"error reading restart file; N2 differs\n") ;
    fail( FAIL_MPI_BASIC ,0) ; 
  }
  fscanf(fp, "%d", &idum  );
  if(idum != N3) {
    fprintf(stderr,"error reading restart file; N3 differs\n") ;
    fail( FAIL_MPI_BASIC,0 ) ; 
  }

  fscanf(fp, "%lf", &(GridLength[0]));
  fscanf(fp, "%d" , &nstep    );
  fscanf(fp, "%lf", &gam      );
  fscanf(fp, "%lf", &cour     );
  OUT_LOOP  fscanf(fp, "%d", &(N_out[i]) ); 
  fscanf(fp, "%d" , &rdump_cnt);
  fscanf(fp, "%lf", &(dx[0])    );

  /*************************************************************
	  READ the body of the restart file: 
  *************************************************************/
  RDUMP_LOOP  PLOOP fscanf(fp, "%lf", &(p[i][j][k][l])); 

  TRACE_END;

  return ;
}

