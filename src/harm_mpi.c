
#include "decs.h"
#include "harm_mpi.h"

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************
   M P I     R O U T I N E S 
*********************************************************************************/
/********************************************************************************/
/********************************************************************************
   -- Physical points are those points we want to output.  They are usually 
      those cells that are not ghost cells.  Boundary cells are usually ghost 
      cells, but we try to differentiate them as much as we can since we 
      want these routines to be as general as possible.  So, physical cells 
      are non-boundary cells.  We do not want to output boundary cells 
      since they are shared amongst other grids and are therefore redundant. 
   -- The "global domain" is the collection of physical cells over all subdomains. 
      It includes no boundary cells and represents all the cells that are dumped.
   -- The term "local" refers to aspects belonging to a CPU-specific quantity 
      in a particular subdomain.  That is, each CPU has its own local domain. 
      The union of all physical cells from all local domains is the global 
      domain. 
********************************************************************************/
/********************************************************************************/

// arrays: (physical points are ones to be outputted) */
int     global_index;       /* Global index of domain's 1st physical cell if you were to make the global 3D array 1D */
int     bc_range_send[NDIM][2][NDIM][2];  /* Contains index ranges used when sending BCs */
int     bc_range_recv[NDIM][2][NDIM][2];  /* Contains index ranges used when receiving BCs */
int     bc_npts[NDIM][2];         /* Number of points to transfer (back/forth) at a give face */
static double *bc_data_send[NDIM][2];   /* pointers to data arrays used to send BC data */
static double *bc_data_recv[NDIM][2];   /* pointers to data arrays used to receive BC data */
static int    *bc_pflag_send[NDIM][2];  /* pointers to pflag mpi arrays used to send BC data */
static int    *bc_pflag_recv[NDIM][2];  /* pointers to pflag mpi arrays used to receive BC data */
int     mpi_chunk_size;          /* Number of points in global output */ 
double *mpi_comm_chunk;          /* pointer to MPI work array for global output */ 

//mpi-specific stuff:
#if( USEMPI ) 
int    procnamelen;
char myidtxt[MAXFILENAME], processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status mpi_status;
static MPI_Request request_out, request_in, req_bc_out[2*NDIM], req_bc_in[2*NDIM];
#endif

void read_chunk( int pid, int n_chunk, double *f_chunk, double *f_glob); 
void setup_mpibc_info( void );
void setup_mpibc_info_sym( void );
void  get_global_ijk( int i, int j, int k, int  *i_glob, int *j_glob, int *k_glob) ;
void  get_global_ijk_pid( int pid, int i, int j, int k, int *i_glob, int *j_glob, int *k_glob) ;
int get_global_index( int i, int j, int k );
int get_global_index_pid( int pid, int i, int j, int k );
void get_cpupos_pid( int pid, int tmpcpupos[NDIM] ) ;
static void myargs(int argc, char *argv[]);

#if( USEMPI ) 
int init_MPI(int argc, char *argv[]);
#endif 

static void set_pids(void);

/*********************************************************************************/
/*********************************************************************************
  setup_mpi():
 -------------
     -- collects a given grid function to the master node from all 
        the cells from every domain; 
     -- this means that the master node must have enough memory to allocate 
        a double array with as many elements as there are total points in the 
        simulation (n_cells_global) 
     -- the collection is done using a simple loop through the cpus;
*********************************************************************************/
void setup_mpi(int argc, char *argv[])
{
  int i ; 

#if( USEMPI ) 
  init_MPI(argc, argv);  
#else 
  for( i = 0; i < NDIM; i++ ) {  ncpux[i] = 1;  }
  numprocs = 1;
  myid = 0;
  set_pids();
#endif 

  myargs(argc,argv); // set constants from command line arguments (layout of global grid's segmentation)
  set_mpi_grid();    // these parameters are used for non-mpi operations as well
  set_mpi_misc();


  return;
}

/*********************************************************************************/
/*********************************************************************************
  get_global_gfunc(): 
  -----------------
     -- collects a given grid function to the master node from all 
        the cells from every domain; 
     -- this means that the master node must have enough memory to allocate 
        a double array with as many elements as there are total points in the 
        simulation (n_cells_global);  this allocation must take place before 
        this call;  the reason why this function does not allocate it, is because 
        this function will invariably be called frequently and so the lack of 
        an "alloc()" call eliminates many alloc's and dealloc's; 
     -- the collection is done using a simple loop through the cpus;
     -- the pointer is to a 1-dimensional array ordered as the 
         "p[][][][]" array and can be indexed: 
            f[ k + totalsize[3]*(j + i*totalsize[2]) ]

           i.e.  "k" runs fastest, then "j", then "i"

     -- note that the slaves use "mpi_comm_chunk()" as the sending array and 
        the master uses it as the receiving ; 
*********************************************************************************/
void get_global_gfunc( int ip, double *gfunc ) 
{
  int i, j, k ; 
  int tag, pid;

  TRACE_BEG;

#if(!USEMPI)
  write_chunk(ip, mpi_comm_chunk) ;
  read_chunk(myid, mpi_chunk_size, mpi_comm_chunk, gfunc); 

#else 

  /********************************************************************************
    Do master specific operations : 
  ********************************************************************************/
  if( myid == io_pid )  {
    /*  Read the master's data directly (sidestep MPI send/recv process):  */
    write_chunk(ip,  mpi_comm_chunk) ;
    read_chunk(myid, mpi_chunk_size, mpi_comm_chunk, gfunc); 
  }

  /********************************************************************************
    Get Loop over all slaves since we have already done the master:
  ********************************************************************************/
  for( pid = 1; pid < numprocs ; pid++ ) { 
    tag = pid + BASE_TAG_DUMP; 

    // Do the recieve and send requests first from the master:
    if( myid == io_pid ) { 
      //      fprintf(stdout,"get_global_gfunc(): proc = %d  \n", myid); fflush(stdout);
      MPI_Irecv( mpi_comm_chunk, 
		 mpi_chunk_size, 
		 MPI_DOUBLE, 
		 pid , 
		 tag, 
		 MPI_COMM_WORLD, 
		 &request_in );
    }

    // if it is this cpu's turn in the loop, then send its stuff: 
    if( pid == myid ) { 
      //      fprintf(stdout,"get_global_gfunc(): proc = %d \n", myid); fflush(stdout);
      write_chunk(ip, mpi_comm_chunk) ;
      MPI_Isend( mpi_comm_chunk, 
		 mpi_chunk_size, 
		 MPI_DOUBLE, 
		 io_pid ,  
		 tag, 
		 MPI_COMM_WORLD, 
		 &request_out );
    }

    // Then wait for them to complete, read in the data after we have received it :
    if( myid == io_pid ) {
      //      fprintf(stdout,"get_global_gfunc(): proc = %d \n", myid); fflush(stderr);
      MY_MPI_Wait( &request_in, &mpi_status);
      read_chunk(pid, mpi_chunk_size, mpi_comm_chunk, gfunc); 
    }
    if( pid == myid ) { 
      MY_MPI_Wait( &request_out, &mpi_status);
    }

  }

  return;

#endif 

  TRACE_END;
  
}


/*********************************************************************************/
/*********************************************************************************
  set_mpi_grid(): 
  ----------
     -- setup the domain decomposition grid parameter
     -- assigns each grid to a cpu 
*********************************************************************************/
void set_mpi_grid(void)
{
  int i;

  TRACE_BEG;

  set_global_domain_variables();

  fprintf(stdout," percpusize:  ");
  fprintf(stdout, " NM1=%d NM2=%d NM3=%d ", NM1,NM2,NM3); 

  /* Dimensions of global physical domain */
  totalsize[0] = 1;
  totalsize[1] = ncpux[1] * NM1; 
  totalsize[2] = ncpux[2] * NM2; 
  totalsize[3] = ncpux[3] * NM3; 

  n_cells_glob = 1;
  for(i = 0 ; i < NDIM ; i++ ) {     
    n_cells_glob *=  totalsize[i];   /* Num. of physical points in entire evolution over all CPUs */
  }

  /* Set indices (or order in each direction) of the CPU's domain within the domain decomposition */
  get_cpupos_pid( myid, cpupos ); 

  /* Set the x1,x2,x3 indices of local domain's 1st physical cell position in the global domain: */
  get_global_ijk_pid( myid, 0,0,0, &(globalpos[1]), &(globalpos[2]), &(globalpos[3]) ); 
  globalpos[0] = 0;

  /* set the number of the local domain's first cell in the global enumeration  
     (order of first cell in the global array)      */
  global_index = get_global_index_pid( myid, 0, 0, 0 );


  /* Set the number of independent spatial dimensions */ 
  n_spatial_dims = 0; 
  if( totalsize[1] > 1 )  n_spatial_dims++; 
  if( totalsize[2] > 1 )  n_spatial_dims++; 
  if( totalsize[3] > 1 )  n_spatial_dims++; 

  /* Dump domain information : */
  fprintf(stdout, "\n################# MPI DOMAIN INFO ################################\n");
  fprintf(stdout,"MPI-DOM: pid=%d : ncpux[0-3]     = %4d %4d %4d %4d \n", 
	  myid,ncpux[0],ncpux[1],ncpux[2],ncpux[3]);
  fprintf(stdout,"MPI-DOM: pid=%d : cpupos[0-3]    = %4d %4d %4d %4d \n", 
	  myid,cpupos[0],cpupos[1],cpupos[2],cpupos[3]);
  fprintf(stdout,"MPI-DOM: pid=%d : totalsize[0-3] = %4d %4d %4d %4d \n", 
	  myid,totalsize[0], totalsize[1], totalsize[2], totalsize[3]); 
  fprintf(stdout,"MPI-DOM: pid=%d : globalpos[0-3] = %4d %4d %4d %4d \n", 
	  myid,globalpos[0], globalpos[1], globalpos[2], globalpos[3]); 
  fprintf(stdout,"MPI-DOM: pid=%d : n_cells_glob   = %ld \n", myid,n_cells_glob ); 
  fprintf(stdout,"MPI-DOM: pid=%d : global_index   = %d \n", myid,global_index ); 
  fprintf(stdout,"MPI-DOM: pid=%d : numprocs       = %d \n", myid,numprocs ); 
  fprintf(stdout,"MPI-DOM: pid=%d : n_spatial_dims = %d \n", myid,n_spatial_dims ); 

  fprintf(stdout, "###################################################################\n");
  fflush(stdout);



  /* Set boundary condition information (e.g. neighbors, etc.) : */
  setup_mpibc_info();

  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  set_mpi_misc():
  ----------
    -- set misc. MPI variables that are indep. of the grid:
    -- even used for dump purposes; 
*********************************************************************************/
void set_mpi_misc(void)
{

  TRACE_BEG;

  // this is the work array for data transfer:
  mpi_chunk_size = NM1*NM2*NM3;  /* If you change this, then change LOOP bounds in [read,write]_chunk() */

  mpi_comm_chunk = (double *) calloc( mpi_chunk_size , sizeof(double) );

  if( mpi_comm_chunk == NULL ) { 
    fprintf(stderr,"set_mpi_misc(): Cannot allocate mpi_comm_chunk[%d] !!\n", 
	    mpi_chunk_size);
    fflush(stderr);
    fail( FAIL_MPI_BASIC,0 );
  }

  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  myexit():
  ----------
     -- calls final MPI routines before exiting simulation;
*********************************************************************************/
void myexit( int ret )
{
  int i,j,k;
  extern void free_global_arrays(void);

  TRACE_BEG;

  free_global_arrays();

#if(USEMPI)

  /* Free dynamically allocated arrays : */
  if( mpi_comm_chunk != NULL ) { free( mpi_comm_chunk ); }
  for( i = 1; i < NDIM; i++ )  for( j = 0; j < 2; j++ ) { 
    if( bc_data_send[i][j] != NULL ) {  free( bc_data_send[i][j] ); }
    if( bc_data_recv[i][j] != NULL ) {  free( bc_data_recv[i][j] ); }
    if( bc_pflag_send[i][j] != NULL ) { free( bc_pflag_send[i][j] ); }
    if( bc_pflag_recv[i][j] != NULL ) { free( bc_pflag_recv[i][j] ); }
  }

  /* Finish up MPI  */

  /* Do a full abort if there is an error -- it could be that there is an MPI problem on only one processor: */
  if( ret ) {   MPI_Abort(MPI_COMM_WORLD,ret);  }

  MPI_Barrier(MPI_COMM_WORLD);        /*  Required!  */
  MPI_Finalize();
#endif


  TRACE_END;

  exit( ret ) ; 

  return;
}


/*********************************************************************************/
/*********************************************************************************
  exit_status():
  ----------
     -- makes sure that all processes are running without failures;
     -- processes will reach here in one of three ways: 
         1) at the end of the full timestep procedure (see call in main())
         2) via some fatal error during the timestep procedure (see call in diag())
         3) before a mpi communication to ensure that all siblings are on track (see calls in harm_mpi.c)

     -- consequently, all the processes that fail in a timestep will be able to 
        report that it has failed;  
*********************************************************************************/
void exit_status( void ) 
{
  short int fail_max, sendbuf; 
  static short int exiting = 0 ; 

#if( FAST_AND_FURIOUS ) 
  return;
#endif

  TRACE_BEG;

#if( USEMPI ) 
  /* Make sure we let mpi-related exiting procedures to continue  */
  if( exiting ) { return; } 

  fail_max = failed ;

  if( failed ) { 
    fprintf(stdout,"exit_status: first-fail pid=%d   failed = %d   \n",myid, failed);
    fflush(stdout);
  }

  sendbuf = fail_max;
  MPI_Allreduce(&sendbuf, &fail_max, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD);

  if( failed || fail_max ) { 
    exiting = 1; 
  }

  /* Only call fail() if we have not done so already */ 
  /* Note that failed is only set when a fatal error occurs, so fail_max should always 
      be a fatal error state.  */
  if( fail_max && (!failed) ) { 
    sibling_failed = 1; 
    fprintf(stdout,"exit_status: sibling-failed  pid=%d   failed = %d    fail_max = %d \n",myid, failed,fail_max); 
    fflush(stdout);
    fail(fail_max,0);
  }
#endif 

  TRACE_END;

  return;
}


/*********************************************************************************/
/*********************************************************************************
  get_pid():
  ----------
     --  returns the processor's id number (order in group) that is located
         at a given cpupos[]
     -- returns with BAD_PID if the position is not valid ; 
*********************************************************************************/
int get_pid( int tmpcpupos[NDIM] ) 
{ 
  int i ; 

  TRACE_BEG;

  /* First, check to see if this is a valid position: */
  SDLOOP1 { 
    if( (tmpcpupos[i] < 0) || (tmpcpupos[i] >= ncpux[i]) ) { 
      return( BAD_PID ); 
    }
  }

  TRACE_END;

  return( tmpcpupos[1]  +  ncpux[1] * (tmpcpupos[2] + tmpcpupos[3]*ncpux[2]) );
}

/*********************************************************************************/
/*********************************************************************************
  setup_mpibc_info():
  ------------------
     -- used with and without MPI;
     -- setup data structures that hold neighbor information so that 
        we know which CPUs to send/recv boundary information to/from;
     -- if neighbor equals "BAD_PID", then we consider it a physical boundary 
        and let the unigrid part handle the BC;
     -- bc_pid[DIM][DNUP]  (DN=0, UP=1, DIM=0,1,2,3)
     -- this routine assumes only block ordering and does not take into 
        account an symmetries that may align edges of domains not normally aligned 
        in a cartesian layout; (e.g. phi-symmetry still needs inter-domain comm.);
*********************************************************************************/
void setup_mpibc_info( void )
{
  int i,j,k,l, npts, nbr_pos[NDIM]; 

  TRACE_BEG;

  /* Start assuming that all BCs are physical: */
  DLOOP1 { bc_pid[i][BCDN] = bc_pid[i][BCUP] = BC_PHYS; } 

#if( USEMPI )

  /* Use pid position on global domain to determine neighbors */  
  DLOOP1 nbr_pos[i] = cpupos[i]; 

  /* DN boundaries: */
  SDLOOP1 { 
    nbr_pos[i]--;   
    bc_pid[i][BCDN] = get_pid(nbr_pos);  
    nbr_pos[i]++; 
  }

  /* UP boundaries: */
  SDLOOP1 { 
    nbr_pos[i]++;   
    bc_pid[i][BCUP] = get_pid(nbr_pos);  
    nbr_pos[i]--; 
  }

  /*************************************************************************************
     Set MPI boundary information specific to inherent topologies of the grid 
     (e.g. symmetries, phi-symmetric boundary conditions)
  *************************************************************************************/
  setup_mpibc_info_sym();


  /*************************************************************************************
     Setup index ranges than indicate which cells along each boundary are transferred 
         bc_range_send[BC_DIM][BC_UP/DN][INDEX_DIM][INDEX_BEG/END]
     In setting bc_range_send/recv[], assume that all boundaries are MPI boundaries. 
  *************************************************************************************/
  for(i=0;i<NDIM;i++) for(j=0;j<2;j++) for(k=0;k<NDIM;k++) for(l=0;l<2;l++) {
    bc_range_send[i][j][k][l] = bc_range_recv[i][j][k][l] = 0;
  }

  /* x3-face: */
  for(j=0;j<2;j++) {
    bc_range_send[3][j][1][IBEG] = bc_range_recv[3][j][1][IBEG] = 0;
    bc_range_send[3][j][2][IBEG] = bc_range_recv[3][j][2][IBEG] = 0;
    bc_range_send[3][j][1][IEND] = bc_range_recv[3][j][1][IEND] = NM1_TOT - 1;
    bc_range_send[3][j][2][IEND] = bc_range_recv[3][j][2][IEND] = NM2_TOT - 1;
  }
  bc_range_send[3][BCDN][3][IBEG] = NM3S               ;     
  bc_range_send[3][BCDN][3][IEND] = NM3S + NM3_BND - 1 ;
  bc_range_send[3][BCUP][3][IBEG] = NM3E - NM3_BND + 1 ;     
  bc_range_send[3][BCUP][3][IEND] = NM3E               ;

  bc_range_recv[3][BCDN][3][IBEG] = bc_range_send[3][BCDN][3][IBEG] - NM3_BND  ;     
  bc_range_recv[3][BCDN][3][IEND] = bc_range_send[3][BCDN][3][IEND] - NM3_BND  ;     
  bc_range_recv[3][BCUP][3][IBEG] = bc_range_send[3][BCUP][3][IBEG] + NM3_BND ;
  bc_range_recv[3][BCUP][3][IEND] = bc_range_send[3][BCUP][3][IEND] + NM3_BND ;


  /* x2-face: */
  for(j=0;j<2;j++) {
    bc_range_send[2][j][1][IBEG] = bc_range_recv[2][j][1][IBEG] = 0;
    bc_range_send[2][j][3][IBEG] = bc_range_recv[2][j][3][IBEG] = 0;
    bc_range_send[2][j][1][IEND] = bc_range_recv[2][j][1][IEND] = NM1_TOT - 1;
    bc_range_send[2][j][3][IEND] = bc_range_recv[2][j][3][IEND] = NM3_TOT - 1;
  }
  bc_range_send[2][BCDN][2][IBEG] = NM2S               ;     
  bc_range_send[2][BCDN][2][IEND] = NM2S + NM2_BND - 1 ;
  bc_range_send[2][BCUP][2][IBEG] = NM2E - NM2_BND + 1 ;     
  bc_range_send[2][BCUP][2][IEND] = NM2E               ;

  bc_range_recv[2][BCDN][2][IBEG] = bc_range_send[2][BCDN][2][IBEG] - NM2_BND  ;     
  bc_range_recv[2][BCDN][2][IEND] = bc_range_send[2][BCDN][2][IEND] - NM2_BND  ;     
  bc_range_recv[2][BCUP][2][IBEG] = bc_range_send[2][BCUP][2][IBEG] + NM2_BND ;
  bc_range_recv[2][BCUP][2][IEND] = bc_range_send[2][BCUP][2][IEND] + NM2_BND ;


  /* x1-face: */
  for(j=0;j<2;j++) {
    bc_range_send[1][j][2][IBEG] = bc_range_recv[1][j][2][IBEG] = 0;
    bc_range_send[1][j][3][IBEG] = bc_range_recv[1][j][3][IBEG] = 0;
    bc_range_send[1][j][2][IEND] = bc_range_recv[1][j][2][IEND] = NM2_TOT - 1;
    bc_range_send[1][j][3][IEND] = bc_range_recv[1][j][3][IEND] = NM3_TOT - 1;
  }
  bc_range_send[1][BCDN][1][IBEG] = NM1S               ;     
  bc_range_send[1][BCDN][1][IEND] = NM1S + NM1_BND - 1 ;
  bc_range_send[1][BCUP][1][IBEG] = NM1E - NM1_BND + 1 ;     
  bc_range_send[1][BCUP][1][IEND] = NM1E               ;

  bc_range_recv[1][BCDN][1][IBEG] = bc_range_send[1][BCDN][1][IBEG] - NM1_BND  ;     
  bc_range_recv[1][BCDN][1][IEND] = bc_range_send[1][BCDN][1][IEND] - NM1_BND  ;     
  bc_range_recv[1][BCUP][1][IBEG] = bc_range_send[1][BCUP][1][IBEG] + NM1_BND ;
  bc_range_recv[1][BCUP][1][IEND] = bc_range_send[1][BCUP][1][IEND] + NM1_BND ;


  /*************************************************************************************
    Using the above values of bc_range_send/recv[][][], allocate the send/recv 
    arrays for future use in bounds_mpi(): 
    -- Assumes that we are sending as many points as we are receiving (again, 
       assuming for these purposes that all boundaries are MPI boundaries)
          bc_npts[BC_DIM][I_UP/DN]
    -- we want to allocate an array for each send and receive since we want to 
       do this via nonblocking calls so we cannot re-use an array for each send/recv;
    -- this should not be a significant memory burden since it is only boundary data ;
    -- make bc array for pflag[] too (it needs only be 1/NP as large);
  *************************************************************************************/
  for( i = 0; i < NDIM; i++ )  for( j = 0; j < 2; j++ ) {  bc_npts[i][j] = 0; }
  
  for( i = 1; i < NDIM; i++ )  for( j = 0; j < 2; j++ ) { 
    bc_npts[i][j] = NP; 
    for( k = 1; k < NDIM; k++ ) { 
      bc_npts[i][j]  *=  bc_range_recv[i][j][k][IEND] - bc_range_recv[i][j][k][IBEG] + 1; 
    }

    bc_data_send[i][j] = (double *) calloc( bc_npts[i][j] , sizeof(double) );    
    if( bc_data_send[i][j] == NULL ) { 
      fprintf(stderr,
	      "setup_mpibc_info(): Cannot allocate bc_data_send[%d][%d] of length %d \n",
	      i,j,bc_npts[i][j]);
      fflush(stderr);
      fail( FAIL_MPI_BASIC,0 );
    }
    
    bc_data_recv[i][j] = (double *) calloc( bc_npts[i][j] , sizeof(double) );    
    if( bc_data_recv[i][j] == NULL ) { 
      fprintf(stderr,
	      "setup_mpibc_info(): Cannot allocate bc_data_recv[%d][%d] of length %d \n",
	      i,j,bc_npts[i][j]);
      fflush(stderr);
      fail( FAIL_MPI_BASIC,0 );
    }

    /* pflag allocations : */
    bc_pflag_send[i][j] = (int *) calloc( (bc_npts[i][j]/NP) , sizeof(int) );    
    if( bc_pflag_send[i][j] == NULL ) { 
      fprintf(stderr,
	      "setup_mpibc_info(): Cannot allocate bc_pflag_send[%d][%d] of length %d \n",
	      i,j,(bc_npts[i][j]/NP));
      fflush(stderr);
      fail( FAIL_MPI_BASIC,0 );
    }
    
    bc_pflag_recv[i][j] = (int *) calloc( (bc_npts[i][j]/NP) , sizeof(int) );    
    if( bc_pflag_recv[i][j] == NULL ) { 
      fprintf(stderr,
	      "setup_mpibc_info(): Cannot allocate bc_pflag_recv[%d][%d] of length %d \n",
	      i,j,(bc_npts[i][j]/NP));
      fflush(stderr);
      fail( FAIL_MPI_BASIC,0 );
    }
    
  }

#endif 


  /* Dump bc type information : */ 
  fprintf(stdout,"\n################# MPI BC INFO ######################################\n");
  DLOOP1 { 
    fprintf(stdout,"MPI-BC: pid=%d : DN :  dir=%d  :  bctype = %d \n", 
	    myid, i, bc_pid[i][BCDN]);
    fprintf(stdout,"MPI-BC: pid=%d : UP :  dir=%d  :  bctype = %d \n", 
	    myid, i, bc_pid[i][BCUP]);
  }
  fprintf(stdout,"###################################################################\n");
  fflush(stdout);
	 

  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  get_global_ijk():
  ----------------
     -- return with a cell's GLOBAL indices that correspond given 
        the LOCAL ones:
*********************************************************************************/
void  get_global_ijk( int i, int j, int k, int  *i_glob, int *j_glob, int *k_glob) 
{

  //  TRACE_BEG;

  *i_glob = i + globalpos[1] ;
  *j_glob = j + globalpos[2] ;
  *k_glob = k + globalpos[3] ;

  //  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  get_global_ijk_pid():
  ----------------
     -- return with a cell's GLOBAL indices that correspond given 
        the LOCAL ones;  the cell belong to the domain of processor "pid"
*********************************************************************************/
void  get_global_ijk_pid( int pid, int i, int j, int k, int *i_glob, int *j_glob, int *k_glob) 
{
  int tmppos[NDIM];

  TRACE_BEG;

  get_cpupos_pid( pid, tmppos );
  
  *i_glob = i + tmppos[1] * NM1; 
  *j_glob = j + tmppos[2] * NM2; 
  *k_glob = k + tmppos[3] * NM3; 
 
  TRACE_END;

  return;
}


/*********************************************************************************/
/*********************************************************************************
  get_global_index():
  ----------------
     -- return with the GLOBAL total index (enumerated order 
        in the global array, not including ghosts) of a cell from the local 
        processor that correspond given the LOCAL indices;
     -- we set it so we loop fastest in k (x3), then j (x2), then i (x1)
*********************************************************************************/
int get_global_index( int i, int j, int k )
{
  int ig,jg,kg;

  TRACE_BEG;
  
  get_global_ijk(i,j,k,&ig,&jg,&kg);

  TRACE_END;

  return( kg + totalsize[3] * ( jg  +  totalsize[2] * ig )  );
}


/*********************************************************************************/
/*********************************************************************************
  get_global_index_pid():
  ----------------
     -- return with the GLOBAL total index (enumerated order 
        in the global array, not including ghosts) of a cell from processor "pid"
        that correspond given its LOCAL indices;
     -- we set it so we loop fastest in k (x3), then j (x2), then i (x1)
*********************************************************************************/
int get_global_index_pid( int pid, int i, int j, int k )
{
  int ig,jg,kg;

  TRACE_BEG;
  
  get_global_ijk_pid(pid,i,j,k,&ig,&jg,&kg);

  TRACE_END;

  return( kg + totalsize[3] * ( jg  +  totalsize[2] * ig )  );
}

/*********************************************************************************/
/*********************************************************************************
  get_cpupos_pid():
  ----------
     -- Given the processor id number, 
        return the MPI position of the point's subgrid in the global domain 
     -- subdomains are distributed fastest in x1, and slowest in x3 since x1 is 
        always used 
*********************************************************************************/
void get_cpupos_pid( int pid, int tmpcpupos[NDIM] ) 
{ 
  TRACE_BEG;
  tmpcpupos[0] = 0;
  tmpcpupos[1] = ( pid                       ) % ncpux[1] ;
  tmpcpupos[2] = ( pid /  ncpux[1]           ) % ncpux[2] ; 
  tmpcpupos[3] = ( pid / (ncpux[1]*ncpux[2]) )            ; 
  TRACE_END;
  return;
}


/*********************************************************************************/
/*********************************************************************************
  write_chunk(): 
  ----------
       -- copy elements from whereever (based on ivar) to the mpi array for sending;
       -- we assume that NG=N_BND (w/ our use of "LOOP" here)
       -- if we want to output that is not always stored, then the user needs to 
          specific it here; 
*********************************************************************************/
#define IJCON  (NP+4)
//#define IJCON2 (NP+4+4)
#define IFAIL  (NP+4+4)

void write_chunk( int ivar, double *f_chunk )
{
  int i,j,k,l, ifail;
  double gamma; 
  struct of_geom  *geom;

  TRACE_BEG;

  /* Assuming here that (e.g.) N1E-N1S+1 = NM1 */
  l = 0;

  switch( ivar ) { 

    /* Primitive variables : */
  case RHO :  case UU  :  case U1  :  case U2  :  case U3  :  case B1  :  case B2  :  case B3  : 
    LOOP { f_chunk[l++] = p[i][j][k][ivar] ; } 
    break;


    /* b^2  : */
  case NP : 
    LOOP {
      get_geometry(i,j,k,CENT,ncurr,geom);
      f_chunk[l++] = bsq_calc(p[i][j][k], geom);
    }
    break;


    /* Lorentz gamma factor : */    
  case (NP+1) : 
    LOOP { 
      get_geometry(i,j,k,CENT,ncurr,geom);
      gamma_calc(p[i][j][k], geom, &gamma);
      f_chunk[l++] = gamma ; 
    }
    break;


    /* divergence of B^i using Flux-CT stencil (should be 0 to roundoff w/ our method) */
  case (NP+2) : 
    switch( n_spatial_dims ) { 
    case 1 : 
      LOOP { 
	if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS)) {
	  f_chunk[l++] = 0. ; 
	}
	else { f_chunk[l++] = divb_calc(i,j,k); }
      }
      break;
    case 2 : 
      if( totalsize[3] == 1 ) { 
	LOOP { 
	  if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	      (j==N2S && bc_pid[2][BCDN]==BC_PHYS)   ) { 
	    f_chunk[l++] = 0. ; 
	  }
	  else { f_chunk[l++] = divb_calc(i,j,k); }
	}
      }
      else if( totalsize[2] == 1 ) { 
	LOOP { 
	  if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	      (k==N3S && bc_pid[3][BCDN]==BC_PHYS)   ) { 
	    f_chunk[l++] = 0. ; 
	  }
	  else { f_chunk[l++] = divb_calc(i,j,k); }
	}
      }
      else { 
	LOOP { 
	  if( (k==N3S && bc_pid[3][BCDN]==BC_PHYS) || 
	      (j==N2S && bc_pid[2][BCDN]==BC_PHYS)   ) { 
	    f_chunk[l++] = 0. ; 
	  }
	  else { f_chunk[l++] = divb_calc(i,j,k); }
	}
      }
      break;
    case 3 :
      LOOP { 
	if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	    (j==N2S && bc_pid[2][BCDN]==BC_PHYS) || 
	    (k==N3S && bc_pid[3][BCDN]==BC_PHYS) ) { 
	  f_chunk[l++] = 0. ; 
	}
	else { f_chunk[l++] = divb_calc(i,j,k); }
      }
      break;
    }
    break; 
    

    /* divergence of B^i using cell centered stencil (may not be 0 to roundoff w/ our method) */
  case (NP+3) : 
    switch( n_spatial_dims ) { 
    case 1 : 
      LOOP { 
	if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	    (i==N1E && bc_pid[1][BCUP]==BC_PHYS)  ) {
	  f_chunk[l++] = 0. ; 
	}
	else { f_chunk[l++] = divb_cen_calc(i,j,k); }
      }
      break;
    case 2 : 
      if( totalsize[3] == 1 ) { 
	LOOP { 
	  if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	      (i==N1E && bc_pid[1][BCUP]==BC_PHYS) || 
	      (j==N2S && bc_pid[2][BCDN]==BC_PHYS) || 
	      (j==N2E && bc_pid[2][BCUP]==BC_PHYS) ) { 
	    f_chunk[l++] = 0. ; 
	  }
	  else { f_chunk[l++] = divb_cen_calc(i,j,k); }
	}
      }
      else if( totalsize[2] == 1 ) { 
	LOOP { 
	  if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	      (i==N1E && bc_pid[1][BCUP]==BC_PHYS) || 
	      (k==N3S && bc_pid[3][BCDN]==BC_PHYS) || 
	      (k==N3E && bc_pid[3][BCUP]==BC_PHYS) ) { 
	    f_chunk[l++] = 0. ; 
	  }
	  else { f_chunk[l++] = divb_cen_calc(i,j,k); }
	}
      }
      else { 
	LOOP { 
	  if( (k==N3S && bc_pid[3][BCDN]==BC_PHYS) || 
	      (k==N3E && bc_pid[3][BCUP]==BC_PHYS) || 
	      (j==N2S && bc_pid[2][BCDN]==BC_PHYS) || 
	      (j==N2E && bc_pid[2][BCUP]==BC_PHYS) ) { 
	    f_chunk[l++] = 0. ; 
	  }
	  else { f_chunk[l++] = divb_cen_calc(i,j,k); }
	}
      }
      break;
    case 3 :
      LOOP { 
	if( (i==N1S && bc_pid[1][BCDN]==BC_PHYS) || 
	    (i==N1E && bc_pid[1][BCUP]==BC_PHYS) || 
	    (j==N2S && bc_pid[2][BCDN]==BC_PHYS) || 
	    (j==N2E && bc_pid[2][BCUP]==BC_PHYS) || 
	    (k==N3S && bc_pid[3][BCDN]==BC_PHYS) || 
	    (k==N3E && bc_pid[3][BCUP]==BC_PHYS) ) { 
	  f_chunk[l++] = 0. ; 
	}
	else { f_chunk[l++] = divb_cen_calc(i,j,k); }
      }
      break;
    }
    break;


    /* jcon1 dumps : */
  case (IJCON) : case (IJCON+1) :  case (IJCON+2) :  case (IJCON+3) : 
    ifail  = ivar - IJCON;
#if( CALC_CURRENT ) 
    LOOP {
      f_chunk[l++] = jcon[i][j][k][ifail];
    }
#endif 
    break;


    /* jcon2 dumps : */
//  case (IJCON2) : case (IJCON2+1) :  case (IJCON2+2) :  case (IJCON2+3) : 
//    ifail  = ivar - IJCON2;
//    LOOP {
//      f_chunk[l++] = jcon2[i][j][k][ifail];
//    }
//    break;
//

    /* nfail dumps : */
  case (IFAIL) : case (IFAIL+1) :  case (IFAIL+2) :  case (IFAIL+3) :  case (IFAIL+4) : 
    ifail  = ivar - IFAIL;
    LOOP {
      f_chunk[l++] = (double) nfail[ifail][i][j][k];
    }
    LOOP {
      nfail[ifail][i][j][k] = 0;
    }
    break;


  default : 
    fprintf(stderr,"write_chunk(): Function %d is not implemented !! \n", ivar ) ; 
    fflush(stderr);
    fail( FAIL_MPI_BASIC,0 ) ; 
    break;

  }

  TRACE_END;
    
  return;
}

#undef IJCON  
#undef IJCON2 
#undef IFAIL  

/*********************************************************************************/
/*********************************************************************************
  read_chunk():
  ----------
     -- copy elements from the mpi array to the output array:
     -- responsible for mapping the local indexing to the global indexing;
     -- the "global array" and chunk contain only physical cells (i.e. no ghosts);
     -- assumes that we are looping fastest over x3, then x2, then x1 (k,j,i)
     -- assumes there are no boundary/ghost cells in the global domain;
*********************************************************************************/
void read_chunk( int pid,  int n_chunk, double *f_chunk, double *f_glob)
{
  int n, i_glob;
  int i,j,k, n_j,n_i;

  TRACE_BEG;

  if( n_chunk != NM1*NM2*NM3 ) { 
    fprintf(stderr,"read_chunk(): Bad value of n_chunk :  %d \n", n_chunk) ; 
    fflush(stderr);
    fail( FAIL_MPI_BASIC,0 ) ; 
  }
  
  n_i = totalsize[3] * (totalsize[2] - NM2 + 1) - NM3;
  n_j = totalsize[3] - NM3;
  
  /* Start off at this domain's starting point in the global array :  */
  i_glob = get_global_index_pid( pid, 0,0,0);

  /*******************************************************************************
    Temporarily use "fortran" index notation to eliminate extra additions (i.e. 
     compare to NM3 instead of (NM3-1);   this should work for all combinations 
     of NM1,NM2,NM3 (i.e. 1D, 2D,3D problems). 
  ********************************************************************************/
  i = j = k = 1 ; 
  for( n = 0 ; n < n_chunk ; n++ ) { 
    f_glob[i_glob++] = f_chunk[n]; 
    if( k >= NM3 ) { 
      if( j >= NM2 ) { 	j = 1;  i_glob += n_i;  }
      else { 	        j++  ;  i_glob += n_j;  }
      k = 1;
    }
    else { k++; }
  }

  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_min():
  ----------
     -- finds the minimum value among all domains of a scalar quantity;
     -- routine is called with a pointer to the domain's minimum value; 
     -- on exit, pointer is to the global minimum value;
*********************************************************************************/
void mpi_global_min(double *minval)
{
  double sendbuf;
  
  TRACE_BEG;

#if(USEMPI)
  sendbuf = *minval;
  exit_status();
  MPI_Allreduce(&sendbuf, minval, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_max():
  ----------
     -- finds the maximum value among all domains of a scalar quantity;
     -- routine is called with a pointer to the domain's maximum value; 
     -- on exit, pointer is to the global maximum value;
*********************************************************************************/
void mpi_global_max(double *maxval)
{
  double sendbuf;

  TRACE_BEG;
  
#if(USEMPI)
  sendbuf = *maxval;
  exit_status();
  MPI_Allreduce(&sendbuf, maxval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_vmax():
  ----------
     -- finds the maximum value among all domains of a vector quantity;
     -- routine is called with a pointer to the domain's maximum value; 
     -- on exit, pointer is to the global maximum values;
*********************************************************************************/
void mpi_global_vmax(int n, double *maxval)
{
  int i;
  double *sendbuf;

  TRACE_BEG;
  
#if(USEMPI)
  sendbuf = (double *) calloc( n, sizeof(double) );    
  if( sendbuf == NULL ) { 
    fprintf(stderr,"mpi_global_vamx(): Cannot allocate sendbuf !! \n"); 
    fflush(stderr);
    fail(FAIL_MPI_BASIC,0);
  }
  exit_status();
  for(i=0;i<n;i++) sendbuf[i] = maxval[i];
  MPI_Allreduce(sendbuf, maxval, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  free(sendbuf);
#endif
  
  TRACE_END;
  return;
 
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_imax():
  ----------
     -- finds the maximum value among all domains of a scalar quantity;
     -- routine is called with a pointer to the domain's maximum value; 
     -- on exit, pointer is to the global maximum value;
*********************************************************************************/
void mpi_global_imax(int *maxval)
{
  int sendbuf;
  
  TRACE_BEG;
#if(USEMPI)
  sendbuf = *maxval;
  exit_status();
  MPI_Allreduce(&sendbuf, maxval, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_imax2():
  ----------
     -- finds the maximum value among all domains of of two integers ;
     -- routine is called with a pointer to the domain's maximum value; 
     -- on exit, pointer is to the global maximum value;
*********************************************************************************/
void mpi_global_imax2(int *maxval1,int *maxval2)
{
  int sendbuf[2],writebuf[2];

  TRACE_BEG;
  
#if(USEMPI)
  sendbuf[0] = *maxval1;
  sendbuf[1] = *maxval2;
  exit_status();
  MPI_Allreduce(sendbuf, writebuf, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  *maxval1 = writebuf[0];
  *maxval2 = writebuf[1];
#endif
  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_sum():
  ----------
     -- sums the scalars on all domains and returns with the sum ; 
     -- on exit, pointer is to the global sum of the value ;
*********************************************************************************/
void mpi_global_sum(double *sumval)
{
  double sendbuf;
  TRACE_BEG;
  
#if(USEMPI)
  sendbuf = *sumval;
  exit_status();
  MPI_Allreduce(&sendbuf, sumval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_sum_to_head():
  ----------
     -- sums a vector (per element) amongst the processors and returns 
        the result to the head (aka master) node
     -- on exit, pointer is to the global sum of the value ;
*********************************************************************************/
void mpi_sum_to_head(double *v_send, double *v_recv, int n)
{
  int i; 
  
  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Reduce(v_send, v_recv, n, MPI_DOUBLE, MPI_SUM, io_pid, MPI_COMM_WORLD);
#else 
  for(i=0;i<n;i++) { v_recv[i] = v_send[i]; } 
#endif
  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_ulong_sum_to_head():
  ----------
     -- sums a vector (per element) amongst the processors and returns 
        the result to the head (aka master) node
     -- on exit, pointer is to the global sum of the value ;
*********************************************************************************/
void mpi_ulong_sum_to_head(unsigned long int *v_send, unsigned long int *v_recv, int n)
{
  int i; 
  TRACE_BEG;
  
#if(USEMPI)
  exit_status();
  MPI_Reduce(v_send, v_recv, n, MPI_UNSIGNED_LONG, MPI_SUM, io_pid, MPI_COMM_WORLD);
#else 
  for(i=0;i<n;i++) { v_recv[i] = v_send[i]; } 
#endif
  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  mpi_global_avg():
  ----------
     -- computes the average of scalars on all domains and returns with the average ; 
     -- on exit, pointer is to the global average of the value ;
*********************************************************************************/
void mpi_global_avg(double *avgval)
{
  double sendbuf;
  TRACE_BEG;
  
  mpi_global_sum(avgval); 
  *avgval /= numprocs; 
  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  sync_val():
  ----------
     -- make sure cpu has same value as master node; 
*********************************************************************************/
void sync_val(double *val)
{
 
  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Bcast(val, 1, MPI_DOUBLE, io_pid, MPI_COMM_WORLD);
#endif
  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  sync_val_from_rank():
  ----------
     -- make sure cpu has same value as a specified rank;
*********************************************************************************/
void sync_val_from_rank(double *val, int rank)
{
 
  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Bcast(val, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
#endif
  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  sync_int_from_rank():
  ----------
     -- make sure cpu has same value as master node; 
*********************************************************************************/
void sync_int_from_rank(int *val, int rank)
{
 
  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Bcast(val, 1, MPI_INT, rank, MPI_COMM_WORLD);
#endif
  TRACE_END;

  return;
}

/*********************************************************************************/
/*********************************************************************************
  sync_vect():
  ----------
     -- make sure cpu has same value as master node; 
*********************************************************************************/
void sync_vect(double *val, int n )
{
 
  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Bcast(val, n, MPI_DOUBLE, io_pid, MPI_COMM_WORLD);
#endif
  TRACE_END;
  
  return;
}

/*********************************************************************************/
/*********************************************************************************
  sync_vect_from_rank():
  ----------
     -- make sure cpu has same value as master node; 
*********************************************************************************/
void sync_vect_from_rank(double *val, int n, int rank )
{
 
  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Bcast(val, n, MPI_DOUBLE, rank, MPI_COMM_WORLD);
#endif
  TRACE_END;
  
  return;
}

/*********************************************************************************/
/*********************************************************************************
  sync_int_vect_from_rank():
  ----------
     -- make sure cpu has same value as master node; 
*********************************************************************************/
void sync_int_vect_from_rank(int *val, int n, int rank )
{

  TRACE_BEG;
#if(USEMPI)
  exit_status();
  MPI_Bcast(val, n, MPI_INT, rank, MPI_COMM_WORLD);
#endif
  TRACE_END;
  
  return;
}

/*********************************************************************************/
/*********************************************************************************
  collect_mask_values():
  ----------
     -- returns with unchanged array if not an MPI run; 
     -- if an MPI run, then returns with an array containing the values "vals[i]"
        for which mask[i] was set to be non-zero;
     -- collects it to rank = "rank"; 
*********************************************************************************/
void collect_mask_values(double *val, int *mask, int n, int rank )
{
 
  TRACE_BEG;
#if(USEMPI)
  int i ; 
  of_double_int *sendv;
  of_double_int *recv;
  MPI_Op myOp;
  void reduce_mask_operator( of_double_int *in, of_double_int *inout, int *len, MPI_Datatype *dtype) ;

  /* First load up structure array: */
  ALLOC_ARRAY(sendv,n);
  //  if( myid == rank ) { ALLOC_ARRAY( recv,n); }
  ALLOC_ARRAY( recv,n);

  for(i = 0 ; i < n ; i++) { 
    sendv[i].val  = val[i];
    sendv[i].mask = mask[i];
  }
  
  /* Define MPI operation */
  exit_status();
  /* MPI_Op_create((MPI_User_function *) reduce_mask_operator, 0, &myOp); */
  MPI_Op_create((MPI_User_function *) reduce_mask_operator, 1, &myOp);
  MPI_Reduce(sendv, recv, n, MPI_DOUBLE_INT, myOp, rank, MPI_COMM_WORLD);

  if( myid == rank ) { 
    for(i = 0 ; i < n ; i++) { 
       val[i] = recv[i].val ;
      mask[i] = recv[i].mask;
    }
    //    DEALLOC_ARRAY( recv,n); 
  }

  DEALLOC_ARRAY( recv,n); 
  DEALLOC_ARRAY(sendv,n);
  MPI_Op_free(&myOp); 

#endif
  TRACE_END;
  
  return;
}

/****************************************************************************************/
/*****************************************************************************************
  check_boundary_pflag()
  --------------
    -- responsible for setting non-zero status for variables "boundary_phys_pflag"
       and "boundary_mpi_pflag"  when this routine is called with indices for a cell 
        that is within a physical or mpi boundary communication;
    -- note that we do not expect there to be many instances where boundary_[mpi,phys]_pflag
       will be set, so we should optimize assuming that most cases will not 
*****************************************************************************************/
void check_boundary_pflag( int i, int j, int k ) 
{

  TRACE_BEG;

  if( boundary_mpi_pflag && boundary_phys_pflag ) { 
    TRACE_END;
    return; 
  }

  if( i < NM1S+NM1_BND ) { 
    if( bc_pid[1][BCDN] == BC_PHYS ) { boundary_phys_pflag = 1; }
    else                             { boundary_mpi_pflag  = 1; }
  }
  else if( i > NM1E-NM1_BND ) { 
    if( bc_pid[1][BCUP] == BC_PHYS ) { boundary_phys_pflag = 1; }
    else                             { boundary_mpi_pflag  = 1; }
  }

  if( j < NM2S+NM2_BND ) { 
    if( bc_pid[2][BCDN] == BC_PHYS ) { boundary_phys_pflag = 1; }
    else                             { boundary_mpi_pflag  = 1; }
  }
  else if( j > NM2E-NM2_BND ) { 
    if( bc_pid[2][BCUP] == BC_PHYS ) { boundary_phys_pflag = 1; }
    else                             { boundary_mpi_pflag  = 1; }
  }

  if( k < NM3S+NM3_BND ) { 
    if( bc_pid[3][BCDN] == BC_PHYS ) { boundary_phys_pflag = 1; }
    else                             { boundary_mpi_pflag  = 1; }
  }
  else if( k > NM3E-NM3_BND ) { 
    if( bc_pid[3][BCUP] == BC_PHYS ) { boundary_phys_pflag = 1; }
    else                             { boundary_mpi_pflag  = 1; }
  }

  TRACE_END;
  return;
}


/**********************************************************************************

  The following are MPI-only routines, not to be used with non-mpi runs:

**********************************************************************************/


/*********************************************************************************/
/*********************************************************************************
  myargs(): 
  ----------
     -- parse command-line arguments and assign global variables that 
        define the layout of the domain decomposition on the global grid
     -- set the number of CPUs used in each direction;
*********************************************************************************/
static void myargs(int argc, char *argv[])
{
  int i;

  TRACE_BEG;
  for( i = 0; i < NDIM; i++ ) { 
    ncpux[i] = 1;
  }

  if(argc<7){
    if(myid==master_pid){
      fprintf(stderr,"proc: %04d : Bad command line: argc: %d needed=%d, please specify:\n",
	      myid,(argc-1),3);
      fprintf(stderr,"proc: %04d : mpirun <mpirunoptions> <progname>  ncpux[1] ncpux[2] ncpux[3] N1_glob N2_glob N3_glob \n",
	      myid);
    }
    fail( FAIL_BASIC,0 ) ; 
  }

  /* Set the number of cpus used in each direction : */
  for( i = 1; i <= 3; i++ ) { 
    ncpux[i] = atoi(argv[i]);
  }
  N1_glob = atoi(argv[i++]);
  N2_glob = atoi(argv[i++]);
  N3_glob = atoi(argv[i++]);

#if(!USEMPI)
  for( i = 1; i <= 3; i++ ) { 
    if( ncpux[i] > 1 ) { 
      fprintf(stderr,"myargs(): For serial job, invalid number of processors:  ncpux[%3d] = %d \n",i,ncpux[i]); fflush(stderr);
      fail(FAIL_BASIC,0);
    }
  }
#endif

  fprintf(stdout, "myid=[%10d]:   ncpux[1-3]      = [ %6d , %6d , %6d ]\n",myid,ncpux[1],ncpux[2],ncpux[3]); 
  fprintf(stdout, "myid=[%10d]:       N[1-3]_glob = [ %6d , %6d , %6d ]\n",myid,N1_glob,N2_glob,N3_glob);
  fflush(stdout);

  TRACE_END;
  return;
}


#if(USEMPI)

/*********************************************************************************/
/*********************************************************************************
  init_MPI():
  ----------
     -- initialize MPI run; 
     -- set global MPI parametesr, CPU's identity, etc. 
*********************************************************************************/
// initialize stuff: 
int init_MPI(int argc, char *argv[])
{
  int i;

  TRACE_BEG;
  fprintf(stdout, "Begin: init_MPI\n");
  fflush(stdout);

  /* not expecting MPI command line arguments so this is a formality; */
  fprintf(stdout, "MPI_Init....\n"); fflush(stdout);
  MPI_Init(&argc, &argv);   

  fprintf(stdout, "MPI_Comm_size....\n"); fflush(stdout);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);   // get the number of cpus in our group

  fprintf(stdout, "MPI_Comm_rank....\n"); fflush(stdout);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);    // get this cpu's id number

  fprintf(stdout, "MPI_Get_processor_name...\n"); fflush(stdout);
  sprintf(myidtxt, CPUTXT, myid);   
  MPI_Get_processor_name(processor_name, &procnamelen);   // get cpu name (e.g. hostname)

  fprintf(stdout, "myargs....\n"); fflush(stdout);
  fprintf(stdout, "MYID = %d \n",myid); fflush(stdout);
  // currently INIT provides args to rest of processes

  if (MAXCPUS < numprocs) {
    fprintf(stdout, "Must increase MAXCPUS in global.h, %d is too many\n", numprocs);
    fail( FAIL_MPI_BASIC,0 ) ; 
  }

  set_pids(); 

  if( myid == printer_pid ) { 
    fprintf(stdout, "numprocs=%d \n", numprocs) ;
    fprintf(stdout, "proc: %s on %s\n", myidtxt, processor_name);
    fprintf(stdout, "end: init_MPI\n");
    fflush(stdout);
  }

  TRACE_END;
  return (0);
}


/*********************************************************************************/
/*********************************************************************************
  setup_mpibc_info_sym():
  ------------------
     -- sets special MPI-BC information that is is not obvious and is (usually)
        due to an imposed symmetry;
     -- uses macro "BC_TYPE_CHOICE" to handle choices; 
*********************************************************************************/
void setup_mpibc_info_sym( void )
{
  int i, nbr_pos[NDIM]; 

  TRACE_BEG;

#if( BC_TYPE_CHOICE == BC_CARTESIAN_PERIODIC )
  DLOOP1 nbr_pos[i] = cpupos[i]; 
  SDLOOP1 { 
    /* Down condition: */
    if( cpupos[i] == 0 ) { 
      nbr_pos[i] = ncpux[i]-1;
      bc_pid[i][BCDN] = get_pid(nbr_pos);  
      nbr_pos[i] = cpupos[i];
    }
    /* Up condition: */
    if( cpupos[i] == (ncpux[i]-1) ) { 
      nbr_pos[i] = 0;
      bc_pid[i][BCUP] = get_pid(nbr_pos);  
      nbr_pos[i] = cpupos[i];
    }
  }

#elif(TOP_TYPE_CHOICE == TOP_SPHERICAL || TOP_TYPE_CHOICE == TOP_CYLINDRICAL )
  DLOOP1 nbr_pos[i] = cpupos[i]; 
  i = 3;  /* Assume this is the phi index */

  /* Do not use MPI to send boundary information if we do not have to : */
  if( ncpux[i] == 1 ) { 
    bc_pid[i][BCDN] = bc_pid[i][BCUP] = BC_PHYS ; 
  }
  else { 
    /* Down condition: */
    if( cpupos[i] == 0 ) { 
      nbr_pos[i] = ncpux[i]-1;
      bc_pid[i][BCDN] = get_pid(nbr_pos);  
      nbr_pos[i] = cpupos[i];
    }
    /* Up condition: */
    if( cpupos[i] == (ncpux[i]-1) ) { 
      nbr_pos[i] = 0;
      bc_pid[i][BCUP] = get_pid(nbr_pos);  
      nbr_pos[i] = cpupos[i];
    }
  }

#endif 

  TRACE_END;
  return; 

}


/*********************************************************************************/
/*********************************************************************************
  bounds_mpi()
  ------------------
     -- responsible for the CPUs sending and receiving of boundary information; 
     -- in order to make sure that boundary cells have been set correctly,
        this routine is called after bounds(); 
     -- all cells (even ghosts) in the boundaries are sent;
     -- array indexing always follows the convention that x3(k) runs 
        fastest, then x2(j), then x1(i).  
     -- since the face transfer have the same dimensions as this cell's dimensions, 
          then we know what to expect;
     -- we loop over the faces and do  pack->send->recv per face;
        this allows the CPU to receive data while it is send and packing 
        others;
     -- the waiting is done after all the send/recv's have been initiated;
     -- since all communication is done between adjacent faces, we need to 
        communicate corners/edges (areas from "kiddie-corner" cells) 
        indirectly through the adjacent cell;  in order to do this, we 
        need to complete all the pack/send/recv/unpack operations 
        in a given direction before proceeding to the next direction;
     -- note that we are unpacking in a set order (not for any face at 
        any time)
*********************************************************************************/
void bounds_mpi( double ****prim_arr )
{
  int dir, side;
  int i,j,k, l, ipk, tag_send, tag_recv, req_id; 
  int ibeg, iend, jbeg, jend, kbeg, kend;


  TRACE_BEG;

  exit_status();  /* Make sure that all processes are alive before we wait */

  /*****************************************************************************
    PACK/SEND/RECV LOOP :
  *****************************************************************************/
  /* Loop over the directions (dir) times the two sides (side=up,dn): */ 
  for( dir = 1; dir < NDIM; dir++ ) {

    /* Need to pack/send/recv all out before we wait */
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 
      tag_send  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * myid;
      tag_recv  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * bc_pid[dir][side] - side + (!side);

      //      fprintf(stdout,"bounds_mpi() %d :  %d %d \n", myid,dir,side); fflush(stdout);

      /* "Pack" data to contiguous space for MPI transfer:  */
      ibeg = bc_range_send[dir][side][1][IBEG];   iend = bc_range_send[dir][side][1][IEND]; 
      jbeg = bc_range_send[dir][side][2][IBEG];   jend = bc_range_send[dir][side][2][IEND]; 
      kbeg = bc_range_send[dir][side][3][IBEG];   kend = bc_range_send[dir][side][3][IEND]; 
      ipk = 0;
      for(i=ibeg; i<=iend; i++)  for(j=jbeg; j<=jend; j++) for(k=kbeg; k<=kend; k++) {
	PLOOP bc_data_send[dir][side][ipk++] = prim_arr[i][j][k][l]; 
      }

      //      fprintf(stdout,"bounds_mpi() %d 1:  %d %d \n", myid,dir,side); fflush(stdout);

      //-check      if( bc_npts[dir][side] != ipk ) { 
      //-check	       fprintf(stderr,"bc_npts[%d][%d] mismatch with ipk on send %d : %d %d \n", 
      //-check		ir, side, myid, bc_npts[dir][side], ipk );
      //-check	       fflush(stderr);
      //-check      }

      /* Send data to neighbor :  */
      MPI_Isend( bc_data_send[dir][side], 
		 bc_npts[dir][side],
		 MPI_DOUBLE, 
		 bc_pid[dir][side],
		 tag_send,
		 MPI_COMM_WORLD, 
		 &(req_bc_out[req_id]) );

      //      fprintf(stdout,"bounds_mpi() %d 2:  %d %d \n", myid,dir,side); fflush(stdout);

      /* Receive data from same neighbor at same face (swapping faces, ala "Face-off"): */
      MPI_Irecv( bc_data_recv[dir][side], 
		 bc_npts[dir][side],
		 MPI_DOUBLE, 
		 bc_pid[dir][side],
		 tag_recv,
		 MPI_COMM_WORLD, 
		 &(req_bc_in[req_id]) );
    }


    /*****************************************************************************
       WAIT-RECV/UNPACK LOOP :
    *****************************************************************************/
    /* Now wait and unpack along this direction before continuing on to the next */
    for( side = 0; side < 2; side++ )  if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      //      fprintf(stdout,"bounds_mpi() %d 4:  %d %d \n", myid,dir,side); fflush(stdout);

      /* Wait for the data to arrive so that we can unpack it:  */
      MY_MPI_Wait( &(req_bc_in[req_id]), &mpi_status);

      /* Unpack the boundary into the pid's local grid array : */
      ibeg = bc_range_recv[dir][side][1][IBEG];   iend = bc_range_recv[dir][side][1][IEND]; 
      jbeg = bc_range_recv[dir][side][2][IBEG];   jend = bc_range_recv[dir][side][2][IEND]; 
      kbeg = bc_range_recv[dir][side][3][IBEG];   kend = bc_range_recv[dir][side][3][IEND]; 
      ipk = 0;
      for(i=ibeg; i<=iend; i++)  for(j=jbeg; j<=jend; j++) for(k=kbeg; k<=kend; k++) {
	PLOOP prim_arr[i][j][k][l] = bc_data_recv[dir][side][ipk++] ; 
      }

      //-check      if( bc_npts[dir][side] != ipk ) { 
      //-check	fprintf(stderr,"bc_npts[%d][%d] mismatch with ipk on recv %d : %d %d \n", 
      //-check		dir, side, myid, bc_npts[dir][side], ipk );
      //-check	fflush(stderr);
      //-check      }
    }

  }  /* end of for(dir..) and for(side...)  loops  */

  /*****************************************************************************
    SEND-WAIT LOOP : 
  *****************************************************************************/
  for( dir = 1; dir < NDIM; dir++ ) {
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      /* Wait for the data to be sent so that we can leave this routine:  */
      MY_MPI_Wait( &(req_bc_out[req_id]), &mpi_status);

    }
  }

  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;

  return; 
}

/*********************************************************************************/
/*********************************************************************************
  bounds_pflag_mpi()
  ------------------
     -- should be the same as bounds_mpi() except that we are sending 
        pflag[] data instead of grid function data ; 
     -- note that we are sending integers here; 
     -- we are passing (1/NP) as much data as in bounds_mpi();
*********************************************************************************/
void bounds_pflag_mpi( int ***pflag_arr )
{
  int dir, side;
  int i,j,k, l, ipk, tag_send, tag_recv, req_id; 
  int ibeg, iend, jbeg, jend, kbeg, kend;

  TRACE_BEG;

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    PACK/SEND/RECV LOOP :
  *****************************************************************************/
  /* Loop over the directions (dir) times the two sides (side=up,dn): */ 
  for( dir = 1; dir < NDIM; dir++ ) {

    /* Need to pack/send/recv all out before we wait */
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 
      tag_send  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * myid;
      tag_recv  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * bc_pid[dir][side] - side + (!side);

      /* "Pack" data to contiguous space for MPI transfer:  */
      ibeg = bc_range_send[dir][side][1][IBEG];   iend = bc_range_send[dir][side][1][IEND]; 
      jbeg = bc_range_send[dir][side][2][IBEG];   jend = bc_range_send[dir][side][2][IEND]; 
      kbeg = bc_range_send[dir][side][3][IBEG];   kend = bc_range_send[dir][side][3][IEND]; 
      ipk = 0;
      for(i=ibeg; i<=iend; i++)  for(j=jbeg; j<=jend; j++) for(k=kbeg; k<=kend; k++) {
	bc_pflag_send[dir][side][ipk++] = pflag_arr[i][j][k]; 
      }

      //-check      if( (bc_npts[dir][side]/NP) != ipk ) { 
      //-check	fprintf(stderr,"bc_npts[%d][%d] pflag mismatch with ipk on send %d : %d %d \n", 
      //-check		dir, side, myid, (bc_npts[dir][side]/NP), ipk );
      //-check	fflush(stderr);
      //-check      }

      /* Send data to neighbor :  */
      MPI_Isend( bc_pflag_send[dir][side], 
		 (bc_npts[dir][side]/NP),
		 MPI_INT, 
		 bc_pid[dir][side],
		 tag_send,
		 MPI_COMM_WORLD, 
		 &(req_bc_out[req_id]) );


      /* Receive data from same neighbor at same face (swapping faces, ala "Face-off"): */
      MPI_Irecv( bc_pflag_recv[dir][side], 
		 (bc_npts[dir][side]/NP),
		 MPI_INT, 
		 bc_pid[dir][side],
		 tag_recv,
		 MPI_COMM_WORLD, 
		 &(req_bc_in[req_id]) );

    }  /* end "side" loop */

    /*****************************************************************************
       WAIT-RECV/UNPACK LOOP :
    *****************************************************************************/
    /* Now wait and unpack along this direction before continuing on to the next */
    for( side = 0; side < 2; side++ )  if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      /* Wait for the data to arrive so that we can unpack it:  */
      MY_MPI_Wait( &(req_bc_in[req_id]), &mpi_status);

      /* Unpack the boundary into the pid's local grid array : */
      ibeg = bc_range_recv[dir][side][1][IBEG];   iend = bc_range_recv[dir][side][1][IEND]; 
      jbeg = bc_range_recv[dir][side][2][IBEG];   jend = bc_range_recv[dir][side][2][IEND]; 
      kbeg = bc_range_recv[dir][side][3][IBEG];   kend = bc_range_recv[dir][side][3][IEND]; 
      ipk = 0;
      for(i=ibeg; i<=iend; i++)  for(j=jbeg; j<=jend; j++) for(k=kbeg; k<=kend; k++) {
	pflag_arr[i][j][k] = bc_pflag_recv[dir][side][ipk++] ; 
      }

      //-check      if( (bc_npts[dir][side]/NP) != ipk ) { 
      //-check  	fprintf(stderr,"bc_npts[%d][%d] pflag mismatch with ipk on recv %d : %d %d \n", 
      //-check		dir, side, myid, (bc_npts[dir][side]/NP), ipk );
      //-check  	fflush(stderr);
      //-check      }
    }

  }  /* end of "dir" loop  */

    
  //  fprintf(stdout,"Send-wait  bounds_pflag_mpi() %d \n", myid); fflush(stdout);

  /*****************************************************************************
    SEND-WAIT LOOP : 
  *****************************************************************************/
  for( dir = 1; dir < NDIM; dir++ ) {
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      /* Wait for the data to be sent so that we can leave this routine:  */
      MY_MPI_Wait( &(req_bc_out[req_id]), &mpi_status);

    }
  }

  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;

  return; 
}

/*********************************************************************************/
/*********************************************************************************
  bounds_mpi_fast()
  ------------------
     -- supposedly a faster version of bounds_mpi()
     -- responsible for the CPUs sending and receiving of boundary information; 
     -- in order to make sure that boundary cells have been set correctly,
        this routine is called after bounds(); 
     -- all cells (even ghosts) in the boundaries are sent;
     -- array indexing always follows the convention that x3(k) runs 
        fastest, then x2(j), then x1(i).  
     -- since the face transfer have the same dimensions as this cell's dimensions, 
          then we know what to expect;
     -- we loop over the faces and do  pack->send->recv per face;
        this allows the CPU to receive data while it is sending and packing 
        others;
     -- all the pack->send->recv commands for all directions are executed before 
        starting the wait->unpack process;
     -- the order in which we unpack sets the boundaries in the correct order;
*********************************************************************************/
void bounds_mpi_fast( double ****prim_arr )
{
  int dir, side;
  int i,j,k, l, ipk, tag_send, tag_recv, req_id; 

  TRACE_BEG;

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    PACK/SEND/RECV LOOP :
  *****************************************************************************/
  /* Loop over the directions (dir) times the two sides (side=up,dn): */ 
  for( dir = 1; dir < NDIM; dir++ ) {

    /* Need to pack/send/recv all out before we wait */
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 
      tag_send  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * myid;
      tag_recv  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * bc_pid[dir][side] - side + (!side);

      //      fprintf(stdout,"bounds_mpi() %d :  %d %d \n", myid,dir,side); fflush(stdout);

      /* "Pack" data to contiguous space for MPI transfer:  */
      ipk = 0;
      for(i=bc_range_send[dir][side][1][IBEG]; i<=bc_range_send[dir][side][1][IEND]; i++)  
	for(j=bc_range_send[dir][side][2][IBEG]; j<=bc_range_send[dir][side][2][IEND]; j++) 
	  for(k=bc_range_send[dir][side][3][IBEG]; k<=bc_range_send[dir][side][3][IEND]; k++) {
	    PLOOP bc_data_send[dir][side][ipk++] = prim_arr[i][j][k][l]; 
	  }


      //      fprintf(stdout,"bounds_mpi() %d 1:  %d %d \n", myid,dir,side); fflush(stdout);

      //-check      if( bc_npts[dir][side] != ipk ) { 
      //-check	       fprintf(stderr,"bc_npts[%d][%d] mismatch with ipk on send %d : %d %d \n", 
      //-check		ir, side, myid, bc_npts[dir][side], ipk );
      //-check	       fflush(stderr);
      //-check      }

      /* Send data to neighbor :  */
      MPI_Isend( bc_data_send[dir][side], 
		 bc_npts[dir][side],
		 MPI_DOUBLE, 
		 bc_pid[dir][side],
		 tag_send,
		 MPI_COMM_WORLD, 
		 &(req_bc_out[req_id]) );

      //      fprintf(stdout,"bounds_mpi() %d 2:  %d %d \n", myid,dir,side); fflush(stdout);

      /* Receive data from same neighbor at same face (swapping faces, ala "Face-off"): */
      MPI_Irecv( bc_data_recv[dir][side], 
		 bc_npts[dir][side],
		 MPI_DOUBLE, 
		 bc_pid[dir][side],
		 tag_recv,
		 MPI_COMM_WORLD, 
		 &(req_bc_in[req_id]) );

    } /* end of for(side...) loop   */
  }  /* end of for(dir..) loop  */


  /*****************************************************************************
       WAIT-RECV/UNPACK LOOP :
  *****************************************************************************/
  for( dir = 1; dir < NDIM; dir++ ) {
    /* Now wait and unpack along this direction before continuing on to the next */
    for( side = 0; side < 2; side++ )  if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      //      fprintf(stdout,"bounds_mpi() %d 4:  %d %d \n", myid,dir,side); fflush(stdout);

      /* Wait for the data to arrive so that we can unpack it:  */
      MY_MPI_Wait( &(req_bc_in[req_id]), &mpi_status);

      /* Unpack the boundary into the pid's local grid array : */
      ipk = 0;
      for(i=bc_range_recv[dir][side][1][IBEG]; i<=bc_range_recv[dir][side][1][IEND]; i++)  
	for(j=bc_range_recv[dir][side][2][IBEG]; j<=bc_range_recv[dir][side][2][IEND]; j++) 
	  for(k=bc_range_recv[dir][side][3][IBEG]; k<=bc_range_recv[dir][side][3][IEND]; k++) {
	    PLOOP prim_arr[i][j][k][l] = bc_data_recv[dir][side][ipk++] ; 
	  }
      //-check      if( bc_npts[dir][side] != ipk ) { 
      //-check	fprintf(stderr,"bc_npts[%d][%d] mismatch with ipk on recv %d : %d %d \n", 
      //-check		dir, side, myid, bc_npts[dir][side], ipk );
      //-check	fflush(stderr);
      //-check      }

    } /* end of for(side...) loop   */
  }  /* end of for(dir..) loop  */

    
  /*****************************************************************************
    SEND-WAIT LOOP : 
  *****************************************************************************/
  for( dir = 1; dir < NDIM; dir++ ) {
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      /* Wait for the data to be sent so that we can leave this routine:  */
      MY_MPI_Wait( &(req_bc_out[req_id]), &mpi_status);

    }
  }

  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;

  return; 
}

/*********************************************************************************/
/*********************************************************************************
  bounds_pflag_mpi_fast()
  ------------------
     -- a faster version of bounds_pflag_mpi();
     -- should be the same as bounds_mpi() except that we are sending 
        pflag[] data instead of grid function data ; 
     -- note that we are sending integers here; 
     -- we are passing (1/NP) as much data as in bounds_mpi();
     -- we loop over the faces and do  pack->send->recv per face;
        this allows the CPU to receive data while it is sending and packing 
        others;
     -- all the pack->send->recv commands for all directions are executed before 
        starting the wait->unpack process;
     -- the order in which we unpack sets the boundaries in the correct order;
*********************************************************************************/
void bounds_pflag_mpi_fast( int ***pflag_arr )
{
  int dir, side;
  int i,j,k, l, ipk, tag_send, tag_recv, req_id; 

  TRACE_BEG;

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    PACK/SEND/RECV LOOP :
  *****************************************************************************/
  /* Loop over the directions (dir) times the two sides (side=up,dn): */ 
  for( dir = 1; dir < NDIM; dir++ ) {

    /* Need to pack/send/recv all out before we wait */
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 
      tag_send  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * myid;
      tag_recv  =  req_id  +  BASE_TAG_BOUNDS + 2 * 3 * bc_pid[dir][side] - side + (!side);

      /* "Pack" data to contiguous space for MPI transfer:  */
      ipk = 0;
      for(i=bc_range_send[dir][side][1][IBEG]; i<=bc_range_send[dir][side][1][IEND]; i++)  
	for(j=bc_range_send[dir][side][2][IBEG]; j<=bc_range_send[dir][side][2][IEND]; j++) 
	  for(k=bc_range_send[dir][side][3][IBEG]; k<=bc_range_send[dir][side][3][IEND]; k++) {
	    bc_pflag_send[dir][side][ipk++] = pflag_arr[i][j][k]; 
	  }

      //-check      if( (bc_npts[dir][side]/NP) != ipk ) { 
      //-check	fprintf(stderr,"bc_npts[%d][%d] pflag mismatch with ipk on send %d : %d %d \n", 
      //-check		dir, side, myid, (bc_npts[dir][side]/NP), ipk );
      //-check	fflush(stderr);
      //-check      }

      /* Send data to neighbor :  */
      MPI_Isend( bc_pflag_send[dir][side], 
		 (bc_npts[dir][side]/NP),
		 MPI_INT, 
		 bc_pid[dir][side],
		 tag_send,
		 MPI_COMM_WORLD, 
		 &(req_bc_out[req_id]) );


      /* Receive data from same neighbor at same face (swapping faces, ala "Face-off"): */
      MPI_Irecv( bc_pflag_recv[dir][side], 
		 (bc_npts[dir][side]/NP),
		 MPI_INT, 
		 bc_pid[dir][side],
		 tag_recv,
		 MPI_COMM_WORLD, 
		 &(req_bc_in[req_id]) );

    } /* end of for(side...) loop   */
  }  /* end of for(dir..) loop  */


    /*****************************************************************************
       WAIT-RECV/UNPACK LOOP :
    *****************************************************************************/
    /* Now wait and unpack along this direction before continuing on to the next */
  for( dir = 1; dir < NDIM; dir++ ) {
    for( side = 0; side < 2; side++ )  if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      /* Wait for the data to arrive so that we can unpack it:  */
      MY_MPI_Wait( &(req_bc_in[req_id]), &mpi_status);

      /* Unpack the boundary into the pid's local grid array : */
      ipk = 0;
      for(i=bc_range_recv[dir][side][1][IBEG]; i<=bc_range_recv[dir][side][1][IEND]; i++)  
	for(j=bc_range_recv[dir][side][2][IBEG]; j<=bc_range_recv[dir][side][2][IEND]; j++) 
	  for(k=bc_range_recv[dir][side][3][IBEG]; k<=bc_range_recv[dir][side][3][IEND]; k++) {
	    pflag_arr[i][j][k] = bc_pflag_recv[dir][side][ipk++] ; 
	  }

      //-check      if( (bc_npts[dir][side]/NP) != ipk ) { 
      //-check  	fprintf(stderr,"bc_npts[%d][%d] pflag mismatch with ipk on recv %d : %d %d \n", 
      //-check		dir, side, myid, (bc_npts[dir][side]/NP), ipk );
      //-check  	fflush(stderr);
      //-check      }

    } /* end of for(side...) loop   */
  }  /* end of for(dir..) loop  */

  //  fprintf(stdout,"Send-wait  bounds_pflag_mpi() %d \n", myid); fflush(stdout);

  /*****************************************************************************
    SEND-WAIT LOOP : 
  *****************************************************************************/
  for( dir = 1; dir < NDIM; dir++ ) {
    for( side = 0; side < 2; side++ ) if( bc_pid[dir][side] != BC_PHYS ) { 
      req_id = (dir-1)*2  + side  ; 

      /* Wait for the data to be sent so that we can leave this routine:  */
      MY_MPI_Wait( &(req_bc_out[req_id]), &mpi_status);

    }
  }

  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;

  return; 
}

/*********************************************************************************/
/*********************************************************************************
  history_mpi():
  ------------------
     -- handles summing history arrays together; 
     -- collect the data on the  ncpux[2]==ncpux[3]==1   nodes; 
*********************************************************************************/
void history_mpi( double *hist_arr[N_HISTORY][N_HIST_TYPES] )
{
  int i,j,k,l, pid_send, tag, pid_recv, i_mpi,ii,jj,kk;
  int tmpcpupos[NDIM];
  
  static MPI_Request req_hist_in;

  TRACE_BEG;

  if( (ncpux[2]*ncpux[3]) <= 1 ) { return; } 

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    Find the processor that is receiving data : 
  *****************************************************************************/
  tag = BASE_TAG_HIST;
  tmpcpupos[0] = 0;  tmpcpupos[1] = cpupos[1];  tmpcpupos[2] = 0;  tmpcpupos[3] = 0;
  pid_recv = get_pid(tmpcpupos);

  /*****************************************************************************
    Now  SEND/RECV  loop through all domains in the x2,x3 directions : 
  *****************************************************************************/
  i_mpi = 0;

  for(j=0; j<ncpux[2]; j++) for(k=0; k<ncpux[3]; k++)  if( !(j==0 && k==0) )  { 
    tmpcpupos[2] = j;   tmpcpupos[3] = k;
    pid_send = get_pid(tmpcpupos);

    /* Send stuff if we are a sender  : */
    if( myid == pid_send ) { 
      
      l = 0; 
      for(ii=0;ii<N_HISTORY;ii++) for(jj=0;jj<N_HIST_TYPES;jj++) for(kk=0;kk<N1;kk++) { hist_send[l++] = hist_arr[ii][jj][kk];  }

      MPI_Isend( hist_send,
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_recv,
		 tag,
		 MPI_COMM_WORLD, 
		 &request_out    );

    }

    //    fprintf(stdout,"history_mpi(): here2 %d-%d  pid=%d  nstep=%d \n",j,k,myid,nstep);fflush(stdout);

    /* Receive this pid_send's data */
    if( myid == pid_recv ) { 
      MPI_Irecv( history_buffer[i_mpi], 
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_send,
		 tag,
		 MPI_COMM_WORLD, 
		 &req_hist_in );

      /* Wait for the recv to finish: */
      MY_MPI_Wait( &req_hist_in, &mpi_status);

      /* Unpack and add it to the total so far: */
      l = 0; 
      for(ii=0;ii<N_HISTORY;ii++) for(jj=0;jj<N_HIST_TYPES;jj++) for(kk=0;kk<N1;kk++) { hist_arr[ii][jj][kk] += history_buffer[i_mpi][l++]; }

    }

    //    fprintf(stdout,"history_mpi(): here3 %d-%d  pid=%d  nstep=%d \n",j,k,myid,nstep);fflush(stdout);

    /* Wait for the send to finish: */
    if( myid == pid_send ) { 
      MY_MPI_Wait( &request_out, &mpi_status);
    }

    //    fprintf(stdout,"history_mpi(): here4 %d-%d  pid=%d  nstep=%d \n",j,k,myid,nstep);fflush(stdout);

    i_mpi++;
  }
  
  //  MPI_Barrier(MPI_COMM_WORLD);


  TRACE_END;
    
  return;

}

/*********************************************************************************/
/*********************************************************************************
  history_mpi2():
  ------------------
     -- collects all the domain rows into one global row on the master node;
*********************************************************************************/
void history_mpi2( double *hist_arr[N_HISTORY][N_HIST_TYPES], double *hist_out[N_HISTORY][N_HIST_TYPES] )
{
  int i,j,k,l,n,ipos, pid_send, tag, pid_recv, i_mpi,ii,jj,kk;
  int tmpcpupos[NDIM];
  
  int n_mpi;
  static MPI_Request req_hist_in2;

  TRACE_BEG;

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    Find the processor that is receiving data : 
  *****************************************************************************/
  tag = BASE_TAG_HIST2;
  tmpcpupos[0] = 0;  tmpcpupos[1] = 0;  tmpcpupos[2] = 0;  tmpcpupos[3] = 0;
  pid_recv = 0;

  /* First copy over the master's values to the output array:   */
  if( myid == pid_recv ) { 
    for(i=0; i<N_HISTORY; i++)  for(j=0;j<N_HIST_TYPES;j++) for(k=0;k<N1;k++) { 
      hist_out[i][j][k] = hist_arr[i][j][k]; 
    }
  }
  if( ncpux[1] <= 1 ) {  return;  } 

  /*****************************************************************************
    Setup initial values
  *****************************************************************************/
  n_mpi = ncpux[1] - 1  ;

  //  fprintf(stdout,"history_mpi2(): here1  pid=%d  nstep=%d \n",myid,nstep);fflush(stdout);
  /*****************************************************************************
    Now  SEND/RECV  loop through all cpupos[1]==0 domains  : 
  *****************************************************************************/
  i_mpi = 0;

  for(i=1; i<ncpux[1]; i++)  { 
    tmpcpupos[1] = i; 
    pid_send = get_pid(tmpcpupos);

    /* Send stuff if we are a sender  : */
    if( myid == pid_send ) { 

      l = 0; 
      for(ii=0;ii<N_HISTORY;ii++) for(jj=0;jj<N_HIST_TYPES;jj++) for(kk=0;kk<N1;kk++) { hist_send2[l++] = hist_arr[ii][jj][kk];  }

      MPI_Isend( hist_send2,
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_recv,
		 tag,
		 MPI_COMM_WORLD, 
		 &request_out    );

    }

    //    fprintf(stdout,"history_mpi2(): here2 %d  pid=%d  nstep=%d \n",i,myid,nstep);fflush(stdout);

    /* Receive this pid_send's data */
    if( myid == pid_recv ) { 
      MPI_Irecv( history_buffer2[i_mpi], 
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_send,
		 tag,
		 MPI_COMM_WORLD, 
		 &req_hist_in2 );

      /* Wait for the recv to finish: */
      MY_MPI_Wait( &req_hist_in2, &mpi_status);

    }

    //    fprintf(stdout,"history_mpi2(): here3 %d  pid=%d  nstep=%d \n",i,myid,nstep);fflush(stdout);

    if( myid == pid_send ) { 
      MY_MPI_Wait( &request_out, &mpi_status);
    }

    //    fprintf(stdout,"history_mpi2(): here4 %d  pid=%d  nstep=%d \n",i,myid,nstep);fflush(stdout);

    i_mpi++;
  }
  

  /*****************************************************************************
    Copy data from buffers to master nodes array: 
  *****************************************************************************/
  if( myid == pid_recv ) { 
    for( i_mpi=0; i_mpi<n_mpi; i_mpi++ ) { 
      /* Unpack and add it to the total so far: */
      l=0; 
      ipos = (i_mpi+1)*N1;
      for(i=0; i<N_HISTORY; i++)  for(j=0;j<N_HIST_TYPES;j++) for(k=0;k<N1;k++) { 
	hist_out[i][j][k+ipos] = history_buffer2[i_mpi][l++];
      }
    }

  }
    
  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;

  return;

}

/*********************************************************************************/
/*********************************************************************************
  history_mpi_fast():
  ------------------
     -- fast version of history_mpi();
     -- handles summing history arrays together; 
     -- collect the data on the  ncpux[2]==ncpux[3]==1   nodes; 
*********************************************************************************/
#define N_MPI_MAX (200) 
void history_mpi_fast( double hist_arr[N_HISTORY][N_HIST_TYPES][N1] )
{
  int i,j,k,l, pid_send, tag, pid_recv, i_mpi, is_sender=0;
  int tmpcpupos[NDIM];
  
  static int first_time = 1; 
  static int n_mpi;
  static MPI_Request req_hist_in[200], req_hist_out;

  TRACE_BEG;

  if( (ncpux[2]*ncpux[3]) <= 1 ) { return; } 

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    Find the processor that is receiving data : 
  *****************************************************************************/
  tag = BASE_TAG_HIST;
  tmpcpupos[0] = 0;  tmpcpupos[1] = cpupos[1];  tmpcpupos[2] = 0;  tmpcpupos[3] = 0;
  pid_recv = get_pid(tmpcpupos);


  /*****************************************************************************
    Setup initial values and allocations  
  *****************************************************************************/
  if( first_time ) { 
    n_mpi = ncpux[2] * ncpux[3] - 1  ;
    if( n_mpi > N_MPI_MAX ) { 
      fprintf(stderr,"history_mpi_fast(): n_mpi > N_MPI_MAX :  %d %d \n",n_mpi,N_MPI_MAX);
      fflush(stderr);
      fail(FAIL_MPI_BASIC,0);
    }

    /* Allocate the buffer used for recv communication : */
    history_buffer = (double **) calloc( n_mpi, sizeof(double*) );    
    if( history_buffer == NULL ) { 
      fprintf(stderr,
	      "history_mpi(): Cannot allocate history_buffer of length %d \n",
	      n_mpi);
      fflush(stderr);
      fail( FAIL_MPI_BASIC,0 );
    }
    for(i_mpi=0; i_mpi<n_mpi; i_mpi++ ) { 
      history_buffer[i_mpi] = (double *) calloc( N_HIST_POINTS, sizeof(double) );    
      if( history_buffer[i_mpi] == NULL ) { 
	fprintf(stderr,
		"history_mpi(): Cannot allocate history_buffer[%d] of length %d \n",
		i_mpi, N_HIST_POINTS);
	fflush(stderr);
	fail( FAIL_MPI_BASIC,0 );
      }
    }
    first_time = 0;
  }


  /*****************************************************************************
    Now  SEND/RECV  loop through all domains in the x2,x3 directions : 
  *****************************************************************************/
  i_mpi = 0;
  for(j=0; j<ncpux[2]; j++) for(k=0; k<ncpux[3]; k++)  if( !(j==0 && k==0) )  { 
    tmpcpupos[2] = j;   tmpcpupos[3] = k;
    pid_send = get_pid(tmpcpupos);

    /* Send stuff if we are a sender  : */
    if( myid == pid_send ) { 
      is_sender = 1;
      for(l=0; l<N_HIST_POINTS; l++) { hist_send[l] = hist_arr[0][0][l]; }
      MPI_Isend( hist_send,
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_recv,
		 tag,
		 MPI_COMM_WORLD, 
		 &req_hist_out );
    }

    /* Receive this pid_send's data */
    if( myid == pid_recv ) { 
      MPI_Irecv( history_buffer[i_mpi], 
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_send,
		 tag,
		 MPI_COMM_WORLD, 
		 &(req_hist_in[i_mpi]) );
    }
    i_mpi++;
  }
  
  /*****************************************************************************
    Now  WAIT for recvs and unpack 
  *****************************************************************************/
  i_mpi = 0; 
  if( myid == pid_recv ) { 
    for(j=0; j<ncpux[2]; j++) for(k=0; k<ncpux[3]; k++)  if( !(j==0 && k==0) )  { 
      tmpcpupos[2] = j;   tmpcpupos[3] = k;
      pid_send = get_pid(tmpcpupos);

      /* Wait for the recv to finish: */
      MY_MPI_Wait( &(req_hist_in[i_mpi]), &mpi_status);

      /* Unpack and add it to the total so far: */
      for(l=0; l<N_HIST_POINTS; l++) { hist_arr[0][0][l] += history_buffer[i_mpi][l]; }

      i_mpi++;
    }
  }
  
  /*****************************************************************************
    Now  WAIT for send 
  *****************************************************************************/
  if( is_sender ) { 
    MY_MPI_Wait( &req_hist_out, &mpi_status);
  }
    
  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;

  return;

}

/*********************************************************************************/
/*********************************************************************************
  history_mpi2_fast():
  ------------------
     -- fast version of history_mpi2();
     -- collects all the domain rows into one global row on the master node;
*********************************************************************************/
void history_mpi2_fast( double hist_arr[N_HISTORY][N_HIST_TYPES][N1], double *hist_out[N_HISTORY][N_HIST_TYPES] )
{
  int i,j,k,l,n,ipos, pid_send, tag, pid_recv, i_mpi, is_sender=0;
  int tmpcpupos[NDIM];
  
  static int first_time = 1; 
  static int n_mpi;
  static MPI_Request req_hist_in[N_MPI_MAX], req_hist_out;

  TRACE_BEG;

  if( ncpux[1] <= 1 ) {  return;  } 

  exit_status();  /* Make sure that all processes are alive */

  /*****************************************************************************
    Find the processor that is receiving data : 
  *****************************************************************************/
  tag = BASE_TAG_HIST2;
  tmpcpupos[0] = 0;  tmpcpupos[1] = 0;  tmpcpupos[2] = 0;  tmpcpupos[3] = 0;
  pid_recv = 0;

  /* First copy over the master's values to the output array:   */
  if( myid == pid_recv ) { 
    for(i=0; i<N_HISTORY; i++)  for(j=0;j<N_HIST_TYPES;j++) for(k=0;k<N1;k++) { 
      hist_out[i][j][k] = hist_arr[i][j][k]; 
    }
  }

  /*****************************************************************************
    Setup initial values and allocations  
  *****************************************************************************/
  if( first_time ) { 

    n_mpi = ncpux[1] - 1  ;
    if( n_mpi > N_MPI_MAX ) { 
      fprintf(stderr,"history_mpi2_fast(): n_mpi > N_MPI_MAX :  %d %d \n",n_mpi,N_MPI_MAX);
      fflush(stderr);
      fail(FAIL_MPI_BASIC,0);
    }

    /* Allocate the buffer used for communication : */
    history_buffer2 = (double **) calloc( n_mpi, sizeof(double*) );    
    if( history_buffer2 == NULL ) { 
      fprintf(stderr,
	      "history_mpi2(): Cannot allocate history_buffer2 of length %d \n",
	      n_mpi);
      fflush(stderr);
      fail( FAIL_MPI_BASIC,0 );
    }
    for(i_mpi=0; i_mpi<n_mpi; i_mpi++ ) { 
      history_buffer2[i_mpi] = (double *) calloc( N_HIST_POINTS, sizeof(double) );    
      if( history_buffer2[i_mpi] == NULL ) { 
	fprintf(stderr,
		"history_mpi2(): Cannot allocate history_buffer2[%d] of length %d \n",
		i_mpi, N_HIST_POINTS);
	fflush(stderr);
	fail( FAIL_MPI_BASIC,0 );
      }
    }

    first_time = 0;
  }

  /*****************************************************************************
    Now  SEND/RECV  loop through all cpupos[1]==0 domains  : 
  *****************************************************************************/
  i_mpi = 0;
  for(i=1; i<ncpux[1]; i++)  { 
    tmpcpupos[1] = i; 
    pid_send = get_pid(tmpcpupos);

    /* Send stuff if we are a sender  : */
    if( myid == pid_send ) { 
      is_sender = 1;
      for(l=0; l<N_HIST_POINTS; l++) { hist_send2[l] = hist_arr[0][0][l]; }
      MPI_Isend( hist_send2,
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_recv,
		 tag,
		 MPI_COMM_WORLD, 
		 &req_hist_out    );
    }

    /* Receive this pid_send's data */
    if( myid == pid_recv ) { 
      MPI_Irecv( history_buffer2[i_mpi], 
		 N_HIST_POINTS,
		 MPI_DOUBLE, 
		 pid_send,
		 tag,
		 MPI_COMM_WORLD, 
		 &(req_hist_in[i_mpi]) );
    }
    i_mpi++;
  }
  

  /*****************************************************************************
    Copy data from buffers to master nodes array: 
  *****************************************************************************/
  if( myid == pid_recv ) { 
    for( i_mpi=0; i_mpi<n_mpi; i_mpi++ ) { 
      MY_MPI_Wait( &(req_hist_in[i_mpi]), &mpi_status);
      /* Unpack and add it to the total so far: */
      l=0; 
      ipos = (i_mpi+1)*N1;
      for(i=0; i<N_HISTORY; i++)  for(j=0;j<N_HIST_TYPES;j++) for(k=0;k<N1;k++) { 
	hist_out[i][j][k+ipos] = history_buffer2[i_mpi][l++];
      }
    }
  }

  if( is_sender ) { 
    MY_MPI_Wait( &req_hist_out, &mpi_status);
  }
    
  //  MPI_Barrier(MPI_COMM_WORLD);

  TRACE_END;
  return;
}

#undef N_MPI_MAX


/*********************************************************************************/
/*********************************************************************************
  surface_mpi():
  ------------------
     -- handles summing surface arrays together; 
     -- collect the data on the  ncpux[2]==0   nodes; 
*********************************************************************************/
void surface_mpi( double *surf_arr[N_SURFACE][N_SURF_TYPES] )
{
  int i,j,k,l, pid_send, tag, pid_recv,ii,jj,kk;
  int tmpcpupos[NDIM];
  
  static MPI_Request req_surf_in;

  TRACE_BEG;

  exit_status();  /* Make sure that all processes are alive */

  if( ncpux[2] <= 1 ) { return; } 

  /*****************************************************************************
    Find the processor that is receiving data : 
  *****************************************************************************/
  tag = BASE_TAG_SURF;
  tmpcpupos[0] = 0;  tmpcpupos[1] = cpupos[1];  tmpcpupos[2] = 0;  tmpcpupos[3] = cpupos[3];
  pid_recv = get_pid(tmpcpupos);

  //  fprintf(stdout,"surface_mpi(): here1  pid=%d  nstep=%d \n",myid,nstep);fflush(stdout);

  /*****************************************************************************
    Now  SEND/RECV  loop through all domains in the x2,x3 directions : 
  *****************************************************************************/
  const int n13 = N1*N3;

  for(j=1; j<ncpux[2]; j++)  { 
    tmpcpupos[2] = j; 
    pid_send = get_pid(tmpcpupos);

    /* Send stuff if we are a sender  : */
    if( myid == pid_send ) { 
      
      l = 0; 
      for(ii=0;ii<N_SURFACE;ii++) for(jj=0;jj<N_SURF_TYPES;jj++) for(kk=0;kk<n13;kk++) {  surf_send[l++] = surf_arr[ii][jj][kk];   }

      MPI_Isend( surf_send,
		 N_SURF_POINTS,
		 MPI_DOUBLE, 
		 pid_recv,
		 tag,
		 MPI_COMM_WORLD, 
		 &request_out    );

      MY_MPI_Wait( &request_out, &mpi_status);

    }
    //    fprintf(stdout,"surface_mpi(): here2 %d-%d  pid=%d  nstep=%d \n",j,k,myid,nstep);fflush(stdout);

    /* Receive this pid_send's data */
    if( myid == pid_recv ) { 
      MPI_Irecv( surf_send,
		 N_SURF_POINTS,
		 MPI_DOUBLE, 
		 pid_send,
		 tag,
		 MPI_COMM_WORLD, 
		 &req_surf_in );

      /* Wait for the recv to finish: */
      MY_MPI_Wait( &req_surf_in, &mpi_status);

      /* Unpack and add it to the total so far: */
      l = 0; 
      for(ii=0;ii<N_SURFACE;ii++) for(jj=0;jj<N_SURF_TYPES;jj++) for(kk=0;kk<n13;kk++) {  surf_arr[ii][jj][kk] += surf_send[l++];   }
    }

    //    fprintf(stdout,"surface_mpi(): here3 %d-%d  pid=%d  nstep=%d \n",j,k,myid,nstep);fflush(stdout);

  }
  
  //  MPI_Barrier(MPI_COMM_WORLD);


  TRACE_END;
    
  return;

}

/*********************************************************************************/
/*********************************************************************************
  surface_mpi_faster():
  ------------------
     -- optimized ?   version of surface_mpi_faster();
     -- handles summing surface arrays together; 
     -- collect the data on the  ncpux[2]==0   nodes; 
*********************************************************************************/
void surface_mpi_faster( double *surf_arr[N_SURFACE][N_SURF_TYPES] )
{
  int i,j,k,l, pid_recv,ii,jj,kk;
  int tmpcpupos[NDIM];

  const int n13 = N1*N3;

  int *pids_surf;
  MPI_Group group_world, group_surf;
  static MPI_Comm  comm_surf;

  static  double *surf_recv;
  static int local_first_time = 1; 
  
  TRACE_BEG;

  if( ncpux[2] <= 1 ) { return; } 

  /*****************************************************************************
    Find the processor that is receiving data : 
  *****************************************************************************/
  tmpcpupos[0] = 0;  tmpcpupos[1] = cpupos[1];  tmpcpupos[2] = 0;  tmpcpupos[3] = cpupos[3];
  pid_recv = get_pid(tmpcpupos);


  /*****************************************************************************
    Organize the MPI Communication groups: 
  *****************************************************************************/
  if( local_first_time ) { 

    if( pid_recv == myid ) {  ALLOC_ARRAY(surf_recv ,N_SURF_POINTS);   }

    ALLOC_ARRAY(pids_surf ,ncpux[2]);     
    for(i = 0 ; i < ncpux[2] ; i++) { 
      tmpcpupos[2] = i;  
      pids_surf[i] = get_pid(tmpcpupos);
    }
    
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, ncpux[2], pids_surf, &group_surf);
    MPI_Comm_create(MPI_COMM_WORLD, group_surf, &comm_surf);

    DEALLOC_ARRAY(pids_surf ,ncpux[2]);     

    local_first_time = 0 ; 
  }


  /*****************************************************************************
    Load up the data to transfer : 
  *****************************************************************************/
  l = 0; 
  for(ii=0;ii<N_SURFACE;ii++) for(jj=0;jj<N_SURF_TYPES;jj++) for(kk=0;kk<n13;kk++) {  surf_send[l++] = surf_arr[ii][jj][kk];   }

  if( N_SURF_POINTS != l ) { 
    fprintf(stdout,"%s(): N_SURF_POINTS miscount :  %d %d \n", __func__,N_SURF_POINTS,l); 
    fflush(stdout);  fail(FAIL_BASIC,0);
  }

  /*****************************************************************************
    Transfer to "head" of the local communication group: 
  *****************************************************************************/
  exit_status();  
  MPI_Reduce(surf_send,surf_recv,N_SURF_POINTS,MPI_DOUBLE,MPI_SUM,pid_recv,comm_surf);

  /*****************************************************************************
     Unpack the already summed data: 
  *****************************************************************************/
  if( myid == pid_recv )  { 
    l = 0; 
    for(ii=0;ii<N_SURFACE;ii++) for(jj=0;jj<N_SURF_TYPES;jj++) for(kk=0;kk<n13;kk++) {  surf_arr[ii][jj][kk] = surf_recv[l++];   }
  }


  TRACE_END;
    
  return;

}

/***********************************************************************************/
/***********************************************************************************
  reduce_mask_operator():
  -------------------------
   -- function to be used in a MPI_Reduce operation that "returns" with the double 
         value if the mask is set, else returns with the first of the two values; 
   --  in other words,  given two ranks of "of_double_int" values, set1 and set2:

   --   result  = set1.val   if ( set1.mask==1  or (set1.mask==0&&set2.mask==0) else 
                = set2.val  

   -- This operation is commutative because of the constraints we put on its arguments 
       and its results, though does not appear to be so in its definition. 
       Our constraint is that we assume that there is only one unique value for "val" 
        for non-zero vales of "mask"; any numerical differences in practice are 
        assumed to be anomalies or round-off error which we will ignore;

   -- style follows example in open-mpi.org documenation; 
***********************************************************************************/
void reduce_mask_operator( of_double_int *in, of_double_int *inout, int *len, MPI_Datatype *dtype) 
{
  int i; 
  
  for( i = 0 ; i < *len; i++ ) { 

    /* see above: 

       inout->mask = in->mask  if ( in->mask==1  or  ( in->mask==0 && inout->mask==0 ) ) else
                     inout->mask

    inout->mask = ( in->mask || ( !in->mask && !inout->mask ) )  ?
      in->mask : inout->mask ;

    inout->val  = ( in->mask || ( !in->mask && !inout->mask ) )  ?
      in->val  : inout->val ;
    */

    /* if( in->mask || ( !in->mask && !inout->mask ) ) { */
    if( in->mask ) {
      inout->mask = in->mask ;
      inout->val  = in->val;
    }
    else {
      inout->mask = inout->mask ;
      inout->val  = inout->val;
    } 


    /* Advance pointers: */
    in++;
    inout++;  
  }

  return;
}


#endif   /*     #if(USEMPI)    */

/*********************************************************************************/
/*********************************************************************************
  set_pids();
  -----------------
     -- assigns ranks to tasks by assigning values to the various *_pid variables;

  HERE-HERE-HERE  We should put a routine here to implemented a more sophisticated assignment method: 

*********************************************************************************/
void set_pids(void)
{
  TRACE_BEG ; 

#if( USEMPI ) 
#else 
#endif 

  printer_pid  = 0;     /* id of rank that reports messages; */
       io_pid  = 0;     /* id of rank that gathers (if using GATHER_IO) and writes data to disk; */
    master_pid = 0;     /* id of rank that controls the entire job (usually equals "0"); */
  special1_pid = 0;     /* id of rank that performs special computations of type 1 */
  ener_out_pid = 0;     /* id of rank that write *.ener.out file   */
  traj_out_pid = 0;     /* id of rank that write *.ener.out file   */
    timers_pid = 0;     /* id of rank that is responsible for gathering and writing timing data */
  
  out_pid[OUT_HISTORY    ] = io_pid; 
  out_pid[OUT_ASCII      ] = io_pid; 
  out_pid[OUT_SDF        ] = io_pid; 
  out_pid[OUT_HDF5       ] = io_pid; 
  out_pid[OUT_IMAGE      ] = io_pid; 
  out_pid[OUT_STAT       ] = out_pid[OUT_HDF5];  /* we would have to arrange things different for these two variables to be different */
  out_pid[OUT_STAT2      ] = out_pid[OUT_HDF5];  /* we would have to arrange things different for these two variables to be different */
  out_pid[OUT_RADFLUX    ] = out_pid[OUT_HDF5];  /* we would have to arrange things different for these two variables to be different */
  out_pid[OUT_PHOTOSPHERE] = out_pid[OUT_HDF5];  /* we would have to arrange things different for these two variables to be different */
  out_pid[OUT_MIN_DT     ] = io_pid; 
  out_pid[OUT_SURFACE    ] = out_pid[OUT_HDF5];  /* we would have to arrange things different for these two variables to be different */

  TRACE_END ; 
  return;
}

/*********************************************************************************
   THE FOLLOWING ROUTINES ARE NOT CURRENTLY USED, BUT MAY BE USEFUL FOR MORE 
   SOPHISTICATED MPI APPLICATIONS:
*********************************************************************************/

/*********************************************************************************/
/*********************************************************************************
  get_local_ijk():
  -----------------
     -- return with the subgrid's LOCAL indices that correspond to 
        the given GLOBAL ones:
*********************************************************************************/
void  get_local_ijk( int i, int j, int k, int  *i_loc, int *j_loc, int *k_loc) 
{

  TRACE_BEG;

  *i_loc = (i % NM1);
  *j_loc = (j % NM2);
  *k_loc = (k % NM3);

  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  get_cpupos_ijk():
  ----------
     -- Given the X1,X2,X3 (i,j,k) indices of a point in the global domain, 
        return the MPI position of the point's subgrid 
*********************************************************************************/
void get_cpupos_ijk( int i, int j , int k, int tmpcpupos[NDIM] ) 
{ 

  TRACE_BEG;

  tmpcpupos[0] = 0;
  tmpcpupos[1] = i / NM1;
  tmpcpupos[2] = j / NM2;
  tmpcpupos[3] = k / NM3;

  TRACE_END;
  return;
}

/*********************************************************************************/
/*********************************************************************************
  get_pid_ijk():
  ----------
     --  returns the processor's id number (order in group) that contains the 
         given physical point specified by the set of global indices:
*********************************************************************************/
int get_pid_ijk( int i, int j , int k ) 
{ 
  int tmpcpupos[NDIM];
  
  TRACE_BEG;

  get_cpupos_ijk( i, j, k, tmpcpupos );

  TRACE_END;

  return( get_pid( tmpcpupos ) ); 
}

