/***********************************************************************************************/
/***********************************************************************************************

                                    C O R O N A 

------------------------------------------------------------------------------------------------

  FILE CONTAINING ROUTINES FOR CORONA PHYSICS, E.G., CALCULATING THE CORONA'S COOLING FUNCTION 
        
***********************************************************************************************/

#include "decs.h"
#include "units.h"
#if(USEMPI==1)
#include "mpi.h"
#endif


#define MAX_OPTICAL_DEPTH  (1.)
#define SMOOTH_FLUX (1) 
#define INTERP_FLUX (1) 

extern double cooling_func_hr_disk( int i, int j, int k, double *ph );
extern int get_pid( int tmpcpupos[NDIM] ) ;

/*********************************************************************************************/
/*********************************************************************************************
  setup_corona_cooling()
  --------------------
   -- support routine for the orona (inverse Compton) cooling function; 
   -- one must call thls  cooling_func_corona() is called locally; 
   -- this routine does the following:
         1) sets various global variables used from the coronal cooling function; 
         2) integrates the absorption coefficient (density) to find the photosphere; 
         3) integrates the thermal cooling rate over photosphere to get the radiation density in corona;

    -- the integrals require performing a local integral, sending the integral to a "server" node to 
        collect all the results, sum the sub-integrals to get the global integrals, and then communicate 
        back to the "client" nodes; 

    -- note that the opacity/optical-depth needs to be in CGS units, but all other calculations 
       can be in geometrized units; 
**********************************************************************************************/
void  setup_corona_cooling( double ****prim)
{

  enum equator_pos_type { Top_Eq, Over_Eq, Bottom_Eq } ; 
  enum equator_pos_type equator_pos; 

  int i,j,k,l;
  double dth, lumtmp;
  double *prim_loc;
  struct of_geom  *geom;
  struct of_coord *coords;
  struct of_state q ;

  int even_grid;

  /* Constants : */
  double  eta                = 0.06;
  double  mdotNUM            = 3.e-4; 
  double  mdot_frac          = 1.e-2;
  double  gtheta             = 4.;
  double  n_elec_per_H       = 1.2; 
  double  mOBS               = 10.*C_mSUN;; 
  double  mNUM               = 1.; 
  double  rOBS               = C_G*mOBS/(C_c2);; 
  double  rNUM               = mNUM; 
  double  L_edd              = C_4pi * C_G * mOBS * C_c / C_kappaT;
  double  mdot_edd           = L_edd / ( C_c2 * eta ); 
  double  mdotOBS            = mdot_frac * mdot_edd; 
  double  mdot2              = mdotOBS / (mdot_edd*eta); 
  double  c1_corona          = ( 16. * C_pi * C_mp * n_elec_per_H * mdot2 /  (C_me * mdotNUM)  );
  double  c2_corona          = ( C_mp*C_mp / (C_me*C_me*15.*gtheta) );
  double  c3_corona          = ( 4. * C_mp * 1.4 / (C_me * 2.3 ) ); 
  double  lengthScale        = rOBS / rNUM;     
  double  rhoScale           = mdotOBS / ( mdotNUM * lengthScale * lengthScale * (C_c) ); 
  double  tau_scale          = C_kappaT * rhoScale * lengthScale ;


  TRACE_BEG;

  static int local_first_time = 1; 
  
  if( local_first_time ) { 
    fprintf(stdout,"cooling_func_corona():  eta           =  %26.16e \n",eta           ); 
    fprintf(stdout,"cooling_func_corona():  mdotNUM       =  %26.16e \n",mdotNUM       ); 
    fprintf(stdout,"cooling_func_corona():  mdot_frac     =  %26.16e \n",mdot_frac     ); 
    fprintf(stdout,"cooling_func_corona():  gtheta        =  %26.16e \n",gtheta        ); 
    fprintf(stdout,"cooling_func_corona():  n_elec_per_H  =  %26.16e \n",n_elec_per_H  ); 
    fprintf(stdout,"cooling_func_corona():  mOBS          =  %26.16e \n",mOBS          ); 
    fprintf(stdout,"cooling_func_corona():  mNUM          =  %26.16e \n",mNUM          ); 
    fprintf(stdout,"cooling_func_corona():  rOBS          =  %26.16e \n",rOBS          ); 
    fprintf(stdout,"cooling_func_corona():  rNUM          =  %26.16e \n",rNUM          ); 
    fprintf(stdout,"cooling_func_corona():  L_edd         =  %26.16e \n",L_edd         ); 
    fprintf(stdout,"cooling_func_corona():  mdot_edd      =  %26.16e \n",mdot_edd      ); 
    fprintf(stdout,"cooling_func_corona():  mdotOBS       =  %26.16e \n",mdotOBS       ); 
    fprintf(stdout,"cooling_func_corona():  mdot2         =  %26.16e \n",mdot2         ); 
    fprintf(stdout,"cooling_func_corona():  c1_corona     =  %26.16e \n",c1_corona     ); 
    fprintf(stdout,"cooling_func_corona():  c2_corona     =  %26.16e \n",c2_corona     ); 
    fprintf(stdout,"cooling_func_corona():  c3_corona     =  %26.16e \n",c3_corona     ); 
    fprintf(stdout,"cooling_func_corona():  lengthScale   =  %26.16e \n",lengthScale   ); 
    fprintf(stdout,"cooling_func_corona():  rhoScale      =  %26.16e \n",rhoScale      ); 
    fprintf(stdout,"cooling_func_corona():  tau_scale     =  %26.16e \n",tau_scale     ); 
    fflush(stdout);
    local_first_time = 0 ; 
  }

#if( USE_COOLING_FUNCTION != 4 )
  fprintf(stderr,"cooling_func_corona(): Bad value USE_COOLING_FUNCTION  should not be herer  %d   !! \n",USE_COOLING_FUNCTION); 
  fflush(stderr);    fail(FAIL_BASIC,0);
#endif

  /***************************************************************************************************
    Determine proximity to equator and position in the communication chain:
  ***************************************************************************************************/
  /* Always assume that jglob=globalsize[2]/2 is the equator:  */
  if( (totalsize[2] % 2) == 0 ) {  even_grid = 1;  }
  else                          {  even_grid = 0;  } 

  l = totalsize[2]/2;

  int jmid0 = l - globalpos[2] + N2S ;   /* the middle middle */
  int jmid1 = l - globalpos[2] + N2S ;   /* index closest to middle from the left  */
  int jmid2 = l - globalpos[2] + N2S ;   /* index closest to middle from the right */

  if( even_grid ) { jmid1 -= 1; }
    
#if(USEMPI) 
  int pid_recv_top = BC_PHYS;
  int pid_recv_bot = BC_PHYS;
  int pid_send_top = BC_PHYS;
  int pid_send_bot = BC_PHYS;
  int nbr_pos[NDIM]; 

  /* collector's pid is either the equator containing domain or the "bottom" equator-bordering domain: */
  DLOOP1 { nbr_pos[i] = cpupos[i]; }
  nbr_pos[2] = (l / N2);
  int collector_pid = get_pid(nbr_pos);

  /* we integrate from jmid2 to jend (forward), and from jmid1 to jbeg (backward) */
  if( jmid2 < N2S ) { 
    equator_pos = Bottom_Eq ;  /* after the equator */
    jmid1 = jmid2 = N2S;
    pid_recv_bot = bc_pid[2][BCUP];
    pid_send_bot = bc_pid[2][BCDN];
  }
  else if( jmid1 > N2E ) { 
    equator_pos = Top_Eq;   /* before the equator */
    jmid1 = jmid2 = N2E;
    pid_recv_top = bc_pid[2][BCDN];
    pid_send_top = bc_pid[2][BCUP];
  }
  else { 
    equator_pos = Over_Eq ;  /* intersecting the equator */
    pid_recv_top = bc_pid[2][BCDN];
    pid_recv_bot = bc_pid[2][BCUP];
    pid_send_top = pid_send_bot = BC_PHYS;  /* assume that we don't need to send anything */
    if( even_grid ) { 
      /* only need to send anything if we are the top equator domain bordering the equator: */
      if( jmid1 == N2E ) { 
	pid_send_top = collector_pid; 
	pid_recv_bot = BC_PHYS;
	if( collector_pid == myid ) { 
	  fprintf(stderr,"setup_corona_cooling(): should not be here 4341 :  %d \n", collector_pid); 
	  fflush(stderr);   fail(FAIL_BASIC,0);
	}
      }
    }
  }
#else 
  equator_pos = Over_Eq; 
#endif

  if( equator_pos == Over_Eq ) { 
    if( !even_grid ) { 
      jmid1 -= 1;
      jmid2 += 1;
    }
  }
    
  /***************************************************************************************************
    If not on the equator, start waiting for data from domains closer to the equator:
  ***************************************************************************************************/
#if( USEMPI ) 
  MPI_Status mpi_status;
  static MPI_Request req_in_top1, req_in_top2, req_in_bot1, req_in_bot2;
  static MPI_Request req_out_top1, req_out_top2, req_out_bot1, req_out_bot2;
  static MPI_Request req_in_tot1; 
  static MPI_Request req_in_tot2; 
  static MPI_Request req_out_tot1[100]; 
  static MPI_Request req_out_tot2[100]; 
  if( ncpux[2] > 100 ) { 
    fprintf(stderr,"setup_corona_cooling(): need to increase size of req_out_tot1 and req_out_tot2 :  %d \n", ncpux[2]); 
    fflush(stderr);   fail(FAIL_BASIC,0);
  }

  int npts  = N1TOT*N3TOT;
  int npts2 = 2*N1TOT*N3TOT;
  int tag_top1 = 1*numprocs;
  int tag_top2 = 2*numprocs;
  int tag_bot1 = 3*numprocs;
  int tag_bot2 = 4*numprocs;
  int tag_tot1 = 5*numprocs;
  int tag_tot2 = 6*numprocs;
  int *pids_out; 

  /* Receive data neighbor if we need it: */
  exit_status();

#if(TRACE_CALLS) 
  trace_message("MPI_Irecv-1",0);
#endif
  if( pid_recv_top >= 0 ) { 
    MPI_Irecv( irecv_buffer_top, npts,  MPI_INT   , pid_recv_top, tag_top1, MPI_COMM_WORLD, &req_in_top1);
    MPI_Irecv(  recv_buffer_top, npts2, MPI_DOUBLE, pid_recv_top, tag_top2, MPI_COMM_WORLD, &req_in_top2);
  }
  if( pid_recv_bot >= 0 ) { 
    MPI_Irecv( irecv_buffer_bot, npts,  MPI_INT   , pid_recv_bot, tag_bot1, MPI_COMM_WORLD, &req_in_bot1);
    MPI_Irecv(  recv_buffer_bot, npts2, MPI_DOUBLE, pid_recv_bot, tag_bot2, MPI_COMM_WORLD, &req_in_bot2);
  }
  if( myid != collector_pid ) { 
    MPI_Irecv( irecv_buffer_tot, npts2, MPI_INT   , collector_pid, tag_tot1, MPI_COMM_WORLD, &req_in_tot1);
    MPI_Irecv(  recv_buffer_tot, npts,  MPI_DOUBLE, collector_pid, tag_tot2, MPI_COMM_WORLD, &req_in_tot2);
  }
#if(TRACE_CALLS) 
  trace_message("MPI_Irecv-1",1);
#endif

#endif

  /***************************************************************************************************
    Calculate dtau and thermal emissivity over volume:
      -- need to cover all x1-x3 domain so that we can perform boxcar smoothing of the emissivity 
         in the corona to reduce the large cell-to-cell gradients that arise when there is no cooling 
         in a cell;
  ***************************************************************************************************/
  int bound; 
  N1ALL_LOOP  N2_LOOP  N3ALL_LOOP { 
    //HERE-HERE   make sure all this is needed
    prim_loc = prim[i][j][k];

    get_coord(i,j,k,CENT,ncurr,coords);      
    get_geometry(i,j,k,CENT,ncurr,geom); 
    get_state(prim_loc,geom,&q) ;  

    bound = (  (q.ucov[0]*( q.p + prim_loc[UU] + prim_loc[RHO] + q.bsq ))  >  (-prim_loc[RHO]) ) ?  1 : 0 ; 

    //    dth = sqrt(geom->gcov[TH][TH]) * dx[TH] * lengthScale; 
    dth = sqrt(geom->gcov[TH][TH]) * dx[TH];
    dtau[i][j][k] = tau_scale * dth * prim_loc[RHO]; 

    lumtmp = cooling_func_hr_disk(i,j,k,prim_loc); 
    lum_thermal[i][j][k] = lumtmp * bound; 
    dflux[i][j][k]       = dth * lumtmp;
  }

  /***************************************************************************************************
    If at the equator: 
             1) integrate from the equator outward to find optical depth of local volume, tau(r,ph); 
             2) record theta(r,ph) coordinate where optical depth exceeds  MAX_OPTICAL_DEPTH; 
             3) communicate theta(r,ph) and tau(r,ph) to relevant neighbors; 

    If not that at the equator: 
             1) wait for neighbor nearer to the equator to give you theta(r,ph), tau(r,ph); 
             2) using those as boundary data, continue integrating outward to fill in theta(r,ph) further; 
             3) communicate 
         
  ***************************************************************************************************/

  /***************************************************************************************************
      Initialize integration arrays: 
  ***************************************************************************************************/

  N1ALL_LOOP N3ALL_LOOP { 
     j_photo_bot[i][k]  =  j_photo_top[i][k] = -1 ; 
         tau_bot[i][k]  =      tau_top[i][k] =  0.;
        flux_bot[i][k] =      flux_top[i][k] =  0.;
  }


#if( USEMPI ) 
  /* Wait for the data to arrive so that we can unpack it:  */
#if(TRACE_CALLS) 
  trace_message("MPI_Wait-1",0);
#endif
  if( pid_recv_top >= 0 ) { 
    MPI_Wait( &req_in_top1, &mpi_status);
    MPI_Wait( &req_in_top2, &mpi_status);
    l = 0; 
    N1ALL_LOOP  N3ALL_LOOP { j_photo_top[i][k] = irecv_buffer_top[l++];   }
    l = 0; 
    N1ALL_LOOP  N3ALL_LOOP {     tau_top[i][k] =  recv_buffer_top[l++];   }
    N1ALL_LOOP  N3ALL_LOOP {    flux_top[i][k] =  recv_buffer_top[l++];   }
  }
#if(TRACE_CALLS) 
  trace_message("MPI_Wait-1-2",0);
#endif

  if( pid_recv_bot >= 0 ) { 
    MPI_Wait( &req_in_bot1, &mpi_status);
    MPI_Wait( &req_in_bot2, &mpi_status);
    l = 0; 
    N1ALL_LOOP  N3ALL_LOOP { j_photo_bot[i][k] = irecv_buffer_bot[l++];   }
    l = 0; 
    N1ALL_LOOP  N3ALL_LOOP {     tau_bot[i][k] =  recv_buffer_bot[l++];   }
    N1ALL_LOOP  N3ALL_LOOP {    flux_bot[i][k] =  recv_buffer_bot[l++];   }
  }
#if(TRACE_CALLS) 
  trace_message("MPI_Wait-1",1);
#endif

#endif

  /***************************************************************************************************
     Perform integrations :
  ***************************************************************************************************/
  double tau_tmp; 
  int jbeg1;

  /* Integration from top to equator:  */
  N1ALL_LOOP  N3ALL_LOOP { 
    if( j_photo_top[i][k] < 0 ) {   /* if we have already hit the photosphere we can stop integrating */
      tau_tmp = tau_top[i][k]; 
      for( j=N2S; j<=jmid1; j++)  { 
	tau_tmp += dtau[i][j][k] ; 
	if( tau_tmp >= MAX_OPTICAL_DEPTH ) { 
	  j_photo_top[i][k] = j - N2S + globalpos[2] ; 
	  break; 
	}
      }
      tau_top[i][k]     = tau_tmp; 
    }
    /* Note that there is not an "else" here because we want to catch the continuation of the previous loop, too */ 
    if( j_photo_top[i][k] >= 0 ) { 
      jbeg1 = j_photo_top[i][k] - globalpos[2] + N2S;
      jbeg1 = MAX(jbeg1,N2S);
      for( j=jbeg1; j<=jmid1; j++)  { 
	flux_top[i][k] += dflux[i][j][k];
      }
    }
  }
  
  /* Integration from bottom to equator:  */
  N1ALL_LOOP  N3ALL_LOOP { 
    if( j_photo_bot[i][k] < 0 ) {   /* if we have already hit the photosphere we can stop integrating */
      tau_tmp = tau_bot[i][k]; 
      for( j=N2E; j>=jmid2; j--)  { 
	tau_tmp += dtau[i][j][k] ; 
	if( tau_tmp >= MAX_OPTICAL_DEPTH ) { 
	  j_photo_bot[i][k] = j - N2S + globalpos[2] ; 
	  break; 
	}
      }
      tau_bot[i][k]     = tau_tmp; 
    }
    /* Note that there is not an "else" here because we want to catch the continuation of the previous loop, too */ 
    if( j_photo_bot[i][k] >= 0 ) { 
      jbeg1 = j_photo_bot[i][k] - globalpos[2] + N2S;
      jbeg1 = MIN(jbeg1,N2E);
      for( j=jbeg1; j>=jmid2; j--)  { 
	flux_bot[i][k] += dflux[i][j][k];
      }
    }
  }

  /* Include contribution from equator cell: */
  if( !even_grid ) { 
    if( equator_pos == Over_Eq ) { 
      j = jmid0;

      N1ALL_LOOP  N3ALL_LOOP { 
	/*  Top: */
	if( j_photo_top[i][k] < 0 ) { 
	  tau_top[i][k] += 0.5*dtau[i][j][k] ; 
	  if( tau_top[i][k] >= MAX_OPTICAL_DEPTH ) { 
	    j_photo_top[i][k] = j - N2S + globalpos[2] ; 
	  }
	}
	/* If there is a photosphere, then the middle is always in it: */
	if( j_photo_top[i][k] >= 0 ) { 
	  flux_top[i][k] += 0.5*dflux[i][j][k];
	}

	/*  Bottom: */
	if( j_photo_bot[i][k] < 0 ) { 
	  tau_bot[i][k] += 0.5*dtau[i][j][k] ; 
	  if( tau_bot[i][k] >= MAX_OPTICAL_DEPTH ) { 
	    j_photo_bot[i][k] = j - N2S + globalpos[2] ; 
	  }
	}
	/* If there is a photosphere, then the middle is always in it: */
	if( j_photo_bot[i][k] >= 0 ) { 
	  flux_bot[i][k] += 0.5*dflux[i][j][k];
	}
      }
    }
  }
  
  /***************************************************************************************************
    Once we are integrated to the local domain, send on results to the neighbor: 
  ***************************************************************************************************/
#if( USEMPI ) 
#if(TRACE_CALLS) 
  trace_message("MPI-Isend-1",0);
#endif
  if( pid_send_top >= 0 ) { 
      l = 0; 
      N1ALL_LOOP  N3ALL_LOOP { isend_buffer_top[l++] = j_photo_top[i][k];  }
      l = 0; 
      N1ALL_LOOP  N3ALL_LOOP {  send_buffer_top[l++] =  tau_top[i][k];     }
      N1ALL_LOOP  N3ALL_LOOP {  send_buffer_top[l++] = flux_top[i][k];     }
      MPI_Isend( isend_buffer_top, npts,  MPI_INT   , pid_send_top, tag_top1, MPI_COMM_WORLD, &req_out_top1);
      MPI_Isend(  send_buffer_top, npts2, MPI_DOUBLE, pid_send_top, tag_top2, MPI_COMM_WORLD, &req_out_top2);
  }
  if( pid_send_bot >= 0 ) { 
      l = 0; 
      N1ALL_LOOP  N3ALL_LOOP { isend_buffer_bot[l++] = j_photo_bot[i][k];  }
      l = 0; 
      N1ALL_LOOP  N3ALL_LOOP {  send_buffer_bot[l++] =  tau_bot[i][k];     }
      N1ALL_LOOP  N3ALL_LOOP {  send_buffer_bot[l++] = flux_bot[i][k];     }
      MPI_Isend( isend_buffer_bot, npts,  MPI_INT   , pid_send_bot, tag_bot1, MPI_COMM_WORLD, &req_out_bot1);
      MPI_Isend(  send_buffer_bot, npts2, MPI_DOUBLE, pid_send_bot, tag_bot2, MPI_COMM_WORLD, &req_out_bot2);
  }
#if(TRACE_CALLS) 
  trace_message("MPI-Isend-1",1);
#endif

  /***************************************************************************************************
    If we are the collector, then we have already collected neighboring data and integrated 
     to the equator.  Now we need to sum for total disk flux and then communicate this and the 
     location of the photosphere to all the receivers:
  ***************************************************************************************************/
  if( myid == collector_pid ) { 

    N1ALL_LOOP  N3ALL_LOOP { flux_disk[i][k] = flux_top[i][k] + flux_bot[i][k];  }

    l = 0 ; 
    N1ALL_LOOP  N3ALL_LOOP {  send_buffer_tot[l++] =  flux_disk[i][k];     }
    l = 0 ; 
    N1ALL_LOOP  N3ALL_LOOP { isend_buffer_tot[l++] = j_photo_top[i][k];    }
    N1ALL_LOOP  N3ALL_LOOP { isend_buffer_tot[l++] = j_photo_bot[i][k];    }

//      int ktmp2 = 32 - globalpos[3];  /* here the number is the "global iphi" to print */
//      int phtmp2 = ktmp2 + N3S; 
//      if( (ktmp2 >= 0) && (ktmp2 < N3) ) { 
//        N1_LOOP { 
//  	get_coord(i,N2S,phtmp2,CENT,ncurr,coords);      
//  	fprintf(stdout,"fluxdiskcollector %d  %26.16e %26.16e %26.16e %26.16e \n",globalpos[1]+i,t,coords->r,coords->x[3],flux_disk[i][phtmp2]);
//        }
//        fflush(stdout); 
//      }

    /* Setup collector for receiving all integration results */
    ALLOC_ARRAY(pids_out, ncpux[2]); 

    /* make list of processes to send data to : */
    DLOOP1 { nbr_pos[i] = cpupos[i]; }
    for( j=0; j < ncpux[2]; j++ ) { 
      nbr_pos[2] = j; 
      pids_out[j] = get_pid(nbr_pos);
    }

#if(TRACE_CALLS) 
    trace_message("MPI-Isend-collector-2",0);
#endif
    for( j=0; j < ncpux[2]; j++ )  if( myid != pids_out[j] )   { 
      MPI_Isend( isend_buffer_tot, npts2, MPI_INT   , pids_out[j], tag_tot1, MPI_COMM_WORLD, &req_out_tot1[j]);
      MPI_Isend(  send_buffer_tot, npts , MPI_DOUBLE, pids_out[j], tag_tot2, MPI_COMM_WORLD, &req_out_tot2[j]);
    }
#if(TRACE_CALLS) 
    trace_message("MPI-Isend-collector-2",1);
#endif
  }
  else { 
#if(TRACE_CALLS) 
    trace_message("MPI-Wait-2",0);
#endif
    MPI_Wait( &req_in_tot1, &mpi_status);
    l = 0 ; 
    N1ALL_LOOP  N3ALL_LOOP { j_photo_top[i][k] = irecv_buffer_tot[l++] ; }
    N1ALL_LOOP  N3ALL_LOOP { j_photo_bot[i][k] = irecv_buffer_tot[l++] ; }

    MPI_Wait( &req_in_tot2, &mpi_status);
    l = 0 ; 
    N1ALL_LOOP  N3ALL_LOOP {  flux_disk[i][k] = recv_buffer_tot[l++];    }
#if(TRACE_CALLS) 
    trace_message("MPI-Wait-2",1);
#endif
  }

#endif

  /***************************************************************************************************
    Perform box car smoothing over neighboring cells to smooth out the flux: 
  ***************************************************************************************************/
#if( SMOOTH_FLUX )
  N1ALL_LOOP  N3ALL_LOOP {  flux_bot[i][k] = flux_disk[i][k]; }
  lumtmp = 1./15.;
  for(i=2; i<N1TOT-2; i++)  for(k=2; k<N3TOT-2; k++)  { 
      flux_bot[i][k] = 
	3.*flux_disk[i  ][k  ] + 
	flux_disk[i-2][k  ] + 
	flux_disk[i+2][k  ] + 
	flux_disk[i  ][k-2] + 
	flux_disk[i  ][k+2] +
	flux_disk[i-1][k-1] + 
	flux_disk[i-1][k  ] + 
	flux_disk[i-1][k+1] + 
	flux_disk[i  ][k-1] + 
	flux_disk[i  ][k+1] +
	flux_disk[i+1][k-1] + 
	flux_disk[i+1][k  ] + 
	flux_disk[i+1][k+1] ;
    }
  for(i=1  ;i<N1TOT-1;i++)  for(k=1;k<N3TOT-1;k++)  {  flux_disk[i][k] = flux_bot[i][k] * lumtmp;  }
#endif

  /***************************************************************************************************
    Interpolate disk's flux over regions without any flux, else we get patchy, discrete coronal cooling:
  ***************************************************************************************************/
#if( INTERP_FLUX )

#endif

  //  /* adjust photosphere indices to local domain: */
//   N1ALL_LOOP  N3ALL_LOOP if( j_photo_top[i][k] >= 0 ) {  j_photo_top[i][k] -= jtmp; }
//   N1ALL_LOOP  N3ALL_LOOP if( j_photo_bot[i][k] >= 0 ) {  j_photo_bot[i][k] -= jtmp; }

  /* Equally split the disk's flux between the two hemispheres: */
  N1ALL_LOOP  N3ALL_LOOP { flux_disk[i][k] *= 0.5 ; }

  /***************************************************************************************************
    Now each processor should have all the necessary data to calculate the cooling functions : 
  ***************************************************************************************************/
  int jtmp = globalpos[2] - N2S; 
  int jglob;
  double pgas, rho, Te, urad, cool0, cool1, cool2; 

  l = 0 ; 
  LOOP { 
    cool2 = 0.;
    lumtmp = lum_thermal[i][j][k];
    cool0 = cool1 = lumtmp;
    
    jglob = j + jtmp; 

    /* We use coronal cooling function outside the photosphere : */
    if( (j_photo_top[i][k] < 0) || (jglob < j_photo_top[i][k]) || (jglob > j_photo_bot[i][k]) ) { 
      urad = flux_disk[i][k]; 
      if( urad > 0. ) { 
	prim_loc = prim[i][j][k];
	rho   = prim_loc[RHO]; 
	pgas  = (gam-1.)*prim_loc[UU]; 
	Te    =  pgas / ( 2.*rho + c2_corona*urad ) ;
	cool0 = cool2 = c1_corona * rho * Te * urad * ( 1. + c3_corona*Te ); 
      }
    }
    coolflux[0][l] = cool0;  /* final  */
    coolflux[1][l] = cool1;  /* disk   */
    coolflux[2][l] = cool2;  /* corona */
    l++;
  }

  int ktmp = 32 - globalpos[3];  /* here the number is the "global iphi" to print */
  int phtmp = ktmp + N3S; 
  if( (ktmp >= 0) && (ktmp < N3) ) { 
    N1_LOOP { 
      get_coord(i,N2S,phtmp,CENT,ncurr,coords);      
      fprintf(stdout,"flux_disk %d  %d  %d %d  %26.16e %26.16e %26.16e %26.16e \n",
	      globalpos[1]+i,cpupos[2],j_photo_top[i][phtmp],j_photo_bot[i][phtmp],
	      t,coords->r,coords->x[3],flux_disk[i][phtmp]);
    }
    fflush(stdout); 
  }

#if( USEMPI )
  /* Wait for all outstanding MPI communications : */
#if(TRACE_CALLS) 
  trace_message("MPI-Wait-3",0);
#endif
  if( pid_send_top >= 0 ) { 
    MPI_Wait( &req_out_top1, &mpi_status);
    MPI_Wait( &req_out_top2, &mpi_status);
  }
  if( pid_send_bot >= 0 ) { 
    MPI_Wait( &req_out_bot1, &mpi_status);
    MPI_Wait( &req_out_bot2, &mpi_status);
  }
  if( myid == collector_pid )   { 
    for( j=0; j < ncpux[2]; j++ ) if( myid != pids_out[j] ) { 
      MPI_Wait( &req_out_tot1[j], &mpi_status);
      MPI_Wait( &req_out_tot2[j], &mpi_status);
    }
    DEALLOC_ARRAY(pids_out, ncpux[2]); 
  }
#if(TRACE_CALLS) 
  trace_message("MPI-Wait-3",1);
#endif
#endif


  TRACE_END;

  return;
}
