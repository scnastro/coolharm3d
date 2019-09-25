/***********************************************************************************************/
/***********************************************************************************************

                                    C O R O N A 2

------------------------------------------------------------------------------------------------

  FILE CONTAINING ROUTINES FOR CORONA PHYSICS, E.G., CALCULATING THE CORONA'S COOLING FUNCTION 

-- look for 
      "//--OPTIMIZE" 
    for places to optimize
        
***********************************************************************************************/

#include "hdf5.h"
#include "hdf5_hl.h"

#include "decs.h"
#include "units.h"
#if(USEMPI==1)
#include "mpi.h"
#else
  wont-compile
#endif

#define TESTING  (0)
#define VERBOSE  (0) 
#define MAX_OPTICAL_DEPTH  (1.)
#define SMOOTH_FLUX (0)
#define FRAC_COARSEN_PHI (8)
#define FRAC_COARSEN_RAD (8)

#define PHI_WEDGE_LOOP    for( iwedge=0; iwedge < n_phi_wedges; iwedge++ )

// definitions needed for coolfunc table usage
#define NUM_A  (301)
#define NUM_B  (301)
#define NUM_C  (71)
#define CF_FILE "/mnt/c/scratch/sciteam/kinch/a0run/src-new6-4/coolfunc_table_ch.h5"

/* variable for the grid that holds the photospheric data: */
int  ip, kp, np, ifunc, iwedge; 
int  n_phi_wedges;
int  actual_frac_coarsen_tot, actual_frac_coarsen_rad, actual_frac_coarsen_phi;

double coarsen_weight, dphi_wedges; 
double  dxp1_photo, dxp3_photo;
double *cos_ph_photo, *sin_ph_photo;
double *cos_th_photo_bot, *sin_th_photo_bot, *cos_th_photo_top, *sin_th_photo_top;
double *cos_th_photo, *sin_th_photo;
double *data_photo;

extern double cooling_func_hr_disk( int i, int j, int k, double *ph );
extern int get_pid( int tmpcpupos[NDIM] ) ;

void  create_photosphere_grid(void);

/*********************************************************************************************/
/*********************************************************************************************
  setup_corona_cooling2()
  --------------------
   -- support routine for the corona (inverse Compton) cooling function; 
   -- one must call this routine before  cooling_func_corona() is called locally; 
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
void  setup_corona_cooling2( double ****prim, double Dt)
{

  enum equator_pos_type { Top_Eq, Over_Eq, Bottom_Eq } ; 
  enum equator_pos_type equator_pos; 

  int i,j,k,l;
  double dth, lumtmp;
  double th_top, th_bot, dA_cell;
  double *prim_loc;
  struct of_geom  *geom;
  struct of_coord *coords;
  struct of_state q ;

  int even_grid;
  char strout[200];

  /* Constants : */
  double  eta                = 0.06;
  double  mOBS               = 10.*C_mSUN;; 
  double  mdotNUM            = 3.e-4; 
  double  mdot_frac          = 1.e-2;
  double  gtheta             = 4.;
  double  coulomb_log        = 25.;
  double  n_elec_per_H       = 1.2; 

  double  mNUM               = 1.; 
  double  rOBS               = C_G*mOBS/(C_c2);; 
  double  rNUM               = mNUM; 
  double  L_edd              = C_4pi * C_G * mOBS * C_c / C_kappaT;
  double  mdot_edd           = L_edd / ( C_c2 * eta ); 
  double  mdotOBS            = mdot_frac * mdot_edd; 
  double  mdot2              = mdotOBS / (mdot_edd*eta);
  double  c1_corona          = ( 16. * C_pi * C_mp * n_elec_per_H * mdot2 /  (C_me * mdotNUM)  );
  double  c2_corona          = ( (8./3.)*C_mp*C_mp / (C_me*C_me*coulomb_log*gtheta) );
  double  c3_corona          = ( 4. * C_mp * 1.4 / (C_me * 2.3 ) );
  double  c4_corona          = 0.01691 * pow(mOBS/C_mSUN, -0.25) * pow(mdot2, 0.25) * pow(mdotNUM, -0.25);
  double  lengthScale        = rOBS / rNUM;     
  double  rhoScale           = mdotOBS / ( mdotNUM * lengthScale * lengthScale * (C_c) ); 
  double  tau_scale          = C_kappaT * rhoScale * lengthScale ;
  double  mp_o_me            = C_mp / C_me; 


  /* constant factor in expression for U_rad, in which we are integrating over the 
     surface of the photosphere.  See notes. */
  double  area_const1        = 0.5*pow( (0.75*dx[1]*dx[2]*dx[3]/M_PI), (2./3.)) / C_c; 

  TRACE_BEG;

  static int local_first_time = 1;
  
  if( local_first_time ) { 

    //--now in  alloc_grid_arrays() in decomp.c  -->    create_photosphere_grid();
    
    if( myid == printer_pid ) {
      fprintf(stdout,"%s():  eta             =  %26.16e \n",__func__, eta           ); 
      fprintf(stdout,"%s():  mdotNUM         =  %26.16e \n",__func__, mdotNUM       ); 
      fprintf(stdout,"%s():  mdot_frac       =  %26.16e \n",__func__, mdot_frac     ); 
      fprintf(stdout,"%s():  gtheta          =  %26.16e \n",__func__, gtheta        ); 
      fprintf(stdout,"%s():  n_elec_per_H    =  %26.16e \n",__func__, n_elec_per_H  ); 
      fprintf(stdout,"%s():  mOBS            =  %26.16e \n",__func__, mOBS          ); 
      fprintf(stdout,"%s():  mNUM            =  %26.16e \n",__func__, mNUM          ); 
      fprintf(stdout,"%s():  rOBS            =  %26.16e \n",__func__, rOBS          ); 
      fprintf(stdout,"%s():  rNUM            =  %26.16e \n",__func__, rNUM          ); 
      fprintf(stdout,"%s():  L_edd           =  %26.16e \n",__func__, L_edd         ); 
      fprintf(stdout,"%s():  mdot_edd        =  %26.16e \n",__func__, mdot_edd      ); 
      fprintf(stdout,"%s():  mdotOBS         =  %26.16e \n",__func__, mdotOBS       ); 
      fprintf(stdout,"%s():  mdot2           =  %26.16e \n",__func__, mdot2         ); 
      fprintf(stdout,"%s():  c1_corona       =  %26.16e \n",__func__, c1_corona     ); 
      fprintf(stdout,"%s():  c2_corona       =  %26.16e \n",__func__, c2_corona     ); 
      fprintf(stdout,"%s():  c3_corona       =  %26.16e \n",__func__, c3_corona     ); 
      fprintf(stdout,"%s():  lengthScale     =  %26.16e \n",__func__, lengthScale   ); 
      fprintf(stdout,"%s():  rhoScale        =  %26.16e \n",__func__, rhoScale      ); 
      fprintf(stdout,"%s():  tau_scale       =  %26.16e \n",__func__, tau_scale     ); 
      fprintf(stdout,"%s():  coarsen_weight  =  %26.16e \n",__func__, coarsen_weight); 
      fflush(stdout);
    }
    
    local_first_time = 0 ; 
  }

#if( USE_COOLING_FUNCTION != 5 )
  fprintf(stderr,"cooling_func_corona2(): Bad value USE_COOLING_FUNCTION  should not be herer  %d   !! \n",USE_COOLING_FUNCTION); 
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
#define MAX_N2CPUS (100)
  static MPI_Request req_out_tot1[MAX_N2CPUS]; 
  static MPI_Request req_out_tot2[MAX_N2CPUS]; 
  if( ncpux[2] > MAX_N2CPUS ) { 
    fprintf(stderr,"setup_corona_cooling(): need to increase size of req_out_tot1 and req_out_tot2 :  %d \n", ncpux[2]); 
    fflush(stderr);   fail(FAIL_BASIC,0);
  }

  int npts  = N1TOT*N3TOT;    //--OPTIMIZE  :  why do we have to pass the ghost zone data ??
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

#if(TRACE_CALLS||VERBOSE) 
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
    //    MPI_Irecv(  recv_buffer_tot, npts,  MPI_DOUBLE, collector_pid, tag_tot2, MPI_COMM_WORLD, &req_in_tot2);
  }
#if(TRACE_CALLS||VERBOSE) 
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
    prim_loc = prim[i][j][k];
    get_geometry(i,j,k,CENT,ncurr,geom); 
    get_state(prim_loc,geom,&q) ;

#if( TESTING )
    bound = 1.;
#else
    bound = (  (q.ucov[0]*( q.p + prim_loc[UU] + prim_loc[RHO] + q.bsq ))  >  (-prim_loc[RHO]) ) ?  1 : 0 ;
#endif

    //    dth = sqrt(geom->gcov[TH][TH]) * dx[TH] * lengthScale; 
    dth = sqrt(geom->gcov[TH][TH]) * dx[TH];
    dtau[i][j][k] = tau_scale * dth * prim_loc[RHO]; 

#if( TESTING )
    lumtmp = 0.;
    //    if( prim_loc[RHO] > 1. ) {     lumtmp = 1.; }
    get_coord(i,j,k,CENT,ncurr,coords);
    double r = coords->x[RR]; 
    double    x = sqrt(r);
    double x0 = sqrt(r_isco);
    double x1 = 2.*cos(  (acos(a) - M_PI) / 3. );
    double x2 = 2.*cos(  (acos(a) + M_PI) / 3. );
    double x3 = -2.*cos(  (acos(a) ) / 3. );
    
    double f = (1.5/( x*x*(x*x*x - 3.*x + 2.*a))) * ( 
						     x - x0 - 1.5*a*log(x/x0)  
						     - 3.*log((x-x1)/(x0-x1))*(x1-a)*(x1-a)/(x1*(x1-x2)*(x1-x3))
						  - 3.*log((x-x2)/(x0-x2))*(x2-a)*(x2-a)/(x2*(x2-x1)*(x2-x3))
						     - 3.*log((x-x3)/(x0-x3))*(x3-a)*(x3-a)/(x3*(x3-x1)*(x3-x2))
						  );
     double Rr = (2.*r*r/3.) * f;
     //double mdot = mdotOBS*C_c2/mdotNUM;
     double mdot = mdotNUM;
     double qnt = 3. * mdot * Rr / ( 4. * M_PI * r*r*r );
     if( prim_loc[RHO] > 1.e-6 ) {   lumtmp = qnt / (0.1*r); }
     
#else 
    lumtmp = cooling_func_hr_disk(i,j,k,prim_loc);
#endif
    
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
#if(TRACE_CALLS||VERBOSE) 
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
#if(TRACE_CALLS||VERBOSE) 
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
#if(TRACE_CALLS||VERBOSE) 
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
#if(TRACE_CALLS||VERBOSE) 
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
#if(TRACE_CALLS||VERBOSE) 
  trace_message("MPI-Isend-1",1);
#endif

  /***************************************************************************************************
    If we are the collector, then we have already collected neighboring data and integrated 
     to the equator.  Now we need to sum for total disk flux and then communicate this and the 
     location of the photosphere to all the receivers:
  ***************************************************************************************************/
  if( myid == collector_pid ) { 

    //    N1ALL_LOOP  N3ALL_LOOP { flux_disk[i][k] = flux_top[i][k] + flux_bot[i][k];  }

//    l = 0 ; 
//    N1ALL_LOOP  N3ALL_LOOP {  send_buffer_tot[l++] =  flux_disk[i][k];     }
    l = 0 ; 
    N1ALL_LOOP  N3ALL_LOOP { isend_buffer_tot[l++] = j_photo_top[i][k];    }
    N1ALL_LOOP  N3ALL_LOOP { isend_buffer_tot[l++] = j_photo_bot[i][k];    }


    /* Setup collector for receiving all integration results */
    ALLOC_ARRAY(pids_out, ncpux[2]); 

    /* make list of processes to send data to : */
    DLOOP1 { nbr_pos[i] = cpupos[i]; }
    for( j=0; j < ncpux[2]; j++ ) { 
      nbr_pos[2] = j; 
      pids_out[j] = get_pid(nbr_pos);
    }

#if(TRACE_CALLS||VERBOSE) 
    trace_message("MPI-Isend-collector-2",0);
#endif
    for( j=0; j < ncpux[2]; j++ )  if( myid != pids_out[j] )   { 
      MPI_Isend( isend_buffer_tot, npts2, MPI_INT   , pids_out[j], tag_tot1, MPI_COMM_WORLD, &req_out_tot1[j]);
      //      MPI_Isend(  send_buffer_tot, npts , MPI_DOUBLE, pids_out[j], tag_tot2, MPI_COMM_WORLD, &req_out_tot2[j]);
    }
#if(TRACE_CALLS||VERBOSE) 
    trace_message("MPI-Isend-collector-2",1);
#endif
  }
  else { 
#if(TRACE_CALLS||VERBOSE) 
    trace_message("MPI-Wait-2",0);
#endif
    MPI_Wait( &req_in_tot1, &mpi_status);
    l = 0 ; 
    N1ALL_LOOP  N3ALL_LOOP { j_photo_top[i][k] = irecv_buffer_tot[l++] ; }
    N1ALL_LOOP  N3ALL_LOOP { j_photo_bot[i][k] = irecv_buffer_tot[l++] ; }

//    MPI_Wait( &req_in_tot2, &mpi_status);
//    l = 0 ; 
//    N1ALL_LOOP  N3ALL_LOOP {  flux_disk[i][k] = recv_buffer_tot[l++];    }
#if(TRACE_CALLS||VERBOSE) 
    trace_message("MPI-Wait-2",1);
#endif
  }

#endif  /* USEMPI */


  /*******************************************************************
    If we are the collector, then we have already collected
     neighboring data and integrated to the equator.  Now we need to
     sum for total disk flux, interpolate/coarsen data to photosphere
     grid, and then global reduce the data (including photosphere
     height) while sharing it to all the processors:
  ********************************************************************/
  double xtmp[NDIM],xptmp[NDIM];
  
  PHOTODATA_LOOP {  data_photo[ip] = 0. ;       }
  PHOTODATA_LOOP {  final_data_photo[ip] = 0. ; }
  
  if( myid == collector_pid ) {

    /* Since we are summing to a global photospheric grid, we will
       ignore ghost cell values since those values reside on other
       domains: */
    
    N1_LOOP {
      ip = (int) ((i+globalpos[1]-N1S) / actual_frac_coarsen_rad);   /* radial location in coarsened photosphere grid */

      N3_LOOP {
	
	kp = (int) ((k+globalpos[3]-N3S) / actual_frac_coarsen_phi); /* azimuthal location in coarsened photosphere grid */
	np = ip*n3_photo + kp;

	get_coord(i,N2S,k,CENT,ncurr,coords);

	j = j_photo_bot[i][k];
	if( j >= 0 ) {
	  th_bot = startx[2] + (j + NG + 0.5)*dx[2];
	  xptmp[0] = coords->xp[0]; 
	  xptmp[1] = coords->xp[1]; 
	  xptmp[2] = coords->xp[2]; 
	  xptmp[3] = coords->xp[3]; 
	  xtmp[ 0] = coords->x[ 0]; 
	  xtmp[ 1] = coords->x[ 1]; 
	  xtmp[ 2] = coords->x[ 2]; 
	  xtmp[ 3] = coords->x[ 3];
	  xptmp[TH] = th_bot;
	  x_of_xp(xtmp, xptmp);
	  th_bot = xtmp[TH];
	}
	else {
	  th_bot = 0.5*M_PI;
	}

	j = j_photo_top[i][k];
	if( j >= 0 ) {
	  th_top = startx[2] + (j + NG + 0.5)*dx[2];
	  xptmp[0] = coords->xp[0]; 
	  xptmp[1] = coords->xp[1]; 
	  xptmp[2] = coords->xp[2]; 
	  xptmp[3] = coords->xp[3]; 
	  xtmp[ 0] = coords->x[ 0]; 
	  xtmp[ 1] = coords->x[ 1]; 
	  xtmp[ 2] = coords->x[ 2]; 
	  xtmp[ 3] = coords->x[ 3];
	  xptmp[TH] = th_top;
	  x_of_xp(xtmp, xptmp);
	  th_top = xtmp[TH]; 
	}
	else {
	  th_top = 0.5*M_PI;
	}

#if( 1 && TESTING  )
	fprintf(stdout,"flux-disk %26.16e %26.16e %26.16e %26.16e \n", coords->r, 0.5*(flux_top[i][k] + flux_bot[i][k])/(0.1*coords->r),tau_bot[i][k],tau_top[i][k]); fflush(stdout);
#endif
#if( 0 && TESTING  )
	fprintf(stdout,"flux-disk %5d %5d : %5d %5d  %26.16e %26.16e %26.16e %26.16e \n", i,k, j_photo_top[i][k], j_photo_bot[i][k], th_top, th_bot, coords->r, coords->x[PH]); fflush(stdout);
#endif
	dA_cell = coords->r * coords->dx_dxp[1][1] * dx[1] * dx[3];
	data_photo[N_PHOTO_FUNCS*np+PHOTO_FLUX  ] += dA_cell * 0.5*(flux_top[i][k] + flux_bot[i][k]);
	//	data_photo[N_PHOTO_FUNCS*np+PHOTO_FLUX  ] += dA_cell ;
	data_photo[N_PHOTO_FUNCS*np+PHOTO_TH_BOT] += dA_cell * th_bot;
	data_photo[N_PHOTO_FUNCS*np+PHOTO_TH_TOP] += dA_cell * th_top;

#if( 0 && TESTING  )
	fprintf(stdout,"flux-disk %5d %5d : %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", i,k, coords->r, coords->x[PH], dA_cell, coords->dx_dxp[1][1],  dx[1], dx[3] ); fflush(stdout);
#endif
	
      }
    }
    //-HERE HERE   -- we need to geometrically weight the fluxes  since the flux itself is  not geometrically weighted
    //    PHOTODATA_LOOP {  data_photo[ip] *= coarsen_weight ; }
  }
  
#if(TRACE_CALLS||VERBOSE) 
  trace_message("MPI_Allreduce",0);
#endif
 
#if(USEMPI)    
    MPI_Allreduce(data_photo, final_data_photo, n_data_photo, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    PHOTODATA_LOOP {  final_data_photo[ip] = data_photo[ip]; }
#endif

#if(TRACE_CALLS||VERBOSE) 
  trace_message("MPI_Allreduce",1);
#endif

  /*******************************************************************
    Calculate auxiliary data based on the photospheric grid data 
    collected so far. 
  ********************************************************************/
  np = 0 ;
  PHOTO_N1_LOOP PHOTO_N3_LOOP {
    th_bot = final_data_photo[N_PHOTO_FUNCS*np+PHOTO_TH_BOT] / dA_photo[ip];  /* Calculate area average of the poloidal location of photospheric cell */
    th_top = final_data_photo[N_PHOTO_FUNCS*np+PHOTO_TH_TOP] / dA_photo[ip];  /* Calculate area average of the poloidal location of photospheric cell */
    cos_th_photo_bot[np] = cos(th_bot); 
    sin_th_photo_bot[np] = sin(th_bot); 
    cos_th_photo_top[np] = cos(th_top); 
    sin_th_photo_top[np] = sin(th_top);
#if( TESTING )
    fprintf(stdout,"flux-disk-area-mpi %4d %4d : %12.5e %12.5e %12.5e \n", ip, kp, th_bot, th_top, dA_photo[ip]); fflush(stdout); 
#endif
    np++;
  }

  
  /*******************************************************************
    Now we have the photosphere data, we have to use it if we are in
    the corona, and use the local cooling function if we are in the
    disk:

    NOTE: we still need to use j_photo_[bot,top] to identify
          corona/disk cells because otherwise we would see cooling
          patterns in the disk flux on the coarsening interval (i.e.
          discrete jumps in the disk component every "FRAC_COARSEN_*"
          cell.
           
  ********************************************************************/

  /*****************************************************************************************
    Now use the coarsened data to calculate the local cooling
         function, either via coronal cooling or thermal;
  ******************************************************************************************/
  int jtmp = globalpos[2] - N2S; 
  int jglob;
  double pgas, rho, Te, urad, cool0, cool1, cool2, phi0, cos_phi0, sin_phi0;
  double r_loc, cos_ph_loc, sin_ph_loc, cos_th_loc, sin_th_loc, G_sign, cos_dphi, rdotr, dist_sq, dist, dOmega, G_pm, cos_Psi, rp_over_r, r_over_R, r_over_R_sq, inv_dist_sq;
  double kT_e, Theta_e; 

#if( TESTING )
  double xloc, yloc, zloc;
  double xphoto, yphoto, zphoto, dist_tmp;
#endif

  // coolfunc table read-in (should almost certainly be done somewhere else)
  hid_t file_id, dataset_id;
  herr_t status;

  double A_grid[NUM_A];
  double B_grid[NUM_B];
  double C_grid[NUM_C];
  double th_e_table[NUM_A*NUM_B*NUM_C];

  double mean_E;
  int n_bad;
  double a, b, dt;
  double ucon[NDIM];

#define N_READ_MAX_TRIES (100)

  file_id = H5Fopen(CF_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
  while (file_id < 0) { file_id = H5Fopen(CF_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);  }

  n_bad = 0; 
  while( (H5LTread_dataset( file_id, "/A_grid" ,  H5T_NATIVE_DOUBLE, A_grid ) < 0) &&  (n_bad < N_READ_MAX_TRIES) ) { n_bad++; }
  if( n_bad >= N_READ_MAX_TRIES ) { 
    fprintf(stdout,"%s() %d %d %d : line %d : too many read attempts \n",__func__,myid,nstep,n_substep, __LINE__); 
    fflush(stdout); fail(FAIL_HDF,0);
  }
    
  n_bad = 0; 
  while( (H5LTread_dataset( file_id, "/B_grid" ,  H5T_NATIVE_DOUBLE, B_grid ) < 0) &&  (n_bad < N_READ_MAX_TRIES) ) { n_bad++; }
  if( n_bad >= N_READ_MAX_TRIES ) { 
    fprintf(stdout,"%s() %d %d %d : line %d : too many read attempts \n",__func__,myid,nstep,n_substep, __LINE__); 
    fflush(stdout); fail(FAIL_HDF,0);
  }

  n_bad = 0; 
  while( (H5LTread_dataset( file_id, "/C_grid" ,  H5T_NATIVE_DOUBLE, C_grid ) < 0) &&  (n_bad < N_READ_MAX_TRIES) ) { n_bad++; }
  if( n_bad >= N_READ_MAX_TRIES ) { 
    fprintf(stdout,"%s() %d %d %d : line %d : too many read attempts \n",__func__,myid,nstep,n_substep, __LINE__); 
    fflush(stdout); fail(FAIL_HDF,0);
  }

  n_bad = 0; 
  while( (H5LTread_dataset( file_id, "/th_e" ,  H5T_NATIVE_DOUBLE, th_e_table ) < 0) &&  (n_bad < N_READ_MAX_TRIES) ) { n_bad++; }
  if( n_bad >= N_READ_MAX_TRIES ) { 
    fprintf(stdout,"%s() %d %d %d : line %d : too many read attempts \n",__func__,myid,nstep,n_substep, __LINE__); 
    fflush(stdout); fail(FAIL_HDF,0);
  }

  while (H5Fclose(file_id) < 0) { ; } 

#undef N_READ_MAX_TRIES 
    
    // end of coolfunc table read-in block

  l = 0 ; 
  LOOP { 
    urad = cool2 = 0.;
    mean_E = 0.0;
    lumtmp = lum_thermal[i][j][k];
    cool0 = cool1 = lumtmp;
    prim_loc = prim[i][j][k];
    kT_e = (gam-1.)*prim_loc[UU] / prim_loc[RHO];   /* for diagnostic output reasons  */

    jglob = j + jtmp; 

    /* We use coronal cooling function outside the photosphere : */
    if( (j_photo_top[i][k] < 0) || (jglob < j_photo_top[i][k]) || (jglob > j_photo_bot[i][k]) ) {

      /* Integrate thermal emission over the photosphere in the appropriate hemisphere : */
      get_coord(i,j,k,CENT,ncurr,coords);
      cos_ph_loc = cos(coords->x[PH]);   /* maybe optimize this? */
      sin_ph_loc = sin(coords->x[PH]);
      cos_th_loc = cos(coords->x[TH]); 
      sin_th_loc = sin(coords->x[TH]); 
      r_loc = coords->r;


      if( jglob < totalsize[2]/2 )  {  /* top */
	cos_th_photo = cos_th_photo_top;
	sin_th_photo = sin_th_photo_top;
	G_sign = -1.;
#if( 0 && TESTING )
	fprintf(stdout,"angles-top %5d %5d %5d - %5d :  %12.5e %12.5e %12.5e \n", i,j,k,jglob,r_loc, coords->x[TH], coords->x[PH] ); fflush(stdout);
#endif
      }
      else {  /* bottom */
	cos_th_photo = cos_th_photo_bot;
	sin_th_photo = sin_th_photo_bot;
	G_sign = 1.;
#if( 0 && TESTING )
	fprintf(stdout,"angles-bot %5d %5d %5d - %5d :  %12.5e %12.5e %12.5e \n", i,j,k,jglob,r_loc, coords->x[TH], coords->x[PH] ); fflush(stdout);
#endif
      }
      
      PHI_WEDGE_LOOP { 
	np = 0;
	phi0 = iwedge * dphi_wedges;
	cos_phi0 = cos(phi0);  //--OPTIMIZE
	sin_phi0 = sin(phi0);  //--OPTIMIZE

	PHOTO_N1_LOOP PHOTO_N3_LOOP {
	  //	  cos_dphi = ( cos_ph_loc * cos_ph_photo[kp] + sin_ph_loc * sin_ph_photo[kp] ) ; 
	  cos_dphi = ( cos_ph_loc * (cos_ph_photo[kp] * cos_phi0  - sin_ph_photo[kp] * sin_phi0)   +
		       sin_ph_loc * (sin_ph_photo[kp] * cos_phi0  + cos_ph_photo[kp] * sin_phi0) ) ; //--ok

	  cos_Psi =  sin_th_loc * sin_th_photo[np] * cos_dphi +  cos_th_loc * cos_th_photo[np] ;
	  rdotr = r_loc * r_photo[ip] * cos_Psi; //--ok
	  rp_over_r = r_photo[ip] / r_loc; 
	  r_over_R_sq = 1./( 1. + rp_over_r * ( rp_over_r - 2.*cos_Psi ) );
	  r_over_R = sqrt(r_over_R_sq);
	  
	  //	  dist_sq = r_loc * r_loc  +  r_photo[ip] * r_photo[ip] - 2.*rdotr ;    //--ok
	  //	  dist = sqrt(dist_sq);

	  inv_dist_sq = r_over_R_sq/(r_loc*r_loc);
	  
	  G_pm  =  sin_th_loc * cos_th_photo[np] * cos_dphi  -  cos_th_loc * sin_th_photo[np] ;
	  G_pm  *= G_sign ;
	  if( G_pm < 0. ) { G_pm = 0.; }
	  //	  urad += G_pm * final_data_photo[N_PHOTO_FUNCS*np + PHOTO_FLUX] * pow(dist_sq, -1.5);  /* area factor already included in the MPI integration  */

	  //	  dOmega = 2.*M_PI * ( 1. - 1./sqrt( 1. + dA_photo[ip]/(M_PI*dist_sq)  ) ); 
	  dOmega = 2.*M_PI * ( 1. - 1./sqrt( 1. + dA_photo[ip]*inv_dist_sq/(M_PI)  ) ); 
	  //	  urad += G_pm * final_data_photo[N_PHOTO_FUNCS*np + PHOTO_FLUX] * dOmega / (dist*dA_photo[ip]) ;  /* area factor already included in the MPI integration  */

	  
	  urad   += G_pm * final_data_photo[N_PHOTO_FUNCS*np + PHOTO_FLUX] * dOmega * r_over_R / (dA_photo[ip]) ;  /* area factor already included in the MPI integration  */
      mean_E += pow(final_data_photo[N_PHOTO_FUNCS*np + PHOTO_FLUX], 1.25) * (G_pm * dOmega * r_over_R / (dA_photo[ip]));

#if( 0 &&  TESTING )
	  if( (i == 50) && (jglob==0 || jglob==(totalsize[2]-1)) ) { 
	    fprintf(stdout,"photo-angle-special %5d %5d %5d  %12.5e %12.5e %12.5e %12.5e %12.5e %6d \n", i, j, k, coords->r, coords->x[TH], coords->x[PH], urad, dist_sq, ip); fflush(stdout);
	  }
#endif
	  
#if( 0 &&  TESTING )
	  fprintf(stdout,"photo %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n", r_photo[ip], ph_photo[kp], phi0, dA_photo[ip], final_data_photo[N_PHOTO_FUNCS*np + PHOTO_FLUX], final_data_photo[N_PHOTO_FUNCS*np + PHOTO_FLUX]/(dA_photo[ip]*r_photo[ip]*0.1)); fflush(stdout);
#endif
#if( 0 && TESTING )
	  xphoto = r_photo[ip] * sin_th_photo[np] * cos(phi0 + ph_photo[kp]);
	  yphoto = r_photo[ip] * sin_th_photo[np] * sin(phi0 + ph_photo[kp]);
	  zphoto = r_photo[ip] * cos_th_photo[np];
	  xloc   = r_loc * sin_th_loc * cos_ph_loc ;
	  yloc   = r_loc * sin_th_loc * sin_ph_loc ;
	  zloc   = r_loc * cos_th_loc              ;
	  dist_tmp =
	    (xphoto-xloc)*(xphoto-xloc) +
	    (yphoto-yloc)*(yphoto-yloc) + 
	    (zphoto-zloc)*(zphoto-zloc) ;

	  double rd_tmp1 = REL_DIFF_FUNC(dist_sq, dist_tmp);
	  if( rd_tmp1 > 1e-11 ) { 
	    fprintf(stdout,"distances sq = %26.16e  %26.16e    rd=[ %26.16e ] \n", dist_sq, dist_tmp,rd_tmp1 );
	  }
	  
	  fflush(stdout); 
#endif
	  
	  np++;
	}
      }
      
      //      urad *= coords->r ;

#if( 0 &&  TESTING )
      if( i == 50 ) { 
	fprintf(stdout,"photo-angle %5d %5d %5d  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %6d \n", i, j, k, coords->r, coords->x[TH], coords->x[PH], urad, G_pm, G_sign,np); fflush(stdout);
      }
#endif

      if( urad > 0. ) {
	prim_loc = prim[i][j][k];
	rho      = prim_loc[RHO]; 
	pgas     = (gam-1.)*prim_loc[UU];
	
	//	printf("%d %d\n", A_ndx, B_ndx);

//	Theta_e = (1.0/((A_grid[A_ndx+1] - A_grid[A_ndx])*(B_grid[B_ndx+1] - B_grid[B_ndx]))) * (th_e_table[A_ndx*NUM_B + B_ndx]*(A_grid[A_ndx+1] - A)*(B_grid[B_ndx+1] - B) + th_e_table[(A_ndx+1)*NUM_B + B_ndx]*(A - A_grid[A_ndx])*(B_grid[B_ndx+1] - B) + th_e_table[A_ndx*NUM_B + B_ndx+1]*(A_grid[A_ndx+1] - A)*(B - B_grid[B_ndx]) + th_e_table[(A_ndx+1)*NUM_B + B_ndx+1]*(A - A_grid[A_ndx])*(B - B_grid[B_ndx]));

    Theta_e = trilin_interp(pgas/rho, urad/rho, 3.832 * ((c4_corona * mean_E)/urad), A_grid, B_grid, C_grid, th_e_table);

//  Theta_e = (1.0/(1.0 + n_elec_per_H)) * A * (C_mp/C_me);
	
	//	printf("%e\n", 511.0 * Theta_e);
	kT_e    = (C_me/C_mp) * Theta_e;

	// "k T_e" of electrons in code units, or in other words it is  "p_e_code / rho_code"  where p_e_code is the electron partial pressure in code units
	//          kT_e    =  pgas / ( 2.*rho + c2_corona*urad ) ;
	//          Theta_e  = mp_o_me * kT_e;

	//--orig cool0 = cool2 = c1_corona * rho * Te * urad * ( 1. + c3_corona*Te );
//  cool0 = cool2 = c1_corona * rho * kT_e * urad * ( 1. + 4.*Theta_e );
    cool0 = cool2 = c1_corona * rho * kT_e * urad * ( 1. + 4.*Theta_e ) - 0.25 * (C_me/C_mp) * c1_corona * rho * urad * 3.832 * ((c4_corona * mean_E)/urad);

    /*
    get_geometry(i,j,k,CENT,ncurr,geom);
    ucon_calc(prim_loc, geom, ucon);
    dt = Dt/ucon[0];
    a  = (1.0/(1.0 + n_elec_per_H)) * c1_corona * (gam - 1.0) * urad;
    b  = ((4.0/(1.0 + n_elec_per_H)) * (gam - 1.0) * (C_mp/C_me) * prim_loc[UU])/rho;
    cool0 = cool2 = (prim_loc[UU]/dt) * (1.0 - 1.0/((1.0 + b)*exp(a*dt) - b));
    */
      }
    }

    //  gtheta = 4.
    //  const = 2.666666*((1837.)^2)/(40.*gtheta)
    // tetp = 1./(1. + const*urad/rho)
    // tpte = 1/tetp = (1. + const*urad/rho)
    // c2_corona  =      = ( C_mp*C_mp / (C_me*C_me*15.*gtheta) );
    // 2.6666  =  2 + 2/3 =  (6+2)/3 = 8/3
    //  8/3/40 = 1/(3*5) = 1/15
    //
    // 1837 = mp/me
    //
    // c2_corona is ok ,  = const = c2_corona
    //

    // tempion = ((gam-1)*uu/rho)/(1. + tetp)
    // tempion = ((gam-1)*uu/rho)/(1. + 1/tpte)
    // tempelec = tetp*tempion
    // tempelec = tetp*((gam-1)*uu/rho)/(1. + 1/tpte)
    // tempelec = (1/tpte)*((gam-1)*uu/rho)/(1. + 1/tpte)
    // tempelec = ((gam-1)*uu/rho)/(tpte + 1)
    // tempelec = ((gam-1)*uu)/(rho*tpte + rho)
    // tempelec = ((gam-1)*uu)/(const*urad + rho  + rho)
    // tempelec = ((gam-1)*uu)/(const*urad + 2*rho)
    //
    // tempelec = Te 

    //        theta = tempelec/511.
    //        pe = 0.666*uu
    //         icconst = 16.*!PI*(mdot/mdotcode)*1837.*1.2
    //   coolcor = icconst*pe*urad*(1. + 4*theta)
    //
    //   c1_corona = iccons  is OK!!
    // 
    //   pe  ?=  rho * Te
    //   urad = urad  (by assumption)
    //   4 * theta  ?=   c3_corona * Te  = 4 * (1.4*1837/2.3) Te
    //    theta ?= (1.4*1837/2.3) Te
    //   tempelec/511  ?= 1118 Te
    // !!!!!!! HERE 
    //   
    //   c1_corona  = ( 16. * C_pi * C_mp * n_elec_per_H * mdot2 /  (C_me * mdotNUM)  );
    //    c3_corona =  ( 4. * C_mp * 1.4 / (C_me * 2.3 ) ); 

    

    //  p V = N k T  ->  p =  n k T = rho k T / mp
    // k T =  (p/rho) * mp

      //  pCGS    = EOS1(uCGS) ;         /* pCGS is the electrons' contribution to the pressure */
      //  *TeCGS  = pCGS/(nrhoCGS*C_k) ;                  /* TeCGS is the electron temperature */
      //   nrhoCGS = rho       * units[RTU_nrhoScale] ; 
    //   rhoScale  = rhoScale 
    //    nrhoScale = rhoScale / ((C_me) + (C_mp)) ;
    //   pressureScale = energyDensityScale
      //     units[RTU_energyDensityScale] = units[RTU_rhoScale] * (C_c2) ;
      //
      // TeCGS = p * energyDensityScale / ( rho * nrhoScale * C_k)
      // TeCGS = p * rhoScale * C_c2 / ( rho * rhoScale * C_k/(mp+me))
      // TeCGS = p * C_c2 (mp+me) / ( rho * C_k)
      // TeCGS = (p /  rho)    * C_c2 (mp+me) / C_k

    double Te_eV = kT_e * C_mp * C_c2 / C_eV; 
#if( TESTING )    
    fprintf(stdout,"Te factor = %26.16e  %26.16e %26.16e %26.16e %26.16e %26.16e \n", (C_mp * C_c2 / C_eV), kT_e, Te_eV, pgas/rho, urad); fflush(stdout); 
#endif
    
    //    coolflux[0][l] = cool0;  /* final  */
    //    coolflux[1][l] = cool1;  /* disk   */
    //    coolflux[2][l] = cool2;  /* corona */
    //    coolflux[3][l] = urad;  /* corona */
    coolflux[0][l] = cool0;  /* final  */
    coolflux[1][l] = cool1;  /* disk   */
    coolflux[2][l] = cool2;  /* corona */
    coolflux[3][l] = urad;  /* corona */
#if( USE_COOLING_FUNCTION == 5 )
    radfunc_Te_eV[l] = Te_eV;
#endif
    
    l++;
  }
  
#if(0)  
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
#endif
    

#if( USEMPI )
  /* Wait for all outstanding MPI communications : */
#if(TRACE_CALLS||VERBOSE) 
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
      //      MPI_Wait( &req_out_tot2[j], &mpi_status);
    }
    DEALLOC_ARRAY(pids_out, ncpux[2]); 
  }
#if(TRACE_CALLS||VERBOSE) 
  trace_message("MPI-Wait-3",1);
#endif
#endif


  TRACE_END;

  return;
}

/*********************************************************************************************/
/*********************************************************************************************
  create_photosphere_grid():
  --------------------

   -- allocates and assigns the coordinate array that represents where
      we will store photospheric data;

   -- we use an even coarsening of the global domain so that we can conserve flux; 

   -- not only is flux_disk[] averaged within each cell, but the
      height of the top and bottom locations too;

**********************************************************************************************/
void  create_photosphere_grid(void)
{
  char strout[200];
  double ftmp;
  
  TRACE_BEG;

#if( (COORD_TYPE_CHOICE != COORD_IDENTITY ) && \
     (COORD_TYPE_CHOICE != COORD_DIAGONAL ) && \
     (COORD_TYPE_CHOICE != COORD_DIAGONAL2) && \
     (COORD_TYPE_CHOICE != COORD_DIAGONAL3)    )

  ERROR_MESSAGE("COORD_TYPE_CHOICE IS INVALID");
  fail(FAIL_BASIC,0);
#endif
  

  if( (totalsize[1] % FRAC_COARSEN_RAD) != 0 ) {
    sprintf(strout,"total radial cell count must be divisible by the Coarsening factor: %d %d ",totalsize[1], FRAC_COARSEN_RAD); 
    ERROR_MESSAGE(strout);
    fail(FAIL_BASIC,0);
  }
  
  if( (totalsize[3] % FRAC_COARSEN_PHI) != 0  && (totalsize[3] > 1) ) {
    sprintf(strout,"total azimuthal cell count must be divisible by the Coarsening factor: %d %d ",totalsize[3], FRAC_COARSEN_PHI); 
    ERROR_MESSAGE(strout);
    fail(FAIL_BASIC,0);
  }

  /***************************************************************************************************
    Figure out if the simulation covers the full azimuthal extent or not: 
       -- this affects how the phi-integration happens over the photosphere;
       -- its easiest to think of two azimuthal loops:
            Loop1 (n3_photo)     : over a coarsened or exact representation of the phi-grid;
            Loop2 (n_phi_wedges) : over phi wedges, which may be equal to 2*pi or smaller in extent. 
  ***************************************************************************************************/
  if( IS_DIFFERENT_DOUBLE( GridLength[3] , 2.*M_PI ) ) {

    /* If we are using a faction of the full extent: */
    n_phi_wedges = (int) (2.*M_PI / GridLength[3]  + 0.5);
    ftmp = n_phi_wedges * GridLength[3]; 
    if( IS_DIFFERENT_DOUBLE( ftmp , 2.*M_PI ) ) {
      sprintf(strout,"The azimuthal wedge needs to be an integer divisor of 2 pi : %26.16e  ", GridLength[3]); 
      ERROR_MESSAGE(strout);
      fail(FAIL_BASIC,0);
    }
  }
  else {
    n_phi_wedges = 1 ; 
  }

  dphi_wedges = 2.*M_PI / n_phi_wedges;

  n1_photo = totalsize[1] / FRAC_COARSEN_RAD;
  n3_photo = totalsize[3] / FRAC_COARSEN_PHI;

  if( n1_photo <= 0 ) { n1_photo = 1; }
  if( n3_photo <= 0 ) { n3_photo = 1; }

  if( (totalsize[1] % n1_photo) != 0 ) {
    sprintf(strout,"total radial cell count must be divisible by photo cell count: %d %d ",totalsize[3], n1_photo);
    ERROR_MESSAGE(strout);
    fail(FAIL_BASIC,0);
  }
  
  if( (totalsize[3] % n3_photo) != 0  ) {
    sprintf(strout,"total azimuthal cell count must be divisible by photo cell count: %d %d ",totalsize[3], n3_photo);
    ERROR_MESSAGE(strout);
    fail(FAIL_BASIC,0);
  }

  fprintf(stderr, "totalsize[1] = %d, n1_photo = %d\n", totalsize[1], n1_photo); fflush(stderr);
  actual_frac_coarsen_rad = totalsize[1] / n1_photo;
  actual_frac_coarsen_phi = totalsize[3] / n3_photo;
  actual_frac_coarsen_tot = actual_frac_coarsen_rad * actual_frac_coarsen_phi;

  coarsen_weight =  1./( (double) (actual_frac_coarsen_tot) );
  
  n_photo = n1_photo * n3_photo;
  n_data_photo = n_photo * N_PHOTO_FUNCS;

  dxp1_photo = dx[1] * actual_frac_coarsen_rad;
  dxp3_photo = dx[3] * actual_frac_coarsen_phi;


  if( myid == printer_pid ) {
    fprintf(stdout,"%s(): n_phi_wedges =  %3d      \n", __func__, n_phi_wedges );
    fprintf(stdout,"%s(): dphi_wedges  =  %26.16e  \n", __func__, dphi_wedges );
    fprintf(stdout,"%s(): n1_photo     =  %3d      \n", __func__, n1_photo    );
    fprintf(stdout,"%s(): n3_photo     =  %3d      \n", __func__, n3_photo );
    fprintf(stdout,"%s(): dxp1_photo   =  %26.16e  \n", __func__, dxp1_photo);
    fprintf(stdout,"%s(): dxp3_photo   =  %26.16e  \n", __func__, dxp3_photo );
    fflush(stdout);
  }


  //  double darea_num = dxp1_photo * dxp3_photo / n_phi_wedges; 
  double darea_num = dxp1_photo * dxp3_photo;

  ALLOC_ARRAY(xp1_photo,n1_photo);
  ALLOC_ARRAY(xp3_photo,n3_photo);
  ALLOC_ARRAY( r_photo,n1_photo);
  ph_photo = xp3_photo;    /* they are the same the assumed coordinates */
  ALLOC_ARRAY(dA_photo,n1_photo);
  ALLOC_ARRAY(cos_ph_photo,n3_photo);
  ALLOC_ARRAY(sin_ph_photo,n3_photo);
  ALLOC_ARRAY(cos_th_photo_top,n_photo);
  ALLOC_ARRAY(sin_th_photo_top,n_photo);
  ALLOC_ARRAY(cos_th_photo_bot,n_photo);
  ALLOC_ARRAY(sin_th_photo_bot,n_photo);
  ALLOC_ARRAY(data_photo,n_data_photo);
  ALLOC_ARRAY(final_data_photo,n_data_photo);

  PHOTO_N1_LOOP { xp1_photo[ip]    = startx[1] + (NG)*dx[1] + (0.5 + ((double) ip))*dxp1_photo; }
  PHOTO_N3_LOOP { xp3_photo[kp]    = startx[3] + (NG)*dx[3] + (0.5 + ((double) kp))*dxp3_photo; }

  PHOTO_N1_LOOP { r_photo[ip]      = R0 + exp(xp1_photo[ip]);   }

  /*  dA_photo =  r * (r-R0) * dxp1_photo * dxp3_photo / n_phi_wedges; 
               =  r *  (dr/dxp1) * dxp1  * (2*pi/(N3*n_phi_wedges))
   */
  PHOTO_N1_LOOP { dA_photo[ip]     = r_photo[ip]  * (r_photo[ip] - R0) * darea_num; }
  PHOTO_N3_LOOP { cos_ph_photo[kp] = cos(ph_photo[kp]);     }
  PHOTO_N3_LOOP { sin_ph_photo[kp] = sin(ph_photo[kp]);     }

#if( VERBOSE )
  double area_out, area_sol; 
  if( myid == printer_pid ) {
    PHOTO_N1_LOOP {  fprintf(stdout,"photo-r  %10d  %26.16e  %26.16e  %26.16e \n", ip, xp1_photo[ip], r_photo[ip], dA_photo[ip]);   }
    fflush(stdout); 

    PHOTO_N3_LOOP {  fprintf(stdout,"photo-phi %10d  %26.16e  %26.16e  %26.16e \n", kp, xp3_photo[kp], cos_ph_photo[kp], sin_ph_photo[kp]);   }
    fflush(stdout);

    area_out = 0.;
    PHOTO_N1_LOOP PHOTO_N3_LOOP PHI_WEDGE_LOOP {
      area_out += dA_photo[ip]; 
    }
    area_sol = M_PI*(Rout*Rout-6.*6.);
    fprintf(stdout,"area-test = %26.16e  %26.16e   %26.16e  \n", area_out, area_sol, REL_DIFF_FUNC(area_out,area_sol)); 
    fflush(stdout);

  }
#endif

  /*   At the z-axis,    \int  r dr / ( r^2 + 2*a*r + c )   =  0.5 * log( 2*a*r + c + r^2) - a arctan( (a+r)/sqrt(c-a^2) ) / sqrt(c-a^2)  

       -- here  
                c = z^2 
                a = z*cos(theta) = z*cos(pi/2 +/- h_o_r) =  z*sin(+/- h_o_r) =  +/- z*sin(h_o_r) 
   */ 

  
  TRACE_END;

  return;
}


/*********************************************************************************************/
/*********************************************************************************************
  free_photosphere_arrays():
  --------------------
   -- deallocates all the arrays used by the corona2 routines;
**********************************************************************************************/
void  free_photosphere_arrays(void)
{
  TRACE_BEG;

  DEALLOC_ARRAY(xp1_photo,n1_photo);
  DEALLOC_ARRAY(xp3_photo,n3_photo);
  DEALLOC_ARRAY( r_photo,n1_photo);
  //  ph_photo = xp3_photo;    /* they are the same the assumed coordinates */
  DEALLOC_ARRAY(dA_photo,n1_photo);
  DEALLOC_ARRAY(cos_ph_photo,n3_photo);
  DEALLOC_ARRAY(sin_ph_photo,n3_photo);
  DEALLOC_ARRAY(cos_th_photo_top,n_photo);
  DEALLOC_ARRAY(sin_th_photo_top,n_photo);
  DEALLOC_ARRAY(cos_th_photo_bot,n_photo);
  DEALLOC_ARRAY(sin_th_photo_bot,n_photo);
  DEALLOC_ARRAY(data_photo,n_data_photo);
  DEALLOC_ARRAY(final_data_photo,n_data_photo);
  
  TRACE_END;
  
  return;
}

double trilin_interp(double A, double B, double C, double *A_grid, double *B_grid, double *C_grid, double *th_e_table) {
    int A_ndx, B_ndx, C_ndx;
    double x0, x1, y0, y1, z0, z1, c000, c001, c010, c011, c100, c101, c110, c111, a0, a1, a2, a3, a4, a5, a6, a7;

	A_ndx = (int)(log(A/A_grid[0])/log(A_grid[1]/A_grid[0])); 
	B_ndx = (int)(log(B/B_grid[0])/log(B_grid[1]/B_grid[0]));
	C_ndx = (int)(log(C/C_grid[0])/log(C_grid[1]/C_grid[0]));
	
	if (A_ndx < 0) {
	  A_ndx = 0;
	  printf("A_ndx < 0\n");
	}
	else if (A_ndx >= NUM_A-2) {
	  A_ndx = NUM_A-2;
	  printf("A_ndx > NUM_A-2: %e\n", log10(A));
	}
	if (B_ndx < 0) {
	  B_ndx = 0;
	  printf("B_ndx < 0\n");
	}
	else if (B_ndx >= NUM_B-2) {
	  B_ndx = NUM_B-2;
	  printf("B_ndx > NUM_B-2\n");
	}
	if (C_ndx < 0) {
	  C_ndx = 0;
	  printf("C_ndx < 0\n");
	}
	else if (C_ndx >= NUM_C-2) {
	  C_ndx = NUM_C-2;
	  printf("C_ndx > NUM_C-2\n");
	}

    x0 = A_grid[A_ndx];
    x1 = A_grid[A_ndx+1];
    y0 = B_grid[B_ndx];
    y1 = B_grid[B_ndx+1];
    z0 = C_grid[C_ndx];
    z1 = C_grid[C_ndx+1];

    c000 = th_e_table[(A_ndx)*(NUM_B*NUM_C)   + (B_ndx)*NUM_C   + (C_ndx)];
    c001 = th_e_table[(A_ndx)*(NUM_B*NUM_C)   + (B_ndx)*NUM_C   + (C_ndx+1)];
    c010 = th_e_table[(A_ndx)*(NUM_B*NUM_C)   + (B_ndx+1)*NUM_C + (C_ndx)];
    c011 = th_e_table[(A_ndx)*(NUM_B*NUM_C)   + (B_ndx+1)*NUM_C + (C_ndx+1)];
    c100 = th_e_table[(A_ndx+1)*(NUM_B*NUM_C) + (B_ndx)*NUM_C   + (C_ndx)];
    c101 = th_e_table[(A_ndx+1)*(NUM_B*NUM_C) + (B_ndx)*NUM_C   + (C_ndx+1)];
    c110 = th_e_table[(A_ndx+1)*(NUM_B*NUM_C) + (B_ndx+1)*NUM_C + (C_ndx)];
    c111 = th_e_table[(A_ndx+1)*(NUM_B*NUM_C) + (B_ndx+1)*NUM_C + (C_ndx+1)];

    a0 = -(c000 * x1 * y1 * z1) + (c001 * x1 * y1 * z0) + (c010 * x1 * y0 * z1) - (c011 * x1 * y0 * z0) + (c100 * x0 * y1 * z1) - (c101 * x0 * y1 * z0) - (c110 * x0 * y0 * z1) + (c111 * x0 * y0 * z0);
    a1 = (c000 * y1 * z1) - (c001 * y1 * z0) - (c010 * y0 * z1) + (c011 * y0 * z0) - (c100 * y1 * z1) + (c101 * y1 * z0) + (c110 * y0 * z1) - (c111 * y0 * z0);
    a2 = (c000 * x1 * z1) - (c001 * x1 * z0) - (c010 * x1 * z1) + (c011 * x1 * z0) - (c100 * x0 * z1) + (c101 * x0 * z0) + (c110 * x0 * z1) - (c111 * x0 * z0);
    a3 = (c000 * x1 * y1) - (c001 * x1 * y1) - (c010 * x1 * y0) + (c011 * x1 * y0) - (c100 * x0 * y1) + (c101 * x0 * y1) + (c110 * x0 * y0) - (c111 * x0 * y0);
    a4 = -(c000 * z1) + (c001 * z0) + (c010 * z1) - (c011 * z0) + (c100 * z1) - (c101 * z0) - (c110 * z1) + (c111 * z0);
    a5 = -(c000 * y1) + (c001 * y1) + (c010 * y0) - (c011 * y0) + (c100 * y1) - (c101 * y1) - (c110 * y0) + (c111 * y0);
    a6 = -(c000 * x1) + (c001 * x1) + (c010 * x1) - (c011 * x1) + (c100 * x0) - (c101 * x0) - (c110 * x0) + (c111 * x0);
    a7 = c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111;

    return (a0 + (a1 * A) + (a2 * B) + (a3 * C) + (a4 * A * B) + (a5 * A * C) + (a6 * B * C) + (a7 * A * B * C))/((x0 - x1)*(y0 - y1)*(z0 - z1));
}

