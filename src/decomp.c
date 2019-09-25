/*************************************************************************************/
/******   DOMAIN DECOMPOSITION SUBROUTINES   *****************************************/
/*************************************************************************************/

#include "decs.h"
#include "harm_mpi.h"   
#include "recon.h"


/**********************************************************************************************/
/**********************************************************************************************
  set_global_domain_variables():
 ------------------
   -- sets global variables describing the shape/extent of the domain; 
   -- this should be called early in the process as we have to allocate global arrays with 
       these variables;
   -- these should already be set:  N1_glob, N2_glob, N3_glob;

**********************************************************************************************/
void set_global_domain_variables(void) 
{

  
  N1TOT_glob =  N1_glob+2*NG;
  N2TOT_glob =  N2_glob+2*NG;
  N3TOT_glob =  N3_glob+2*NG;
  NTOT_glob  =  N1TOT_glob*N2TOT_glob*N3TOT_glob;
  N1E_glob   =  N1_glob+NG-1;
  N2E_glob   =  N2_glob+NG-1;
  N3E_glob   =  N3_glob+NG-1;

  /* Determine the max. dimension of the data arrays : */
  if( N1TOT_glob >= N2TOT_glob ) { 
    if( N1TOT_glob >= N3TOT_glob ) {  MAX_NTOT_glob = N1TOT_glob;   }
    else                           {  MAX_NTOT_glob = N3TOT_glob;   }
  }
  else { 
    if( N2TOT_glob >= N3TOT_glob ) {  MAX_NTOT_glob = N2TOT_glob;   }
    else                           {  MAX_NTOT_glob = N3TOT_glob;   }
  }


#if( MAKE_HISTORY && MAKE_HDF5 ) 
  N_HIST_POINTS_glob = N_HISTORY*N_HIST_TYPES*N1_glob;
  N_SURF_POINTS_glob = N_SURFACE*N_SURF_TYPES*N1_glob*N3_glob;
#else 
  N_HIST_POINTS_glob = N_HISTORY*N_HIST_TYPES*n_r_bins;
  N_SURF_POINTS_glob = N_SURFACE*N_SURF_TYPES*n_r_bins*n_phi_bins;
#endif

  N_COORD_glob = N1TOT_glob*N2TOT_glob*N3TOT_glob*NPOS;
  
  N1_CONN_glob = N1_glob+2*NDC;
  N2_CONN_glob = N2_glob+2*NDC;
  N3_CONN_glob = N3_glob+2*NDC;

#if(METRIC_DIM==3) 
  N_GEOM_glob  =  N1TOT_glob*N2TOT_glob*N3TOT_glob*NPOS; 
  N_CONN_glob  =  N1_CONN_glob*N2_CONN_glob*N3_CONN_glob;
#elif(METRIC_DIM==2) 
  N_GEOM_glob  =  N1TOT_glob*N2TOT_glob*NPOS;
  N_CONN_glob  =  N1_CONN_glob*N2_CONN_glob;
#elif(METRIC_DIM==1) 
  N_GEOM_glob  =  N1TOT_glob*NPOS;
  N_CONN_glob  =  N1_CONN_glob;
#else
  N_GEOM_glob  = 1;
  N_CONN_glob  = 1;
#endif

  CONN2_END1_glob = N1E_glob+NDC;
  CONN2_END2_glob = N2E_glob+NDC;
  CONN2_END3_glob = N3E_glob+NDC;

  NM1_glob =  N1_glob;
  NM2_glob =  N2_glob;
  NM3_glob =  N3_glob;

  NM1E_glob = NM1_glob + NM1_BND - 1;
  NM2E_glob = NM2_glob + NM2_BND - 1;
  NM3E_glob = NM3_glob + NM3_BND - 1;

  NM1_TOT_glob = NM1_glob + 2*NM1_BND;
  NM2_TOT_glob = NM2_glob + 2*NM2_BND;
  NM3_TOT_glob = NM3_glob + 2*NM3_BND;

  DEGENERATE_glob = 0;

  if( N1_glob == 1 )  { 
    N1_R_S_glob = N1S; 
    N1_R_E_glob = N1S; 
    DEGENERATE_glob = 1;
  }
  else {
    N1_R_S_glob = N1S-BUF; 
    N1_R_E_glob = N1E+BUF; 
  }

  if( N2 == 1 ) {
    N2_R_S_glob = N2S;
    N2_R_E_glob = N2S;
    DEGENERATE_glob = 1;
  }
  else { 
    N2_R_S_glob = N2S-BUF;
    N2_R_E_glob = N2E+BUF;
  }

  if( N3 == 1 ) {
    N3_R_S_glob  = N3S;
    N3_R_E_glob  = N3S;
    DEGENERATE_glob = 1;
  }
  else { 
    N3_R_S_glob = N3S-BUF;
    N3_R_E_glob = N3E+BUF;
  }

  NM1_UP_BEG_glob = NM1E-NM1_BND+1;
  NM2_UP_BEG_glob = NM2E-NM2_BND+1;
  NM3_UP_BEG_glob = NM3E-NM3_BND+1;

  do_recon_dir[0] = (N1!=1);
  do_recon_dir[1] = (N2!=1);
  do_recon_dir[2] = (N3!=1);

  NCELLS_glob = N1_glob * N2_glob * N3_glob;

  return;
}


/**********************************************************************************************/
/**********************************************************************************************
  alloc_grid_arrays():
 ------------------
   -- allocates all the global grid functions given the subdomain dimensions;

**********************************************************************************************/
void alloc_grid_arrays(void) 
{

  unsigned long int i,j,k,l,n,d;
  unsigned int ndim_m1 = (NDIM-1);
  unsigned int npos_m1 = (NPOS-1);

/* 123  MHD  functions per cell:     */
  ALLOC_4D_ARRAY(p    ,N1TOT,N2TOT,N3TOT,NP);
  ALLOC_4D_ARRAY(ph   ,N1TOT,N2TOT,N3TOT,NP);
  ALLOC_4D_ARRAY(p_old,N1TOT,N2TOT,N3TOT,NP);

  ALLOC_5D_ARRAY(p_L,N1TOT,N2TOT,N3TOT,ndim_m1,NP);
  ALLOC_5D_ARRAY(p_R,N1TOT,N2TOT,N3TOT,ndim_m1,NP);
  ALLOC_5D_ARRAY(F  ,N1TOT,N2TOT,N3TOT,ndim_m1,NP);

#if(USE_PRESSURE_FLUX_FIX)
  ALLOC_4D_ARRAY(F2 ,N1TOT,N2TOT,N3TOT,ndim_m1);
#endif

  for(n=0;n<N0TOT;n++) { ALLOC_4D_ARRAY(U_gf[n],N1TOT,N2TOT,N3TOT,NP); }

  ALLOC_3D_ARRAY(p_gamma,N1TOT,N2TOT,N3TOT);

  ALLOC_5D_ARRAY(c_gf  ,N1TOT,N2TOT,N3TOT,ndim_m1,2);
  ALLOC_5D_ARRAY(ucon_L,N1TOT,N2TOT,N3TOT,ndim_m1,2);
  ALLOC_5D_ARRAY(ucon_R,N1TOT,N2TOT,N3TOT,ndim_m1,2);

  for(n=0;n<N0_GEOM ;n++) { ALLOC_ARRAY(geom_arr[n] ,N_GEOM ) ; }
  for(n=0;n<N0_COORD;n++) { ALLOC_ARRAY(coord_arr[n],N_COORD); }

  ALLOC_4D_ARRAY(conn,N_CONN,NDIM,NDIM,NDIM);

  ALLOC_3D_ARRAY(pflag,N1TOT,N2TOT,N3TOT);

  for(n=0;n<N_NFAIL;n++) { ALLOC_3D_ARRAY(nfail[n],N1TOT,N2TOT,N3TOT); }

#if( USE_LOCAL_RECON_TYPE ) 
  ALLOC_4D_ARRAY(recon_type,N1TOT,N2TOT,N3TOT,ndim_m1);
#endif

#if(USE_MASK)
  for(n=0;n<N0_GEOM;n++) {  ALLOC_3D_ARRAY(evol_mask[n],N1TOT,N2TOT,N3TOT); }
#endif

#if( RESCALE_REGULARIZE )
  ALLOC_5D_ARRAY(  regularize_gf,N1TOT,N2TOT,N3TOT,npos_m1,ndim_m1);
  ALLOC_5D_ARRAY(unregularize_gf,N1TOT,N2TOT,N3TOT,npos_m1,ndim_m1);
#endif

#if( CALC_CURRENT ) 
  ALLOC_4D_ARRAY(jcon ,N1TOT,N2TOT,N3TOT,NDIM);
  ALLOC_4D_ARRAY(jcon2,N1TOT,N2TOT,N3TOT,NDIM);

  for(n=0;n<N0TOT;n++) {   ALLOC_5D_ARRAY(faraday[n],N1TOT,N2TOT,N3TOT,NDIM,NDIM); }
#endif

#if( USE_KINEMATIC_VISCOSITY )
  ALLOC_4D_ARRAY(visc_source,N1TOT,N2TOT,N3TOT,NDIM);
  for(n=0;n<3;n++) {   ALLOC_3D_ARRAY(visc_funcs[n],N1TOT,N2TOT,N3TOT); }
#endif

#if( USE_COOLING_FUNCTION == 4 ||  USE_COOLING_FUNCTION == 5 )
  ALLOC_3D_ARRAY(lum_thermal,N1TOT,N2TOT,N3TOT);
  ALLOC_3D_ARRAY(      dflux,N1TOT,N2TOT,N3TOT);
  ALLOC_3D_ARRAY(       dtau,N1TOT,N2TOT,N3TOT);
  ALLOC_2D_ARRAY(j_photo_bot,N1TOT,N3TOT);
  ALLOC_2D_ARRAY(j_photo_top,N1TOT,N3TOT);
  ALLOC_2D_ARRAY(    tau_bot,N1TOT,N3TOT);
  ALLOC_2D_ARRAY(    tau_top,N1TOT,N3TOT);
  ALLOC_2D_ARRAY(   flux_top,N1TOT,N3TOT);
  ALLOC_2D_ARRAY(   flux_bot,N1TOT,N3TOT);
  ALLOC_2D_ARRAY(  flux_disk,N1TOT,N3TOT);
  
# if( USEMPI ) 
  i = N1TOT*N3TOT;
  n = 2*i; 
  ALLOC_ARRAY( recv_buffer_top,n);
  ALLOC_ARRAY( recv_buffer_bot,n);
  ALLOC_ARRAY( send_buffer_top,n);
  ALLOC_ARRAY( send_buffer_bot,n);
  ALLOC_ARRAY( send_buffer_tot,i);
  ALLOC_ARRAY( recv_buffer_tot,i);
  ALLOC_ARRAY(irecv_buffer_top,i);
  ALLOC_ARRAY(irecv_buffer_bot,i);
  ALLOC_ARRAY(isend_buffer_top,i);
  ALLOC_ARRAY(isend_buffer_bot,i);
  ALLOC_ARRAY(isend_buffer_tot,n);
  ALLOC_ARRAY(irecv_buffer_tot,n);
# endif
#endif
#if( USE_COOLING_FUNCTION == 5 )
  ALLOC_ARRAY(radfunc_Te_eV,NCELLS); 
#endif


  ALLOC_ARRAY(Katm   ,N1TOT);
  ALLOC_ARRAY(nu_visc,N1TOT);

#if( MAKE_HISTORY && MAKE_HDF5 ) 
  for(i=0;i<N_HISTORY;i++) for(j=0;j<N_HIST_TYPES;j++) { ALLOC_ARRAY(hist_data[i][j]    ,N1          ); }
  for(i=0;i<N_HISTORY;i++) for(j=0;j<N_HIST_TYPES;j++) { ALLOC_ARRAY(hist_data_out[i][j],totalsize[1]); }

  /* Allocate the buffer used for recv communication for history_mpi() : */
  int n_mpi = ncpux[2] * ncpux[3] - 1  ;
  if( n_mpi > 0 ) {  ALLOC_2D_ARRAY(history_buffer,n_mpi,N_HIST_POINTS);  }

  /* Allocate the buffer used for recv communication for history_mpi2() : */
  n_mpi = ncpux[1] - 1  ;
  if( n_mpi > 0 ) {  ALLOC_2D_ARRAY(history_buffer2,n_mpi,N_HIST_POINTS);  }
#endif

#if( MAKE_HISTORY_GEN && MAKE_HDF5 ) 
  /* Allocate the buffer used for recv communication for history_mpi() : */
  int n_mpi = ncpux[2] * ncpux[3] - 1  ;
  if( n_mpi > 0 ) {  ALLOC_2D_ARRAY(history_buffer,n_mpi,N_HIST_POINTS);  }

  /* Allocate the buffer used for recv communication for history_mpi2() : */
  n_mpi = ncpux[1] - 1  ;
  if( n_mpi > 0 ) {  ALLOC_2D_ARRAY(history_buffer2,n_mpi,N_HIST_POINTS);  }
#endif

#if( MAKE_SURFACE && MAKE_HDF5 ) 
  n = N1*N3;
  for(i=0;i<N_SURFACE;i++) for(j=0;j<N_SURF_TYPES;j++) { ALLOC_ARRAY(surf_data[i][j],n); }
#endif

  FACE_LOOP {    ALLOC_3D_ARRAY(emf[d],N1TOT,N2TOT,N3TOT); }

  ALLOC_ARRAY( p_vect,MAX_NTOT);
  ALLOC_ARRAY(dq_vect,MAX_NTOT);
  ALLOC_ARRAY(pL_vect,MAX_NTOT);
  ALLOC_ARRAY(pR_vect,MAX_NTOT);

#if( RESCALE_R )
  for(i=0;i<6;i++) { ALLOC_ARRAY(    rscale[i],N1TOT); }
  for(i=0;i<6;i++) { ALLOC_ARRAY(inv_rscale[i],N1TOT); }
#endif


#if( MAKE_SDF )
  ALLOC_ARRAY(f_sdf ,n_cells_glob);
  ALLOC_ARRAY(f_sdf2,n_cells_glob);
#endif

#if( MAKE_HDF5 )
  for(i=0;i<N_HDF_FUNCS;i++) { ALLOC_ARRAY(f_hdf[i],NCELLS);  }
#endif


#if( USEMPI ) 
  ALLOC_ARRAY(hist_send ,N_HIST_POINTS); 
  ALLOC_ARRAY(hist_send2,N_HIST_POINTS); 
  ALLOC_ARRAY(surf_send ,N_SURF_POINTS); 
#endif


#if( DUMP_ALL_STAT )
  PLOOP  { ALLOC_ARRAY(U_avg[l],NCELLS); }
#endif

#if( MAKE_STAT && DUMP_ALL_STAT ) 
  int g; 
  PLOOP { ALLOC_ARRAY(U_out[l],NCELLS); }
  PLOOP { d = 0;  LOOP U_out[l][d++] = 0.; }

  PLOOP { ALLOC_ARRAY(U_pre[l],NCELLS); }
  PLOOP { d = 0;  LOOP U_pre[l][d++] = 0.; }

  NPH_LOOP for(g=0; g<N_STAT; g++)  { 
    ALLOC_ARRAY(dU_stat[l][g],NCELLS);
    d = 0;  LOOP dU_stat[l][g][d++] = 0. ; 
  }
#endif

#if( MAKE_STAT2 )
  int g; 
  FACE_LOOP  PLOOP { 
    ALLOC_ARRAY(FL_stat2[d][l],NCELLS);
    g = 0;  LOOP FL_stat2[d][l][g++] = 0.;
  }
  FACE_LOOP  PLOOP { 
    ALLOC_ARRAY(FR_stat2[d][l],NCELLS);
    g = 0;  LOOP FR_stat2[d][l][g++] = 0.;
  }
  FACE_LOOP  PLOOP { 
    ALLOC_ARRAY(UL_stat2[d][l],NCELLS);
    g = 0;  LOOP UL_stat2[d][l][g++] = 0.;
  }
  FACE_LOOP  PLOOP { 
    ALLOC_ARRAY(UR_stat2[d][l],NCELLS);
    g = 0;  LOOP UR_stat2[d][l][g++] = 0.;
  }
  FACE_LOOP  { 
    ALLOC_ARRAY(ctop_stat2[d],NCELLS);
    g = 0;  LOOP ctop_stat2[d][g++] = 0.;
  }
  for(l=0;l<2;l++) { 
    ALLOC_ARRAY(S_stat2[l],NCELLS);
    g = 0;  LOOP S_stat2[l][g++] = 0.;
  }
#endif

#if( MAKE_RADFLUX ) 
  for(l=0;l<NDIM;l++) { ALLOC_ARRAY(coolflux[l],NCELLS); }
  for(l=0;l<NDIM;l++) for(i=0; i<NCELLS; i++) {  coolflux[l][i] = 0.; }
#endif

#if( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL )  
  for(l=0; l<3; l++) { 
    coord_params = &(coord_params_all[l]);
    ALLOC_2D_ARRAY( coord_params->square_xp1_3,N1TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_xp1_4,N1TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->dsquare_xp1_3,N1TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->dsquare_xp1_4,N1TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_per_xp2_1,N2TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_per_xp2_2,N2TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_per_xp3_1,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_per_xp3_2,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_per_xp3_3,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_per_xp3_4,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_int_per_xp2_1,N2TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_int_per_xp2_2,N2TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_int_per_xp3_1,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_int_per_xp3_2,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_int_per_xp3_3,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->square_int_per_xp3_4,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->dsquare_per_xp2_1,N2TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->dsquare_per_xp2_2,N2TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->dsquare_per_xp3_3,N3TOT,NPOS);
    ALLOC_2D_ARRAY( coord_params->dsquare_per_xp3_4,N3TOT,NPOS);
    ALLOC_3D_ARRAY( coord_params->r_of_yt_,N1TOT,NPOS,3);
    ALLOC_3D_ARRAY( coord_params->drdy_of_yt_,N1TOT,NPOS,3);
  }
#endif

#if( COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN ) 
  for(l=0; l<3; l++) { 
    coord_params = &(coord_params_all[l]);
    ALLOC_2D_ARRAY(  coord_params->square_per_xp2_1      N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp2_2      N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp2_3      N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp2_4      N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->dsquare_per_xp2_3     N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->dsquare_per_xp2_4     N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_int_per_xp2_1  N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_int_per_xp2_2  N2TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp1_1      N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp1_2      N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp1_3      N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_per_xp1_4      N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->dsquare_per_xp1_3     N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->dsquare_per_xp1_4     N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_int_per_xp1_1  N1TOT,NPOS) ; 
    ALLOC_2D_ARRAY(  coord_params->square_int_per_xp1_2  N1TOT,NPOS) ; 
  }
#endif

#if( OUTPUT_MAX_VCHAR )
  ALLOC_4D_ARRAY( max_vchar ,N1TOT,N2TOT,N3TOT,ndim_m1);
#endif

  
  return;
}
