
#include "decs.h"


#if( MAKE_SURFACE && MAKE_HDF5 ) 
#include <hdf5.h>

#if( USEMPI )
#include "mpi.h"
extern void surface_mpi( double *surf_arr[N_SURFACE][N_SURF_TYPES] );
#endif 

/* Rank of surface datasets (time versus X1): */
#define SURF_RANK (2) 

/* Use chunks for transfer or not: */ 
#define USE_CHUNKS (1) 

/* Extents of memory and file spaces */
static hsize_t surf_memdims[SURF_RANK],surf_filedims[SURF_RANK];
static hsize_t offset[SURF_RANK];
static hsize_t  count[SURF_RANK]; 


/* Name of surface functions :  */
static char *S_names[N_SURFACE] = { 
  "S_bsq"  , 
  "S_pavg" , 
  "S_rhoav", 
  "S_ddot" , 
  "S_hUz"  , 
  "S_hUt"  , 
  "S_hUx"  , 
  "S_lrho" , 
  "S_Txza" , 
  "S_Txta" , 
  "S_Txxa" , 
  "S_Ttta" , 
  "S_Txzb" , 
  "S_Txtb" , 
  "S_Txxb" , 
  "S_Tttb" , 
  "S_Txzc" , 
  "S_Txtc" , 
  "S_Txxc" , 
  "S_Tttc" , 
  "S_vol"  
#if( CALC_CURRENT )
 ,"S_Jx"  ,      
  "S_Jy"  ,      
  "S_Jz"  ,      
  "S_Jt"  ,      
  "S_Jsq"  
#endif
#if( MAKE_RADFLUX )
 ,"S_Lut"  ,      
  "S_Lut2" ,      
  "S_Lup"        
#endif
#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
 ,"S_dTdr_1",     
  "S_dTdr_2",     
  "S_dTdr_3",     
  "S_dTdr_4",     
  "S_dTdr_5",     
  "S_dTdr_6"     
#endif
};


enum surf_func_index  { S_bsq   ,S_pavg  ,S_rhoav ,S_ddot  ,S_hUz   ,S_hUt   ,S_hUx   ,S_lrho  
			,S_Txza  ,S_Txta  ,S_Txxa  ,S_Ttta  
			,S_Txzb  ,S_Txtb  ,S_Txxb  ,S_Tttb  
			,S_Txzc  ,S_Txtc  ,S_Txxc  ,S_Tttc  
			,S_vol 	
#if( CALC_CURRENT )
			,S_Jx, S_Jy, S_Jz, S_Jt, S_Jsq
#endif
#if( MAKE_RADFLUX )
			,S_Lut, S_Lut2, S_Lup
#endif
#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
			,S_dTdr_1,S_dTdr_2,S_dTdr_3,S_dTdr_4,S_dTdr_5,S_dTdr_6
#endif
			,S_last                    /* used to verify length of surf_func_index */
};

/* Names of hdf groups under which surfaces live */
static char *surf_dirnames[] = { "Unbound" , "Bound", "Jet" };  

/* Name of surface file */
static char surf_filename[200];   


/* Routines : */
extern void write_dump_header(hid_t file_id, char *hdf_name, int all_procs_dump ) ;

static void write_surface_data( void ) ;
static void myH5_write_surffunc( hid_t loc_id, char *name, double *value );
void myH5_write_array( hid_t loc_id, char *name, hid_t type_id, int size, void *value ) ;


/**********************************************************************************************/
/**********************************************************************************************
  dump_surface_gen(): 
 -------------
   -- "generalized" version of dump_surface() that bins data on a static radial and azimuthal 
          coordinate system;
       -- most useful for non-spherical or dynamic coordinate systems;

   -- like dump_history() but does not integrate in phi, so outputs data in  r,phi,t  
      with functions reduced in theta 

   -- The below quanities use HARM notation.  Any translation between HARM 
       definitions and Hawley's definitions should be done at the IDL/analysis 
       level.

   -- we do not use primtoflux() to calculate stress terms since we want to separate
       the parts out ;  also, we want to calculate the components in KS (or r,phi) coordinates not 
       code coordinates;

   -- The choice in nomenclature w.r.t. "x,y,z" was done to be consistent with 
      Beckwith's nomenclature.

   -- note that the surf_data[][][] array is ordered that way so that it is easier to 
        pass it to  surface_mpi() ;


LEGEND: 

Group : Name  : Quantity      (x=r, y=theta , z=phi)
------------------------------------------------
/  : 
        t    : time of this surface slice
       dt    : value of timestep 

    /Unbound :
    /Bound   :
    /Jet     :
           (the below reside in each of the subgroups "Unbound","Bound","Jet");
               bsq     : ||b^2||
               pavg    : pressure
               rhoav   : \rho (proper rest-mass density)
               ddot    : \rho u^r (rest-mass flux)
               hUt     : \rho h u^t (enthalpy times time four-velocity)
               hUx     : \rho h u^r (enthalpy times time four-velocity)
               hUz     : \rho h u^\phi (enthalpy times azimuthal four-velocity)
               lrho    : l \rho (total angular momentum) (where l = -u_\phi / u_t  )
               Ttta    : \rho h u^t u_t (conserved energy density)
               Txta    : \rho h u^r u_t (energy flux of matter)
               Txxa    : \rho h u^r u_r (momentum flux of matter)
               Txza    : \rho h u^r u_\phi (angular momentum flux of matter)
               Tttb    : ||b^2|| u^t u_t  
               Txtb    : ||b^2|| u^r u_t  
               Txxb    : ||b^2|| u^r u_r
               Txzb    : ||b^2|| u^r u_\phi  
               Tttc    : -b^t b_t
               Txtc    : -b^r b_t
               Txxc    : -b^r b_r
               Txzc    : -b^r b_\phi  
               Jt      : J^t  current
               Jx      : J^r  current
               Jy      : J^\theta current
               Jz      : J^\phi  current
               Jsq     : J^\mu J_\mu  current density
               Lut     : L u_t  (L = fluid-frame cooling rate)
               Lut2    : L u_t / u^t 
               Lup     : L u_\phi
               dTdr_1  : rho   u^a u_b \Gamma^b_{a \phi}
               dTdr_2  : rho h u^a u_b \Gamma^b_{a \phi}
               dTdr_3  : bsq   u^a u_b \Gamma^b_{a \phi}
               dTdr_4  : p             \Gamma^a_{a \phi}
               dTdr_5  : bsq/2         \Gamma^a_{a \phi}
               dTdr_6  :     - b^a b_b \Gamma^b_{a \phi}
               vol     : \sqrt{-g} \Delta t \Delta r d\theta d\phi 



  -- hdf details: 
      -- the datasets in the surface file have dimensions [totalsize1,totalsize3]
      -- even though they are written frequently, we will write an individual file per 
         timestep.   Otherwise the files will be too large as there are many function dumped. 
      -- hence, the order of operations is the following:
            if new: 
                -- create file, create dataset structures
            if old: 
                -- open file, open dataset, get dimensions,
                   write the data, and close things ;

**********************************************************************************************/
//void dump_surface_gen( void )
void dump_surface( void )
{
  int i,j,k,l,id,itype,ii,jj,ic,kk,pos;
  int set_type[N_SURF_TYPES];

  double *ptmp, dV_prop;
  double bcon_ks[NDIM],bcov_ks[NDIM],ucon_ks[NDIM],ucov_ks[NDIM];
  double bsq, enthalpy, rho, jcon_i[NDIM],jcon_ks[NDIM], jsq;
  double coolrate;
  double conn_trace, uu_dot_conn, bb_dot_conn, conntmp[NDIM][NDIM];
  struct of_geom  *geom; 
  struct of_state q;
  struct of_coord *coords;
  struct of_coord *coords_m1;
  struct of_coord *coords_p1;
  struct of_coord *coords_p2;
  static const unsigned short int shift[NDIM][NDIM] = {{ 1,0,0,0},
						       { 0,1,0,0},
						       { 0,0,1,0},
						       { 0,0,0,1}};

  double ***surf_data_binned;
  double *r_beg, *r_end, *r_bin_mid; 
  double *phi_beg, *phi_end, *phi_bin_mid; 
  double weight_loc, weight_loc2; 
  double phi1,phi2;
  ulint icell, ibin, ibins_b, ibins_e, kbin, kbins_b, kbins_e, ifunc, n_all_bins;
  usint *used_bins;

  static double surf_data_loc[N_SURFACE];

  TRACE_BEG;

  if( N_SURFACE != S_last ) { 
    fprintf(stderr,"dump_surface_gen():  Invalid values of N_SURFACE and S_last:  %d %d \n",N_SURFACE,S_last); 
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  /*   fprintf(stdout,"Doing surface dump (%d) %d ... \n", myid,N_out[OUT_SURFACE]);   fflush(stdout);  */

  /*********************************************************************************** 
     First reset the surface data so we can begin tallying : 
  ***********************************************************************************/ 
  n_all_bins = n_r_bins * n_phi_bins; 

  if( myid == out_pid[OUT_SURFACE] ) { 
    for(i=0;i<N_SURFACE;i++) for(j=0;j<N_SURF_TYPES;j++) { ALLOC_ARRAY(surf_data[i][j]   ,n_all_bins); }
    for(i=0;i<N_SURFACE;i++) for(j=0;j<N_SURF_TYPES;j++)  for(k=0; k<n_all_bins; k++) { surf_data[i][j][k] = 0.; }
  }

  /*********************************************************************************** 
     Setup bin arrays and weights:
  ***********************************************************************************/ 
  ALLOC_ARRAY(r_beg,NCELLS); 
  ALLOC_ARRAY(r_end,NCELLS); 
  ALLOC_ARRAY(phi_beg,NCELLS); 
  ALLOC_ARRAY(phi_end,NCELLS); 

  icell = 0; 

  LOOP { 
#if( TOP_TYPE_CHOICE == TOP_CARTESIAN )
    get_coord(i  ,j  ,k  ,CORN,ncurr,coords_m1);
    get_coord(i+1,j+1,k+1,CORN,ncurr,coords_p1);
    r_beg[icell] = coords_m1->r; 
    r_end[icell] = coords_p1->r; 

    phi1 = atan2( coords_m1->x[YY], coords_m1->x[XX] );
    phi2 = atan2( coords_p1->x[YY], coords_p1->x[XX] );
    if( phi1 < 0. ) { phi1 += 2.*M_PI; }
    if( phi2 < 0. ) { phi2 += 2.*M_PI; }
#else
    get_coord(i  ,j,k,FACE1,ncurr,coords_m1);
    get_coord(i+1,j,k,FACE1,ncurr,coords_p1);
    r_beg[icell] = coords_m1->r; 
    r_end[icell] = coords_p1->r; 

    get_coord(i,j,k  ,FACE3,ncurr,coords_m1);
    get_coord(i,j,k+1,FACE3,ncurr,coords_p1);
    phi1 = coords_m1->x[PH];
    phi2 = coords_p1->x[PH];
#endif

    phi_beg[icell] = phi1; 
    phi_end[icell] = phi2; 

    icell++;
  }

  set_bin_weights_r(  NCELLS,  r_beg,  r_end,&bin_weights_r  ,&ibins_beg,&ibins_end); 
  set_bin_weights_phi(NCELLS,phi_beg,phi_end,&bin_weights_phi,&kbins_beg,&kbins_end); 

  for(icell=0;icell<NCELLS;icell++) {   r_beg[icell] =   r_end[icell] -   r_beg[icell]; }   /* reassign r_beg to be dr_cell */
  for(icell=0;icell<NCELLS;icell++) { phi_beg[icell] = phi_end[icell] - phi_beg[icell]; }   /* reassign phi_beg to be dphi_cell */

  FREE(  r_end);
  FREE(phi_end);

  ALLOC_ARRAY(used_bins,n_all_bins); 
  for(ibin=0;ibin<n_all_bins;ibin++) { used_bins[ibin] = 0; } 
  for(icell=0;icell<NCELLS;icell++) { 
    ibins_b = ibins_beg[icell];     
    ibins_e = ibins_end[icell]; 
    kbins_b = kbins_beg[icell];     
    kbins_e = kbins_end[icell]; 
    for(ibin=ibins_b;ibin<=ibins_e;ibin++) { 
      i = ibin*n_phi_bins + kbins_b;
      for(kbin=kbins_b;kbin<=kbins_e;kbin++) { used_bins[i++] = 1;  }
    }
  }

  i = 0; 
  for(ibin=0;ibin<n_r_bins;ibin++) { 
    fprintf(stdout,"used");
    for(kbin=0;kbin<n_phi_bins;kbin++) {  fprintf(stdout,"%hu",used_bins[i]);    }
    fprintf(stdout,"\n");
  }
  fflush(stdout);
  

  ALLOC_ARRAY(surf_data_binned,n_all_bins); 
  for(ibin=0;ibin<n_all_bins;ibin++) { 
    if( used_bins[ibin] ) { 
      ALLOC_2D_ARRAY(surf_data_binned[ibin],N_SURF_TYPES,N_SURFACE); 
      for( itype = 0; itype < N_SURF_TYPES ; itype++ )  for(ifunc=0; ifunc<N_SURFACE; ifunc++) { 
	  surf_data_binned[ibin][itype][ifunc] = 0. ; 
	}
    }
    else { 
      surf_data_binned[ibin] = NULL;
    }
  }


  /*********************************************************************************** 
     Begin the loop over all space, integrating along x2 and x3 directions : 
  ***********************************************************************************/ 
  icell = 0; 

  LOOP  {
    get_geometry(i,j,k,CENT,ncurr,geom);  
    get_coord(   i,j,k,CENT,ncurr,coords);
    dV_prop = dV * geom->g  / (r_beg[icell] * phi_beg[icell]); 

    /* Get local copies of fluid state: */
    ptmp = p[i][j][k];
    get_state(  ptmp, geom, &q );

    rho = ptmp[RHO];
    bsq = q.bsq ;
    enthalpy = rho + gam*ptmp[UU];

#if( CALC_CURRENT )
    /* Calculate the current density :  */
    for(l=0; l<NDIM; l++)  { jcon_i[l] = jcon[i][j][k][l]; }
	
    jsq = 0.;
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {
	jsq += geom->gcov[ii][jj] * jcon_i[ii] * jcon_i[jj] ; 
      }
#endif
	
    /* Transform the 4-vectors from code coordinates to KS coordinates : 
       -- We do not use the transform_rank*co*() routines since it would 
       add considerable additional operations to the calculation : 
    */
    for(ii=0;ii<NDIM;ii++) {bcon_ks[ii] = bcov_ks[ii] = jcon_ks[ii] = ucon_ks[ii] = ucov_ks[ii] = 0.;}
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {bcon_ks[ii] += coords->dx_dxp[ii][jj]*q.bcon[jj];}
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {ucon_ks[ii] += coords->dx_dxp[ii][jj]*q.ucon[jj];}
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {bcov_ks[ii] += coords->dxp_dx[jj][ii]*q.bcov[jj];}
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {ucov_ks[ii] += coords->dxp_dx[jj][ii]*q.ucov[jj];}
#if( CALC_CURRENT )
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {jcon_ks[ii] += coords->dx_dxp[ii][jj]*jcon_i[jj];}
#endif

#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
	/* dxp_dxp_dphi[i][j] =  d^2 x^{\hat{i}} / ( dx^{\hat{j}} dphi )   =  d/dx^{\hat{j}}  dxp_dx[i][PH] 	*/
	/* Calculate time derivatives using  2nd-order vertex-centered forwards differencing  since only future 
	   coordinate data is available when this routine is called */
	jj = 0 ; 
	pos = CENT;
	get_coord(i,j,k,pos,n_mid,coords_p1);
	get_coord(i,j,k,pos,n_end,coords_p2);
	for(ii=0;ii<NDIM;ii++) { 
	    conntmp[ii][jj]  = (-coords_p2->dxp_dx[ii][PH] + 4.*coords_p1->dxp_dx[ii][PH] - 3.*coords->dxp_dx[ii][PH]) / dx[0];
	  }

	/* Calculate spatial derivatives using standard centered 2nd-order finite differences using  FACE-centered data: */
	for(jj=1;jj<NDIM;jj++)  { 
	  pos = FACE1 + jj - 1;
	  get_coord(i             ,j             ,k             ,pos,ncurr,coords_m1);
	  get_coord(i+shift[1][jj],j+shift[2][jj],k+shift[3][jj],pos,ncurr,coords_p1);
	  for(ii=0;ii<NDIM;ii++) { 
	      conntmp[ii][jj] = (coords_p1->dxp_dx[ii][PH] - coords_m1->dxp_dx[ii][PH]) * invdx[jj]; 
	    }
	}

	/* conntmp[nu][mu] is the substitute of \Gamma^{\hat{\nu}}_{\hat{\mu} \phi} for general warped coordinates;  
	   See the coords.tex  notes.   It should be equal to 
	   d/dx^{\hat{\mu}} dx^{\hat{\nu}}/dx^{\phi} + dx^{\hat{\kappa}}/dx^{\phi} \Gamma^{\hat{\nu}}_{\hat{\mu} \hat{\kappa}}
	*/
	ic = CONN_ID(i,j,k);
	for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
	      conntmp[ii][jj] += coords->dxp_dx[kk][PH] * conn[ic][ii][jj][kk];
	  }

	conn_trace = uu_dot_conn = bb_dot_conn = 0.;
	for(ii=0;ii<NDIM;ii++) { conn_trace  += conntmp[ii][ii]; }
	for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) { uu_dot_conn += conntmp[ii][jj] * q.ucon[jj] * q.ucov[ii]; }
	for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) { bb_dot_conn += conntmp[ii][jj] * q.bcon[jj] * q.bcov[ii]; }
#endif


    /* Determine the cell's type (e.g. bound,unbound,jet) : 
       -- make sure the order of set_type[] corresponds to order of surf_dirnames[] */
    set_type[UNBOUND] = ((enthalpy+bsq)*ucov_ks[TT]/rho < -1.);           /* Unbound */
    set_type[  BOUND] = !set_type[UNBOUND];                         /* Bound   */
    set_type[    JET] = set_type[UNBOUND]  && (ucon_ks[RR] > 0.) ;  /* Jet     */

    /* Add to the surface integration : */
    for( itype = 0; itype < N_SURF_TYPES ; itype++ )  if( set_type[itype] )  { 
	surf_data_loc[S_bsq  ] =  (bsq);
	surf_data_loc[S_pavg ] =  ((gam-1.)*ptmp[UU]); 
	surf_data_loc[S_rhoav] =  (rho); 
	surf_data_loc[S_ddot ] =  (rho      *  ucon_ks[RR]); 
	surf_data_loc[S_hUt  ] =  (enthalpy *  ucon_ks[TT]); 
	surf_data_loc[S_hUx  ] =  (enthalpy *  ucon_ks[RR]); 
	surf_data_loc[S_hUz  ] =  (enthalpy *  ucon_ks[PH]); 
	surf_data_loc[S_lrho ] =  (-rho     * (ucov_ks[PH] / ucov_ks[TT])); 
	surf_data_loc[S_Ttta ] =  (enthalpy *  ucon_ks[TT] * ucov_ks[TT] ); 
	surf_data_loc[S_Txta ] =  (enthalpy *  ucon_ks[RR] * ucov_ks[TT] ); 
	surf_data_loc[S_Txxa ] =  (enthalpy *  ucon_ks[RR] * ucov_ks[RR] ); 
	surf_data_loc[S_Txza ] =  (enthalpy *  ucon_ks[RR] * ucov_ks[PH] ); 
	surf_data_loc[S_Tttb ] =  (bsq      *  ucon_ks[TT] * ucov_ks[TT] ); 
	surf_data_loc[S_Txtb ] =  (bsq      *  ucon_ks[RR] * ucov_ks[TT] ); 
	surf_data_loc[S_Txxb ] =  (bsq      *  ucon_ks[RR] * ucov_ks[RR] ); 
	surf_data_loc[S_Txzb ] =  (bsq      *  ucon_ks[RR] * ucov_ks[PH] ); 
	surf_data_loc[S_Tttc ] =  (           -bcon_ks[TT] * bcov_ks[TT] ); 
	surf_data_loc[S_Txtc ] =  (           -bcon_ks[RR] * bcov_ks[TT] ); 
	surf_data_loc[S_Txxc ] =  (           -bcon_ks[RR] * bcov_ks[RR] ); 
	surf_data_loc[S_Txzc ] =  (           -bcon_ks[RR] * bcov_ks[PH] ); 
	surf_data_loc[S_vol  ] =  1.;

#if( CALC_CURRENT )
	surf_data_loc[S_Jt   ] =  (jcon_ks[TT]); 
	surf_data_loc[S_Jx   ] =  (jcon_ks[RR]); 
	surf_data_loc[S_Jy   ] =  (jcon_ks[TH]); 
	surf_data_loc[S_Jz   ] =  (jcon_ks[PH]); 
	surf_data_loc[S_Jsq  ] =  (jsq); 
#endif


#if( MAKE_RADFLUX )
# if( USE_COOLING_FUNCTION == 1 ) 
	coolrate = cooling_func_hr_disk( i,j,k, ptmp );
# elif( USE_COOLING_FUNCTION == 2 || USE_COOLING_FUNCTION == 3 ) 
	coolrate = cooling_func_isentropic_disk( i, j, k, ptmp, &q ); //Dennis
# elif( USE_COOLING_FUNCTION >= 4 )
	  coolrate = coolflux[0][k-N3S + N3*((j-N2S) + N2*(i-N1S))];
# endif
        surf_data_loc[S_Lut  ] =  (coolrate*ucov_ks[TT]); 
	surf_data_loc[S_Lut2 ] =  (coolrate*ucov_ks[TT]/ucon_ks[TT]); 
	surf_data_loc[S_Lup  ] =  (coolrate*ucov_ks[PH]);
#endif

#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
	surf_data_loc[S_dTdr_1] = (rho                * uu_dot_conn);
	surf_data_loc[S_dTdr_2] = (enthalpy           * uu_dot_conn);
	surf_data_loc[S_dTdr_3] = (bsq                * uu_dot_conn);
	surf_data_loc[S_dTdr_4] = ((gam-1.)*ptmp[UU]  * conn_trace );
	surf_data_loc[S_dTdr_5] = (0.5*bsq            * conn_trace );
	surf_data_loc[S_dTdr_6] = (-bb_dot_conn);
#endif

	ibins_b = ibins_beg[icell]; 	
	ibins_e = ibins_end[icell]; 
	kbins_b = kbins_beg[icell]; 	
	kbins_e = kbins_end[icell]; 
	jj = 0; 
	for(ibin=ibins_b; ibin<=ibins_e; ibin++) {
	  weight_loc = dV_prop * bin_weights_r[icell][jj];
	  kk = 0; 
	  ii = ibin*n_phi_bins + kbins_b;
	  for(kbin=kbins_b; kbin<=kbins_e; kbin++) {
	    weight_loc2 = weight_loc * bin_weights_phi[icell][kk];
	    for(ifunc=0;ifunc<N_SURFACE;ifunc++) {
	      /* If there is a seg fault here, it is because we have not set used_bin[] correctly above */
	      surf_data_binned[ii][itype][ifunc] += weight_loc2*surf_data_loc[ifunc] ; 
	    }
	    ii++;
	    kk++;
	  }
	  jj++;
	}

      }
    icell++;
  }

  if( bin_weights_r   != NULL ) {  for(icell=0;icell<NCELLS;icell++) { FREE(bin_weights_r[icell]  ); }  }
  if( bin_weights_phi != NULL ) {  for(icell=0;icell<NCELLS;icell++) { FREE(bin_weights_phi[icell]); }  }
  FREE( bin_weights_r    ); 
  FREE( bin_weights_phi  ); 
  FREE( ibins_beg );   FREE( ibins_end ); 
  FREE( kbins_beg );   FREE( kbins_end ); 
  FREE(  r_beg);
  FREE(phi_beg);

  k = 0; 
  for(ibin=0; ibin<n_r_bins; ibin++) {
    weight_loc = dr_bins[ibin] * dphi_bins;
    for(kbin=0; kbin<n_phi_bins; kbin++) {
      if( used_bins[k] ) { 
	for(j=0; j<N_SURF_TYPES; j++) for(i=0; i<N_SURFACE; i++)  {
	    surf_data_binned[k][j][i] *= weight_loc;
	  }
      }
      k++;
    }
  }

  /*********************************************************************************** 
    If we are running in parallel, then sum over the grids in x2 directions : 
  ***********************************************************************************/ 
#if( USEMPI ) 
  surface_mpi(  surf_data_binned ); 
#else 
  for(k=0; k<n_all_bins; k++) { 
    if( used_bins[k] ) { 
      for(j=0; j<N_SURF_TYPES; j++)  for(i=0; i<N_SURFACE; i++)  {
	  surf_data[i][j][k] = surf_data_binned[k][j][i]; 
	}
    }
  }
#endif 

  for(ibin=0;ibin<n_all_bins;ibin++) { 
    if( surf_data_binned[ibin] != NULL ) { 
      for( itype = 0; itype < N_SURF_TYPES ; itype++ )  { FREE(surf_data_binned[ibin][itype]);  }
      FREE(surf_data_binned[ibin]); 
    }
  }
  FREE(surf_data_binned);
  FREE(used_bins);

  /*********************************************************************************** 
    Write the surface data to the open file: write, and close: 
  ***********************************************************************************/ 
  if( cpupos[2] == 0 ) { 
    write_surface_data();
  }


#if( USEMPI ) 
  //  exit_status();
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 

  TRACE_END;

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 write_surface_data():
 ---------------------------
   -- initializes arrays and file used for all surface files; 
   -- note that based on how globalpos, totalsize and N[1-3] are defined, 
       a serial run will be automatically handled;
*******************************************************************************************/
static void write_surface_data( void ) 
{
  int i,j,k ;
  hid_t file_id, group_id;
  char dataname[200];

  /************************************************************************
    Set various static arrays used in writing to the surface file : 
  *************************************************************************/ 
  /* Setup the extent of each slice of a surface dataset in memory : */
  surf_memdims[     0] = n_r_bins;               surf_memdims[     1] = n_phi_bins;

  for(i=0;i<SURF_RANK;i++) { 
    surf_filedims[i] = surf_memdims[i];
    offset[i] = 0; 
    count[i] = 1;
  }

  /************************************************************************
    Set the filename of the surface file : 
  *************************************************************************/ 
#if( USE_MPI_IO_SURF || (!USEMPI) ) 
  sprintf(surf_filename, "%s/%s.surface.%06d.h5", DIR_out[OUT_SURFACE],RUN_TAG,N_out[OUT_SURFACE]);
#else
  sprintf(surf_filename, "%s/%s.surface.%06d.i%05d.k%05d.h5", DIR_out[OUT_SURFACE],RUN_TAG,N_out[OUT_SURFACE],cpupos[1],cpupos[3]);
#endif

  fprintf(stdout,"Dumping %s   [pid=%d].... \n",surf_filename,myid); fflush(stdout);

  /*********************************************************************************
    Create the file and subgroups in the root group if this is a new run:
  **********************************************************************************/ 
  file_id = H5Fcreate(surf_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if( file_id  < 0 ) { 
    fprintf(stderr,"write_surface_data(): Cannot create file %s \n", surf_filename);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Write header data : */ 
  write_dump_header(file_id,surf_filename,1);

  /* Make the groups or the subdirectories : */
  for(i=0; i<N_SURF_TYPES; i++ ) { 
    group_id = H5Gcreate(file_id, surf_dirnames[i], 0);
    if( group_id < 0 ) { 
      fprintf(stderr,"write_surface_data(): Failed to create group %s  !! \n", surf_dirnames[i]);
      fflush(stderr);
      fail(FAIL_HDF,0);
    }
    H5Gclose(group_id); 
  }

  /* Write the functions : */ 
  for(i=0; i<N_SURFACE; i++)  for(j=0; j<N_SURF_TYPES; j++) { 
      sprintf(dataname,"/%s/%s",surf_dirnames[j],S_names[i]);
      myH5_write_surffunc( file_id, dataname, surf_data[i][j] );
  }

  if( myid == out_pid[OUT_SURFACE] ) { 
    /* Write  "r_bins" : */ 
    sprintf(dataname,"/r");
    myH5_write_array( file_id, dataname, H5T_NATIVE_DOUBLE, n_r_bins, r_bins_mid );

    sprintf(dataname,"/dr");
    myH5_write_array( file_id, dataname, H5T_NATIVE_DOUBLE, n_r_bins, dr_bins );

    /* Write  "phi_bins" : */ 
    sprintf(dataname,"/phi");
    myH5_write_array( file_id, dataname, H5T_NATIVE_DOUBLE, n_phi_bins, phi_bins_mid );

    sprintf(dataname,"/dphi");
    myH5_write_array( file_id, dataname, H5T_NATIVE_DOUBLE, 1, &dphi_bins );
  }

  
  H5Fclose(file_id);

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_surffunc(): 
 ---------------------------
   -- driver routine for writing a grid function within a given group/object
      "loc_id" for the surface data dumps; 
   -- H5T_NATIVE_DOUBLE  is the assumed datatype of all grid functions ;
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
*******************************************************************************************/
void myH5_write_surffunc( hid_t loc_id, char *name, double *value )
{
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;

  //  fprintf(stdout,"myH5_write_surffunc(): BEG    %s  \n ",name); fflush(stdout);

  /************************************************************************************
     Create the data set: 
   ************************************************************************************/
  /* Create the data space for the dataset, always using 3 dimensions since we want 
     filespace to match memory space : */
  filespace_id = H5Screate_simple(SURF_RANK, surf_filedims, NULL); 
  memspace_id  = H5Screate_simple(SURF_RANK, surf_memdims , NULL); 

  /* Set the property list to for creation */ 
  prop_id = H5Pcreate(H5P_DATASET_CREATE);

#if( USE_CHUNKS ) 
  /* Set properties to use chunks, set chunk's extent : */
  H5Pset_chunk(prop_id, SURF_RANK, surf_memdims);
#endif

  /* Create a dataset assuming double var type */
  dataset_id = H5Dcreate(loc_id, name, H5T_NATIVE_DOUBLE, filespace_id, prop_id);
  H5Pclose(prop_id);   
  H5Sclose(filespace_id);


  /************************************************************************************
    Select the hyperslab, or section of space, to which our gridfunction will be written
       -- I am not sure why we need to close and then reget the dataset's dataspace, 
          but that is what they usually do in the example programs in the hdf5 tutorial;
  ************************************************************************************/
  filespace_id = H5Dget_space(dataset_id);  

#if( USE_CHUNKS ) 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, surf_memdims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, surf_memdims, NULL);
#endif 
  
  
  /************************************************************************************
    Write the hyperslab to the dataset : 
  ************************************************************************************/
  /* Setup the properties of the write based upon type of run : */
  prop_id = H5Pcreate(H5P_DATASET_XFER);

#if( USE_MPI_IO ) 
  H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
#endif

  /* Write the dataset using defined dataspace and properties. */
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, prop_id, value);

  /* Close the open object handles : */
  H5Dclose(dataset_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
  H5Pclose(prop_id);

  //  fprintf(stdout,"myH5_write_surffunc(): END    %s  \n ",name); fflush(stdout);

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_array(): 
 ---------------------------
   -- driver routine for writing a general 1d array under group/file id "loc_id";  
      the array is identified by the name "name", and the array
      has type "type_id"; 
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
*******************************************************************************************/
void myH5_write_array( hid_t loc_id, char *name, hid_t type_id, int size, void *value ) 
{
  hid_t   filespace_id, memspace_id, dataset_id, prop_id;
  hsize_t dims[1];

  dims[0] = size; 

  /* Create the dataspace for the dataset, single number format:  */
  memspace_id  = H5Screate_simple(1, dims, NULL); 
  filespace_id = H5Screate_simple(1, dims, NULL); 

  /* Create a dataset, which is like a small dataset  */
  dataset_id = H5Dcreate(loc_id, name, type_id, filespace_id, H5P_DEFAULT);

  /* Write the dataset using defined dataspace and default properties. */
  H5Dwrite(dataset_id, type_id, memspace_id, filespace_id, H5P_DEFAULT, value);

  /* Close all open object handles : */
  H5Dclose(dataset_id);  
  H5Sclose(filespace_id);        
  H5Sclose(memspace_id);         

  return;
}



#else
void dump_surface_gen(void) { return; } 
#endif


