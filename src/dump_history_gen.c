
#include "decs.h"


#if( MAKE_HISTORY && MAKE_HDF5 ) 
#include <hdf5.h>

#if( USEMPI )
#include "mpi.h"
extern void history_mpi_gen( ulint n_used_bins,  usint *used_bins, double ***hist_binned );
extern void history_mpi2_gen( double ***hist_binned, double *hist_data_out[N_HISTORY][N_HIST_TYPES] ); 
#endif 

/* Rank of history datasets (time versus X1): */
#define HIST_RANK (2) 

/* Extents of memory and file spaces */
static hsize_t hist_memdims[HIST_RANK],hist_filedims[HIST_RANK],hist_max_filedims[HIST_RANK];
static hsize_t offset[HIST_RANK];


/* Name of history functions :  (has length  N_HISTORY+2  because of t and dt) */
static char *H_names[N_HISTORY+2] = { 
  "H_bsq"  , 
  "H_pavg" , 
  "H_rhoav", 
  "H_ddot" , 
  "H_hUz"  , 
  "H_hUt"  , 
  "H_hUx"  , 
  "H_lrho" , 
  "H_Txza" , 
  "H_Txta" , 
  "H_Txxa" , 
  "H_Ttta" , 
  "H_Txzb" , 
  "H_Txtb" , 
  "H_Txxb" , 
  "H_Tttb" , 
  "H_Txzc" , 
  "H_Txtc" , 
  "H_Txxc" , 
  "H_Tttc" , 
  "H_Txzd" , 
  "H_Txtd" , 
  "H_Txxd" , 
  "H_Tttd" , 
  "H_vol"  
#if( CALC_CURRENT )
 ,"H_Jx"  ,      
  "H_Jy"  ,      
  "H_Jz"  ,      
  "H_Jt"  ,      
  "H_Jsq"  
#endif
#if( MAKE_RADFLUX )
 ,"H_Lut"  ,      
  "H_Lut2" ,      
  "H_Lup"        
#endif
#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
 ,"H_dTdr_1",     
  "H_dTdr_2",     
  "H_dTdr_3",     
  "H_dTdr_4",     
  "H_dTdr_5",     
  "H_dTdr_6"     
#endif
  ,"H_t"    , "H_dt"  
};

enum hist_func_index  {  H_bsq   ,H_pavg  ,H_rhoav ,H_ddot  ,H_hUz   ,H_hUt   ,H_hUx   ,H_lrho  
			,H_Txza  ,H_Txta  ,H_Txxa  ,H_Ttta  
			,H_Txzb  ,H_Txtb  ,H_Txxb  ,H_Tttb  
			,H_Txzc  ,H_Txtc  ,H_Txxc  ,H_Tttc  
			,H_Txzd  ,H_Txtd  ,H_Txxd  ,H_Tttd  
                        ,H_vol   
#if( CALC_CURRENT )
			,H_Jx, H_Jy, H_Jz, H_Jt, H_Jsq
#endif
#if( MAKE_RADFLUX )
			,H_Lut, H_Lut2, H_Lup
#endif
#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
			,H_dTdr_1,H_dTdr_2,H_dTdr_3,H_dTdr_4,H_dTdr_5,H_dTdr_6
#endif
			,H_last                    /* used to verify length of surf_func_index */
};

/* Names of hdf groups under which histories live */
static char *hist_dirnames[] = { "Unbound" , "Bound", "Jet" };  

/* Name of history file */
static char hist_filename[200];   


/* Routines : */
static void setup_history( void ) ;
static void myH5_write_histfunc(hid_t loc_id, int new_func, char *name, double *value );
static void myH5_write_histscalar(hid_t loc_id,int new_scalar,char *name,hid_t type_id,void *value);

extern hsize_t *myH5_get_simple_dims( hid_t dataspace_id,  int *ndims );

/**********************************************************************************************/
/**********************************************************************************************
  dump_history_gen(): 
 -------------
   -- "generalized" version of dump_history() that bins data on a static radial coordinate system;
       -- most useful for non-spherical or dynamic coordinate systems;

   -- routine that outputs the "history" information;
   -- this is frequent output of usually shell-integrated quantities; 
   -- The below quanities use HARM notation.  Any translation between HARM 
       definitions and Hawley's definitions should be done at the IDL/analysis 
       level.

   -- we do not use primtoflux() to calculate stress terms since we want to separate
       the parts out ;  also, we want to calculate the components in KS coordinates not 
       code coordinates;

   -- The choice in nomenclature w.r.t. "x,y,z" was done to be consistent with 
      Beckwith's nomenclature.

   -- note that the hist_data[][][] array is ordered that way so that it is easier to 
        pass it to  history_mpi() ;



LEGEND: 

Group : Name  : Quantity      (x=r, y=theta , z=phi)
------------------------------------------------
/  : 
        t    : time of this history slice
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
               Tttd    : ||b^2|| u^t u_t  - b^t b_t
               Txtd    : ||b^2|| u^r u_t  - b^r b_t
               Txxd    : ||b^2|| u^r u_r  - b^r b_r
               Txzd    : ||b^2|| u^r u_\phi  - b^r b_\phi
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
      -- the datasets in the history file have dimensions [n_r_bins,n_hist_dump_tot] 
      -- hence, the order of operations is the following:
            if new: 
                -- create file, create dataset structures
            if old: 
                -- open file, open dataset, get dimensions,
                   write the data, and close things ;

**********************************************************************************************/
//void dump_history_gen( void )
void dump_history( void )
{
  int i,j,k,l,id,itype,ii,jj,ic,kk,pos ;
  int set_type[N_HIST_TYPES];

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

  char dataname[200];

  hid_t file_id, prop_id;

  double ***hist_data_binned;
  double *r_beg, *r_end, *r_bin_mid; 
  double weight_loc; 
  ulint icell, ibin, ibins_b, ibins_e, ifunc, n_used_bins;
  usint *used_bins;

  static double hist_data_loc[N_HISTORY];

  static int first_time = 1; 

  TRACE_BEG;

  if( N_HISTORY != H_last ) { 
    fprintf(stderr,"dump_history_gen():  Invalid values of N_HISTORY and H_last:  %d %d \n",N_HISTORY,H_last); 
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  /******************************************************************************
    Initialize function names, hdf5 file, group structure, etc. 
      -- do not set first_time = 0  until the end of the routine; 
  ******************************************************************************/ 
  if( first_time ) {  

    if( myid == out_pid[OUT_HISTORY] ) { 
      setup_history() ;
    }
    if( using_restart ) { first_time = 0; }
  }

  if( myid == printer_pid ) {  fprintf(stdout,"Doing history dump (%d) %d ... \n", myid,N_hist_dump);  fflush(stdout); }

  /*********************************************************************************** 
     First reset the history data so we can begin tallying : 
  ***********************************************************************************/ 
  if( myid == out_pid[OUT_HISTORY] ) { 
    for(i=0;i<N_HISTORY;i++) for(j=0;j<N_HIST_TYPES;j++) { ALLOC_ARRAY(hist_data_out[i][j]   ,n_r_bins); }
    for(i=0;i<N_HISTORY;i++) for(j=0;j<N_HIST_TYPES;j++)  for(k=0; k<n_r_bins; k++) { hist_data_out[i][j][k] = 0.; }
  }

  /*********************************************************************************** 
     Setup bin arrays and weights:
  ***********************************************************************************/ 
  ALLOC_ARRAY(r_beg,NCELLS); 
  ALLOC_ARRAY(r_end,NCELLS); 

  icell = 0; 

  LOOP { 
#if( TOP_TYPE_CHOICE == TOP_CARTESIAN )
    get_coord(i  ,j  ,k  ,CORN,ncurr,coords_m1);
    get_coord(i+1,j+1,k+1,CORN,ncurr,coords_p1);
    r_beg[icell] = coords_m1->r; 
    r_end[icell] = coords_p1->r; 
#else
    get_coord(i  ,j,k,FACE1,ncurr,coords_m1);
    get_coord(i+1,j,k,FACE1,ncurr,coords_p1);
    r_beg[icell] = coords_m1->r; 
    r_end[icell] = coords_p1->r; 
#endif
    icell++;
  }

  set_bin_weights_r(NCELLS,r_beg,r_end,&bin_weights_r,&ibins_beg,&ibins_end); 

  for(icell=0;icell<NCELLS;icell++) { r_beg[icell] = r_end[icell] - r_beg[icell]; }   /* reassign r_beg to be dr_cell */

  FREE(r_end);

  ALLOC_ARRAY(used_bins,n_r_bins); 
  for(ibin=0;ibin<n_r_bins;ibin++) { used_bins[ibin] = 0; } 
  for(icell=0;icell<NCELLS;icell++) { 
    ibins_b = ibins_beg[icell];     
    ibins_e = ibins_end[icell]; 
    for(ibin=ibins_b;ibin<=ibins_e;ibin++) { 
      used_bins[ibin] = 1; 
    }
  }

  n_used_bins = 0;
  for(ibin=0;ibin<n_r_bins;ibin++) if( used_bins[ibin] ) { n_used_bins++ ; }

  ALLOC_ARRAY(hist_data_binned,n_r_bins); 
  for(ibin=0;ibin<n_r_bins;ibin++) { 
    if( used_bins[ibin] ) { 
      ALLOC_2D_ARRAY(hist_data_binned[ibin],N_HIST_TYPES,N_HISTORY); 
      for( itype = 0; itype < N_HIST_TYPES ; itype++ )  for(ifunc=0; ifunc<N_HISTORY; ifunc++) { 
	  hist_data_binned[ibin][itype][ifunc] = 0. ; 
	}
    }
    else { 
      hist_data_binned[ibin] = NULL;
    }
  }

  /*********************************************************************************** 
     Begin the loop over all space, integrating along x2 and x3 directions : 
  ***********************************************************************************/
  icell = 0; 

  LOOP { 
    get_geometry(i  ,j,k,CENT,ncurr,geom);  
    get_coord(   i  ,j,k,CENT,ncurr,coords);
    
    dV_prop = dV * geom->g / r_beg[icell] ; 

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
       -- make sure the order of set_type[] corresponds to order of hist_dirnames[] */
    set_type[UNBOUND] = ((enthalpy+bsq)*ucov_ks[TT]/rho < -1.);           /* Unbound */
    set_type[  BOUND] = !set_type[UNBOUND];                         /* Bound   */
    set_type[    JET] = set_type[UNBOUND]  && (ucon_ks[RR] > 0.) ;  /* Jet     */

    /* Add to the history integration : */
    for( itype = 0; itype < N_HIST_TYPES ; itype++ )  if( set_type[itype] )  { 
	hist_data_loc[H_bsq   ]  =  (bsq);
	hist_data_loc[H_pavg  ]  =  ((gam-1.)*ptmp[UU]); 
	hist_data_loc[H_rhoav ]  =  (rho); 
	hist_data_loc[H_ddot  ]  =  (rho      *  ucon_ks[RR]); 
	hist_data_loc[H_hUt   ]  =  (enthalpy *  ucon_ks[TT]); 
	hist_data_loc[H_hUx   ]  =  (enthalpy *  ucon_ks[RR]); 
	hist_data_loc[H_hUz   ]  =  (enthalpy *  ucon_ks[PH]); 
	hist_data_loc[H_lrho  ]  =  (-rho     * (ucov_ks[PH] / ucov_ks[TT])); 
	hist_data_loc[H_Ttta  ]  =  (enthalpy *  ucon_ks[TT] * ucov_ks[TT] ); 
	hist_data_loc[H_Txta  ]  =  (enthalpy *  ucon_ks[RR] * ucov_ks[TT] ); 
	hist_data_loc[H_Txxa  ]  =  (enthalpy *  ucon_ks[RR] * ucov_ks[RR] ); 
	hist_data_loc[H_Txza  ]  =  (enthalpy *  ucon_ks[RR] * ucov_ks[PH] ); 
	hist_data_loc[H_Tttb  ]  =  (bsq      *  ucon_ks[TT] * ucov_ks[TT] ); 
	hist_data_loc[H_Txtb  ]  =  (bsq      *  ucon_ks[RR] * ucov_ks[TT] ); 
	hist_data_loc[H_Txxb  ]  =  (bsq      *  ucon_ks[RR] * ucov_ks[RR] ); 
	hist_data_loc[H_Txzb  ]  =  (bsq      *  ucon_ks[RR] * ucov_ks[PH] ); 
	hist_data_loc[H_Tttc  ]  =  (           -bcon_ks[TT] * bcov_ks[TT] ); 
	hist_data_loc[H_Txtc  ]  =  (           -bcon_ks[RR] * bcov_ks[TT] ); 
	hist_data_loc[H_Txxc  ]  =  (           -bcon_ks[RR] * bcov_ks[RR] ); 
	hist_data_loc[H_Txzc  ]  =  (           -bcon_ks[RR] * bcov_ks[PH] ); 
	hist_data_loc[H_Tttd  ]  =  (bsq      *  ucon_ks[TT] * ucov_ks[TT] -bcon_ks[TT] * bcov_ks[TT] ); 
	hist_data_loc[H_Txtd  ]  =  (bsq      *  ucon_ks[RR] * ucov_ks[TT] -bcon_ks[RR] * bcov_ks[TT] ); 
	hist_data_loc[H_Txxd  ]  =  (bsq      *  ucon_ks[RR] * ucov_ks[RR] -bcon_ks[RR] * bcov_ks[RR] ); 
	hist_data_loc[H_Txzd  ]  =  (bsq      *  ucon_ks[RR] * ucov_ks[PH] -bcon_ks[RR] * bcov_ks[PH] ); 
	hist_data_loc[H_vol   ]  =  1. ;
#if( CALC_CURRENT )
	hist_data_loc[H_Jt    ]  =  (jcon_ks[TT]); 
	hist_data_loc[H_Jx    ]  =  (jcon_ks[RR]); 
	hist_data_loc[H_Jy    ]  =  (jcon_ks[TH]); 
	hist_data_loc[H_Jz    ]  =  (jcon_ks[PH]); 
	hist_data_loc[H_Jsq   ]  =  (jsq); 
#endif

#if( MAKE_RADFLUX )
# if( USE_COOLING_FUNCTION == 1 ) 
	coolrate = cooling_func_hr_disk( i,j,k, ptmp );
# elif( USE_COOLING_FUNCTION == 2 || USE_COOLING_FUNCTION == 3 ) 
	coolrate = cooling_func_isentropic_disk( i, j, k, ptmp, &q ); //Dennis
# elif( USE_COOLING_FUNCTION >= 4 )
	  coolrate = coolflux[0][k-N3S + N3*((j-N2S) + N2*(i-N1S))];
# endif
	hist_data_loc[H_Lut   ]  =  (coolrate*ucov_ks[TT]); 
	hist_data_loc[H_Lut2  ]  =  (coolrate*ucov_ks[TT]/ucon_ks[TT]); 
	hist_data_loc[H_Lup   ]  =  (coolrate*ucov_ks[PH]);
#endif

#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
	hist_data_loc[H_dTdr_1]  =  (rho                * uu_dot_conn);
	hist_data_loc[H_dTdr_2]  =  (enthalpy           * uu_dot_conn);
	hist_data_loc[H_dTdr_3]  =  (bsq                * uu_dot_conn);
	hist_data_loc[H_dTdr_4]  =  ((gam-1.)*ptmp[UU]  * conn_trace );
	hist_data_loc[H_dTdr_5]  =  (0.5*bsq            * conn_trace );
	hist_data_loc[H_dTdr_6]  =  (-bb_dot_conn);
#endif
	ibins_b = ibins_beg[icell]; 	
	ibins_e = ibins_end[icell]; 
	jj = 0; 
	for(ibin=ibins_b; ibin<=ibins_e; ibin++) {
	  weight_loc = dV_prop * bin_weights_r[icell][jj];
	  for(ifunc=0;ifunc<N_HISTORY;ifunc++) {
	    hist_data_binned[ibin][itype][ifunc] += weight_loc*hist_data_loc[ifunc] ; 
	  }
	  jj++;
	}
      }
    icell++;
  }

  if( bin_weights_r != NULL ) {  for(icell=0;icell<NCELLS;icell++) { FREE(bin_weights_r[icell]); }  }
  FREE( bin_weights_r  ); 
  FREE( ibins_beg );   FREE( ibins_end ); 
  FREE(r_beg);

  for(k=0; k<n_r_bins; k++) {
    if( used_bins[k] ) { 
      jsq = dr_bins[k]; 
      for(j=0; j<N_HIST_TYPES; j++) for(i=0; i<N_HISTORY; i++)  {
	  hist_data_binned[k][j][i] *= jsq;
	}
    }
  }

  /*********************************************************************************** 
    If we are running in parallel, then sum over the grids in x2,x3 directions : 
  ***********************************************************************************/ 
#if( USEMPI ) 
  history_mpi_gen( n_used_bins, used_bins, hist_data_binned ); 
  history_mpi2_gen( hist_data_binned, hist_data_out ); 
#else 
  for(k=0; k<n_r_bins; k++) { 
    if( used_bins[k] ) { 
      for(j=0; j<N_HIST_TYPES; j++)  for(i=0; i<N_HISTORY; i++)  {
	  hist_data_out[i][j][k] = hist_data_binned[k][j][i]; 
	}
    }
  }
#endif 

  for(ibin=0;ibin<n_r_bins;ibin++) { 
    if( hist_data_binned[ibin] != NULL ) { 
      for( itype = 0; itype < N_HIST_TYPES ; itype++ )  { FREE(hist_data_binned[ibin][itype]);  }
      FREE(hist_data_binned[ibin]); 
    }
  }
  FREE(hist_data_binned);
  FREE(used_bins);

  /*********************************************************************************** 
    Append the data to the history file : open, write, and close: 
  ***********************************************************************************/ 
  if( myid == out_pid[OUT_HISTORY] ) { 

    prop_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fopen(hist_filename, H5F_ACC_RDWR, prop_id);
    if( file_id  < 0 ) { 
      fprintf(stderr,"Cannot open file %s \n", hist_filename);
      fflush(stderr);     fail(FAIL_HDF,0);
    }
    H5Pclose(prop_id);    


    /* Write  "t" : */ 
    myH5_write_histscalar(file_id, first_time, H_names[N_HISTORY],H5T_NATIVE_DOUBLE, &t);

    /* Write  "dt" : */ 
    myH5_write_histscalar(file_id, first_time, H_names[N_HISTORY+1],H5T_NATIVE_DOUBLE, &(dx[0]));

    /* Write  "r_bins" : */ 
    sprintf(dataname,"/r");
    myH5_write_histfunc( file_id, first_time, dataname, r_bins_mid );

    sprintf(dataname,"/dr");
    myH5_write_histfunc( file_id, first_time, dataname, dr_bins );

    /* Write the functions : */ 
    for(i=0; i<N_HISTORY; i++)  for(j=0; j<N_HIST_TYPES; j++) { 
      sprintf(dataname,"/%s/%s",hist_dirnames[j],H_names[i]);
      myH5_write_histfunc( file_id, first_time, dataname, hist_data_out[i][j] );
    }
    
    /* Close hdf objects : */
    H5Fclose(file_id);  
  }

  for(i=0; i<N_HISTORY; i++)  { 
    for(j=0; j<N_HIST_TYPES; j++) { FREE( hist_data_out[i][j] ); }
  }

  /* Make sure that it is no longer our first time */
  first_time = 0; 

  N_hist_dump++;

#if( USEMPI ) 
  //  exit_status();
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 

  TRACE_END;

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 setup_history():
 ---------------------------
   -- initializes arrays and file used for all history file; 
   -- note that a serial run will be automatically handled;
   -- does not use extended datasets;
   -- uses GridLength[0] and DT_out[OUT_HISTORY] to determine the length of the dataset;
         -- uses GridLength[0] / DT_out[OUT_HISTORY]  + 5 
         -- the 5 extra sets are to accomodate for any unforeseen extra calls to dump_history_gen()
              like the OUT_INIT and OUT_FINAL
*******************************************************************************************/
static void setup_history( void ) 
{
  int i,j,k ;
  hid_t file_id, prop_id;
  hid_t group_id[N_HIST_TYPES];

  /************************************************************************
    Set the filename of the history file : 
  *************************************************************************/ 
#if( USE_MPI_IO_HIST ) 
  sprintf(hist_filename, "%s/%s.history_gen.h5", DIR_out[OUT_HISTORY],RUN_TAG);
#else
  sprintf(hist_filename, "%s/%s.history_gen.p%05d.h5", DIR_out[OUT_HISTORY],RUN_TAG,myid);
#endif

  fprintf(stdout,"\n######################################\n");
  fprintf(stdout,"History file name = %s (%d) \n",hist_filename,myid);
  fprintf(stdout,"####################################\n\n");
  fflush(stdout);

  /************************************************************************
    Determine the number of history dumps we may perform: 
        -- adds a few extra so that we do not overflow the allocated limit;
  *************************************************************************/ 
  if( !using_restart ) { 
    N_hist_dump_tot = 5  + ((int)  (GridLength[0] / DT_out[OUT_HISTORY]) );
  }

  /************************************************************************
    Set various static arrays used in writing to the history file : 
  *************************************************************************/ 
  /* Setup the extent of each slice of a history dataset in memory : */
  hist_memdims[    0] = 1;                   hist_memdims[    1] = n_r_bins;
  hist_filedims[   0] = 1            ;       hist_filedims[   1] = n_r_bins;
  hist_max_filedims[0] = N_hist_dump_tot;   hist_max_filedims[1] = hist_filedims[1];

  offset[0] = 0;                       /* We do not really know this yet */
  offset[1] = 0;

  /*********************************************************************************
    Create the file and subgroups in the root group if this is a new run:
  **********************************************************************************/ 
  if( !using_restart ) { 
    prop_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate(hist_filename, H5F_ACC_EXCL, H5P_DEFAULT, prop_id);
    if( file_id  < 0 ) { 
      fprintf(stderr,"Cannot create file %s \n", hist_filename);
      fflush(stderr);     fail(FAIL_HDF,0);
    }
    H5Pclose(prop_id);

    /* Make the groups or the subdirectories : */
    for(i=0; i<N_HIST_TYPES; i++ ) { 
      group_id[i] = H5Gcreate(file_id, hist_dirnames[i], 0);
      if( group_id[i] < 0 ) { 
	fprintf(stderr,"setup_history(): Failed to create group %s  !! \n", hist_dirnames[i]);
	fflush(stderr);
	fail(FAIL_HDF,0);
      }
    }

    /* After making them, close them : */
    for(i=0; i<N_HIST_TYPES; i++ ) { 
      H5Gclose(group_id[i]); 
    }
    H5Fclose(file_id);
  }


  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_histscalar(): 
 ---------------------------
   -- driver routine for writing a scalar in the history file specified by "loc_id"; 
   -- uses datasets of fixed dimensions
   -- this scalar is really a dataset  f[NT]  where NT is the number of timeslices in 
   -- H5T_NATIVE_DOUBLE  is the assumed datatype of all grid functions ;
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between MPI and serial calls; 
       -- makes sure that only the master node writes the scalar data;
   -- if new_scalar != 0 , then assume that this function has not been created yet; 
   -- I cannot believe that it takes this much code to add one number to a data column!!!
*******************************************************************************************/
static void myH5_write_histscalar(hid_t loc_id, int new_scalar, char *name, 
			   hid_t type_id, void *value )
{
  int ndims;
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;
  hsize_t dims[1], maxdims[1],  *fulldims, count[1]; 
  hsize_t scalar_offset[1];

  count[0] = dims[0]    = 1; 
  maxdims[0] = N_hist_dump_tot;

  /*  Either create or open/extend the dataset so that it is ready for writing new data */
  if( new_scalar ) { 
    prop_id      = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(prop_id, 1, dims);
    filespace_id = H5Screate_simple(1, maxdims, NULL); 
    dataset_id = H5Dcreate(loc_id, name, type_id, filespace_id, prop_id);
    H5Pclose(prop_id);  H5Sclose(filespace_id);    
    scalar_offset[0] = 0;   /* start at the beginning */
  }
  else { 
    dataset_id = H5Dopen(loc_id, name);   /* get dataset */
    filespace_id = H5Dget_space(dataset_id);  /* get dataspace of dataset */
    scalar_offset[0] = N_hist_dump;  /* set the offset to the next spot in dataset */
  }

  memspace_id  = H5Screate_simple(1, dims, NULL); 
  filespace_id = H5Dget_space(dataset_id);  /* get new dataspace  */
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, scalar_offset, NULL, count, dims);  

  prop_id = H5Pcreate(H5P_DATASET_XFER);

  H5Dwrite(dataset_id, type_id, memspace_id, filespace_id, prop_id, value);

  /* Close all open object handles : */
  H5Dclose(dataset_id);	 
  H5Sclose(filespace_id);	 
  H5Sclose(memspace_id);	 
  H5Pclose(prop_id);
  
  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_histfunc(): 
 ---------------------------
   -- driver routine for writing a function n the history file specified by "loc_id"; 
   -- this routine uses datasets of fixed length
   -- this function  is really a dataset  f[NT][N1]  where NT is the number of timeslices 
   -- H5T_NATIVE_DOUBLE  is the assumed datatype of all grid functions ;
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between MPI and serial calls; 
       -- checks cpupos[] to see if the cpu should write to history file; 
       -- all processors with  cpupos[2] = cpupos[3] = 0  need to write to history file;
   -- if new_func != 0 , then assume that this function has not been created yet; 
   -- I cannot believe that it takes this much code to add one number to a data column!!!
*******************************************************************************************/
static void myH5_write_histfunc(hid_t loc_id, int new_func, char *name, double *value )
{
  int ndims,i;
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;
  hsize_t *fulldims, count[2];

  for( i=0; i < HIST_RANK; i++ ) {  count[i] = 1; } 

  /*  Either create or open/extend the dataset so that it is ready for writing new data */
  if( new_func ) { 
    prop_id      = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(prop_id, HIST_RANK, hist_memdims);
    filespace_id = H5Screate_simple(HIST_RANK, hist_max_filedims, NULL); 
    dataset_id = H5Dcreate(loc_id, name, H5T_NATIVE_DOUBLE, filespace_id, prop_id);
    H5Pclose(prop_id);  H5Sclose(filespace_id);    
    /* Never touch offset[1], just change offset[0] as array gets longer */
    offset[0] = 0;   /* start at the beginning */
  }
  else { 
    dataset_id = H5Dopen(loc_id, name);   /* get dataset */
    filespace_id = H5Dget_space(dataset_id);  /* get dataspace of dataset */
    fulldims = myH5_get_simple_dims( filespace_id,  &ndims ); /* get extent of dataspace */
    H5Sclose(filespace_id);	 
    if( fulldims[0] != N_hist_dump_tot ) { 
      fprintf(stderr,"myH5_write_histfunc(): dimensions offset, %d %d ", 
	      ((int) fulldims[0]), N_hist_dump_tot); 
      fflush(stderr); H5Dclose(dataset_id);	 
      fail(FAIL_HDF,0);
    }
    offset[0] = N_hist_dump;  /* set the offset to the next spot in dataset */
  }

  memspace_id  = H5Screate_simple(HIST_RANK, hist_memdims, NULL); 
  filespace_id = H5Dget_space(dataset_id);  /* get new dataspace  */
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, hist_memdims);  

  prop_id = H5Pcreate(H5P_DATASET_XFER);

  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, prop_id, value);

  /* Close all open object handles : */
  H5Dclose(dataset_id);	 
  H5Sclose(filespace_id);	 
  H5Sclose(memspace_id);	 
  H5Pclose(prop_id);

  return;
}


#else
void dump_history_gen(void) { return; } 
#endif


