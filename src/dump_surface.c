
#include "decs.h"


#if( MAKE_SURFACE && MAKE_HDF5 ) 
#include <hdf5.h>

#if( USEMPI )
#include "mpi.h"
extern void surface_mpi( double *surf_arr[N_SURFACE][N_SURF_TYPES] );
extern void surface_mpi_faster( double *surf_arr[N_SURFACE][N_SURF_TYPES] );
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


/**********************************************************************************************/
/**********************************************************************************************
  dump_surface(): 
 -------------
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
void dump_surface( void )
{
  int i,j,k,l,id,itype,ii,jj,ic,kk,pos;
  int set_type[N_SURF_TYPES];

  double *ptmp, dV_prop;
  double bcon_ks[NDIM],bcov_ks[NDIM],ucon_ks[NDIM],ucov_ks[NDIM];
  double bsq, enthalpy, rho, jcon_i[NDIM],jcon_ks[NDIM], jsq;
  double coolrate=0.;
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

  TRACE_BEG;

  if( N_SURFACE != S_last ) { 
    fprintf(stderr,"dump_surface():  Invalid values of N_SURFACE and S_last:  %d %d \n",N_SURFACE,S_last); 
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  /*   fprintf(stdout,"Doing surface dump (%d) %d ... \n", myid,N_out[OUT_SURFACE]);   fflush(stdout);  */

  /*********************************************************************************** 
     First reset the surface data so we can begin tallying : 
  ***********************************************************************************/ 
  l = N1*N3;
  for(i=0; i<N_SURFACE; i++) for(j=0; j<N_SURF_TYPES; j++) for(k=0; k<l; k++)  { surf_data[i][j][k] = 0.; }

  /*********************************************************************************** 
     Begin the loop over all space, integrating along x2 and x3 directions : 
  ***********************************************************************************/ 
  LOOP  {
    id = (i - N1S)*N3 + (k-N3S);

    get_geometry(i,j,k,CENT,ncurr,geom);  
    get_coord(   i,j,k,CENT,ncurr,coords);
    dV_prop = dV * geom->g ; 

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
	  surf_data[S_bsq  ][itype][id] +=  dV_prop * (bsq);
	  surf_data[S_pavg ][itype][id] +=  dV_prop * ((gam-1.)*ptmp[UU]); 
	  surf_data[S_rhoav][itype][id] +=  dV_prop * (rho); 
	  surf_data[S_ddot ][itype][id] +=  dV_prop * (rho      *  ucon_ks[RR]); 
	  surf_data[S_hUt  ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[TT]); 
	  surf_data[S_hUx  ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[RR]); 
	  surf_data[S_hUz  ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[PH]); 
	  surf_data[S_lrho ][itype][id] +=  dV_prop * (-rho     * (ucov_ks[PH] / ucov_ks[TT])); 
	  surf_data[S_Ttta ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[TT] * ucov_ks[TT] ); 
	  surf_data[S_Txta ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[RR] * ucov_ks[TT] ); 
	  surf_data[S_Txxa ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[RR] * ucov_ks[RR] ); 
	  surf_data[S_Txza ][itype][id] +=  dV_prop * (enthalpy *  ucon_ks[RR] * ucov_ks[PH] ); 
	  surf_data[S_Tttb ][itype][id] +=  dV_prop * (bsq      *  ucon_ks[TT] * ucov_ks[TT] ); 
	  surf_data[S_Txtb ][itype][id] +=  dV_prop * (bsq      *  ucon_ks[RR] * ucov_ks[TT] ); 
	  surf_data[S_Txxb ][itype][id] +=  dV_prop * (bsq      *  ucon_ks[RR] * ucov_ks[RR] ); 
	  surf_data[S_Txzb ][itype][id] +=  dV_prop * (bsq      *  ucon_ks[RR] * ucov_ks[PH] ); 
	  surf_data[S_Tttc ][itype][id] +=  dV_prop * (           -bcon_ks[TT] * bcov_ks[TT] ); 
	  surf_data[S_Txtc ][itype][id] +=  dV_prop * (           -bcon_ks[RR] * bcov_ks[TT] ); 
	  surf_data[S_Txxc ][itype][id] +=  dV_prop * (           -bcon_ks[RR] * bcov_ks[RR] ); 
	  surf_data[S_Txzc ][itype][id] +=  dV_prop * (           -bcon_ks[RR] * bcov_ks[PH] ); 
	  surf_data[S_vol  ][itype][id] +=  dV_prop ;

#if( CALC_CURRENT )
	surf_data[S_Jt   ][itype][id] +=  dV_prop * (jcon_ks[TT]); 
	surf_data[S_Jx   ][itype][id] +=  dV_prop * (jcon_ks[RR]); 
	surf_data[S_Jy   ][itype][id] +=  dV_prop * (jcon_ks[TH]); 
	surf_data[S_Jz   ][itype][id] +=  dV_prop * (jcon_ks[PH]); 
	surf_data[S_Jsq  ][itype][id] +=  dV_prop * (jsq); 
#endif


#if( MAKE_RADFLUX )
# if( USE_COOLING_FUNCTION == 1 ) 
	coolrate = cooling_func_hr_disk( i,j,k, ptmp );
# elif( USE_COOLING_FUNCTION == 2 || USE_COOLING_FUNCTION == 3 ) 
	coolrate = cooling_func_isentropic_disk( i, j, k, ptmp, &q ); //Dennis
# elif( USE_COOLING_FUNCTION >= 4 )
	  coolrate = coolflux[0][k-N3S + N3*((j-N2S) + N2*(i-N1S))];
# endif
        surf_data[S_Lut  ][itype][id] +=  dV_prop * (coolrate*ucov_ks[TT]); 
	surf_data[S_Lut2 ][itype][id] +=  dV_prop * (coolrate*ucov_ks[TT]/ucon_ks[TT]); 
	surf_data[S_Lup  ][itype][id] +=  dV_prop * (coolrate*ucov_ks[PH]);
#endif

#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
	surf_data[S_dTdr_1][itype][id] +=  dV_prop * (rho                * uu_dot_conn);
	surf_data[S_dTdr_2][itype][id] +=  dV_prop * (enthalpy           * uu_dot_conn);
	surf_data[S_dTdr_3][itype][id] +=  dV_prop * (bsq                * uu_dot_conn);
	surf_data[S_dTdr_4][itype][id] +=  dV_prop * ((gam-1.)*ptmp[UU]  * conn_trace );
	surf_data[S_dTdr_5][itype][id] +=  dV_prop * (0.5*bsq            * conn_trace );
	surf_data[S_dTdr_6][itype][id] +=  dV_prop * (-bb_dot_conn);
#endif

      }
  }

  /*********************************************************************************** 
    If we are running in parallel, then sum over the grids in x2 directions : 
  ***********************************************************************************/ 
#if( USEMPI ) 
  //  surface_mpi_faster(  surf_data ); 
  surface_mpi(  surf_data ); 
#else 
  /* else, we need to copy data over to output arrays : */
  if( totalsize[1] != N1 ) { 
    fprintf(stderr,"dump_surface(): totalsize[1] should equal N1 :  %d %d \n", totalsize[1],N1);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }
  /* No copy necessary as we always use surf_data */
#endif 

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
  surf_memdims[     0] = N1;               surf_memdims[     1] = N3;

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



#else
void dump_surface(void) { return; } 
#endif


