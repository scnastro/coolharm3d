
#include "decs.h"

#if( MAKE_HDF5 ) 
#include <hdf5.h>
#include <hdf5_hl.h>
//#include <H5LT.h>

#if( USEMPI )
#include "mpi.h"
#endif 




/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************
  THIS FILE CONTAINS THE FOLLOWING ROUTINES:

    dump_hdf5()    : writes local gridfunctions to an hdf5 formatted file ; 

**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/

/* Use chunks for transfer or not: */ 
#define USE_CHUNKS (1) 

/* Rank of hdf5 gridfunctions : */ 
#define HDF_RANK       (3) 

#if(METRIC_DIM==0||METRIC_DIM==1)
#  define HDF_GDUMP_RANK (1) 
#elif(METRIC_DIM==2)
#  define HDF_GDUMP_RANK (2) 
#else
#  define HDF_GDUMP_RANK (3) 
#endif


/* global persistent variables : */
/* h5_dims is the same for all n_spatial_dims since we always pass a 3d array to H5Dwrite() */
static hsize_t   h5_memdims[HDF_RANK], h5_filedims[HDF_RANK];
static hsize_t   gdump_memdims[HDF_GDUMP_RANK], gdump_filedims[HDF_GDUMP_RANK];
static hsize_t   count[HDF_RANK], gdump_count[HDF_GDUMP_RANK];
static hsize_t   offset[HDF_RANK], gdump_offset[HDF_GDUMP_RANK];

static char      f_hdf_name[N_HDF_FUNCS][100];
static int       write_to_hdf[N_HDF_FUNCS];

static hsize_t dim_scal[] = {(hsize_t) 1};


/*  hdf5 utility routines : */
void     myH5_read_basic( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
static void     myH5_read_scalar( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
static char    *myH5_get_full_pathname( hid_t loc_id );
static void     set_hdf5_gfuncs(void); 
static void myH5_write_gdump_func( hid_t loc_id, char *name, double *value );

void     setup_hdf5(void); 
void     myH5_write_gfunc_orig( hid_t loc_id, char *name, double *value ) ;
void     myH5_write_gfunc_gather( hid_t loc_id, char *name, double *value ) ;
void     myH5_write_int_gfunc_orig( hid_t loc_id, char *name, int *value );
void     myH5_write_int_gfunc_gather( hid_t loc_id, char *name, int *value );
void     myH5_write_scalar( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
void     myH5_write_scalar2_orig( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
void     myH5_write_scalar2_orig2( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
void     myH5_write_scalar2_gather( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
hid_t    myH5_Fopen(char *name);
hsize_t *myH5_get_simple_dims( hid_t dataspace_id,  int *ndims ); 
void     write_dump_header(hid_t file_id, char *hdf_name, int all_procs_dump ) ;
void     myH5_read_scalar2_orig( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
void     myH5_read_scalar2_gather( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
void     myH5_read_gfunc_orig(hid_t loc_id, char *name, double *value ) ;
void     myH5_read_gfunc_gather(hid_t loc_id, char *name, double *value ) ;
hid_t    myH5_Fcreate_orig(char *name, int overwrite);
hid_t    myH5_Fcreate_gather(char *name, int overwrite);


/**********************************************************************************************/
/**********************************************************************************************
  dump_hdf5(): 
  -------------
      -- wrapper for dump_hdf5_gen() that makes full dumps and not restart files;
      -- necessary for how "dump" routines are called from diag();
**********************************************************************************************/
void dump_hdf5( void ) {  
  void dump_hdf5_gen(int do_restart);
  dump_hdf5_gen(0); 
  return; 
}

/**********************************************************************************************/
/**********************************************************************************************
  dump_hdf5_gen(): 
 -----------
   -- routine that generates dumps of various quantities in HDF5 (binary) format;
   -- can write full dumps (do_restart=0) or restart dumps (do_restart!=0); 
   -- we clobber existing files so that this routine can do restart dumps too; 

   -- full dump file names take the form of : 

        RUN_TAG.<dump_id>.h5

          where  
 
            <dump_id> = time index number

    -- restart dump file names take the form: 

           rdump_N.h5 
 
                 where N = N_out
                

    -- general notes for making hdf5 communications efficient: 
         -- http://www.hdfgroup.org/HDF5/faq/parallelchunk.html   says that the filespace 
             and memory space should have the same shape (3d memspace array goes to 
            3d filespace);
         -- collective i/o  gathers contiguous chunks from processors before writing/reading
             in order to reduce the number of i/o operations with the filesystem (since this 
             is a slow process);  independent i/o allows each processor to handle its 
             i/o communication;   not knowing how this is done, it is unclear which one 
             will be faster (need to test); 
           
     -- use the USE_CHUNKS macro to control whether or not to use chunks in parallel IO;
         -- we assume that the chunk size is the memory of one gridfunction in the 
             local domain (i.e. N1*N2*N3);

**********************************************************************************************/
void dump_hdf5_gen(int do_restart)
{
  int i,j,k,l,ind;
  int ii,jj,n1;
  hid_t file_id;
  char hdf_name[200], headvar_name[200];

  TRACE_BEG;

  /******************************************************************************
    Set up persistent data used for all time by hdf5 routines: 
  ******************************************************************************/ 
  setup_hdf5();    
  

  /* Set the name of this time's hdf file : */
  if( do_restart ) { 
#if( USE_MPI_IO || (!USEMPI) || GATHER_IO )
    if( do_restart > 1 ) { sprintf(hdf_name,"%s/rdump_starttest.h5", DIR_out[OUT_HDF5]); }
    else {                 sprintf(hdf_name,"%s/rdump_%d.h5", DIR_out[OUT_HDF5], (rdump_cnt%2)); }
#else 
    if( do_restart > 1 ) { sprintf(hdf_name,"%s/rdump_starttest.p%05d.h5",DIR_out[OUT_HDF5],myid); }
    else {                 sprintf(hdf_name,"%s/rdump_%d.p%05d.h5",DIR_out[OUT_HDF5],(rdump_cnt%2),myid); }
#endif
  }
  else { 
#if( USE_MPI_IO || (!USEMPI) || GATHER_IO )
    sprintf(hdf_name,"%s/%s.%06d.h5",DIR_out[OUT_HDF5],RUN_TAG,N_out[OUT_HDF5]);
#else 
    sprintf(hdf_name,"%s/%s.%06d.p%05d.h5",DIR_out[OUT_HDF5],RUN_TAG,N_out[OUT_HDF5],myid);
#endif
  }


  if( myid == printer_pid ) {   fprintf(stdout,"Dumping %s .... \n",hdf_name); fflush(stdout); }
  
  /******************************************************************************
    Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Create a new, clobbering file using serial or parallel properties */
  file_id = myH5_Fcreate(hdf_name, do_restart);

  write_dump_header(file_id, hdf_name, (!GATHER_IO) );
  
  /******************************************************************************
    Set and write gridfunctions : 
         -- restart is identical to dump, except that only the primitive functions 
            are dumped to restart files; 
  ******************************************************************************/ 
  if( do_restart ) { 
    /* Primitive variables: */
    ind = 0 ; 
    LOOP { 
      PLOOP { f_hdf[l][ind] = p[i][j][k][l]; } 
      ind++; 
    }

    /* Write only the primitive variables: */
    PLOOP { myH5_write_gfunc(file_id, f_hdf_name[l], f_hdf[l] ); }
  }
  else { 
    /* Calculate or copy over the grid functions to be written : */ 

    set_hdf5_gfuncs(); 

    /* Write the primitive variables to the hdf file : */
    for(i=0; i<N_HDF_FUNCS; i++)  if( write_to_hdf[i] ) { 
      myH5_write_gfunc( file_id, f_hdf_name[i], f_hdf[i] ); 
    }

    /* Dump the conserved variables : */
#if( OUTPUT_CONSERVED_VARS ) 
    ind = 0 ; 
    LOOP { 
      PLOOP { f_hdf[l][ind] = U_gf[0][i][j][k][l]; } 
      ind++; 
    }
    PLOOP { 
      sprintf(headvar_name, "U%1d", l);
      myH5_write_gfunc( file_id, headvar_name, f_hdf[l] ); 
    }
#endif

#if( OUTPUT_COORDS ) 
  struct of_coord *coords;
    ind = 0 ; 
    LOOP { 
      get_coord(i  ,j  ,k  ,CENT,ncurr,coords);
      l = 0 ; 
      f_hdf[l++][ind] = coords->x[1];
      f_hdf[l++][ind] = coords->x[2];
      f_hdf[l++][ind] = coords->x[3];
      f_hdf[l++][ind] = coords->dx_dxp[1][0];
      f_hdf[l++][ind] = coords->dx_dxp[1][1];
      f_hdf[l++][ind] = coords->dx_dxp[1][2];
      f_hdf[l++][ind] = coords->dx_dxp[1][3];
      f_hdf[l++][ind] = coords->dx_dxp[2][0];
      f_hdf[l++][ind] = coords->dx_dxp[2][1];
      f_hdf[l++][ind] = coords->dx_dxp[2][2];
      f_hdf[l++][ind] = coords->dx_dxp[2][3];
      f_hdf[l++][ind] = coords->dx_dxp[3][0];
      f_hdf[l++][ind] = coords->dx_dxp[3][1];
      f_hdf[l++][ind] = coords->dx_dxp[3][2];
      f_hdf[l++][ind] = coords->dx_dxp[3][3];
      ind++; 
    }
    l = 0; 
    myH5_write_gfunc( file_id, "x1"      , f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "x2"      , f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "x3"      , f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp10", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp11", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp12", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp13", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp20", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp21", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp22", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp23", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp30", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp31", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp32", f_hdf[l] ); l++;
    myH5_write_gfunc( file_id, "dx_dxp33", f_hdf[l] ); l++;
#endif

#if( OUTPUT_CONN ) 
    int kk;
    char fname[200];
    set_general_conn( t );
    for(ii=0;ii<NDIM;ii++)  { 
      ind = 0; 
      GDUMP_LOOP {
	l = 0;
	n1 = CONN_ID(i,j,k);
	for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
	    f_hdf[l++][ind] = conn[n1][ii][jj][kk];
	}
	ind++;
      }
      l = 0;
      for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
	  sprintf(fname,"conn%1d%1d%1d",ii,jj,kk);
	  myH5_write_gfunc( file_id, fname, f_hdf[l++] );      // 9 - 72
      }
    }
#endif

#if( OUTPUT_MAX_VCHAR )
    int ndim_m1=NDIM-1;
    char vcharname[200];
    for(ii=0;ii<ndim_m1;ii++)  { 
      ind = 0; 
      LOOP { 
	f_hdf[0][ind++] = max_vchar[i][j][k][ii];
      }
      sprintf(vcharname,"max_vchar%1d",(ii+1));
      myH5_write_gfunc( file_id, vcharname, f_hdf[0] );
    }
#endif

    
  } /* end of "else,  not do_restart " */

 
  /******************************************************************************
    Close objects : 
  ******************************************************************************/ 
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_HDF5] ) { 
#else
  {
#endif
    H5Fclose(file_id);      /* Terminate access to the file. */
  }


#if( USEMPI ) 
  //  exit_status();
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 


  TRACE_END;

  return;
}


/**********************************************************************************************/
/**********************************************************************************************
  write_dump_header():
 -----------
   -- routine that writes the header information for a hdf5 dump file ;

**********************************************************************************************/
void write_dump_header(hid_t file_id, char *hdf_name, int all_procs_dump ) 
{
  int i,n1,n2,n3,ntmp,ng,np;
  double ftmp;
  hid_t header_id, grid_id, macros_id, ioparams_id;
  char headvar_name[200];
  void (*loc_write_scalar2) ( hid_t loc_id, char *name, hid_t type_id, void *value ) ;

  TRACE_BEG;

  if( all_procs_dump )  { 
    loc_write_scalar2 = myH5_write_scalar2_orig2; 
  }
  else { 
    loc_write_scalar2 = myH5_write_scalar2_gather;
  }
   
  if( (myid == out_pid[OUT_HDF5]) || all_procs_dump ) { 
    /* Create a subgroup to which we attach attributes or header information */
    header_id = H5Gcreate(file_id, "Header", 0);
    if( header_id  < 0 ) { 
      fprintf(stderr,"write_dump_header(): Cannot create Header group in file %s \n", hdf_name);
      fflush(stderr);     fail(FAIL_HDF,0);
    }

    /* Create a header subgroup that contains the grid parameters: */
    grid_id = H5Gcreate(header_id, "Grid", 0);
    if( grid_id  < 0 ) { 
      fprintf(stderr,"write_dump_header(): Cannot create Grid group in file %s \n", hdf_name);
      fflush(stderr);     fail(FAIL_HDF,0);
    }

    /* Create a header subgroup that contains the I/O parameters: */
    ioparams_id = H5Gcreate(header_id, "IO", 0);
    if( ioparams_id  < 0 ) { 
      fprintf(stderr,"write_dump_header(): Cannot create IO group in file %s \n", hdf_name);
      fflush(stderr);     fail(FAIL_HDF,0);
    }

    /* Create a header subgroup that contains the compile time macro definitions: */
    macros_id = H5Gcreate(header_id, "Macros", 0);
    if( macros_id  < 0 ) { 
      fprintf(stderr,"write_dump_header(): Cannot create Macros group in file %s \n", hdf_name);
      fflush(stderr);     fail(FAIL_HDF,0);
    }

    /******************************************************************************
    Dump the grid parameters into the HDF file : 
    ******************************************************************************/ 
    n1 = N1;   n2 = N2;  n3 = N3; ng=NG; np = NP;
    loc_write_scalar2(grid_id, "N1"              , H5T_NATIVE_INT   ,  &n1               );
    loc_write_scalar2(grid_id, "N2"              , H5T_NATIVE_INT   ,  &n2               );
    loc_write_scalar2(grid_id, "N3"              , H5T_NATIVE_INT   ,  &n3               );
    loc_write_scalar2(grid_id, "NG"              , H5T_NATIVE_INT   ,  &ng               );
    loc_write_scalar2(grid_id, "NP"              , H5T_NATIVE_INT   ,  &np               );
    loc_write_scalar2(grid_id, "totalsize0"      , H5T_NATIVE_INT   ,  &(totalsize[0])   );
    loc_write_scalar2(grid_id, "totalsize1"      , H5T_NATIVE_INT   ,  &(totalsize[1])   );
    loc_write_scalar2(grid_id, "totalsize2"      , H5T_NATIVE_INT   ,  &(totalsize[2])   );
    loc_write_scalar2(grid_id, "totalsize3"      , H5T_NATIVE_INT   ,  &(totalsize[3])   );
    loc_write_scalar2(grid_id, "cpupos0"         , H5T_NATIVE_INT   ,  &(cpupos[0])      );
    loc_write_scalar2(grid_id, "cpupos1"         , H5T_NATIVE_INT   ,  &(cpupos[1])      );
    loc_write_scalar2(grid_id, "cpupos2"         , H5T_NATIVE_INT   ,  &(cpupos[2])      );
    loc_write_scalar2(grid_id, "cpupos3"         , H5T_NATIVE_INT   ,  &(cpupos[3])      );
    loc_write_scalar2(grid_id, "startx0"         , H5T_NATIVE_DOUBLE,  &(startx[0])      );
    loc_write_scalar2(grid_id, "startx1"         , H5T_NATIVE_DOUBLE,  &(startx[1])      );
    loc_write_scalar2(grid_id, "startx2"         , H5T_NATIVE_DOUBLE,  &(startx[2])      );
    loc_write_scalar2(grid_id, "startx3"         , H5T_NATIVE_DOUBLE,  &(startx[3])      );
    loc_write_scalar2(grid_id, "dx0"             , H5T_NATIVE_DOUBLE,  &(dx[0])          );
    loc_write_scalar2(grid_id, "dx1"             , H5T_NATIVE_DOUBLE,  &(dx[1])          );
    loc_write_scalar2(grid_id, "dx2"             , H5T_NATIVE_DOUBLE,  &(dx[2])          );
    loc_write_scalar2(grid_id, "dx3"             , H5T_NATIVE_DOUBLE,  &(dx[3])          );
    loc_write_scalar2(grid_id, "gridlength0"     , H5T_NATIVE_DOUBLE,  &(GridLength[0])  );
    loc_write_scalar2(grid_id, "gridlength1"     , H5T_NATIVE_DOUBLE,  &(GridLength[1])  );
    loc_write_scalar2(grid_id, "gridlength2"     , H5T_NATIVE_DOUBLE,  &(GridLength[2])  );
    loc_write_scalar2(grid_id, "gridlength3"     , H5T_NATIVE_DOUBLE,  &(GridLength[3])  );
    loc_write_scalar2(grid_id, "t"               , H5T_NATIVE_DOUBLE,  &t                );
    loc_write_scalar2(grid_id, "cour"            , H5T_NATIVE_DOUBLE,  &cour    	        );
    loc_write_scalar2(grid_id, "gam"             , H5T_NATIVE_DOUBLE,  &gam     	        );
    loc_write_scalar2(grid_id, "a"               , H5T_NATIVE_DOUBLE,  &a	        );
    loc_write_scalar2(grid_id, "h_slope"         , H5T_NATIVE_DOUBLE,  &h_slope	        );
    loc_write_scalar2(grid_id, "th_cutout"       , H5T_NATIVE_DOUBLE,  &th_cutout        );
    loc_write_scalar2(grid_id, "th_beg"          , H5T_NATIVE_DOUBLE,  &th_beg  	        );
    loc_write_scalar2(grid_id, "th_end"          , H5T_NATIVE_DOUBLE,  &th_end	        );
    loc_write_scalar2(grid_id, "X1_slope"        , H5T_NATIVE_DOUBLE,  &X1_slope         );
    loc_write_scalar2(grid_id, "X1_0"            , H5T_NATIVE_DOUBLE,  &X1_0	        );
    loc_write_scalar2(grid_id, "R0"              , H5T_NATIVE_DOUBLE,  &R0	        );
    loc_write_scalar2(grid_id, "Rin"             , H5T_NATIVE_DOUBLE,  &Rin	        );
    loc_write_scalar2(grid_id, "Rout"            , H5T_NATIVE_DOUBLE,  &Rout	        );
    loc_write_scalar2(grid_id, "r_isco"          , H5T_NATIVE_DOUBLE,  &r_isco	        );
    loc_write_scalar2(grid_id, "r_horizon"       , H5T_NATIVE_DOUBLE,  &r_horizon        );
    loc_write_scalar2(grid_id, "initial_bbh_separation"  , H5T_NATIVE_DOUBLE   ,  &initial_bbh_separation );
    loc_write_scalar2(grid_id, "m_bh1"           , H5T_NATIVE_DOUBLE   ,  &m_bh1);
    loc_write_scalar2(grid_id, "m_bh2"           , H5T_NATIVE_DOUBLE   ,  &m_bh2);
    loc_write_scalar2(grid_id, "m_bh_tot"        , H5T_NATIVE_DOUBLE   ,  &m_bh_tot);
    loc_write_scalar2(grid_id, "M"               , H5T_NATIVE_DOUBLE   ,  &M);
    loc_write_scalar2(grid_id, "bh1_traj0"       , H5T_NATIVE_DOUBLE   ,  &(bh1_traj[ncurr][0]));
    loc_write_scalar2(grid_id, "bh1_traj1"       , H5T_NATIVE_DOUBLE   ,  &(bh1_traj[ncurr][1]));
    loc_write_scalar2(grid_id, "bh1_traj2"       , H5T_NATIVE_DOUBLE   ,  &(bh1_traj[ncurr][2]));
    loc_write_scalar2(grid_id, "bh1_traj3"       , H5T_NATIVE_DOUBLE   ,  &(bh1_traj[ncurr][3]));
    loc_write_scalar2(grid_id, "bh2_traj0"       , H5T_NATIVE_DOUBLE   ,  &(bh2_traj[ncurr][0]));
    loc_write_scalar2(grid_id, "bh2_traj1"       , H5T_NATIVE_DOUBLE   ,  &(bh2_traj[ncurr][1]));
    loc_write_scalar2(grid_id, "bh2_traj2"       , H5T_NATIVE_DOUBLE   ,  &(bh2_traj[ncurr][2]));
    loc_write_scalar2(grid_id, "bh2_traj3"       , H5T_NATIVE_DOUBLE   ,  &(bh2_traj[ncurr][3]));

    loc_write_scalar2(grid_id, "n_within_horizon", H5T_NATIVE_INT   ,  &n_within_horizon );
    loc_write_scalar2(grid_id, "diag3_exponent"  , H5T_NATIVE_INT   ,  &diag3_exponent   );

#if(COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
    loc_write_scalar2(grid_id, "warped_spherical.delta_x1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_x1));
    loc_write_scalar2(grid_id, "warped_spherical.delta_x2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_x2));
    loc_write_scalar2(grid_id, "warped_spherical.delta_x3",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_x3));
    loc_write_scalar2(grid_id, "warped_spherical.delta_x4",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_x4));
    loc_write_scalar2(grid_id, "warped_spherical.delta_y1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_y1));
    loc_write_scalar2(grid_id, "warped_spherical.delta_y2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_y2));
    loc_write_scalar2(grid_id, "warped_spherical.delta_y3",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_y3));
    loc_write_scalar2(grid_id, "warped_spherical.delta_y4",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_y4));
    loc_write_scalar2(grid_id, "warped_spherical.delta_z1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_z1));
    loc_write_scalar2(grid_id, "warped_spherical.delta_z2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->delta_z2));
    loc_write_scalar2(grid_id, "warped_spherical.a_x10",            H5T_NATIVE_DOUBLE   ,  &(coord_params->a_x10));
    loc_write_scalar2(grid_id, "warped_spherical.a_x20",            H5T_NATIVE_DOUBLE   ,  &(coord_params->a_x20));
    loc_write_scalar2(grid_id, "warped_spherical.a_z10",            H5T_NATIVE_DOUBLE   ,  &(coord_params->a_z10));
    loc_write_scalar2(grid_id, "warped_spherical.h_x1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_x1));
    loc_write_scalar2(grid_id, "warped_spherical.h_x2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_x2));
    loc_write_scalar2(grid_id, "warped_spherical.h_x3",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_x3));
    loc_write_scalar2(grid_id, "warped_spherical.h_x4",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_x4));
    loc_write_scalar2(grid_id, "warped_spherical.h_y1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_y1));
    loc_write_scalar2(grid_id, "warped_spherical.h_y2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_y2));
    loc_write_scalar2(grid_id, "warped_spherical.h_y3",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_y3));
    loc_write_scalar2(grid_id, "warped_spherical.h_y4",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_y4));
    loc_write_scalar2(grid_id, "warped_spherical.h_z1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_z1));
    loc_write_scalar2(grid_id, "warped_spherical.h_z2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->h_z2));
    loc_write_scalar2(grid_id, "warped_spherical.s_1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->s_[1]));
    loc_write_scalar2(grid_id, "warped_spherical.s_2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->s_[2]));
    loc_write_scalar2(grid_id, "warped_spherical.s_3",            H5T_NATIVE_DOUBLE   ,  &(coord_params->s_[0]));
    loc_write_scalar2(grid_id, "warped_spherical.br_1",            H5T_NATIVE_DOUBLE   ,  &(coord_params->br_[1]));
    loc_write_scalar2(grid_id, "warped_spherical.br_2",            H5T_NATIVE_DOUBLE   ,  &(coord_params->br_[2]));
    loc_write_scalar2(grid_id, "warped_spherical.br_3",            H5T_NATIVE_DOUBLE   ,  &(coord_params->br_[0]));
#endif

#if(COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD)
    loc_write_scalar2(grid_id, "diagonal3_dyn_rad.f0_r",             H5T_NATIVE_DOUBLE   ,  &(coord_params->f0_r            ));
    loc_write_scalar2(grid_id, "diagonal3_dyn_rad.fa_r",             H5T_NATIVE_DOUBLE   ,  &(coord_params->fa_r            ));
    loc_write_scalar2(grid_id, "diagonal3_dyn_rad.upsilon_r",        H5T_NATIVE_DOUBLE   ,  &(coord_params->upsilon_r       ));
    loc_write_scalar2(grid_id, "diagonal3_dyn_rad.Rin_of_t" ,        H5T_NATIVE_DOUBLE   ,  &(coord_params->Rin_of_t        ));
    loc_write_scalar2(grid_id, "diagonal3_dyn_rad.dln_Rin_of_t_dt" , H5T_NATIVE_DOUBLE   ,  &(coord_params->dln_Rin_of_t_dt ));
#endif

    /******************************************************************************
    Dump macros: 
    ******************************************************************************/ 

    /* ints */
    n1=USEMPI      	        ;loc_write_scalar2(macros_id,"USEMPI"      	       ,H5T_NATIVE_INT,&n1);
    n1=USE_MPI_IO  	        ;loc_write_scalar2(macros_id,"USE_MPI_IO"  	       ,H5T_NATIVE_INT,&n1);
    n1=TRACE_CALLS 	        ;loc_write_scalar2(macros_id,"TRACE_CALLS" 	       ,H5T_NATIVE_INT,&n1);
    n1=RUNTIME_SECONDS	        ;loc_write_scalar2(macros_id,"RUNTIME_SECONDS"	       ,H5T_NATIVE_INT,&n1);
    n1=METRIC_TYPE_CHOICE         ;loc_write_scalar2(macros_id,"METRIC_TYPE_CHOICE"       ,H5T_NATIVE_INT,&n1);
    n1=METRIC_DYNAMIC_TYPE_CHOICE ;loc_write_scalar2(macros_id,"METRIC_DYNAMIC_TYPE_CHOICE",H5T_NATIVE_INT,&n1);
    n1=TOP_TYPE_CHOICE 	        ;loc_write_scalar2(macros_id,"TOP_TYPE_CHOICE" 	  ,H5T_NATIVE_INT,&n1);
    n1=BC_TYPE_CHOICE             ;loc_write_scalar2(macros_id,"BC_TYPE_CHOICE"           ,H5T_NATIVE_INT,&n1);
    n1=COORD_TYPE_CHOICE	        ;loc_write_scalar2(macros_id,"COORD_TYPE_CHOICE"	 ,H5T_NATIVE_INT,&n1);
    n1=RECON_TYPE_CHOICE	        ;loc_write_scalar2(macros_id,"RECON_TYPE_CHOICE"	  ,H5T_NATIVE_INT,&n1);
    n1=FLUX_TYPE_CHOICE 	        ;loc_write_scalar2(macros_id,"FLUX_TYPE_CHOICE" 	  ,H5T_NATIVE_INT,&n1);
    n1=USE_LOCAL_RECON_TYPE       ;loc_write_scalar2(macros_id,"USE_LOCAL_RECON_TYPE"     ,H5T_NATIVE_INT,&n1);
    n1=CALC_CURRENT	        ;loc_write_scalar2(macros_id,"CALC_CURRENT"	       ,H5T_NATIVE_INT,&n1);
    n1=FIXUP_TREE 	        ;loc_write_scalar2(macros_id,"FIXUP_TREE" 	       ,H5T_NATIVE_INT,&n1);
    n1=USE_SIMPLE_EOS 	        ;loc_write_scalar2(macros_id,"USE_SIMPLE_EOS" 	       ,H5T_NATIVE_INT,&n1);
    n1=ALLOW_NEGATIVE_DENSITIES   ;loc_write_scalar2(macros_id,"ALLOW_NEGATIVE_DENSITIES" ,H5T_NATIVE_INT,&n1);
    n1=USE_GAMMA_CEILING_X1DN_BC  ;loc_write_scalar2(macros_id,"USE_GAMMA_CEILING_X1DN_BC",H5T_NATIVE_INT,&n1);
    n1=USE_ENTROPY_EQ 	        ;loc_write_scalar2(macros_id,"USE_ENTROPY_EQ" 	       ,H5T_NATIVE_INT,&n1);
    n1=KEEP_CONSERVED_VARS        ;loc_write_scalar2(macros_id,"KEEP_CONSERVED_VARS"      ,H5T_NATIVE_INT,&n1);
    n1=FIX_CUTOUT_PRESSURE        ;loc_write_scalar2(macros_id,"FIX_CUTOUT_PRESSURE"      ,H5T_NATIVE_INT,&n1);
    n1=FIX_CUTOUT_GAMMA           ;loc_write_scalar2(macros_id,"FIX_CUTOUT_GAMMA"         ,H5T_NATIVE_INT,&n1);
    n1=USE_TEMPERATURE_MAX        ;loc_write_scalar2(macros_id,"USE_TEMPERATURE_MAX"      ,H5T_NATIVE_INT,&n1);
    n1=OUTPUT_CONSERVED_VARS      ;loc_write_scalar2(macros_id,"OUTPUT_CONSERVED_VARS"    ,H5T_NATIVE_INT,&n1);
    n1=CHECK_INVERSION	        ;loc_write_scalar2(macros_id,"CHECK_INVERSION"	       ,H5T_NATIVE_INT,&n1);
    n1=CHECK_PRIM_CHANGE 	        ;loc_write_scalar2(macros_id,"CHECK_PRIM_CHANGE"        ,H5T_NATIVE_INT,&n1);
    n1=RESCALE_B                  ;loc_write_scalar2(macros_id,"RESCALE_B"                ,H5T_NATIVE_INT,&n1);
    n1=RESCALE_R                  ;loc_write_scalar2(macros_id,"RESCALE_R"                ,H5T_NATIVE_INT,&n1);
    n1=FLOOR_TYPE_CHOICE	        ;loc_write_scalar2(macros_id,"FLOOR_TYPE_CHOICE"        ,H5T_NATIVE_INT,&n1);
    n1=COORDSINGFIX 	        ;loc_write_scalar2(macros_id,"COORDSINGFIX"             ,H5T_NATIVE_INT,&n1);
    n1=USE_COOLING_FUNCTION       ;loc_write_scalar2(macros_id,"USE_COOLING_FUNCTION"     ,H5T_NATIVE_INT,&n1);
    n1=USE_LIGHT_SPEED            ;loc_write_scalar2(macros_id,"USE_LIGHT_SPEED"          ,H5T_NATIVE_INT,&n1);
    n1=DUMP_ALL_STAT 	        ;loc_write_scalar2(macros_id,"DUMP_ALL_STAT"            ,H5T_NATIVE_INT,&n1);
    n1=N_HIST_TYPES  	        ;loc_write_scalar2(macros_id,"N_HIST_TYPES"             ,H5T_NATIVE_INT,&n1);
    n1=N_HIST_POINTS              ;loc_write_scalar2(macros_id,"N_HIST_POINTS"            ,H5T_NATIVE_INT,&n1);
    n1=N_SURF_TYPES  	        ;loc_write_scalar2(macros_id,"N_SURF_TYPES"             ,H5T_NATIVE_INT,&n1);
    n1=N_SURF_POINTS              ;loc_write_scalar2(macros_id,"N_SURF_POINTS"            ,H5T_NATIVE_INT,&n1);

    n1=DYNAMIC_SPACETIME          ;loc_write_scalar2(macros_id,"DYNAMIC_SPACETIME"        ,H5T_NATIVE_INT,&n1);
    n1=DYNAMIC_COORDINATES        ;loc_write_scalar2(macros_id,"DYNAMIC_COORDINATES"      ,H5T_NATIVE_INT,&n1);
    n1=CONN_METHOD_CHOICE         ;loc_write_scalar2(macros_id,"CONN_METHOD_CHOICE"       ,H5T_NATIVE_INT,&n1);

    /* doubles :  */
    ftmp=BETA_MIN          ;loc_write_scalar2(macros_id,"BETA_MIN"          ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=TEMPERATURE_MAX   ;loc_write_scalar2(macros_id,"TEMPERATURE_MAX"   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=INVERSION_THRESH  ;loc_write_scalar2(macros_id,"INVERSION_THRESH"  ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=PRIM_CHANGE_THRESH;loc_write_scalar2(macros_id,"PRIM_CHANGE_THRESH",H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=GAMMAMAX    	 ;loc_write_scalar2(macros_id,"GAMMAMAX"    	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=GAMMAMAX2   	 ;loc_write_scalar2(macros_id,"GAMMAMAX2"   	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=RHOMIN      	 ;loc_write_scalar2(macros_id,"RHOMIN"      	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=UUMIN       	 ;loc_write_scalar2(macros_id,"UUMIN"       	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=RHOMINLIMIT 	 ;loc_write_scalar2(macros_id,"RHOMINLIMIT" 	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=UUMINLIMIT  	 ;loc_write_scalar2(macros_id,"UUMINLIMIT"  	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=RHOPOWER    	 ;loc_write_scalar2(macros_id,"RHOPOWER"    	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=UUPOWER    	 ;loc_write_scalar2(macros_id,"UUPOWER"    	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=SMALL	  	 ;loc_write_scalar2(macros_id,"SMALL"	           ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=SAFE      	 ;loc_write_scalar2(macros_id,"SAFE"     	   ,H5T_NATIVE_DOUBLE,&ftmp);
    ftmp=SINGSMALL         ;loc_write_scalar2(macros_id,"SINGSMALL"         ,H5T_NATIVE_DOUBLE,&ftmp);

    /******************************************************************************
    Dump the IO parameters into the HDF header: 
    ******************************************************************************/ 
    loc_write_scalar2(ioparams_id,  "nstep"          , H5T_NATIVE_INT ,  &nstep          );
    loc_write_scalar2(ioparams_id,  "n_restart"      , H5T_NATIVE_INT ,  &n_restart      );
    loc_write_scalar2(ioparams_id,  "N_hist_dump"    , H5T_NATIVE_INT ,  &N_hist_dump    );
    loc_write_scalar2(ioparams_id,  "N_hist_dump_tot", H5T_NATIVE_INT ,  &N_hist_dump_tot);
    ntmp = N_OUT_TYPES;
    loc_write_scalar2(ioparams_id,  "N_OUT_TYPES"    , H5T_NATIVE_INT ,  &ntmp           );

    for(i=0;i<N_OUT_TYPES;i++) { 
      sprintf(headvar_name, "N_out%1d", i);
      loc_write_scalar2(ioparams_id, headvar_name, H5T_NATIVE_INT   ,  &(N_out[i]) );
    }
    for(i=0;i<N_OUT_TYPES;i++) { 
      sprintf(headvar_name, "T_out%1d", i);
      loc_write_scalar2(ioparams_id, headvar_name, H5T_NATIVE_DOUBLE,  &(T_out[i]) );
    }
    for(i=0;i<N_OUT_TYPES;i++) { 
      sprintf(headvar_name, "DT_out%1d",i);
      loc_write_scalar2(ioparams_id, headvar_name, H5T_NATIVE_DOUBLE,  &(DT_out[i]));
    }

#if( COORD_RADIUS_DIM > 1 )
    loc_write_scalar2(ioparams_id,  "n_r_bins"      , H5T_NATIVE_INT ,    &n_r_bins    );
    loc_write_scalar2(ioparams_id,  "n_phi_bins"    , H5T_NATIVE_INT ,    &n_phi_bins  );
    loc_write_scalar2(ioparams_id,  "r0_bins"       , H5T_NATIVE_DOUBLE , &r0_bins     );
    loc_write_scalar2(ioparams_id,  "r_min_bins"    , H5T_NATIVE_DOUBLE , &r_min_bins  );
    loc_write_scalar2(ioparams_id,  "r_max_bins"    , H5T_NATIVE_DOUBLE , &r_max_bins  );
    loc_write_scalar2(ioparams_id,  "xp1_min_bins"  , H5T_NATIVE_DOUBLE , &xp1_min_bins);
    loc_write_scalar2(ioparams_id,  "dxp1_bins"     , H5T_NATIVE_DOUBLE , &dxp1_bins   );
#endif

    H5Gclose(grid_id);      /* Terminate access to the group /Header/Grid. */
    H5Gclose(ioparams_id);  /* Terminate access to the group /Header/IO. */
    H5Gclose(macros_id);  /* Terminate access to the group /Header/IO. */
    H5Gclose(header_id);    /* Terminate access to the group /Header. */

  }


  TRACE_END;

  return;
}


/**********************************************************************************************/
/**********************************************************************************************
  gdump_hdf5():
 -----------
   -- routine that generates an hdf5 file containing coordinate and metric 
      information that never changes over time and is assume here to be 
      symmetric in the x3-direction (except for the x3 coordinate of course);

   -- like dump_hdf5_gen() above;

   -- writes data to a file in the HDF5 directory, filename is: 

        RUN_TAG.gdump.h5


     -- use the USE_CHUNKS macro to control whether or not to use chunks in parallel IO;
         -- we assume that the chunk size is the memory of one gridfunction in the 
             local domain (i.e. N1*N2*N3);

**********************************************************************************************/
void gdump_hdf5( void )
{
  int i,j,k,l,ind,n1,n2,n3,ng,np,ii,jj,kk,ig,jg;
  hid_t file_id, header_id, grid_id, ioparams_id;
  char hdf_name[200], fname[200];

  double dx_dxp_dxp[NDIM][NDIM][NDIM];

  struct of_geom *geom;
  struct of_coord *coords;

#if( NPOS > NDIM ) 
  const int N_f_gdump = NDIM*NDIM*NPOS;
#else
  const int N_f_gdump = NDIM*NDIM*NDIM;
#endif
  double *f_gdump[N_f_gdump];


  TRACE_BEG;

  //  exit_status();  /* Make sure that all processes are alive */

  setup_hdf5();    

  /******************************************************************************
    Set the name of this time's hdf file : 
  ******************************************************************************/ 
#if( USE_MPI_IO || (!USEMPI) || GATHER_IO ) 
  sprintf(hdf_name,"%s/%s.gdump.h5",DIR_out[OUT_HDF5],RUN_TAG);
#else
  sprintf(hdf_name,"%s/%s.gdump.p%05d.h5",DIR_out[OUT_HDF5],RUN_TAG,myid);
#endif
  if( myid == printer_pid ) {   fprintf(stdout,"Dumping %s .... \n",hdf_name); fflush(stdout); }

  /******************************************************************************
    Setup dimensionality and extent of grid functions in hdf5 file:
  ******************************************************************************/ 

  /* Dimensions of data in memory space */
#if(METRIC_DIM>=0)
  gdump_memdims[0] = 1;
#endif
#if(METRIC_DIM>=1)
  gdump_memdims[0] = N1;
#endif
#if(METRIC_DIM>=2)
  gdump_memdims[1] = N2;
#endif
#if(METRIC_DIM>=3)
  gdump_memdims[2] = N3;
#endif


#if( (USE_MPI_IO && USEMPI) || GATHER_IO )
  for(i=0;i<HDF_GDUMP_RANK;i++) { 
    gdump_filedims[i] = (hsize_t) totalsize[i+1];  /* Dimensions of data in file space */
    gdump_offset[i]   = (hsize_t) globalpos[i+1];  /* The starting position in the file to write/read data for this processor : */
  }

#else
  for(i=0;i<HDF_GDUMP_RANK;i++) { 
    gdump_filedims[i] = gdump_memdims[i] ;  /* Dimensions of data in file space */
    gdump_offset[i]   = 0;                  /* The starting position in the file to write/read data for this processor : */
  }
#endif

  /* Array for use with chunks:  */
  for(i=0;i<HDF_GDUMP_RANK;i++) { gdump_count[i] = 1; }


  /* Allocate gridfunction array: */
  for(i=0; i<N_f_gdump; i++) { 
    //-teragrid    f_gdump[i] = (double *) calloc(N1*N2, sizeof(double));
    ALLOC_ARRAY(f_gdump[i],NCELLS);
  }

  /******************************************************************************
    Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Create a new, clobbering file using serial or parallel properties */
  file_id = myH5_Fcreate(hdf_name,0);

  /* Create a subgroup to which we attach attributes or header information */
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_HDF5] ) { 
#else
  {
#endif
    header_id = H5Gcreate(file_id, "Header", 0);
    if( header_id  < 0 ) { 
      fprintf(stderr,"dump_hdf5(): Cannot create Header group in file %s \n", hdf_name);
      fflush(stderr);     fail(FAIL_HDF,0);
    }

    /******************************************************************************
    Dump the grid parameters into the HDF file : 
    ******************************************************************************/ 
    n1 = N1;   n2 = N2;  n3 = N3; ng=NG; np = NP;
    myH5_write_scalar2(header_id, "N1"              , H5T_NATIVE_INT   ,  &n1              );
    myH5_write_scalar2(header_id, "N2"              , H5T_NATIVE_INT   ,  &n2              );
    myH5_write_scalar2(header_id, "N3"              , H5T_NATIVE_INT   ,  &n3              );
    myH5_write_scalar2(header_id, "NG"              , H5T_NATIVE_INT   ,  &ng             );
    myH5_write_scalar2(header_id, "NP"              , H5T_NATIVE_INT   ,  &np             );
    myH5_write_scalar2(header_id, "totalsize0"      , H5T_NATIVE_INT   ,  &(totalsize[0])  );
    myH5_write_scalar2(header_id, "totalsize1"      , H5T_NATIVE_INT   ,  &(totalsize[1])  );
    myH5_write_scalar2(header_id, "totalsize2"      , H5T_NATIVE_INT   ,  &(totalsize[2])  );
    myH5_write_scalar2(header_id, "totalsize3"      , H5T_NATIVE_INT   ,  &(totalsize[3])  );
    myH5_write_scalar2(header_id, "cpupos0"         , H5T_NATIVE_INT   ,  &(cpupos[0])    );
    myH5_write_scalar2(header_id, "cpupos1"         , H5T_NATIVE_INT   ,  &(cpupos[1])    );
    myH5_write_scalar2(header_id, "cpupos2"         , H5T_NATIVE_INT   ,  &(cpupos[2])    );
    myH5_write_scalar2(header_id, "cpupos3"         , H5T_NATIVE_INT   ,  &(cpupos[3])    );
    myH5_write_scalar2(header_id, "startx0"         , H5T_NATIVE_DOUBLE,  &(startx[0])     );
    myH5_write_scalar2(header_id, "startx1"         , H5T_NATIVE_DOUBLE,  &(startx[1])     );
    myH5_write_scalar2(header_id, "startx2"         , H5T_NATIVE_DOUBLE,  &(startx[2])     );
    myH5_write_scalar2(header_id, "startx3"         , H5T_NATIVE_DOUBLE,  &(startx[3])     );
    myH5_write_scalar2(header_id, "dx0"             , H5T_NATIVE_DOUBLE,  &(dx[0])         );
    myH5_write_scalar2(header_id, "dx1"             , H5T_NATIVE_DOUBLE,  &(dx[1])         );
    myH5_write_scalar2(header_id, "dx2"             , H5T_NATIVE_DOUBLE,  &(dx[2])         );
    myH5_write_scalar2(header_id, "dx3"             , H5T_NATIVE_DOUBLE,  &(dx[3])         );
    myH5_write_scalar2(header_id, "gridlength0"     , H5T_NATIVE_DOUBLE,  &(GridLength[0]) );
    myH5_write_scalar2(header_id, "gridlength1"     , H5T_NATIVE_DOUBLE,  &(GridLength[1]) );
    myH5_write_scalar2(header_id, "gridlength2"     , H5T_NATIVE_DOUBLE,  &(GridLength[2]) );
    myH5_write_scalar2(header_id, "gridlength3"     , H5T_NATIVE_DOUBLE,  &(GridLength[3]) );
    myH5_write_scalar2(header_id, "cour"            , H5T_NATIVE_DOUBLE,  &cour     	);
    myH5_write_scalar2(header_id, "gam"             , H5T_NATIVE_DOUBLE,  &gam     	);
    myH5_write_scalar2(header_id, "a"               , H5T_NATIVE_DOUBLE,  &a	    	);
    myH5_write_scalar2(header_id, "h_slope"         , H5T_NATIVE_DOUBLE,  &h_slope	 );
    myH5_write_scalar2(header_id, "th_cutout"       , H5T_NATIVE_DOUBLE,  &th_cutout  	 );
    myH5_write_scalar2(header_id, "th_beg"          , H5T_NATIVE_DOUBLE,  &th_beg  	 );
    myH5_write_scalar2(header_id, "th_end"          , H5T_NATIVE_DOUBLE,  &th_end	         );
    myH5_write_scalar2(header_id, "X1_slope"        , H5T_NATIVE_DOUBLE,  &X1_slope        );
    myH5_write_scalar2(header_id, "X1_0"            , H5T_NATIVE_DOUBLE,  &X1_0	        );
    myH5_write_scalar2(header_id, "R0"              , H5T_NATIVE_DOUBLE,  &R0	        );
    myH5_write_scalar2(header_id, "Rin"             , H5T_NATIVE_DOUBLE,  &Rin	        );
    myH5_write_scalar2(header_id, "Rout"            , H5T_NATIVE_DOUBLE,  &Rout	        );
    myH5_write_scalar2(header_id, "r_isco"          , H5T_NATIVE_DOUBLE,  &r_isco	        );
    myH5_write_scalar2(header_id, "r_horizon"       , H5T_NATIVE_DOUBLE,  &r_horizon       );
    myH5_write_scalar2(header_id, "n_within_horizon", H5T_NATIVE_INT   ,  &n_within_horizon);
    myH5_write_scalar2(header_id, "diag3_exponent"  , H5T_NATIVE_INT   ,  &diag3_exponent   );
    myH5_write_scalar2(header_id, "t"               , H5T_NATIVE_DOUBLE,  &t                );

    H5Gclose(header_id);    /* Terminate access to the group /Header. */

  }

  /******************************************************************************
    Set and write gridfunctions : 
         -- restart is identical to dump, except that only the primitive functions 
            are dumped to restart files; 
  ******************************************************************************/ 
#define DUMP_DX_DXP_DXP (0)

  k = N3S;

  /*   x1,x2  xp1,xp2  dx_dxp[1,2][1,2]  */
  ind = 0; 
  N1_LOOP N2_LOOP N3_LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);

    l = 0; 
    f_gdump[l++][ind] = coords->x[1];          
    f_gdump[l++][ind] = coords->x[2];           
    f_gdump[l++][ind] = coords->x[3];           
    f_gdump[l++][ind] = coords->xp[1];          
    f_gdump[l++][ind] = coords->xp[2];          
    f_gdump[l++][ind] = coords->xp[3];          
    f_gdump[l++][ind] = coords->dx_dxp[RR][1];  
    f_gdump[l++][ind] = coords->dx_dxp[RR][2];  
    f_gdump[l++][ind] = coords->dx_dxp[RR][3];  
    f_gdump[l++][ind] = coords->dx_dxp[TH][1];  
    f_gdump[l++][ind] = coords->dx_dxp[TH][2];  
    f_gdump[l++][ind] = coords->dx_dxp[TH][3];  
    f_gdump[l++][ind] = coords->dx_dxp[PH][1];  
    f_gdump[l++][ind] = coords->dx_dxp[PH][2];  
    f_gdump[l++][ind] = coords->dx_dxp[PH][3];  

# if( DUMP_DX_DXP_DXP )
    dx_dxp_dxp_calc( coords->x, coords->xp, dx_dxp_dxp ) ;
    f_gdump[l++][ind] = dx_dxp_dxp[RR][1][1];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][1][2];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][1][3];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][2][1];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][2][2];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][2][3];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][3][1];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][3][2];
    f_gdump[l++][ind] = dx_dxp_dxp[RR][3][3];

    f_gdump[l++][ind] = dx_dxp_dxp[TH][1][1];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][1][2];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][1][3];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][2][1];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][2][2];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][2][3];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][3][1];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][3][2];
    f_gdump[l++][ind] = dx_dxp_dxp[TH][3][3];

    f_gdump[l++][ind] = dx_dxp_dxp[PH][1][1];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][1][2];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][1][3];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][2][1];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][2][2];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][2][3];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][3][1];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][3][2];
    f_gdump[l++][ind] = dx_dxp_dxp[PH][3][3];
# endif


    ind++;
  }
  l = 0; 

  myH5_write_gfunc( file_id, "x1"     , f_gdump[l++] );   // 1
  myH5_write_gfunc( file_id, "x2"     , f_gdump[l++] );   // 2
  myH5_write_gfunc( file_id, "x3"     , f_gdump[l++] );   // 3
  myH5_write_gfunc( file_id, "xp1"    , f_gdump[l++] );   // 4
  myH5_write_gfunc( file_id, "xp2"    , f_gdump[l++] );   // 5
  myH5_write_gfunc( file_id, "xp3"    , f_gdump[l++] );   // 6
  myH5_write_gfunc( file_id, "dxdxp11", f_gdump[l++] );   // 7
  myH5_write_gfunc( file_id, "dxdxp12", f_gdump[l++] );   // 8
  myH5_write_gfunc( file_id, "dxdxp13", f_gdump[l++] );   // 9
  myH5_write_gfunc( file_id, "dxdxp21", f_gdump[l++] );   // 10
  myH5_write_gfunc( file_id, "dxdxp22", f_gdump[l++] );   // 11
  myH5_write_gfunc( file_id, "dxdxp23", f_gdump[l++] );   // 12
  myH5_write_gfunc( file_id, "dxdxp31", f_gdump[l++] );   // 13
  myH5_write_gfunc( file_id, "dxdxp32", f_gdump[l++] );   // 14
  myH5_write_gfunc( file_id, "dxdxp33", f_gdump[l++] );   // 15


# if( DUMP_DX_DXP_DXP )
  myH5_write_gfunc( file_id, "dxdxpdxp111", f_gdump[l++] );   // 16
  myH5_write_gfunc( file_id, "dxdxpdxp112", f_gdump[l++] );   // 17
  myH5_write_gfunc( file_id, "dxdxpdxp113", f_gdump[l++] );   // 18
  myH5_write_gfunc( file_id, "dxdxpdxp121", f_gdump[l++] );   // 19
  myH5_write_gfunc( file_id, "dxdxpdxp122", f_gdump[l++] );   // 20
  myH5_write_gfunc( file_id, "dxdxpdxp123", f_gdump[l++] );   // 21
  myH5_write_gfunc( file_id, "dxdxpdxp131", f_gdump[l++] );   // 22
  myH5_write_gfunc( file_id, "dxdxpdxp132", f_gdump[l++] );   // 23
  myH5_write_gfunc( file_id, "dxdxpdxp133", f_gdump[l++] );   // 24

  myH5_write_gfunc( file_id, "dxdxpdxp211", f_gdump[l++] );   // 25
  myH5_write_gfunc( file_id, "dxdxpdxp212", f_gdump[l++] );   // 26
  myH5_write_gfunc( file_id, "dxdxpdxp213", f_gdump[l++] );   // 27
  myH5_write_gfunc( file_id, "dxdxpdxp221", f_gdump[l++] );   // 28
  myH5_write_gfunc( file_id, "dxdxpdxp222", f_gdump[l++] );   // 29
  myH5_write_gfunc( file_id, "dxdxpdxp223", f_gdump[l++] );   // 30
  myH5_write_gfunc( file_id, "dxdxpdxp231", f_gdump[l++] );   // 31
  myH5_write_gfunc( file_id, "dxdxpdxp232", f_gdump[l++] );   // 32
  myH5_write_gfunc( file_id, "dxdxpdxp233", f_gdump[l++] );   // 33

  myH5_write_gfunc( file_id, "dxdxpdxp311", f_gdump[l++] );   // 34
  myH5_write_gfunc( file_id, "dxdxpdxp312", f_gdump[l++] );   // 35
  myH5_write_gfunc( file_id, "dxdxpdxp313", f_gdump[l++] );   // 36
  myH5_write_gfunc( file_id, "dxdxpdxp321", f_gdump[l++] );   // 37
  myH5_write_gfunc( file_id, "dxdxpdxp322", f_gdump[l++] );   // 38
  myH5_write_gfunc( file_id, "dxdxpdxp323", f_gdump[l++] );   // 39
  myH5_write_gfunc( file_id, "dxdxpdxp331", f_gdump[l++] );   // 40
  myH5_write_gfunc( file_id, "dxdxpdxp332", f_gdump[l++] );   // 41
  myH5_write_gfunc( file_id, "dxdxpdxp333", f_gdump[l++] );   // 42
# endif

  /*  Connection conn[][][] */
  set_general_conn( t );
  ind = 0; 
  GDUMP_LOOP {
    l = 0;
    n1 = CONN_ID(i,j,k);
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
      f_gdump[l++][ind] = conn[n1][ii][jj][kk];
    }
    ind++;
  }

  l = 0;
  for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
    sprintf(fname,"conn%1d%1d%1d",ii,jj,kk);
    //-teragrid    myH5_write_gdump_func( file_id, fname, f_gdump[l++] );      // 9 - 72
    myH5_write_gfunc( file_id, fname, f_gdump[l++] );      // 9 - 72
  }
    
  /*  Metric contravariant gcon[] */
  ind = 0; 
  GDUMP_LOOP { 
    l = 0;
    for(kk=0;kk<NPOS; kk++) { 
      get_geometry(i,j,k,kk,ncurr,geom); 
      for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {
	f_gdump[l++][ind] = geom->gcon[ii][jj];
      }
    }
    ind++; 
  }
  l = 0;
  for(kk=0;kk<NPOS; kk++) for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {
    sprintf(fname,"gcon%1d%1d%1d",kk,ii,jj);
    //-teragrid    myH5_write_gdump_func( file_id, fname, f_gdump[l++] );      // 73 - 136
    myH5_write_gfunc( file_id, fname, f_gdump[l++] );      // 73 - 136
  }


  /*  Metric covariant gcov[] */
  ind = 0; 
  GDUMP_LOOP { 
    l = 0;
    for(kk=0;kk<NPOS; kk++) { 
      get_geometry(i,j,k,kk,ncurr,geom); 
      for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {
	f_gdump[l++][ind] = geom->gcov[ii][jj];
      }
    }
    ind++; 
  }
  l = 0;
  for(kk=0;kk<NPOS; kk++) for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {
    sprintf(fname,"gcov%1d%1d%1d",kk,ii,jj);
    //-teragrid   myH5_write_gdump_func( file_id, fname, f_gdump[l++] );      // 137 - 200
    myH5_write_gfunc( file_id, fname, f_gdump[l++] );      // 137 - 200
  }


  /*  Sqrt of (-det(g))  gdet[] */
  ind = 0; 
  GDUMP_LOOP { 
    l = 0;
    for(kk=0;kk<NPOS; kk++)  {
      get_geometry(i,j,k,kk,ncurr,geom); 
      f_gdump[l++][ind] = geom->g;
    }
    ind++; 
  }
  l = 0;
  for(kk=0;kk<NPOS; kk++)   {
    sprintf(fname,"gdet%1d",kk);
    //-teragrid    myH5_write_gdump_func( file_id, fname, f_gdump[l++] );      // 201 - 204
    myH5_write_gfunc( file_id, fname, f_gdump[l++] );      // 201 - 204
  }


  /******************************************************************************
    Close objects : 
  ******************************************************************************/ 
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_HDF5] ) {
#else
  {
#endif
    H5Fclose(file_id);      /* Terminate access to the file. */
  }


  /* Free memory */
  for(i=0; i<N_f_gdump; i++) { FREE(f_gdump[i]); }

#if( USEMPI ) 
  //  exit_status();
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 

  TRACE_END;

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 setup_hdf5():
 ---------------------------
   -- initializes arrays used for all grid functions written to hdf5 files;
   -- note that based on how globalpos, totalsize and N[1-3] are defined, 
       a serial run will be automatically handled;
*******************************************************************************************/
void setup_hdf5( void ) 
{
  int i,l;
  static int first_time = 1; 

  if( !first_time ) { 
    return; 
  }
  
  first_time = 0; 

  /* Dimensions of data in memory space */
  h5_memdims[0] = N1;
  h5_memdims[1] = N2;
  h5_memdims[2] = N3;

#if( (USE_MPI_IO && USEMPI) || GATHER_IO )
  /* Dimensions of data in file space */
  h5_filedims[0] = (hsize_t) totalsize[1];
  h5_filedims[1] = (hsize_t) totalsize[2];
  h5_filedims[2] = (hsize_t) totalsize[3];

  /* The starting position in the file to write/read data for this processor : */
  offset[0] = (hsize_t) globalpos[1]; 
  offset[1] = (hsize_t) globalpos[2]; 
  offset[2] = (hsize_t) globalpos[3]; 
#else 
  /* Dimensions of data in file space */
  h5_filedims[0] = N1;
  h5_filedims[1] = N2;
  h5_filedims[2] = N3;

  /* The starting position in the file to write/read data for this processor : */
  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0;
#endif

  /*********************************************************************************
    Array for use with chunks:  
      -- this means that we send one chunk (a full block) of many points, 
           instead of the other way that sends many blocks of one point; 
  **********************************************************************************/
  count[0] = 1;
  count[1] = 1;
  count[2] = 1;


  /*******************************************************************************
      Set the names of the gridfunctions.  
       -- make sure that the order here matches the order of f_hdf[] assignment below:
  ********************************************************************************/
  l = 0;

  /* Set write_to_hdf[] to non-zero  next to name  in order to include that function in hdf dumps */
  sprintf(f_hdf_name[l], "/rho"    );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/uu"     );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/v1"     );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/v2"     );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/v3"     );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/B1"     );  write_to_hdf[l++] = (!HYDRO_ONLY);
  sprintf(f_hdf_name[l], "/B2"     );  write_to_hdf[l++] = (!HYDRO_ONLY);
  sprintf(f_hdf_name[l], "/B3"     );  write_to_hdf[l++] = (!HYDRO_ONLY);
  sprintf(f_hdf_name[l], "/bsq"    );  write_to_hdf[l++] = (!HYDRO_ONLY);
  sprintf(f_hdf_name[l], "/gamma"  );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/divb"   );  write_to_hdf[l++] = (!HYDRO_ONLY);
  sprintf(f_hdf_name[l], "/divbcen");  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov00" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov01" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov02" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov03" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov11" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov12" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov13" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov22" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov23" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gcov33" );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/gdet"   );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/pflag"  );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/evmask" );  write_to_hdf[l++] = (USE_MASK);
//  sprintf(f_hdf_name[l], "/nfail00");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail01");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail02");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail03");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail04");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail05");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail06");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail07");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail08");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail09");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail10");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail11");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail12");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail13");  write_to_hdf[l++] = 0;
//  sprintf(f_hdf_name[l], "/nfail14");  write_to_hdf[l++] = 0;

#if( USE_KINEMATIC_VISCOSITY ) 
  sprintf(f_hdf_name[l], "/visc0" );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/visc1" );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/visc2" );  write_to_hdf[l++] = 1;
  sprintf(f_hdf_name[l], "/visc3" );  write_to_hdf[l++] = 1;
#endif

#if( CALC_CURRENT )     
  sprintf(f_hdf_name[l], "/jcon0"  );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/jcon1"  );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/jcon2"  );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/jcon3"  );  write_to_hdf[l++] = 0;
  sprintf(f_hdf_name[l], "/jsq"    );  write_to_hdf[l++] = 0;
#endif 

  if( l != N_HDF_FUNCS ) { 
    fprintf(stderr,"setup_hdf5(): f_hdf count mismatch : %d  %d \n", l, N_HDF_FUNCS);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_gfunc_orig(): 
 ---------------------------
   -- original version of the routine to write standard gfuncs; 
   -- driver routine for writing a grid function within a given group/object
      "loc_id"; 
   -- the grid function is a simple dataspace of dimensions given by the 
      macros "N1*N2*N3";
   -- H5T_NATIVE_DOUBLE  is the assumed datatype of all grid functions ;
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between MPI and serial calls; 
      -- assumes that all processes write gridfunctions in the typical domain 
         decomposition pattern;
*******************************************************************************************/
void myH5_write_gfunc_orig( hid_t loc_id, char *name, double *value )
{
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;

  TRACE_BEG;

  /************************************************************************************
     Create the data set: 
   ************************************************************************************/
  /* Create the data space for the dataset, always using 3 dimensions since we want 
     filespace to match memory space : */
  filespace_id = H5Screate_simple(HDF_RANK, h5_filedims, NULL); 
  memspace_id  = H5Screate_simple(HDF_RANK, h5_memdims , NULL); 

  /* Set the property list to for creation */ 
  prop_id = H5Pcreate(H5P_DATASET_CREATE);

#if( USE_CHUNKS ) 
  /* Set properties to use chunks, set chunk's extent : */
  H5Pset_chunk(prop_id, HDF_RANK, h5_memdims);
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
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, h5_memdims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, h5_memdims, NULL);
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

  TRACE_END;

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_gfunc_gather(): 
 ---------------------------
   -- unlike the original version of the routine in that it writes the gfuncs after doing 
      a MPI_Gather to collect all the data before writing; 
   -- using (edited) blocks of code from Hotaka Shiokawa's miharm3d routine  "dump_single()" with permission;
   -- driver routine for writing a grid function within a given group/object
      "loc_id"; 
   -- the grid function is a simple dataspace of dimensions given by the 
      macros "N1*N2*N3";
   -- H5T_NATIVE_DOUBLE  is the assumed datatype of all grid functions ;
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between MPI and serial calls; 
      -- assumes that all processes write gridfunctions in the typical domain 
         decomposition pattern;
*******************************************************************************************/
void myH5_write_gfunc_gather( hid_t loc_id, char *name, double *value )
{
  int i,j,k,temp_id;
  ulint n; 
  static ulint *sort; 
  static double *data;
  static double *recv;
  static usint local_first_time = 1; 

  TRACE_BEG;

#if(USEMPI) /* nothing needs to happen for non-mpi jobs */

  if( local_first_time ) { 
    if( myid == out_pid[OUT_HDF5] ) {
      ALLOC_ARRAY(sort,n_cells_glob); 
      ALLOC_ARRAY(data,n_cells_glob); 
      ALLOC_ARRAY(recv,n_cells_glob); 

      /* Construct the array "sort" which contains how the gathered data "recv" should be sorted */
      n = 0;
      for(temp_id = 0; temp_id < numprocs; temp_id++) {
	for(i=0; i<N1; i++) for(j=0; j<N2; j++) for(k=0; k<N3; k++) {
	      sort[n] = get_global_index_pid(temp_id,i,j,k);
	      n++;
	    }
      }
    }
    local_first_time = 0 ; 
  }

  /* Gather data to processor 0 */
  exit_status();
  MPI_Gather(value,NCELLS,MPI_DOUBLE,recv,NCELLS,MPI_DOUBLE,out_pid[OUT_HDF5],MPI_COMM_WORLD);

  /* Sort received data */
  if( myid == out_pid[OUT_HDF5] ) {
    for(n=0; n<n_cells_glob; n++) {
      data[sort[n]] = recv[n];
    }

    /* Write to file */
    H5LTmake_dataset_double(loc_id, name, 3, h5_filedims, data);
  }

  //  MPI_Barrier(MPI_COMM_WORLD);

#endif

  TRACE_END;

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_int_gfunc_orig(): 
 ---------------------------
   -- just like myH5_write_gfunc() except dumps integer arrays;
*******************************************************************************************/
void myH5_write_int_gfunc_orig( hid_t loc_id, char *name, int *value )
{
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;

  /************************************************************************************
     Create the data set: 
   ************************************************************************************/
  /* Create the data space for the dataset, always using 3 dimensions since we want 
     filespace to match memory space : */
  filespace_id = H5Screate_simple(HDF_RANK, h5_filedims, NULL); 
  memspace_id  = H5Screate_simple(HDF_RANK, h5_memdims , NULL); 

  /* Set the property list to for creation */ 
  prop_id = H5Pcreate(H5P_DATASET_CREATE);

#if( USE_CHUNKS ) 
  /* Set properties to use chunks, set chunk's extent : */
  H5Pset_chunk(prop_id, HDF_RANK, h5_memdims);
#endif

  /* Create a dataset assuming double var type */
  dataset_id = H5Dcreate(loc_id, name, H5T_NATIVE_INT, filespace_id, prop_id);
  H5Pclose(prop_id);   
  H5Sclose(filespace_id);


  /************************************************************************************
    Select the hyperslab, or section of space, to which our gridfunction will be written
       -- I am not sure why we need to close and then reget the dataset's dataspace, 
          but that is what they usually do in the example programs in the hdf5 tutorial;
  ************************************************************************************/
  filespace_id = H5Dget_space(dataset_id);  

#if( USE_CHUNKS ) 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, h5_memdims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, h5_memdims, NULL);
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
  H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace_id, filespace_id, prop_id, value);

  /* Close the open object handles : */
  H5Dclose(dataset_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
  H5Pclose(prop_id);

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_int_gfunc_gather(): 
 ---------------------------
   -- just like myH5_write_gfunc() except dumps integer arrays;
*******************************************************************************************/
void myH5_write_int_gfunc_gather( hid_t loc_id, char *name, int *value )
{
  int i,j,k,temp_id;
  ulint n; 
  static ulint *sort; 
  static int *data;
  static int *recv;
  static usint local_first_time = 1; 

  TRACE_BEG;

#if(USEMPI) /* nothing needs to happen for non-mpi jobs */

  if( local_first_time ) { 
    if( myid == out_pid[OUT_HDF5] ) {
      ALLOC_ARRAY(sort,n_cells_glob); 
      ALLOC_ARRAY(data,n_cells_glob); 
      ALLOC_ARRAY(recv,n_cells_glob); 

      /* Construct the array "sort" which contains how the gathered data "recv" should be sorted */
      n = 0;
      for(temp_id = 0; temp_id < numprocs; temp_id++) {
	for(i=0; i<N1; i++) for(j=0; j<N2; j++) for(k=0; k<N3; k++) {
	      sort[n] = get_global_index_pid(temp_id,i,j,k);
	      n++;
	    }
      }
    }
    local_first_time = 0 ; 
  }

  /* Gather data to processor 0 */
  exit_status();
  MPI_Gather(value,NCELLS,MPI_INT,recv,NCELLS,MPI_INT,out_pid[OUT_HDF5],MPI_COMM_WORLD);

  /* Sort received data */
  if( myid == out_pid[OUT_HDF5] ) {
    for(n=0; n<n_cells_glob; n++) {
      data[sort[n]] = recv[n];
    }

    /* Write to file */
    H5LTmake_dataset_int(loc_id, name, 3, h5_filedims, data);
  }

  //  MPI_Barrier(MPI_COMM_WORLD);

#endif

  TRACE_END;

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_gdump_func(): 
 ---------------------------
   -- driver routine for writing a grid function for gdumps within a given group/object
      "loc_id"; 
   -- the grid function is a simple dataspace of dimensions given by the 
      macros "N1*N2";
   -- H5T_NATIVE_DOUBLE  is the assumed datatype of all grid functions ;
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between MPI and serial calls; 
      -- ensures that only those processors at cpupos[3] = 0 dump data;

*******************************************************************************************/
static void myH5_write_gdump_func( hid_t loc_id, char *name, double *value )
{
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;

  /************************************************************************************
     Create the data set: 
   ************************************************************************************/

  /* Create the data space for the dataset, always using 3 dimensions since we want 
     filespace to match memory space : */
  filespace_id = H5Screate_simple(HDF_GDUMP_RANK, gdump_filedims, NULL); 
  memspace_id  = H5Screate_simple(HDF_GDUMP_RANK, gdump_memdims , NULL); 

  /* Set the property list to for creation */ 
  prop_id = H5Pcreate(H5P_DATASET_CREATE);

#if( USE_CHUNKS ) 
  /* Set properties to use chunks, set chunk's extent : */
  H5Pset_chunk(prop_id, HDF_GDUMP_RANK, gdump_memdims);
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
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, gdump_offset, NULL, gdump_count, gdump_memdims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, gdump_offset, NULL, gdump_memdims, NULL);
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

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_scalar(): 
 ---------------------------
   -- driver routine for writing a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- assumes that the master node is the only one that writes scalar data; 
        -- note that all nodes must create dataset, but only one needs to write to it;
*******************************************************************************************/
void myH5_write_scalar( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  hid_t   filespace_id, memspace_id, dataset_id, prop_id;
  hsize_t dims[]={1};

  /* Create the dataspace for the dataset, single number format:  */
  filespace_id = H5Screate_simple(1, dims, NULL); 

  prop_id = H5Pcreate(H5P_DATASET_CREATE);

  /* Create a dataset, which is like a small dataset  */
  dataset_id = H5Dcreate(loc_id, name, type_id, filespace_id, H5P_DEFAULT);
  H5Sclose(filespace_id);	 
  H5Pclose(prop_id);

  /* Create memory and file spaces, removing extents if we are not the master node : */
  memspace_id  = H5Screate_simple(1, dims, NULL); 
  filespace_id = H5Dget_space(dataset_id);
  if( myid != out_pid[OUT_HDF5] ) { 
    H5Sselect_none(memspace_id); 
    H5Sselect_none(filespace_id);
  }

  /* Write the dataset using defined dataspace and default properties. */
  prop_id = H5Pcreate(H5P_DATASET_XFER);

#if( USE_MPI_IO ) 
  H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
#endif 

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
 myH5_rewrite_scalar(): 
 ---------------------------
   -- driver routine for writing an exisitng  scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for creating dataspace, creating dataset, writing to an exisitng,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- assumes that the master node is the only one that writes scalar data; 
        -- note that all nodes must create dataset, but only one needs to write to it;
*******************************************************************************************/
void myH5_rewrite_scalar( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  hid_t   filespace_id, memspace_id, dataset_id;
  hsize_t dims[]={1};

  /* Create a dataset, which is like a small dataset  */
  dataset_id = H5Dopen(loc_id, name); 
  if( dataset_id < 0 ) { 
    fprintf(stderr,"%s(): Cannot open dataset named  %s \n",__func__, name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create memory and file spaces, removing extents if we are not the master node : */
  memspace_id  = H5Screate_simple(1, dims, NULL); 
  filespace_id = H5Dget_space(dataset_id);

  H5Dwrite(dataset_id, type_id, memspace_id, filespace_id, H5P_DEFAULT, value);

  /* Close all open object handles : */
  H5Dclose(dataset_id);	 
  H5Sclose(filespace_id);	 
  H5Sclose(memspace_id);	 

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_scalar2_orig(): 
 ---------------------------
   -- different from myH5_write_scalar() in that it has all processors write the data 
      (which is redundant here of course but does not require the use of "H5Sselect_none()"
       that seems to not work at NCSA);
   -- driver routine for writing a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- assumes that the master node is the only one that writes scalar data; 
        -- note that all nodes must create dataset, but only one needs to write to it;
*******************************************************************************************/
void myH5_write_scalar2_orig( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  hid_t   filespace_id, memspace_id, dataset_id, prop_id;
  static hsize_t scal_memdims[1], scal_filedims[1], scal_count[1] ;
  static hsize_t scal_offset[1];
  static int first_time_loc = 1; 


  /* Setup offsets, dimensions, etc.;  the datasets dump scalar data, so there 
     will be one scalar per processor in the data file;  */
  if( first_time_loc ) { 
    scal_memdims[0]  = (hsize_t) 1;
    scal_count[0]    = (hsize_t) 1;
    first_time_loc = 0; 
#if( USE_MPI_IO ) 
    scal_filedims[0] = (hsize_t) numprocs;
    scal_offset[0]   = (hsize_t) myid; 
#else
    scal_filedims[0] = (hsize_t) 1;
    scal_offset[0]   = (hsize_t) 0;
#endif
  }

  /************************************************************************************
     Create the data set: 
   ************************************************************************************/
  /* Create the data space for the dataset, always using 3 dimensions since we want 
     filespace to match memory space : */
  filespace_id = H5Screate_simple(1, scal_filedims, NULL); 
  memspace_id  = H5Screate_simple(1, scal_memdims , NULL); 

  /* Set the property list to for creation */ 
  prop_id = H5Pcreate(H5P_DATASET_CREATE);

#if( USE_CHUNKS ) 
  /* Set properties to use chunks, set chunk's extent : */
  H5Pset_chunk(prop_id, 1, scal_memdims);
#endif

  /* Create a dataset assuming double var type */
  dataset_id = H5Dcreate(loc_id, name, type_id, filespace_id, prop_id);

  H5Pclose(prop_id);   
  H5Sclose(filespace_id);

  /************************************************************************************
    Select the hyperslab, or section of space, to which our gridfunction will be written
       -- I am not sure why we need to close and then reget the dataset's dataspace, 
          but that is what they usually do in the example programs in the hdf5 tutorial;
  ************************************************************************************/
  filespace_id = H5Dget_space(dataset_id);  

#if( USE_CHUNKS ) 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, scal_offset, NULL, scal_count, scal_memdims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, scal_offset, NULL, scal_memdims, NULL);
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
  H5Dwrite(dataset_id, type_id, memspace_id, filespace_id, prop_id, value);

  /* Close the open object handles : */
  H5Dclose(dataset_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
  H5Pclose(prop_id);

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_scalar2_gather(): 
 --------------------------- 
   -- version that assumes that only master processor write the data;
   -- different from myH5_write_scalar() in that it has all processors write the data 
      (which is redundant here of course but does not require the use of "H5Sselect_none()"
       that seems to not work at NCSA);
   -- driver routine for writing a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- assumes that the master node is the only one that writes scalar data; 
        -- note that all nodes must create dataset, but only one needs to write to it;
*******************************************************************************************/
void myH5_write_scalar2_gather( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  if( myid == out_pid[OUT_HDF5] ) { 
    herr_t  status;
    status = H5LTmake_dataset(loc_id,name,1,dim_scal,type_id,value);
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with making scalar dataset: %s  %d \n",
	      __func__,myid,nstep,n_substep,name,status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }
  }
  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_write_scalar2_orig2(): 
 --------------------------- 
   -- different from myH5_write_scalar() in that it has all processors write the data 
      (which is redundant here of course but does not require the use of "H5Sselect_none()"
       that seems to not work at NCSA);
   -- driver routine for writing a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for creating dataspace, creating dataset, writing dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- assumes that the master node is the only one that writes scalar data; 
        -- note that all nodes must create dataset, but only one needs to write to it;
*******************************************************************************************/
void myH5_write_scalar2_orig2( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  herr_t  status;
  status = H5LTmake_dataset(loc_id,name,1,dim_scal,type_id,value);
  if( status < 0 ) { 
    fprintf(stdout,"%s() %d %d %d : problem with making scalar dataset: %s  %d \n",
	    __func__,myid,nstep,n_substep,name,status); 
    fflush(stdout); fail(FAIL_HDF,0); 
  }
  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_read_scalar(): 
 ---------------------------
   -- driver routine for reading a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for opening dataset, opening dataspace, reading dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- unlike myH5_write_scalar(), all processes read the scalar data (this eliminates
           a broadcast of the data);
*******************************************************************************************/
static void myH5_read_scalar( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  hid_t   filespace_id, memspace_id, dataset_id, prop_id;
  hsize_t dims[]={1};

  /* Open scalar , which is like a small dataset  */
  dataset_id = H5Dopen(loc_id, name); 
  if( dataset_id < 0 ) { 
    fprintf(stderr,"myH5_read_scalar(): Cannot open dataset named  %s \n", name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create memory and file spaces, removing extents if we are not the master node : */
  memspace_id  = H5Screate_simple(1, dims, NULL); 
  filespace_id = H5Dget_space(dataset_id);

  /* Write the dataset using defined dataspace and default properties. */
  prop_id = H5Pcreate(H5P_DATASET_XFER);

#if( USE_MPI_IO ) 
  H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
#endif 

  H5Dread(dataset_id, type_id, memspace_id, filespace_id, prop_id, value);

  /* Close all open object handles : */
  H5Dclose(dataset_id);	 
  H5Sclose(filespace_id);	 
  H5Sclose(memspace_id);	 
  H5Pclose(prop_id);

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_read_scalar2_orig(): 
 ---------------------------
   -- different than myH5_read_scalar() in that it selects just the first element 
      of a potentially larger dataset on disk than a scalar;
   -- driver routine for reading a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for opening dataset, opening dataspace, reading dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- unlike myH5_write_scalar(), all processes read the scalar data (this eliminates
           a broadcast of the data);
*******************************************************************************************/
void myH5_read_scalar2_orig( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  hid_t   filespace_id, memspace_id, dataset_id, prop_id;
  hsize_t dims[]={1}, scal_count[1]={1};
  hsize_t scal_offset[1]={0};

  /* Open scalar , which is like a small dataset  */
  dataset_id = H5Dopen(loc_id, name); 
  if( dataset_id < 0 ) { 
    fprintf(stderr,"%s(): Cannot open dataset named  %s \n",__func__, name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create memory and file spaces, removing extents if we are not the master node : */
  memspace_id  = H5Screate_simple(1, dims, NULL); 
  filespace_id = H5Dget_space(dataset_id);

#if( USE_CHUNKS ) 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, scal_offset, NULL, scal_count, dims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, scal_offset, NULL, dims, NULL);
#endif 


  /* Write the dataset using defined dataspace and default properties. */
  prop_id = H5Pcreate(H5P_DATASET_XFER);

  /* Currently only configured to read a single restart file for MPI runs: */
#if( USE_MPI_IO ) 
  H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
#endif 

  H5Dread(dataset_id, type_id, memspace_id, filespace_id, prop_id, value);

  /* Close all open object handles : */
  H5Dclose(dataset_id);	 
  H5Sclose(filespace_id);	 
  H5Sclose(memspace_id);	 
  H5Pclose(prop_id);

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_read_scalar2_gather(): 
 ---------------------------
   -- version that "broadcasts" the result read in from master process; 
   -- different than myH5_read_scalar() in that it selects just the first element 
      of a potentially larger dataset on disk than a scalar;
   -- driver routine for reading a scalar under given group/file id "loc_id";  
      the scalar is identified by the name "name", and the scalar
      is a single value "value" of type "type_id"; 
   -- responsible for opening dataset, opening dataspace, reading dataset,
      and closing the dataset and the dataspace; 
   -- handles differences between serial and MPI calls; 
        -- unlike myH5_write_scalar(), all processes read the scalar data (this eliminates
           a broadcast of the data);
*******************************************************************************************/
void myH5_read_scalar2_gather( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  if( myid == out_pid[OUT_HDF5] ) { 
    herr_t  status;
    status = H5LTread_dataset( loc_id, name ,  type_id, value );
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with reading dataset: %s  %d \n",__func__,myid,nstep,n_substep,name,status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }
  }

#if(USEMPI) /* nothing needs to happen for non-mpi jobs */
  exit_status();
  if( type_id == H5T_NATIVE_DOUBLE ) { 
    MPI_Bcast(value, 1, MPI_DOUBLE, out_pid[OUT_HDF5], MPI_COMM_WORLD);
  } 
  else if(  type_id == H5T_NATIVE_INT ) { 
    MPI_Bcast(value, 1, MPI_INT, out_pid[OUT_HDF5], MPI_COMM_WORLD);
  }
  else { 
    fprintf(stderr,"%s(): Unknown data type %s   type=%d  \n",__func__, name,((int) type_id));
    fflush(stderr);     fail(FAIL_HDF,0);
  }
#endif

  return;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_read_basic(): 
 ---------------------------
   -- driver routine for reading any simple dataset under given group/file id "loc_id";  
      the dataset is identified by the name "name", and the data 
      is stored under "value" of type "type_id"; 
   -- responsible for assuming a dataspace, opening dataset, reading dataset,
      and closing the dataset and the dataspace; 
*******************************************************************************************/
void myH5_read_basic( hid_t loc_id, char *name, hid_t type_id, void *value ) 
{
  hid_t   space_id, dataset_id;

  /* Open the dataset attached to name : */
  dataset_id = H5Dopen(loc_id, name);

  /* Use the dataspace that comes with the dataset, hopefully it is just a single number: */
  space_id  =  H5Dget_space(dataset_id);  

  /* Read the dataset using the just defined dataspace (which H5S_ALL does)  */
  H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);

  /* Close open object handles : */
  H5Dclose(dataset_id);	 
  H5Sclose(space_id);	 

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_read_gfunc_orig(): 
 ---------------------------
   -- driver routine for reading a gridfunction that is simple dataset under 
      given group/file id "loc_id";  
   -- the dataset is identified by the name "name", and the data 
      is stored under "value" of type "type_id"; the memory space to read i
   -- responsible for assuming a dataspace, opening dataset, reading dataset,
      and closing the dataset and the dataspace; 
   -- assumes same space layout as myH5_write_gfunc();
*******************************************************************************************/
void myH5_read_gfunc_orig(hid_t loc_id, char *name, double *value ) 
{
  hid_t   memspace_id, filespace_id, dataset_id, prop_id;

  /************************************************************************************
     Open the data set: 
   ************************************************************************************/
  /* Open existing dataset :  */
  dataset_id = H5Dopen(loc_id, name); 

  /************************************************************************************
    Select the hyperslab, or section of space, to which our gridfunction will be written
  ************************************************************************************/
  filespace_id = H5Dget_space(dataset_id);  

#if( USE_CHUNKS ) 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, h5_memdims);
#else 
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, h5_memdims, NULL);
#endif 
  
  memspace_id  = H5Screate_simple(HDF_RANK, h5_memdims , NULL); 
  
  /************************************************************************************
    Write the hyperslab to the dataset : 
  ************************************************************************************/
  /* Setup the properties of the write based upon type of run : */
  prop_id = H5Pcreate(H5P_DATASET_XFER);

#if( USE_MPI_IO ) 
  H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
#endif

  /* Write the dataset using defined dataspace and properties. */
  H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, prop_id, value);

  /* Close the open object handles : */
  H5Dclose(dataset_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
  H5Pclose(prop_id);

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_read_gfunc_gather(): 
 ---------------------------
   -- "gather" version of the routine;
   -- driver routine for reading a gridfunction that is simple dataset under 
      given group/file id "loc_id";  
   -- the dataset is identified by the name "name", and the data 
      is stored under "value" of type "type_id"; the memory space to read i
   -- responsible for assuming a dataspace, opening dataset, reading dataset,
      and closing the dataset and the dataspace; 
   -- assumes same space layout as myH5_write_gfunc();
*******************************************************************************************/
void myH5_read_gfunc_gather(hid_t loc_id, char *name, double *value ) 
{
  int i,j,k,temp_id;
  ulint n; 
  static ulint *sort; 
  static double *data;
  static double *sendv;
  static usint local_first_time = 1; 

  TRACE_BEG;

#if( USEMPI ) 
  if( local_first_time ) { 
    if( myid == out_pid[OUT_HDF5] ) {
      ALLOC_ARRAY(sort,n_cells_glob); 
      ALLOC_ARRAY(data,n_cells_glob); 
      ALLOC_ARRAY(sendv,n_cells_glob); 

      /* Construct the array "sort" which contains how the gathered data "sendv" should be sorted */
      n = 0;
      for(temp_id = 0; temp_id < numprocs; temp_id++) {
	for(i=0; i<N1; i++) for(j=0; j<N2; j++) for(k=0; k<N3; k++) {
	      sort[n] = get_global_index_pid(temp_id,i,j,k);
	      n++;
	    }
      }
    }
    local_first_time = 0 ; 
  }

  if( myid == out_pid[OUT_HDF5] ) { 
    herr_t  status;
    status = H5LTread_dataset_double( loc_id, name , data );
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with reading dataset: %s  %d \n",__func__,myid,nstep,n_substep,name,status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }
    for(n=0; n<n_cells_glob; n++) {
      sendv[n] = data[sort[n]];
    }
  }

  exit_status();
  MPI_Scatter( sendv, NCELLS, MPI_DOUBLE, value, NCELLS, MPI_DOUBLE, out_pid[OUT_HDF5], MPI_COMM_WORLD); 
  //  MPI_Barrier(MPI_COMM_WORLD);

#else 

  myH5_read_gfunc_orig(loc_id, name, value );

#endif

  TRACE_END;

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_Fcreate_orig():
 --------------------------
   -- wrapper for creating new hdf5 files; 
   -- handles differences between MPI and serial calls; 
*******************************************************************************************/
hid_t myH5_Fcreate_orig(char *name, int overwrite )
{
  hid_t prop_id, file_id; 

  prop_id = H5Pcreate(H5P_FILE_ACCESS);

#if( USE_MPI_IO ) 
  H5Pset_fapl_mpio(prop_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
  
  if( overwrite ) { file_id = H5Fcreate(name, H5F_ACC_TRUNC,  H5P_DEFAULT, prop_id); }
  else            { file_id = H5Fcreate(name, H5F_ACC_EXCL ,  H5P_DEFAULT, prop_id); }

  if( file_id  < 0 ) { 
    fprintf(stderr,"Cannot create file %s \n", name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  H5Pclose(prop_id);

  return( file_id ); 
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_Fcreate_gather():
 --------------------------
   -- wrapper for creating new hdf5 files; 
   -- handles differences between MPI and serial calls; 
*******************************************************************************************/
hid_t myH5_Fcreate_gather(char *name, int overwrite )
{
  hid_t file_id = ((hid_t) 0); 

  if( myid == out_pid[OUT_HDF5] ) { 
    if( overwrite ) { file_id = H5Fcreate(name, H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT); }
    else            { file_id = H5Fcreate(name, H5F_ACC_EXCL ,  H5P_DEFAULT, H5P_DEFAULT); }
    if( file_id  < 0 ) { 
      fprintf(stderr,"%s(): Cannot create file %s \n", __func__,name);
      fflush(stderr);     fail(FAIL_HDF,0);
    }
  }

  return( file_id ); 
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_Fopen():
 --------------------------
   -- wrapper for opening hdf5 files; 
   -- handles differences between MPI and serial calls; 
*******************************************************************************************/
hid_t myH5_Fopen(char *name)
{
  hid_t prop_id, file_id; 

  prop_id = H5Pcreate(H5P_FILE_ACCESS);

  /* Currently if we open a file we are modifying it, and this needs to be done with MPI-IO 
     when using MPI */
#if( USE_MPI_IO ) 
  H5Pset_fapl_mpio(prop_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
  
  file_id = H5Fopen(name, H5F_ACC_RDWR, prop_id);
  if( file_id  < 0 ) { 
    fprintf(stderr,"Cannot open file %s \n", name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  H5Pclose(prop_id);

  return( file_id ); 

}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_get_full_pathname():
 ------------
   -- determines the full path of the object "loc_id" from the root directory; 
   -- returns with the string; 
*******************************************************************************************/
static char *myH5_get_full_pathname( hid_t loc_id )
{
  size_t size1, size2; 
  char *name;
  
  size1 = sizeof(char)*100;
  name = (char *) malloc(size1);
  size2 = H5Iget_name(loc_id, name, size1 );
  if( size2 > size1 ) { 
    name = (char *) realloc(name,size2);
    size1 = H5Iget_name(loc_id, name, size2 );
  }
  return( name ) ; 
}

/*******************************************************************************************/
/*******************************************************************************************
 myH5_get_simple_dims():
 ------------
   -- given a dataspace object id, returns with a pointer to the dimension 
       extent array and a pointer to the number of dims; 
   -- ignores the maximum allowable dimensions since we are usually using arrays of 
       static size; 
*******************************************************************************************/
hsize_t *myH5_get_simple_dims( hid_t dataspace_id,  int *ndims )
{
  hsize_t *dims, *maxdims;

  /* Find dimensionality of object */
  *ndims  = H5Sget_simple_extent_ndims(dataspace_id); 
  dims    = (hsize_t *) calloc(*ndims,sizeof(hsize_t)); 
  maxdims = (hsize_t *) calloc(*ndims,sizeof(hsize_t));
  *ndims  = H5Sget_simple_extent_dims(dataspace_id, dims, maxdims); 

  free(maxdims);
  return( dims ); 

}

/*******************************************************************************************/
/*******************************************************************************************
 set_hdf5_gfuncs(); 
 ------------------
   -- this routine houses all the information for determining what is written to 
      the hdf5 dump files; 
   -- it is responsible for setting the gridfunction names and their values;
   -- all gridfunctions are stored which allows us to optimize their calculation 
       (eliminate repeated get_geometry() calls, etc); 
*******************************************************************************************/
static void set_hdf5_gfuncs(void)
{
  int i,j,k,l,n, ind,ii,jj; 
  int ih,jh,kh;
  double gamma, prim[NP], jcon_l[NDIM], jsq; 
  struct of_geom  *geom;


  ind = 0; 


  /*******************************************************************************
     Loop over all physical cells, storing output functions along the way. 
  ********************************************************************************/
  LOOP { 
    n = 0 ; 


    /****************************************************************************/
    /****************************************************************************
          NON-EXCISED output:
    *****************************************************************************/
#if( USE_MASK ) 
    //    if( !(evol_mask[ncurr][i][j][k] == MASK_EXCISED) ) { 
    if( 1 ) { 
#else 
      if(1) {
#endif 

	get_geometry(i,j,k,CENT,ncurr,geom);

	/* set local primitives to eliminate memory lookups */
	PLOOP { prim[l] = p[i][j][k][l]; }  

	/* Primitive variables: */
	PLOOP { f_hdf[n++][ind] = prim[l]; }

	/* bsq */
	f_hdf[n++][ind] = bsq_calc(prim, geom);
    
	/* gamma or Lorentz factor */ 
	gamma_calc(prim, geom, &gamma);
	f_hdf[n++][ind] = gamma;


	/* Divergence of B^i using Flux-CT stencil (should be 0 to roundoff w/ our method) 
	   Will be large along lower physical boundaries due to our choice of stencils. 
	*/
	f_hdf[n++][ind] = divb_calc(i,j,k);


	/* Divergence of B^i using cell centered stencil (may not be 0 to roundoff w/ our method)
	   Will be large along lower and upper physical boundaries due to our choice of stencils.
	*/
	f_hdf[n++][ind] = divb_cen_calc(i,j,k);


	/* Metric: */
	f_hdf[n++][ind] = geom->gcov[0][0];
	f_hdf[n++][ind] = geom->gcov[0][1];
	f_hdf[n++][ind] = geom->gcov[0][2];
	f_hdf[n++][ind] = geom->gcov[0][3];
	f_hdf[n++][ind] = geom->gcov[1][1];
	f_hdf[n++][ind] = geom->gcov[1][2];
	f_hdf[n++][ind] = geom->gcov[1][3];
	f_hdf[n++][ind] = geom->gcov[2][2];
	f_hdf[n++][ind] = geom->gcov[2][3];
	f_hdf[n++][ind] = geom->gcov[3][3];
	f_hdf[n++][ind] = geom->g;

	/* Error function */
	f_hdf[n++][ind] = (double) pflag[i][j][k];
    
#if( USE_MASK ) 
	f_hdf[n++][ind] = (double) evol_mask[ncurr][i][j][k];
#else 
	n++;
#endif

//	for(ii=0;ii<N_NFAIL;ii++) {
//	  f_hdf[n++][ind] = (double) nfail[ii][i][j][k];
//	}
//

#if( USE_KINEMATIC_VISCOSITY ) 
	for(ii=0;ii<NDIM;ii++)  { 
	  f_hdf[n++][ind] = visc_source[i][j][k][ii];
	}
#endif

	/* Contravariant current components : */
#if( CALC_CURRENT )     
	for(l=0;l<NDIM;l++) { 
	  f_hdf[n++][ind] = jcon_l[l] = jcon[i][j][k][l]; 
	}

	/* Jsq = Magnitude squared current : */
	jsq = 0.; 
	for(ii=0;ii<NDIM;ii++)  for(jj=0;jj<NDIM;jj++) { 
	    jsq += geom->gcov[ii][jj] * jcon_l[ii] * jcon_l[jj];
	  }
	f_hdf[n++][ind] = jsq; 
#endif
      }
      /****************************************************************************/
      /****************************************************************************
          EXCISED output:
      *****************************************************************************/
      else { 

	PLOOP { prim[l] = p[i][j][k][l]; }  
	PLOOP { f_hdf[n++][ind] = prim[l]; }
	/* bsq */
	f_hdf[n++][ind] = 0.;
	/* gamma or Lorentz factor */ 
	f_hdf[n++][ind] = 1.;
	/* Divergence of B^i using Flux-CT stencil (should be 0 to roundoff w/ our method) 
	   Will be large along lower physical boundaries due to our choice of stencils. 
	*/
	f_hdf[n++][ind] = 0.;
	/* Divergence of B^i using cell centered stencil (may not be 0 to roundoff w/ our method)
	   Will be large along lower and upper physical boundaries due to our choice of stencils.
	*/
	f_hdf[n++][ind] = 0.;
	/* Metric: */
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	f_hdf[n++][ind] = 0.;
	/* Error function */
	f_hdf[n++][ind] = (double) pflag[i][j][k];
#if( USE_MASK ) 
	f_hdf[n++][ind] = (double) evol_mask[ncurr][i][j][k];
#else 
	n++;
#endif
//	for(ii=0;ii<N_NFAIL;ii++) {
//	  f_hdf[n++][ind] = nfail[ii][i][j][k];
//	}

	/* Contravariant current components : */
#if( CALC_CURRENT )     
	for(l=0;l<NDIM;l++) { 
	  f_hdf[n++][ind] = jcon_l[l] = jcon[i][j][k][l]; 
	}
	/* Jsq = Magnitude squared current : */
	jsq = 0.; 
	f_hdf[n++][ind] = jsq; 
#endif

      }

    ind++;
  }


  return;
}


/**********************************************************************************************/
/**********************************************************************************************
  restart_check_hdf5_orig(): 
  -------------------
      -- returns 1 if there is a valid hdf5 restart dump file, 0 otherwise
**********************************************************************************************/
int restart_check_hdf5_orig( char *hdf_name ) 
{
  FILE *fp;

  /* Set the name of the restart file : */
#if( USE_MPI_IO )
  sprintf(hdf_name,"%s/rdump_start.h5", DIR_out[OUT_HDF5]);
#else
  sprintf(hdf_name,"%s/rdump_start.p%05d.h5",DIR_out[OUT_HDF5],myid);
#endif

  /******************************************************************************
    Proceed only if file exists, else tell caller that we are missing it 
  ******************************************************************************/ 
  if( myid == printer_pid ) {   fprintf(stdout,"Trying to read  %s .... \n",hdf_name); fflush(stdout);  }

  fp = fopen(hdf_name,"r") ;
  if( fp == NULL ) { 
    if( myid == printer_pid ) {  fprintf(stdout,"No h5 restart file named %s ... starting anew!!! \n", hdf_name) ; fflush(stdout); }
    return(0) ;
  }
  fclose(fp);

  if( H5Fis_hdf5(hdf_name) <= 0 ) { 
    fprintf(stderr,"Restart file is not in hdf format  %s \n", hdf_name) ; fflush(stderr);
    return(0) ;
  }
  return(1);
}

/**********************************************************************************************/
/**********************************************************************************************
  restart_check_hdf5_gather(): 
  -------------------
      -- returns 1 if there is a valid hdf5 restart dump file, 0 otherwise
**********************************************************************************************/
int restart_check_hdf5_gather( char *hdf_name )
{
  int retval; 

  retval = 0;

  TRACE_BEG;

  if( myid == out_pid[OUT_HDF5] ) {
    FILE *fp;

    /* Set the name of the restart file : */
    sprintf(hdf_name,"%s/rdump_start.h5", DIR_out[OUT_HDF5]);

    /******************************************************************************
    Proceed only if file exists, else tell caller that we are missing it 
    ******************************************************************************/ 
    fprintf(stdout,"Trying to read  %s .... \n",hdf_name); fflush(stdout);

    fp = fopen(hdf_name,"r") ;
    if( fp == NULL ) { 
      fprintf(stderr,"No h5 restart file named %s ... starting anew!!! \n", hdf_name) ; fflush(stderr);
      retval = 0 ; 
    }
    else { 
      fclose(fp);
      retval = 1; 
    }

    if( retval ) { 
      if( H5Fis_hdf5(hdf_name) <= 0 ) { 
	fprintf(stderr,"Restart file is not in hdf format  %s \n", hdf_name) ; fflush(stderr);
	retval = 0 ; 
      }
      else { 
	retval = 1 ; 
      }
    }
  }

  mpi_global_imax(&retval);
  TRACE_END;

  return(retval); 
}


/**********************************************************************************************/
/**********************************************************************************************
  restart_read_hdf5(): 
  -------------------
      -- responsible for reading in gridfunctions and scalar data from a restart 
           file; 
      -- assumes that the restart file is in the directory named in DIR_out[OUT_HDF]
         and is named "rdump_start.h5"
      -- essentially dump_hdf5() except "create" -> "open"  and "write" -> "read";

**********************************************************************************************/
void restart_read_hdf5( char *hdf_name ) 
{
  int i,j,k,l,ind,n1,n2,n3,ntmp,ng,np;
  hid_t file_id, header_id, grid_id, ioparams_id;
  char headvar_name[200];
  FILE *fp;

  file_id     = (hid_t) 0;
  header_id   = (hid_t) 0;
  grid_id     = (hid_t) 0;
  ioparams_id = (hid_t) 0;

#if( GATHER_IO ) 
  if( myid == out_pid[OUT_HDF5] ) { 
#else
   {
#endif
  /******************************************************************************
    Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Open existing file : wrapper handles error checking */ 
  file_id = myH5_Fopen(hdf_name);

  /* Create a subgroup to which we attach attributes or header information */
  header_id = H5Gopen(file_id, "Header");
  if( header_id  < 0 ) { 
    fprintf(stderr,"restart_read_hdf5(): Cannot open Header group in file %s \n", hdf_name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create a header subgroup that contains the grid parameters: */
  grid_id = H5Gopen(header_id, "Grid");
  if( grid_id  < 0 ) { 
    fprintf(stderr,"restart_read_hdf5(): Cannot open Grid group in file %s \n", hdf_name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create a header subgroup that contains the I/O parameters: */
  ioparams_id = H5Gopen(header_id, "IO");
  if( ioparams_id  < 0 ) { 
    fprintf(stderr,"restart_read_hdf5(): Cannot open IO group in file %s \n", hdf_name);
    fflush(stderr);     fail(FAIL_HDF,0);
  }
 }


  /******************************************************************************
    Dump the grid parameters into the HDF file : 
          -- make sure that reading in the restart parameters does not assume 
             a certain domain decomposition
          -- we want to be able to restart anyway we want while just preserving 
             the PHYSICAL extent of the grid, and not how it's subdivided on a cluster;
  ******************************************************************************/ 
//  myH5_read_scalar2(grid_id, "N1"              , H5T_NATIVE_INT   ,  &n1              );
//  myH5_read_scalar2(grid_id, "N2"              , H5T_NATIVE_INT   ,  &n2              );
//  myH5_read_scalar2(grid_id, "N3"              , H5T_NATIVE_INT   ,  &n3              );
//  myH5_read_scalar2(grid_id, "NG"              , H5T_NATIVE_INT   ,  &ng              );
//  myH5_read_scalar2(grid_id, "NP"              , H5T_NATIVE_INT   ,  &np              );
//  myH5_read_scalar2(grid_id, "totalsize0"      , H5T_NATIVE_INT   ,  &(totalsize[0])  );
//  myH5_read_scalar2(grid_id, "totalsize1"      , H5T_NATIVE_INT   ,  &(totalsize[1])  );
//  myH5_read_scalar2(grid_id, "totalsize2"      , H5T_NATIVE_INT   ,  &(totalsize[2])  );
//  myH5_read_scalar2(grid_id, "totalsize3"      , H5T_NATIVE_INT   ,  &(totalsize[3])  );
//  myH5_read_scalar2(grid_id, "cpupos0"         , H5T_NATIVE_INT   ,  &(cpupos[0])     );
//  myH5_read_scalar2(grid_id, "cpupos1"         , H5T_NATIVE_INT   ,  &(cpupos[1])     );
//  myH5_read_scalar2(grid_id, "cpupos2"         , H5T_NATIVE_INT   ,  &(cpupos[2])     );
//  myH5_read_scalar2(grid_id, "cpupos3"         , H5T_NATIVE_INT   ,  &(cpupos[3])     );
//  myH5_read_scalar2(grid_id, "startx0"         , H5T_NATIVE_DOUBLE,  &(startx[0])     );
//  myH5_read_scalar2(grid_id, "startx1"         , H5T_NATIVE_DOUBLE,  &(startx[1])     );
//  myH5_read_scalar2(grid_id, "startx2"         , H5T_NATIVE_DOUBLE,  &(startx[2])     );
//  myH5_read_scalar2(grid_id, "startx3"         , H5T_NATIVE_DOUBLE,  &(startx[3])     );
  myH5_read_scalar2(grid_id, "dx0"             , H5T_NATIVE_DOUBLE,  &(dx[0])         );
//  myH5_read_scalar2(grid_id, "dx1"             , H5T_NATIVE_DOUBLE,  &(dx[1])         );
//  myH5_read_scalar2(grid_id, "dx2"             , H5T_NATIVE_DOUBLE,  &(dx[2])         );
//  myH5_read_scalar2(grid_id, "dx3"             , H5T_NATIVE_DOUBLE,  &(dx[3])         );
//  myH5_read_scalar2(grid_id, "gridlength0"     , H5T_NATIVE_DOUBLE,  &(GridLength[0]) );
//  myH5_read_scalar2(grid_id, "gridlength1"     , H5T_NATIVE_DOUBLE,  &(GridLength[1]) );
//  myH5_read_scalar2(grid_id, "gridlength2"     , H5T_NATIVE_DOUBLE,  &(GridLength[2]) );
//  myH5_read_scalar2(grid_id, "gridlength3"     , H5T_NATIVE_DOUBLE,  &(GridLength[3]) );
  myH5_read_scalar2(grid_id, "t"               , H5T_NATIVE_DOUBLE,  &t               );
//  myH5_read_scalar2(grid_id, "cour"            , H5T_NATIVE_DOUBLE,  &cour     	      );
//  myH5_read_scalar2(grid_id, "gam"             , H5T_NATIVE_DOUBLE,  &gam     	      );
//  myH5_read_scalar2(grid_id, "a"               , H5T_NATIVE_DOUBLE,  &a	    	      );
//  myH5_read_scalar2(grid_id, "h_slope"         , H5T_NATIVE_DOUBLE,  &h_slope	      );
//  myH5_read_scalar2(grid_id, "th_cutout"       , H5T_NATIVE_DOUBLE,  &th_cutout       );
//  myH5_read_scalar2(grid_id, "th_beg"          , H5T_NATIVE_DOUBLE,  &th_beg  	      );
//  myH5_read_scalar2(grid_id, "th_end"          , H5T_NATIVE_DOUBLE,  &th_end	      );
//  myH5_read_scalar2(grid_id, "X1_slope"        , H5T_NATIVE_DOUBLE,  &X1_slope        );
//  myH5_read_scalar2(grid_id, "X1_0"            , H5T_NATIVE_DOUBLE,  &X1_0	      );
//  myH5_read_scalar2(grid_id, "R0"              , H5T_NATIVE_DOUBLE,  &R0	      );
//  myH5_read_scalar2(grid_id, "Rin"             , H5T_NATIVE_DOUBLE,  &Rin	      );
//  myH5_read_scalar2(grid_id, "Rout"            , H5T_NATIVE_DOUBLE,  &Rout	      );
//  myH5_read_scalar2(grid_id, "r_isco"          , H5T_NATIVE_DOUBLE,  &r_isco	      );
//  myH5_read_scalar2(grid_id, "r_horizon"       , H5T_NATIVE_DOUBLE,  &r_horizon       );
//  myH5_read_scalar2(grid_id, "n_within_horizon", H5T_NATIVE_INT   ,  &n_within_horizon);
//  myH5_read_scalar2(grid_id, "diag3_exponent"  , H5T_NATIVE_INT   ,  &diag3_exponent);

  /* Calculate derived global variables that are not in the dumps : */
//  dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */
//  th_length = th_end - th_beg; 

  /******************************************************************************
    Dump the IO parameters into the HDF header: 
  ******************************************************************************/ 
  myH5_read_scalar2(ioparams_id,  "nstep"          , H5T_NATIVE_INT ,  &nstep          );
  myH5_read_scalar2(ioparams_id,  "n_restart"      , H5T_NATIVE_INT ,  &n_restart      );
  myH5_read_scalar2(ioparams_id,  "N_hist_dump"    , H5T_NATIVE_INT ,  &N_hist_dump    );
  myH5_read_scalar2(ioparams_id,  "N_hist_dump_tot", H5T_NATIVE_INT ,  &N_hist_dump_tot);
  myH5_read_scalar2(ioparams_id,  "N_OUT_TYPES"    , H5T_NATIVE_INT ,  &ntmp           );

  if( ntmp < N_OUT_TYPES ) { 
    fprintf(stderr,"restart_read_hdf5(): N_OUT_TYPES > restart value, using default values: %d %d\n",
	    N_OUT_TYPES,ntmp);
    fflush(stderr); 
  }
  
  if( ntmp > N_OUT_TYPES ) { 
    fprintf(stderr,"restart_read_hdf5(): N_OUT_TYPES < restart value, exiting!! : %d %d\n",
	    N_OUT_TYPES,ntmp);
    fflush(stderr);    fail(FAIL_HDF,0); 
  }

  /* It is ok if a restart file lacks one of the following, we'll just use the default value: */
  
  for(i=0;i<ntmp;i++) { 
    sprintf(headvar_name, "N_out%1d", i);
    myH5_read_scalar2(ioparams_id, headvar_name, H5T_NATIVE_INT   ,  &(N_out[i]) );
  }
  for(i=0;i<ntmp;i++) { 
    sprintf(headvar_name, "T_out%1d", i);
    myH5_read_scalar2(ioparams_id, headvar_name, H5T_NATIVE_DOUBLE,  &(T_out[i]) );
  }
  for(i=0;i<ntmp;i++) { 
    sprintf(headvar_name, "DT_out%1d",i);
    myH5_read_scalar2(ioparams_id, headvar_name, H5T_NATIVE_DOUBLE,  &(DT_out[i]));
  }
  
#if( COORD_RADIUS_DIM > 1 )
  myH5_read_scalar2(ioparams_id,  "n_r_bins"      , H5T_NATIVE_INT ,    &n_r_bins    );
  myH5_read_scalar2(ioparams_id,  "n_phi_bins"    , H5T_NATIVE_INT ,    &n_phi_bins  );
  myH5_read_scalar2(ioparams_id,  "r0_bins"       , H5T_NATIVE_DOUBLE , &r0_bins     );
  myH5_read_scalar2(ioparams_id,  "r_min_bins"    , H5T_NATIVE_DOUBLE , &r_min_bins  );
  myH5_read_scalar2(ioparams_id,  "r_max_bins"    , H5T_NATIVE_DOUBLE , &r_max_bins  );
  myH5_read_scalar2(ioparams_id,  "xp1_min_bins"  , H5T_NATIVE_DOUBLE , &xp1_min_bins);
  myH5_read_scalar2(ioparams_id,  "dxp1_bins"     , H5T_NATIVE_DOUBLE , &dxp1_bins   );
  set_bin_arrays();
#endif
  
  /******************************************************************************
     On restarts, we need to initialize things here after the scalars have been 
     set since some of them are used to define the extents of the grid functions:
  ******************************************************************************/ 
  setup_hdf5();    


  /******************************************************************************
    Set and read gridfunctions : 
         -- we only need the primitive variables
         -- we assume that they are the first grid functions in the list;
  ******************************************************************************/ 
  /* Read only the primitive variables: */
  PLOOP { myH5_read_gfunc(file_id, f_hdf_name[l], f_hdf[l] ); }

  /* Copy over to the primitive variable array: */
  ind = 0;
  LOOP { 
    PLOOP p[i][j][k][l] = f_hdf[l][ind]; 
    ind++; 
  }

#if(RUN_FROM_INTERP_DATA)

/*******************************************************************************
    if reading from interpolated data we have to transform vector components to
    the numerical coordinates
 *******************************************************************************/

  double *rho_tmp, *uu_tmp ;
  double *v1_tmp, *v2_tmp, *v3_tmp;
  double *B1_tmp, *B2_tmp, *B3_tmp;

  struct of_coord *coords;
  double Bcon[NDIM], vcon[NDIM];

  ALLOC_ARRAY(rho_tmp, NCELLS);
  ALLOC_ARRAY(uu_tmp , NCELLS);
  ALLOC_ARRAY(v1_tmp , NCELLS);
  ALLOC_ARRAY(v2_tmp , NCELLS);
  ALLOC_ARRAY(v3_tmp , NCELLS);
  ALLOC_ARRAY(B1_tmp , NCELLS);
  ALLOC_ARRAY(B2_tmp , NCELLS);
  ALLOC_ARRAY(B3_tmp , NCELLS);

  myH5_read_gfunc(file_id, "/rho_dest", rho_tmp ) ;
  myH5_read_gfunc(file_id, "/uu_dest" , uu_tmp ) ;
  myH5_read_gfunc(file_id, "/v1_dest" , v1_tmp ) ;
  myH5_read_gfunc(file_id, "/v2_dest" , v2_tmp ) ;
  myH5_read_gfunc(file_id, "/v3_dest" , v3_tmp ) ;
  myH5_read_gfunc(file_id, "/B1_dest" , B1_tmp ) ;
  myH5_read_gfunc(file_id, "/B2_dest" , B2_tmp ) ;
  myH5_read_gfunc(file_id, "/B3_dest" , B3_tmp ) ;

  ind = 0 ;
  LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);	

    p[i][j][k][RHO] = rho_tmp[ind]  ;
    p[i][j][k][UU]  =  uu_tmp[ind]  ;

    vcon[0] = 0 ;
    vcon[1] = v1_tmp[ind] ;
    vcon[2] = v2_tmp[ind] ;
    vcon[3] = v3_tmp[ind] ;
    Bcon[0] = 0 ;
    Bcon[1] = B1_tmp[ind] ;
    Bcon[2] = B2_tmp[ind] ;
    Bcon[3] = B3_tmp[ind] ;

    transform_rank1con2(coords->dxp_dx,vcon);
    transform_rank1con2(coords->dxp_dx,Bcon);

    p[i][j][k][U1]  =  vcon[1]  ;
    p[i][j][k][U2]  =  vcon[2]  ;
    p[i][j][k][U3]  =  vcon[3]  ;
    p[i][j][k][B1]  =  Bcon[1]  ;
    p[i][j][k][B2]  =  Bcon[2]  ;
    p[i][j][k][B3]  =  Bcon[3]  ;
    ind++;
  }

  DEALLOC_ARRAY(rho_tmp, NCELLS);
  DEALLOC_ARRAY(uu_tmp , NCELLS);
  DEALLOC_ARRAY(v1_tmp , NCELLS);
  DEALLOC_ARRAY(v2_tmp , NCELLS);
  DEALLOC_ARRAY(v3_tmp , NCELLS);
  DEALLOC_ARRAY(B1_tmp , NCELLS);
  DEALLOC_ARRAY(B2_tmp , NCELLS);
  DEALLOC_ARRAY(B3_tmp , NCELLS);

#endif
 
    
  /******************************************************************************
    Close objects : 
  ******************************************************************************/ 
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_HDF5] ) { 
#else
   {
#endif
  H5Gclose(grid_id);      /* Terminate access to the group /Header/Grid. */
  H5Gclose(ioparams_id);  /* Terminate access to the group /Header/IO. */
  H5Gclose(header_id);    /* Terminate access to the group /Header. */
  H5Fclose(file_id);      /* Terminate access to the file. */
 }

  /******************************************************************************
    Dump a restart file to test if we have successfully read in data: 
  ******************************************************************************/ 
  dump_hdf5_gen(2); 

  return;
}


/* 
   dump intermediate monopole cleaner results
 */
void  dump_monopole_cleaner(char *h5filename)
{
  int i,j,k,l, ind ;
  hid_t out_id, header_id, grid_id ;

  struct of_coord *coords;

  double *B1_tmp, *B2_tmp, *B3_tmp, *divb_tmp  ;
  double *x1_tmp, *x2_tmp, *x3_tmp ;

  ALLOC_ARRAY(B1_tmp  ,NCELLS); 
  ALLOC_ARRAY(B2_tmp  ,NCELLS); 
  ALLOC_ARRAY(B3_tmp  ,NCELLS); 
  ALLOC_ARRAY(divb_tmp,NCELLS); 
  ALLOC_ARRAY(x1_tmp  ,NCELLS); 
  ALLOC_ARRAY(x2_tmp  ,NCELLS); 
  ALLOC_ARRAY(x3_tmp  ,NCELLS); 


  ind = 0 ;
  LOOP { 

    get_coord(i,j,k,CENT,ncurr,coords);

    x1_tmp[ind] = coords->x[1]; 
    x2_tmp[ind] = coords->x[2]; 
    x3_tmp[ind] = coords->x[3]; 
 
    B1_tmp[ind] = p[i][j][k][B1] ;
    B2_tmp[ind] = p[i][j][k][B2] ;
    B3_tmp[ind] = p[i][j][k][B3] ;

    divb_tmp[ind] = divb_calc(i,j,k);

    ind++ ;
  }



  out_id = myH5_Fcreate(h5filename, 0);

#if( GATHER_IO )
  if( myid == out_pid[OUT_HDF5] ) {
#else
{
#endif
  header_id = H5Gcreate(out_id, "Header", 0);
  grid_id = H5Gcreate(header_id, "Grid", 0);
 }

  myH5_write_scalar2(grid_id, "t", H5T_NATIVE_DOUBLE,  &t );

  myH5_write_gfunc( out_id, "B1", B1_tmp );
  myH5_write_gfunc( out_id, "B2", B2_tmp );
  myH5_write_gfunc( out_id, "B3", B3_tmp );

  myH5_write_gfunc( out_id, "divb", divb_tmp );

  myH5_write_gfunc( out_id, "x1", x1_tmp );
  myH5_write_gfunc( out_id, "x2", x2_tmp );
  myH5_write_gfunc( out_id, "x3", x3_tmp );


#if( GATHER_IO )
  if( myid == out_pid[OUT_HDF5] ) {
#else
{
#endif
  H5Gclose(grid_id);    /* Terminate access to the group /Grid. */
  H5Gclose(header_id);    /* Terminate access to the group /Header. */
  H5Fclose(out_id) ;
  }



  DEALLOC_ARRAY(B1_tmp,NCELLS); 
  DEALLOC_ARRAY(B2_tmp,NCELLS); 
  DEALLOC_ARRAY(B3_tmp,NCELLS); 
  DEALLOC_ARRAY(divb_tmp,NCELLS); 
  DEALLOC_ARRAY(x1_tmp,NCELLS); 
  DEALLOC_ARRAY(x2_tmp,NCELLS); 
  DEALLOC_ARRAY(x3_tmp,NCELLS); 

  return ;
}


/*******************************************************************************************/
/*******************************************************************************************
 myH5_get_dims():
 ------------
   -- given an object name, returns with a pointer to the dimension 
       extent array and a pointer to the number of dims; 
   -- ignores the maximum allowable dimensions since we are usually using arrays of 
       static size; 
*******************************************************************************************/
hsize_t *myH5_get_dims( hid_t loc_id, char *name,  int *ndims )
{
  hsize_t *dims, *maxdims;
  hid_t space_id, set_id;

#if( GATHER_IO ) 
  if( myid == out_pid[OUT_HDF5] ) { 
#else
    {
#endif
      set_id   = H5Dopen(loc_id, name ); /* Open data set so that we can probe it */
      if( set_id < 0 ) { 
	fprintf(stderr,"myH5_get_dims():  H5Dopen() 2 failed for  %s \n",name);
	fflush(stderr); return(NULL);
      }
      space_id = H5Dget_space(set_id);   /* Get dataspace to get dimensionality info */
      if( space_id < 0 ) { 
	fprintf(stderr,"myH5_get_dims(): get_space() failed for  %s \n",name);
	fflush(stderr); return(NULL);
      }

      /* Find dimensionality of object */
      *ndims  = H5Sget_simple_extent_ndims(space_id); 
      dims    = (hsize_t *) calloc(*ndims,sizeof(hsize_t)); 
      maxdims = (hsize_t *) calloc(*ndims,sizeof(hsize_t));
      *ndims  = H5Sget_simple_extent_dims(space_id, dims, maxdims); 
  
      H5Sclose(space_id);
      H5Dclose(set_id);

      free(maxdims);
    }

#if( GATHER_IO )
    int *dims_tmp, i;
    sync_int_from_rank(ndims,out_pid[OUT_HDF5]);
    ALLOC_ARRAY(dims_tmp,(*ndims));
    if( myid == out_pid[OUT_HDF5] ) { 
      for(i=0;i<(*ndims);i++) {  dims_tmp[i] = (int) dims[i]; }
    }
    else { 
      ALLOC_ARRAY(dims,(*ndims));
    }
    sync_int_vect_from_rank(dims_tmp,(*ndims),out_pid[OUT_HDF5]);
    for(i=0;i<(*ndims);i++) {  dims[i] = (hsize_t) dims_tmp[i]; }
    DEALLOC_ARRAY(dims_tmp,(*ndims));
#endif

  return( dims ); 

}


#undef USE_CHUNKS 
#undef HDF_RANK 
#undef N_HDF_FUNCS 

#else 
void dump_hdf5(void) { return; } 
#endif

  
