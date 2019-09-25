
#include "decs.h"

#if( MAKE_RADFLUX ) 
#include <hdf5.h>

#if( USEMPI )
#include "mpi.h"
#endif 

static hid_t file_id;  /* ID variable for hdf5 statistics file */ 
  
extern void     myH5_write_gfunc( hid_t loc_id, char *name, double *value ) ;
extern void     myH5_write_int_gfunc( hid_t loc_id, char *name, int *value );
extern void myH5_write_scalar2( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
extern hid_t myH5_Fcreate(char *name,int overwrite);
extern void     setup_hdf5(void); 

/**********************************************************************************************/
/**********************************************************************************************
  dump_radflux(): 
 -----------
  -- routine that dumps the radiative flux functions 
  -- need a separate routine since we usually want to dump them more frequently than
     the other routines;

  -- code closely follows that in dump_stat.c 

   -- we clobber existing files;

   -- full dump file names take the form of : 

        <RUN_TAG>.RADFLUX.<dump_id>.h5

          where  
 
            <dump_id> = index number of RADFLUX dump

     -- use the USE_CHUNKS macro to control whether or not to use chunks in parallel IO;
         -- we assume that the chunk size is the memory of one gridfunction in the 
             local domain (i.e. N1*N2*N3);

    -- dumps include standard header information;

    -- dumped gridfunctions include : 
           coolfunc[0-3]  :  fcool ucov[0-3]   

              where fcool is the cooling rate in the frame of the fluid

**********************************************************************************************/
void dump_radflux(void) 
{
  int i,j,k,l,d,g,ind;
  int itype;
  char var_name[200];

  char dump_filename[200]; 

  static int local_first_time = 1; 

  void setup_radflux_hdf5(void);

  TRACE_BEG;
  
  //  exit_status();  /* Make sure that all processes are alive */

  /******************************************************************************
     INITIALIZATIONS : 
  ******************************************************************************/ 

  /* Set up persistent data used for all time by all hdf5 routines: */
  setup_hdf5();    
    
  /*****************************************************************************************
   Write the statistical gridfunctions if it is time : 
  *****************************************************************************************/
  setup_radflux_hdf5();

  /* coolfunc[]  */
//  for(l=0;l<NDIM;l++) { 
//    ind = 0;
//    sprintf(var_name,"coolflux%1d",l); 
//    myH5_write_gfunc(file_id, var_name, coolflux[l]);
//  }
  
  /* coolfunc + v^i  */
  myH5_write_gfunc(file_id, "coolfunc", coolflux[0]);

#if( USE_COOLING_FUNCTION >= 4 )
  myH5_write_gfunc(file_id, "coolfunc_disk",  coolflux[1]);
  myH5_write_gfunc(file_id, "coolfunc_corona",coolflux[2]);
#endif

#if( USE_COOLING_FUNCTION == 5 )
  myH5_write_gfunc(file_id, "urad",  coolflux[3]);
  myH5_write_gfunc(file_id, "Te_eV", radfunc_Te_eV);
#endif

  l = 0; 
  LOOP { 
    coolflux[0][l]  = p[i][j][k][RHO];
    coolflux[1][l]  = p[i][j][k][U1];
    coolflux[2][l]  = p[i][j][k][U2];
    coolflux[3][l]  = p[i][j][k][U3];
    l++;
  }

  myH5_write_gfunc(file_id, "rho",coolflux[0]);
  myH5_write_gfunc(file_id, "v1", coolflux[1]);
  myH5_write_gfunc(file_id, "v2", coolflux[2]);
  myH5_write_gfunc(file_id, "v3", coolflux[3]);


#if( GATHER_IO ) 
  if( myid == out_pid[OUT_RADFLUX] ) { 
#else 
  {
#endif
    H5Fclose(file_id);      /* Terminate access to the file. */
  }

  /******************************************************************************
      FINALIZE: CLOSE FILES, ETC. 
  ******************************************************************************/ 

#if( USEMPI ) 
  //  exit_status();  /* Make sure that all processes are alive */
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 

  TRACE_END;

  return;

}



/**********************************************************************************************/
/**********************************************************************************************
  setup_radflux_hdf5():
 ------------------
    -- creates the hdf5 file, writes the header information and closes the header;

**********************************************************************************************/
void setup_radflux_hdf5(void)
{
  int n1,n2,n3;
  char hdf_name[200];  /* Name of current hdf5 statistics file */
  hid_t header_id; 

  /* Set the name of this time's hdf file : */
#if( USE_MPI_IO || (!USEMPI) || GATHER_IO )
  sprintf(hdf_name,"%s/%s.%s.%06d.h5",DIR_out[OUT_RADFLUX],RUN_TAG,"RADFLUX",N_out[OUT_RADFLUX]); 
#else 
  sprintf(hdf_name,"%s/%s.%s.%06d.p%05d.h5",DIR_out[OUT_RADFLUX],RUN_TAG,"RADFLUX",N_out[OUT_RADFLUX],myid); 
#endif
  
  if( myid == printer_pid ) {   fprintf(stdout,"Dumping RADFLUX %s .... \n",hdf_name); fflush(stdout); }
    
  /******************************************************************************
      Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Create a new, clobbering file using serial or parallel properties */
  file_id = myH5_Fcreate(hdf_name,0);
  
  /* Create a subgroup to which we attach attributes or header information */
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_RADFLUX] ) { 
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


#else 
void dump_radflux(void) { return; } 
#endif


