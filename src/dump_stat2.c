
#include "decs.h"

#if( MAKE_STAT2 ) 
#include <hdf5.h>

#if( USEMPI )
#include "mpi.h"
#endif 

static hid_t file_id;  /* ID variable for hdf5 statistics file */ 
  
extern void     myH5_write_gfunc( hid_t loc_id, char *name, double *value ) ;
extern void myH5_write_scalar2( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
extern hid_t myH5_Fcreate(char *name,int overwrite);
extern void     setup_hdf5(void); 

/**********************************************************************************************/
/**********************************************************************************************
  dump_stat2(): 
 -----------
  -- routine that dumps the following functions that are evaluated only at the full dt or 
     corrector step (i.e. the functions use the n+1/2 timestep): 

       PL_d_l   =  left  state of P per direction (d) per equation (l)
       PR_d_l   =  right state of P per direction (d) per equation (l)
       UL_d_l   =  left  state of U per direction (d) per equation (l)
       UR_d_l   =  right state of U per direction (d) per equation (l)
       FL_d_l   =  left  state of F per direction (d) per equation (l)
       FR_d_l   =  right state of F per direction (d) per equation (l)
       ctop_d   =  absolute maximum wave speed per direction (d) 
       S_[2,3]  =  source terms (only S2,S3 are non-zero) at cell center

   -- note that the outermost (per direction) boundary values are excluded from the 
      dumps so that we can use the same dump_hdf.c routines, i.e. so that we are
      still dumping N1xN2xN3 cells of information 
      -- in other words, fluxes resides on the cell faces, which exceed by one 
         in each direction the number of cells (e.g. there are N1 physical cells
         in the x1-direction but there are N1+1 faces)
      -- we chose to neglect the outermost (at largest x[1-3]) faces;
       
   -- note that this dump write a LOT of data, like 1KB per cell, which, for e.g. 
      a 192x192x64 run would be about 2.8GB per file; 
      -- so only dump a few of these;

   -- we clobber existing files;

   -- full dump file names take the form of : 

        RUN_TAG.STAT2.<dump_id>.h5

          where  
 
            <dump_id> = time index number 

     -- use the USE_CHUNKS macro to control whether or not to use chunks in parallel IO;
         -- we assume that the chunk size is the memory of one gridfunction in the 
             local domain (i.e. N1*N2*N3);


    -- full hdf5 dump is done whenever this routine is called; 
             
**********************************************************************************************/
void dump_stat2(void) 
{
  int i,j,k,l,d,g,ind;

  char var_name[200];

  void setup_stat2_hdf5(void) ;

  
  /******************************************************************************
     INITIALIZATIONS : 
  ******************************************************************************/ 

  //  exit_status();  /* Make sure that all processes are alive */

  /* Set up persistent data used for all time by all hdf5 routines: */
  setup_hdf5();    


  /* Create hdf5 file and  print header information for the new hdf5 stat file: */
  setup_stat2_hdf5();   


  /*****************************************************************************************
   Write the STAT2 gridfunctions :
  *****************************************************************************************/
  /*   UL_d_l   =  left  state of P per direction (d) per equation (l)    */
  FACE_LOOP  PLOOP { 
    sprintf(var_name,"UL_%d_%d",d,l);  myH5_write_gfunc(file_id, var_name, UL_stat2[d][l] ); 
  }

  /*   UR_d_l   =  right state of P per direction (d) per equation (l)    */
  FACE_LOOP  PLOOP { 
    sprintf(var_name,"UR_%d_%d",d,l);  myH5_write_gfunc(file_id, var_name, UR_stat2[d][l] ); 
  }

  /*   FL_d_l   =  left  state of P per direction (d) per equation (l)    */
  FACE_LOOP  PLOOP { 
    sprintf(var_name,"FL_%d_%d",d,l);  myH5_write_gfunc(file_id, var_name, FL_stat2[d][l] ); 
  }

  /*   FR_d_l   =  right state of P per direction (d) per equation (l)    */
  FACE_LOOP  PLOOP { 
    sprintf(var_name,"FR_%d_%d",d,l);  myH5_write_gfunc(file_id, var_name, FR_stat2[d][l] ); 
  }

  /*   ctop_d   =  absolute maximum wave speed per direction (d) 	     */
  FACE_LOOP  { 
    sprintf(var_name,"ctop_%d",d);  myH5_write_gfunc(file_id, var_name, ctop_stat2[d] ); 
  }

  /*   S_[2,3]  =  source terms (only S2,S3 are non-zero) at cell center  */
  for(d=0; d<2; d++) { 
    sprintf(var_name,"S_%d",(d+2));  myH5_write_gfunc(file_id, var_name, S_stat2[d] ); 
  }

  /*   PL_d_l   =  left  state of P per direction (d) per equation (l)    */
  FACE_LOOP  PLOOP { 
    ind=0;
    LOOP { UL_stat2[0][0][ind] = p_L[i][j][k][d][l];  ind++; }
    sprintf(var_name,"PL_%d_%d",d,l);  myH5_write_gfunc(file_id, var_name, UL_stat2[0][0] ); 
  }
  /*   PR_d_l   =  right state of P per direction (d) per equation (l)    */
  FACE_LOOP  PLOOP { 
    ind=0;
    LOOP { UL_stat2[0][0][ind] = p_R[i][j][k][d][l];  ind++; }
    sprintf(var_name,"PR_%d_%d",d,l);  myH5_write_gfunc(file_id, var_name, UL_stat2[0][0] ); 
  }

  /******************************************************************************
      FINALIZE: CLOSE FILES, ETC. 
  ******************************************************************************/ 
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_STAT2] ) {
#else
  {
#endif
    H5Fclose(file_id);      /* Terminate access to the file. */
  }

#if( USEMPI ) 
  //  exit_status();
  //  MPI_Barrier(MPI_COMM_WORLD);
#endif 

  return;

}


/**********************************************************************************************/
/**********************************************************************************************
  setup_stat2_hdf5():
 ------------------
    -- creates the hdf5 file, writes the header information and closes the header;

**********************************************************************************************/
void setup_stat2_hdf5(void)
{
  int n1,n2,n3;
  char hdf_name[200];  /* Name of current hdf5 statistics file */
  hid_t header_id; 

  /* Set the name of this time's hdf file : */
#if( USE_MPI_IO  || (!USEMPI) || GATHER_IO )
  sprintf(hdf_name,"%s/%s.%s.%06d.h5",DIR_out[OUT_STAT2],RUN_TAG,"STAT2",N_out[OUT_STAT2]); 
#else 
  sprintf(hdf_name,"%s/%s.%s.%06d.p%05d.h5",DIR_out[OUT_STAT2],RUN_TAG,"STAT2",N_out[OUT_STAT2],myid); 
#endif
  
  if( myid == printer_pid ) {   fprintf(stdout,"Dumping STAT2 %s .... \n",hdf_name); fflush(stdout);  }
    
  /******************************************************************************
      Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Create a new, clobbering file using serial or parallel properties */
  file_id = myH5_Fcreate(hdf_name,0);
  
  /* Create a subgroup to which we attach attributes or header information */
#if( GATHER_IO ) 
  if( myid == out_pid[OUT_STAT2] ) { 
#else
   {
#endif
     header_id = H5Gcreate(file_id, "Header", 0);
     if( header_id  < 0 ) { 
       fprintf(stderr,"setup_stat2_hdf5(): Cannot create Header group in file %s \n", hdf_name);
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
void dump_stat2(void) { return; } 
#endif


