
#include "decs.h"

#if( MAKE_TIMERS && MAKE_HDF5 ) 
#include <hdf5.h>
#include <hdf5_hl.h>

#if( USEMPI )
#include "mpi.h"
#endif 

#define N_TIMER_HDF5_DIMS (3)
hsize_t dims[N_TIMER_HDF5_DIMS];

extern void myH5_rewrite_scalar( hid_t loc_id, char *name, hid_t type_id, void *value ) ;

/************************************************************************************************************************/
/************************************************************************************************************************

  dump_timers(): 
 --------------
    -- writes timing data to disk ; 
    -- uses a simple MPI algorithm: master node collects data, and master node writes data;

************************************************************************************************************************/
/************************************************************************************************************************/
void dump_timers(void) 
{ 
  int i,j,k;

  TRACE_BEG;

  dims[0] = numprocs;
  dims[1] = N_TIMER_TYPES;
  dims[2] = 2;

  /* Package vector for sending to master process:  */
  int npts  = N_TIMER_TYPES*2;
  double *sendv; 
  ALLOC_ARRAY(sendv,npts);
  k = 0 ; 
  for(i=0;i<N_TIMER_TYPES;i++) for(j=0;j<2;j++)  {  sendv[k++] = elapsed_times[i][j];  } 
  
  for(i=0;i<N_TIMER_TYPES;i++) for(j=0;j<2;j++)  {  elapsed_times[i][j]  = 0. ;  } 


  /* Gather up data at master process and reassign meaning of sendv so that we can use the same HDF5 writing calls for 
      both MPI and non-MPI runs : */
#if( USEMPI )
  int nrecv = npts*numprocs;
  double *recv; 

  if( myid==timers_pid )  {  ALLOC_ARRAY(recv,nrecv);  }

  exit_status();
  MPI_Gather(sendv,npts,MPI_DOUBLE,recv,npts,MPI_DOUBLE,timers_pid,MPI_COMM_WORLD);
  
  DEALLOC_ARRAY(sendv,npts); 
  if( myid == timers_pid ) {  sendv = recv; }

#endif


  /* Write the data to the HDF5 file: */
  if( myid == timers_pid ) {  
    double tloc = t; 
    int nstep_loc = nstep;
    herr_t  status;
    hid_t file_id, header_id;
    static char timer_filename[200];   
    static char dataset_name[200];   
    int timer_iter;
    hsize_t dim_scal[1];
    dim_scal[0] = 1; 
    static usint local_first_time = 1; 
    htri_t is_hdf5; 

    sprintf(timer_filename, "%s/%s.timer.h5", DIR_out[OUT_HDF5],RUN_TAG);

    if( local_first_time ) { 
      is_hdf5 = H5Fis_hdf5( timer_filename );
      if( is_hdf5  > 0 ) { 
	fprintf(stdout,"%s(): Using existing timer file  %s \n",__func__,timer_filename); fflush(stdout); 
      }
      if( is_hdf5  == 0 ) { 
	fprintf(stdout,"%s(): Timer file  %s  seems to exist but is not in hdf5 format, exitting!  \n",__func__,timer_filename); 
	fflush(stdout);  fail(FAIL_HDF,0);
      }
      if( is_hdf5  < 0 ) { 
	fprintf(stdout,"%s(): Timer file  %s  does not exist, will create one now!  \n",__func__,timer_filename); fflush(stdout); 

	file_id = H5Fcreate( timer_filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if( file_id < 0 ) { 
	  fprintf(stdout,"%s(): Problem creating timer dump file %s !!\n",
		  __func__,timer_filename); 	fflush(stdout);  fail(FAIL_HDF,0); 
	}

	/* Write Header information : */
	header_id = H5Gcreate1(file_id, "Header", 0);
	timer_iter = nstep / n_timer_interval ; 
	status = H5LTmake_dataset_int(header_id,"timer_iter",1,dim_scal,&timer_iter);
	i = TIMER_TYPE_TOTAL   ;  status = H5LTmake_dataset_int(header_id,"TIMER_TYPE_TOTAL"  ,1,dim_scal,&i);
	i = TIMER_TYPE_ADVANCE ;  status = H5LTmake_dataset_int(header_id,"TIMER_TYPE_ADVANCE",1,dim_scal,&i);
	i = TIMER_TYPE_METRIC  ;  status = H5LTmake_dataset_int(header_id,"TIMER_TYPE_METRIC" ,1,dim_scal,&i);
	i = TIMER_TYPE_BOUNDS  ;  status = H5LTmake_dataset_int(header_id,"TIMER_TYPE_BOUNDS" ,1,dim_scal,&i);
	i = TIMER_TYPE_DIAG    ;  status = H5LTmake_dataset_int(header_id,"TIMER_TYPE_DIAG"   ,1,dim_scal,&i);
	i = TIMER_TYPE_MPI     ;  status = H5LTmake_dataset_int(header_id,"TIMER_TYPE_MPI"    ,1,dim_scal,&i);
	i = N_TIMER_TYPES      ;  status = H5LTmake_dataset_int(header_id,"N_TIMER_TYPES"     ,1,dim_scal,&i);
	i = TIMER_SECONDS      ;  status = H5LTmake_dataset_int(header_id,"TIMER_SECONDS"     ,1,dim_scal,&i);
	i = TIMER_MICROSECONDS ;  status = H5LTmake_dataset_int(header_id,"TIMER_MICROSECONDS",1,dim_scal,&i);
	H5Gclose(header_id);
	H5Fclose(file_id);
      }
    }

    /* From now on assume the file exists ... */
    file_id = H5Fopen( timer_filename, H5F_ACC_RDWR, H5P_DEFAULT );
    if( file_id < 0 ) { 
      fprintf(stdout,"%s():  file %s does not exist or cannot open...\n",__func__,timer_filename); 
      fflush(stdout);
    }
    status = H5LTread_dataset_int(file_id,"/Header/timer_iter",&timer_iter);
    timer_iter++;
    
    /* Write regular timing dataset  : */
    sprintf(dataset_name, "timer-%010d", timer_iter);
    status = H5LTmake_dataset_double(file_id,dataset_name,N_TIMER_HDF5_DIMS, dims, sendv);
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with making dataset: %s  %d \n",__func__,myid,nstep,n_substep,dataset_name,status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }

    status = H5LTset_attribute_double(file_id,dataset_name,"time",&tloc, 1 );
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with setting time attribute:  %d \n",__func__,myid,nstep,n_substep,status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }
    status = H5LTset_attribute_int(file_id,dataset_name,"nstep",&nstep_loc, 1 );
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with setting nstep attribute:  %d \n",__func__,myid,nstep,n_substep,status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }
    myH5_rewrite_scalar(file_id,"/Header/timer_iter",H5T_NATIVE_INT,&timer_iter);
    
    H5Fclose(file_id);

  } /* myid == timers_pid */

  FREE(sendv);

  TRACE_END;
  return; 

} 

#else
void dump_timers(void) { return; } 
#endif


