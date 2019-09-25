
#include "decs.h"


/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************
  THIS FILE CONTAINS THE FOLLOWING ROUTINES:

    dump_ascii()  : writes gridfunctions to an ASCII formatted file ; 
    gdump()       : write metric functions to an ASCII formatted file;
    dump_sdf()    : mpi-gathers and then write gridfunctions to sdf files;

**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/


#if( MAKE_ASCII )
/**********************************************************************************************/
/**********************************************************************************************
  dump_ascii(): 
 --------
   -- routine that dumps the coordinates and primitive variables every "dump step" in 
      ASCII format;
**********************************************************************************************/
void dump_ascii(void)
{
  int i,j,k,l ;
  double X[NDIM];
  char dfnam[50];
  struct of_coord *coords;
  FILE *dump_file;

  TRACE_BEG;
  /* Name and open file : */ 
  fprintf(stdout,"dump_cnt (pid=%d) %d ...",myid,N_out[OUT_ASCII]) ;
  sprintf(dfnam,"%s/dump%04d.%04d",DIR_out[OUT_ASCII],N_out[OUT_ASCII],myid) ;
  dump_file = fopen(dfnam,"w") ;
  if(dump_file==NULL) {
    fprintf(stderr,"error opening dump file\n") ;
    fail( FAIL_BASIC,0 );
  }

  /* Header information : */ 
  fprintf(fp,FMT_INT_OUT,N1           );
  fprintf(fp,FMT_INT_OUT,N2	    );
  fprintf(fp,FMT_INT_OUT,N3	    );
  fprintf(fp,FMT_INT_OUT,totalsize[1] );
  fprintf(fp,FMT_INT_OUT,totalsize[2] );
  fprintf(fp,FMT_INT_OUT,totalsize[3] );
  fprintf(fp,FMT_DBL_OUT,startx[1]    );
  fprintf(fp,FMT_DBL_OUT,startx[2]    );
  fprintf(fp,FMT_DBL_OUT,startx[3]    );
  fprintf(fp,FMT_DBL_OUT,dx[1]	    );
  fprintf(fp,FMT_DBL_OUT,dx[2]	    );
  fprintf(fp,FMT_DBL_OUT,dx[3]	    );
  fprintf(fp,FMT_DBL_OUT,GridLength[1]);
  fprintf(fp,FMT_DBL_OUT,GridLength[2]);
  fprintf(fp,FMT_DBL_OUT,GridLength[3]);
  fprintf(fp,FMT_DBL_OUT,a	    );
  fprintf(fp,FMT_DBL_OUT,h_slope	    );
  fprintf(fp,FMT_DBL_OUT,X1_slope	    );
  fprintf(fp,FMT_DBL_OUT,X1_0	    );
  fprintf(fp,FMT_DBL_OUT,R0	    );
  fprintf(fp,FMT_DBL_OUT,Rin	    );
  fprintf(fp,FMT_DBL_OUT,Rout	    );
  fprintf(fp,FMT_DBL_OUT,r_isco	    );
  fprintf(fp,FMT_DBL_OUT,r_horizon    );
  fprintf(fp,"\n"); 

  /* Body information : */ 
  DUMP_LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);
    fprintf(fp,FMT_DBL_OUT,coords->x[1]);
    fprintf(fp,FMT_DBL_OUT,coords->x[2]);
    fprintf(fp,FMT_DBL_OUT,coords->x[3]);
    PLOOP fprintf(fp,FMT_DBL_OUT,p[i][j][k][l]) ;
    fprintf(fp,"\n") ;
  }

  /* Close the file and leave : */
  fclose(dump_file) ;
  fprintf(stdout,"...done (pid=%d)! \n",myid) ; fflush(stdout);

  TRACE_END;

  return;
}
#else
void dump_ascii(void) { return; } 
#endif


#if( MAKE_GDUMP )
/**********************************************************************************************/
/**********************************************************************************************
  gdump(): 
 --------
   -- routine that dumps the metric functions  in  ASCII format;
   -- made so that we can compare to harm2d easily ; 
   -- done per processor since MPI-ing it too complicated/tedious w/ all the grid functions;
        -- one can easily reconstruct the global file since the global indices are outputted
**********************************************************************************************/
void gdump( void )
{
  int i,j,k,l,n ;
  int ii,jj,kk;
  FILE *fp;
  char dfnam[50] ;
  struct of_geom *geom;
  struct of_coord *coords;


  TRACE_BEG;

#if( MAKE_HDF5 ) 
  extern void gdump_hdf5(void);

  gdump_hdf5();

#else 

  /* Name the file and open it : */
  system("mkdir -p dumps");  /* make sure directory exists */
  sprintf(dfnam,"dumps/gdump-%04d",myid) ;
  fprintf(stdout,"writing to %s ....\n",dfnam) ; fflush(stdout);

  fp = fopen(dfnam,"w") ;
  if(fp==NULL) {
    fprintf(stderr,"gdump(): error opening file = %s \n",dfnam) ;
    fail( FAIL_BASIC,0 );
  }

  /***********************************************************************
    Header information : 
  ***********************************************************************/ 
  fprintf(fp,FMT_INT_OUT,N1         );
  fprintf(fp,FMT_INT_OUT,N2	    );
  fprintf(fp,FMT_INT_OUT,N3	    );
  fprintf(fp,FMT_INT_OUT,totalsize[1]         );
  fprintf(fp,FMT_INT_OUT,totalsize[2]	    );
  fprintf(fp,FMT_INT_OUT,totalsize[3]	    );
  fprintf(fp,FMT_DBL_OUT,startx[1]    );
  fprintf(fp,FMT_DBL_OUT,startx[2]    );
  fprintf(fp,FMT_DBL_OUT,startx[3]    );
  fprintf(fp,FMT_DBL_OUT,dx[1]	    );
  fprintf(fp,FMT_DBL_OUT,dx[2]	    );
  fprintf(fp,FMT_DBL_OUT,dx[3]	    );
  fprintf(fp,FMT_DBL_OUT,GridLength[1]);
  fprintf(fp,FMT_DBL_OUT,GridLength[2]);
  fprintf(fp,FMT_DBL_OUT,GridLength[3]);
  fprintf(fp,FMT_DBL_OUT,a	    );
  fprintf(fp,FMT_DBL_OUT,h_slope	    );
  fprintf(fp,FMT_DBL_OUT,X1_slope	    );
  fprintf(fp,FMT_DBL_OUT,X1_0	    );
  fprintf(fp,FMT_DBL_OUT,R0	    );
  fprintf(fp,FMT_DBL_OUT,Rin	    );
  fprintf(fp,FMT_DBL_OUT,Rout	    );
  fprintf(fp,FMT_DBL_OUT,r_isco	    );
  fprintf(fp,FMT_DBL_OUT,r_horizon    );
  fprintf(fp,FMT_INT_OUT,n_within_horizon  );
  fprintf(fp,"\n"); 


  /***********************************************************************
    Body information: 
      -- only dump out functions over x1,x2  since they are x3-indep. 
  ***********************************************************************/ 
  k = N3S; 

  GDUMP_LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);

    fprintf(fp,FMT_INT_OUT,(i+globalpos[1])); // 1
    fprintf(fp,FMT_INT_OUT,(j+globalpos[2])); // 2 
    fprintf(fp,FMT_INT_OUT,(k+globalpos[3])); // 3

    fprintf(fp,FMT_DBL_OUT,coords->x[1]); //4
    fprintf(fp,FMT_DBL_OUT,coords->x[2]); //5
    fprintf(fp,FMT_DBL_OUT,coords->x[3]); //6

    fprintf(fp,FMT_DBL_OUT,coords->xp[1]); //7
    fprintf(fp,FMT_DBL_OUT,coords->xp[2]); //8
    fprintf(fp,FMT_DBL_OUT,coords->xp[3]); //9

    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[RR][1]); // 10
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[RR][2]); // 11
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[RR][3]); // 12
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[TH][1]); // 13
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[TH][2]); // 14
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[TH][3]); // 15
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[PH][1]); // 16
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[PH][2]); // 17
    fprintf(fp,FMT_DBL_OUT,coords->dx_dxp[PH][3]); // 18

    n = CONN_ID(i,j,k);

#if(USE_STRICT_ARRAY_BOUNDS)

    for(ii=0; ii<NDIM; ii++) for(jj=0; jj<NDIM; jj++) for(kk=0; kk<NDIM; kk++) {
      fprintf(fp,FMT_DBL_OUT,conn[n][ii][jj][kk]); // 19 - 82 
    }

    for(l=0;l<NPOS; l++) { 
      get_geometry(i,j,k,l,ncurr,geom);
      for(ii=0; ii<NDIM; ii++) for(jj=0; jj<NDIM; jj++) {
	    fprintf(fp,FMT_DBL_OUT,geom->gcon[ii][jj]);    // 83 - 162
	  }
    }

    for(l=0;l<NPOS; l++)  { 
      get_geometry(i,j,k,l,ncurr,geom);
      for(ii=0; ii<NDIM; ii++) for(jj=0; jj<NDIM; jj++) {
	fprintf(fp,FMT_DBL_OUT,geom->gcov[ii][jj]);    // 163 - 242
	}
    }

#else  /* #if(USE_STRICT_ARRAY_BOUNDS) */


    for(ii=0; ii<NDIM*NDIM*NDIM; ii++) 
      fprintf(fp,FMT_DBL_OUT,conn[n][0][0][ii]); // 19 - 82 

    for(l=0;l<NPOS; l++) { 
      get_geometry(i,j,k,l,ncurr,geom);
      for(ii=0; ii<NDIM*NDIM; ii++) 
	fprintf(fp,FMT_DBL_OUT,geom->gcon[0][ii]);    // 83 - 162
    }

    for(l=0;l<NPOS; l++)  { 
      get_geometry(i,j,k,l,ncurr,geom);
      for(ii=0; ii<NDIM*NDIM; ii++) 
	fprintf(fp,FMT_DBL_OUT,geom->gcov[0][ii]);    // 163 - 242
    }


#endif



    for(l=0;l<NPOS; l++)  { 
      get_geometry(i,j,k,l,ncurr,geom);
      fprintf(fp,FMT_DBL_OUT,geom->g);     // 243 - 247  
    }

    fprintf(fp,"\n"); 
  }

  fclose(fp) ;

  fprintf(stdout,".... finished writing %s !\n",dfnam); fflush(stdout);

#endif 

  TRACE_END;

  return;
}
#else 
void gdump(void) { return; } 
#endif



#if( MAKE_SDF ) 
/**********************************************************************************************/
/**********************************************************************************************
  dump_sdf(): 
 --------
   -- routine that generates dumps of various quantities in SDF (binary) format;
   -- also responsible for setting up initial SDF-related things like allocating the 
       work array and naming the grid functions that are being dumped;

   -- if you want to add or change the functions dumped here, then one must change the 
       function names "sdfname[]" here and how they are calculated (i.e. the mapping of
       "ivar" to function definition) in  write_chunk() of harm_mpi.c ; 
 
**********************************************************************************************/
#include <sdf.h>
#define N_SDFS   (NP+2*4+N_NFAIL)     /* Number of sdf files to output */ 

void dump_sdf(void)
{
  int i,j,k,l;
  int temp, ind1, ind2, ind3;

  double t_out, gamma, *ptmp;
  double Xl[NDIM], Xh[NDIM];
  double ftmp;

  static int rank=2;
  static int *sdfshape;
  static char sdfname[N_SDFS][100];
  static double *bbox;

  static int first_time = 1; 

  TRACE_BEG;

  --sdf-output-needs-to-be-updated-to-handle-dynamic-domains----wont-compile
	
  /***************************************************************************************
   Allocate the work array and name dumped functions for present and future use: 
  ***************************************************************************************/
  if( (myid == out_pid[OUT_SDF]) && first_time ) { 

    /* Grid function names : */ 
    sprintf( sdfname[ RHO],  "sdf_rho");
    sprintf( sdfname[  UU],  "sdf_uu" );
    sprintf( sdfname[  U1],  "sdf_v1" );
    sprintf( sdfname[  U2],  "sdf_v2" );
    sprintf( sdfname[  U3],  "sdf_v3" );
    sprintf( sdfname[  B1],  "sdf_B1" );
    sprintf( sdfname[  B2],  "sdf_B2" );
    sprintf( sdfname[  B3],  "sdf_B3" );
    sprintf( sdfname[  NP],  "sdf_bsq" );
    sprintf( sdfname[NP+1],  "sdf_gamma" );
    sprintf( sdfname[NP+2],  "sdf_divb" );
    sprintf( sdfname[NP+3],  "sdf_divbcen" );
    sprintf( sdfname[NP+4+0],  "sdf_jcon0" );
    sprintf( sdfname[NP+4+1],  "sdf_jcon1" );
    sprintf( sdfname[NP+4+2],  "sdf_jcon2" );
    sprintf( sdfname[NP+4+3],  "sdf_jcon3" );
//    sprintf( sdfname[NP+4+4+0],  "sdf_jcon20" );
//    sprintf( sdfname[NP+4+4+1],  "sdf_jcon21" );
//    sprintf( sdfname[NP+4+4+2],  "sdf_jcon22" );
//    sprintf( sdfname[NP+4+4+3],  "sdf_jcon23" );
    for( i = 0; i < N_NFAIL; i++ ) { 
      sprintf( sdfname[NP+2*4+i],  "sdf_nfail%02d",i );
    }

    /* Min. coordinates of bounding box (assuming here that the grid is uniform in X1-3 */
    Xl[0] = 0.;
    SDLOOP1 {   Xl[i] = startx[i] + (NG + 0.5)*dx[i] ; } 

    /* Max. coordinates of bounding box */
    Xh[0] = 0.;
    SDLOOP1 {   Xh[i] = startx[i] + (NG + totalsize[i] - 0.5)*dx[i] ; } 

    /**********************************************************************************
        Need to make 2d or 1d sdfs depending on values of totalsize[1-3]
        Below:
          rank       = number of dimensions in sdfs 
          sdfshape[] = number of points along each dimension
    ************************************************************************************/
    rank = 0; 
    if( totalsize[1] > 1 )  rank++;
    if( totalsize[2] > 1 )  rank++;
    if( totalsize[3] > 1 )  rank++; 

    sdfshape = (int *) calloc(rank,sizeof(int));
    bbox     = (double *) calloc((2*rank),sizeof(double));
    
    l = 0 ; 
    if( totalsize[1] > 1 ) { 
      sdfshape[l] = totalsize[1];
      bbox[2*l  ] = Xl[1]; 
      bbox[2*l+1] = Xh[1];
      l++;
    }
    if( totalsize[2] > 1 ) { 
      sdfshape[l] = totalsize[2];
      bbox[2*l  ] = Xl[2]; 
      bbox[2*l+1] = Xh[2];
      l++;
    }
    if( totalsize[3] > 1 ) { 
      sdfshape[l] = totalsize[3];
      bbox[2*l  ] = Xl[3]; 
      bbox[2*l+1] = Xh[3];
      l++;
    }

    fprintf(stdout,"################  SDF INFO ######################################\n");
    fprintf(stdout,"rank = %d \n", rank);
    fprintf(stdout,"sdfshape = ");
    for( l = 0 ; l < rank; l++ ) { fprintf(stdout," %d ",sdfshape[l]); }
    fprintf(stdout,"\n");
    fprintf(stdout,"bbox = ");
    for( l = 0 ; l < rank; l++ ) { fprintf(stdout," %g %g , ",bbox[2*l],bbox[2*l+1]); }
    fprintf(stdout,"\n");
    for( l = 0 ; l < N_SDFS; l++ ) { 
      fprintf(stdout,"sdfname[%d] = %s \n",l,sdfname[l]); 
    }
    fprintf(stdout,"#################################################################\n");

    first_time = 0 ; 
  }


  fprintf(stdout,"dump_sdf(): Dumping sdf count %d  \n",N_out[OUT_SDF]); fflush(stdout);

  /***************************************************************************************
   Now dump all the functions;  (generic, user does not need to change unless you're doing 
     something strange). 
  ***************************************************************************************/
  t_out = t; 

  for( l = 0; l < N_SDFS ; l++ ) { 
    get_global_gfunc( l, f_sdf ) ; 

    /* Transform order so that X1 loops fastest : */ 
    if( myid == out_pid[OUT_SDF] ) { 
      ind1 = 0; 
      for(i=0;i<totalsize[1];i++) for(j=0;j<totalsize[2];j++) for(k=0;k<totalsize[3];k++) {
	ind2 = i + totalsize[1] * ( j + totalsize[2] * k ); 
	f_sdf2[ind2] = f_sdf[ind1++]; 
      }
      temp = gft_out_bbox( sdfname[l],t_out,sdfshape,rank,bbox,f_sdf2);     
    }

  }
 
  TRACE_END;
  return;

}
#else 
void dump_sdf(void) { return; } 
#endif


/*******************************************************************************************/
/*******************************************************************************************
 c_to_fort_array():
 ----------------
   -- change the ordering of the C-style ordered array to fortran-style ordering; 
   -- 
       C style:   row major, outermost index runs fastest, innermost index runs slowest
 Fortran style:   col major, innermost index runs fastest, outermost index runs slowest

   -- just needs the dimensionality of the grid;
   -- arrays must already be allocated to length "npts";
   -- arrays are 1-d arrays, so we need to be a little tricky with the inversion; 
   -- assumes that    0 < ndims < 4 
   -- returns with  non-zero value if the inversion cannot be done ; 
*******************************************************************************************/
int c_to_fort_array( double *f_c, double *f_fort, int npts, int ndims, int *dims )
{
  int i,j,k, ind_c, ind_f; 

  TRACE_BEG;

  ind_c = ind_f = 0; 

  switch( ndims ) { 

    /* Do not need to do anything for one-dimensional grids */
  case 1 : 
    break;
    
  case 2 : 
    for(i=0; i<dims[0]; i++) for(j=0; j<dims[1]; j++) { 
      ind_f = i + dims[0] * j; 
      f_fort[ind_f]  =  f_c[ind_c++]; 
    }
    break;
  
  case 3 : 
    for(i=0; i<dims[0]; i++) for(j=0; j<dims[1]; j++) for(k=0; k<dims[2]; k++) { 
      ind_f = i + dims[0] * ( j  +  dims[1] * k ); 
      f_fort[ind_f]  =  f_c[ind_c++]; 
    }
    break;
  

    /* Only ndims=1,2,3  are implemented so far: */
  default : 
    fprintf(stderr,"c_to_fort_array(): bad value for ndims = %d \n", ndims); 
    fflush(stderr);
    TRACE_END;
    return(1); 
    break;
  }

  TRACE_END;
  return(0); 
}

