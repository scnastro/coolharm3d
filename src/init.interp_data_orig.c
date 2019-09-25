
#include "decs.h"

#if( MAKE_HDF5 ) 
#include <hdf5.h>

#define PRINT_DEBUG_INFO (0)  


/* We will assume the source data is stored in the following file: */
#define FILENAME "rdump_start_other.h5"

static char filename[100]; 
static hid_t file_id, header_id, grid_id, ioparams_id, macros_id;
static int totalsize_src[NDIM], ng_src;;
static double startx_src[NDIM], begx_src[NDIM], dx_src[NDIM], inv_dx_src[NDIM], GridLength_src[NDIM];
static int metric_type_choice_src, metric_dynamic_type_choice_src, top_type_choice_src, coord_type_choice_src;


static void read_parameters(char *filename) ;
static void set_interp_weights( int ****src_indices, double ****weights, int n1e_loc, int n2e_loc, int n3e_loc, int n1, int n2, int n3 );

extern void   myH5_read_scalar2( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
extern void     myH5_read_basic( hid_t loc_id, char *name, hid_t type_id, void *value ) ;

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

   -- This routine reads in data generated in another coordinate system and 
      interpolates it into the given one set by the compile-time macro parameters;  

   -- The source data are assumed to be stored in a regular, merged restart checkpoint file named 
         "FILENAME" (see above macro), using the typical format given in dump_hdf.c and restart.c. 
   
   -- This routine is responsible for figuring out how to convert the source data to the 
       destination coordinate system used in this executable; 

***********************************************************************************/

void init(void)
{
  double ftmp;
  void init_data(void) ;

  TRACE_BEG;

  /* Do some checks on the pre-compiler definitions : */ 
  ftmp = (double) N1TOT;  ftmp *= N2TOT  ;  ftmp *= N3TOT  ;  ftmp *= NP; 
  if( ftmp > ( (double) UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Max. no. of elements too large !!! N = %g \n", ftmp);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  if( n_cells_glob > ( UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Num. of global points too large !!! N = %ld \n", 
	    n_cells_glob);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  init_base();  /* Set global grid parameters and static work arrays */

  init_data();  /* Set MHD grid functions */


  TRACE_END;
  return;
}

/***********************************************************************************/
/***********************************************************************************
  init_base(): 
  ---------
   -- initializes grid structure and constant grid functions; 

***********************************************************************************/
void init_base(void)
{
  int i,j,k ;
  double X[NDIM], x[NDIM];

  void calc_all_geom(void) ;
  void set_special_coord(void);

  TRACE_BEG;

  /* Need to open and read information from the source data file : */ 
  sprintf(filename,"%s",FILENAME); 
  read_parameters(filename);



  /* Set the parameters that define the grid and other constants : */ 
  //	gam = (5./3.) ;
  //	cour = 0.45;

  //	t = 0. ;                     /* Initial time */ 

  /* Coordinate dependent quantities : */ 
  set_special_coord();

  /**************************************************************************
          Length of each dimension : 
  **************************************************************************/
  GridLength[0] = 10. ;         /* Length of X0 dimension (evolution period) */ 

  /*************** SPHERICAL *******************************************/
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )

# if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  GridLength[1] = Rout-Rin  ;   /* Length of X1 dimension */ 
  GridLength[2] = th_length ;   /* Length of X2 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
  GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
  GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
  //	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  GridLength[1] = log((Rout-R0)/(Rin-R0)) ;              /* Length of X1 dimension */ 
  GridLength[2] = xi_diag2[n_diag2_lines]-xi_diag2[0]  ; /* Length of X2 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
  GridLength[1] = log((coord_params->upsilon_r-coord_params->f0_r)/(1.-coord_params->f0_r)) ;              /* Length of X1 dimension */ 
  GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
  //	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
  GridLength[1] = 1.  ;              /* Length of X1 dimension */
  GridLength[2] = 1.  ;                    /* Length of X2 dimension */

# endif
  GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
# if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
  GridLength[3] = 1. ;                 /* Length of X3 dimension */ 
# endif

  /*************** CARTESIAN  *******************************************/
#elif( TOP_TYPE_CHOICE == TOP_CARTESIAN )
# if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  GridLength[1] = 20. ;   /* Length of X1 dimension */ 
  GridLength[2] = 20. ;   /* Length of X1 dimension */ 
  GridLength[3] = 20. ;   /* Length of X1 dimension */ 
	
# else 
  --not-supported-2
# endif
	
#else
    --not-supported-3
#endif

    /**************************************************************************
          Grid discretization scales in each dimension : 
    **************************************************************************/
    SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];
#if(EQUATORIAL_RUN)
  dx[2] = 1.e-10;
#endif

  /**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
  **************************************************************************/
  startx[0] =  0.;                /* Set the Physical Minimum boundary  */

  /*************** SPHERICAL *******************************************/
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
# if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  startx[1] = Rin;                /* Set the Physical Minimum boundary  */
  startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
#  if(EQUATORIAL_RUN)
  startx[2] = th_length*0.5 - 0.5*dx[2] ;
#  else 
  startx[2] = th_beg        ; 
#  endif

# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
  startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
  //	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
  //	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
#  if(EQUATORIAL_RUN)
  startx[2] = 0.5 - 0.5*dx[2] ;
#  else 
  startx[2] = 0.            ; 
#  endif
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
  startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
  startx[1] = 0.;        /* Set the Physical Minimum boundary  */
  startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
  //	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */

# elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
  startx[1] = 1.e-10;        /* Set the Physical Minimum boundary  */
#  if(EQUATORIAL_RUN)
  startx[2] = 0.5-0.5*dx[2] ; 
#  else 
  startx[2] = 0.            ; 
#  endif 

  startx[3] =  0.;                /* Set the Physical Minimum boundary  */
# endif

  /*************** CARTESIAN  *******************************************/
#elif( TOP_TYPE_CHOICE == TOP_CARTESIAN )
# if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  SDLOOP1 startx[i] = -0.5*GridLength[i];
#  if(EQUATORIAL_RUN)
  startx[2] = -0.5*dx[2] ; 
#  endif 

# else 
  --not-supported-2
# endif
	
#else
    --not-supported-3
#endif

    SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

  dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

  /**************************************************************************
          Output frequencies : 
  **************************************************************************/
  DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
  DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
  DT_out[OUT_HISTORY] = 1. ;                	/* history file frequency */
  DT_out[OUT_SURFACE] = 3. ;                	/* surface file frequency */
  DT_out[OUT_HDF5]    = 0.5 ;                     /* hdf frequency */
  DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
  DT_out[OUT_STAT]    = 4.*DT_out[OUT_HDF5];      /* statistics dumps frequency */
  DT_out[OUT_STAT2]   = DT_out[OUT_STAT];         /* statistics dumps frequency */
  DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
  DT_out[OUT_MIN_DT]  = DT_out[OUT_HISTORY];
  n_restart = 1000;                /* number of time steps between restart dumps */
	

  /**************************************************************************
	  Calculate the static metric grid functions 
  **************************************************************************/
  dx[0] = 1.e-4;
  calc_all_geom() ;  
  dt_global_min = find_min_dt();

  //	dx[0] = cour * dt_global_min; /* Discretization size in X0 direction, time step*/

  /**************************************************************************
	  Print out basic grid parameters: 
  **************************************************************************/
  fprintf(stdout,"Length: %10.4e %10.4e %10.4e %10.4e\n", 
	  GridLength[0], GridLength[1], GridLength[2], GridLength[3]);
  fprintf(stdout,"dx    : %10.4e %10.4e %10.4e %10.4e\n", 
	  dx[0], dx[1], dx[2], dx[3]);
  fprintf(stdout,"startx: %10.4e %10.4e %10.4e %10.4e\n", 
	  startx[0], startx[1], startx[2], startx[3]);

  fprintf(stdout,"rhomin,minlimt = %28.18e %28.18e \n", RHOMIN, RHOMINLIMIT);
  fprintf(stdout,"uumin ,minlimt = %28.18e %28.18e \n", UUMIN,  UUMINLIMIT);
	
  fflush(stdout);


  TRACE_END;
  return;

}

/***********************************************************************************/
/***********************************************************************************
  read_parameters(): 
  ---------
   -- opens the source hdf5 file and reads in all the necessary parameters;
   -- follows similar conventions as used in dump_hdf.c and restart.c ;

***********************************************************************************/
static void read_parameters(char *filename) 
{
  TRACE_BEG;

  int i; 
  
  /******************************************************************************
    Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Open existing file : wrapper handles error checking */ 
  file_id = myH5_Fopen(filename);

  /* Create a subgroup to which we attach attributes or header information */
  header_id = H5Gopen(file_id, "Header");
  if( header_id  < 0 ) { 
    fprintf(stderr,"read_parameters(): Cannot open Header group in file %s \n", filename);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create a header subgroup that contains the grid parameters: */
  grid_id = H5Gopen(header_id, "Grid");
  if( grid_id  < 0 ) { 
    fprintf(stderr,"read_parameters(): Cannot open Grid group in file %s \n", filename);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create a header subgroup that contains the I/O parameters: */
  ioparams_id = H5Gopen(header_id, "IO");
  if( ioparams_id  < 0 ) { 
    fprintf(stderr,"read_parameters(): Cannot open IO group in file %s \n", filename);
    fflush(stderr);     fail(FAIL_HDF,0);
  }

  /* Create a header subgroup that contains the macros: */
  macros_id = H5Gopen(header_id, "Macros");
  if( macros_id  < 0 ) { 
    fprintf(stderr,"read_parameters(): Cannot open Macros group in file %s \n", filename);
    fflush(stderr);     fail(FAIL_HDF,0);
  }


  myH5_read_scalar2(grid_id, "totalsize0"      , H5T_NATIVE_INT   ,  &(totalsize_src[0])   );
  myH5_read_scalar2(grid_id, "totalsize1"      , H5T_NATIVE_INT   ,  &(totalsize_src[1])   );
  myH5_read_scalar2(grid_id, "totalsize2"      , H5T_NATIVE_INT   ,  &(totalsize_src[2])   );
  myH5_read_scalar2(grid_id, "totalsize3"      , H5T_NATIVE_INT   ,  &(totalsize_src[3])   );
  myH5_read_scalar2(grid_id, "NG"              , H5T_NATIVE_INT   ,  &ng_src               );
  myH5_read_scalar2(grid_id, "startx0"         , H5T_NATIVE_DOUBLE,  &(startx_src[0])      );
  myH5_read_scalar2(grid_id, "startx1"         , H5T_NATIVE_DOUBLE,  &(startx_src[1])      );
  myH5_read_scalar2(grid_id, "startx2"         , H5T_NATIVE_DOUBLE,  &(startx_src[2])      );
  myH5_read_scalar2(grid_id, "startx3"         , H5T_NATIVE_DOUBLE,  &(startx_src[3])      );
  myH5_read_scalar2(grid_id, "dx0"             , H5T_NATIVE_DOUBLE,  &(dx_src[0])          );
  myH5_read_scalar2(grid_id, "dx1"             , H5T_NATIVE_DOUBLE,  &(dx_src[1])          );
  myH5_read_scalar2(grid_id, "dx2"             , H5T_NATIVE_DOUBLE,  &(dx_src[2])          );
  myH5_read_scalar2(grid_id, "dx3"             , H5T_NATIVE_DOUBLE,  &(dx_src[3])          );
  myH5_read_scalar2(grid_id, "gridlength0"     , H5T_NATIVE_DOUBLE,  &(GridLength_src[0])  );
  myH5_read_scalar2(grid_id, "gridlength1"     , H5T_NATIVE_DOUBLE,  &(GridLength_src[1])  );
  myH5_read_scalar2(grid_id, "gridlength2"     , H5T_NATIVE_DOUBLE,  &(GridLength_src[2])  );
  myH5_read_scalar2(grid_id, "gridlength3"     , H5T_NATIVE_DOUBLE,  &(GridLength_src[3])  );
  myH5_read_scalar2(grid_id, "t"               , H5T_NATIVE_DOUBLE,  &t                );
  myH5_read_scalar2(grid_id, "cour"            , H5T_NATIVE_DOUBLE,  &cour    	        );
  myH5_read_scalar2(grid_id, "gam"             , H5T_NATIVE_DOUBLE,  &gam     	        );
  myH5_read_scalar2(grid_id, "a"               , H5T_NATIVE_DOUBLE,  &a	        );
  myH5_read_scalar2(grid_id, "initial_bbh_separation"  , H5T_NATIVE_DOUBLE   ,  &initial_bbh_separation );
  myH5_read_scalar2(grid_id, "m_bh1"           , H5T_NATIVE_DOUBLE   ,  &m_bh1);
  myH5_read_scalar2(grid_id, "m_bh2"           , H5T_NATIVE_DOUBLE   ,  &m_bh2);
  myH5_read_scalar2(grid_id, "m_bh_tot"        , H5T_NATIVE_DOUBLE   ,  &m_bh_tot);
  myH5_read_scalar2(grid_id, "M"               , H5T_NATIVE_DOUBLE   ,  &M);
  myH5_read_scalar2(macros_id,"METRIC_TYPE_CHOICE"         ,H5T_NATIVE_INT,&metric_type_choice_src        );
  myH5_read_scalar2(macros_id,"METRIC_DYNAMIC_TYPE_CHOICE" ,H5T_NATIVE_INT,&metric_dynamic_type_choice_src);
  myH5_read_scalar2(macros_id,"TOP_TYPE_CHOICE" 	   ,H5T_NATIVE_INT,&top_type_choice_src           );
  myH5_read_scalar2(macros_id,"COORD_TYPE_CHOICE"	   ,H5T_NATIVE_INT,&coord_type_choice_src         );

  DLOOP1 { inv_dx_src[i] = 1./dx_src[i];                             }
  DLOOP1 { begx_src[i]   = startx_src[i] + (ng_src + 0.5)*dx_src[i]; }

  fprintf(stdout,"inv_dx_src = %26.16e %26.16e %26.16e %26.16e\n", inv_dx_src[0],inv_dx_src[1],inv_dx_src[2],inv_dx_src[3]); 
  fprintf(stdout,"begx_src   = %26.16e %26.16e %26.16e %26.16e\n", begx_src[0],begx_src[1],begx_src[2],begx_src[3]); 
  fprintf(stdout,"startx_src = %26.16e %26.16e %26.16e %26.16e\n", startx_src[0],startx_src[1],startx_src[2],startx_src[3]); 
  fflush(stdout);
  

  TRACE_END;
  return;
}

/***********************************************************************************/
/***********************************************************************************
  read_and_interp_all_funcs(): 
  ---------
   -- responsible for reading all grid function data and interpolating it 
      to the local grid and coordinate system; 

***********************************************************************************/
void read_and_interp_all_funcs(void)
{
  TRACE_BEG;
  
#define SOURCE_INTERP(a,b,c)  ( (a[(c[0])])*(b[0]) + (a[(c[1])])*(b[1]) + (a[(c[2])])*(b[2]) + (a[(c[3])])*(b[3]) + (a[(c[4])])*(b[4]) + (a[(c[5])])*(b[5]) + (a[(c[6])])*(b[6]) + (a[(c[7])])*(b[7])  )

  /**************************************************************************************************
     Interpolate rho, uu, v[1-3], B[1-3] to local cell centers:
        -- 1) read in x[1-3] for all use; 
           2) loop over each grid function and interpolate one at a time; 
           3) after interpolating from Cartesian to numerical coordinates, 
               transform components from Cartesian to numerical coordinate system;
  **************************************************************************************************/
  int i,j,k,l;
  int nsrc1  = totalsize_src[1];
  int nsrc2  = totalsize_src[2];
  int nsrc3  = totalsize_src[3];
  int n_data = nsrc1 * nsrc2 * nsrc3;
  double *data_in; 
  double ****interp_weights, *wloc;
  int    ****interp_ind, *iloc;
  
  ALLOC_ARRAY( data_in , n_data );
  ALLOC_4D_ARRAY( interp_weights , N1TOT , N2TOT, N3TOT, 8 ); 
  ALLOC_4D_ARRAY( interp_ind     , N1TOT , N2TOT, N3TOT, 8 ); 


  set_interp_weights(interp_ind, interp_weights,N1E,N2E,N3E,nsrc1,nsrc2,nsrc3);

  /* Now cycle through all the primitives, interpolating one at a time: */ 
  myH5_read_basic( file_id, "rho", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][RHO] = SOURCE_INTERP(data_in,wloc,iloc);  }

  myH5_read_basic( file_id, "uu", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][UU] = SOURCE_INTERP(data_in,wloc,iloc);  }

  myH5_read_basic( file_id, "v1", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][U1] = SOURCE_INTERP(data_in,wloc,iloc);  }

  myH5_read_basic( file_id, "v2", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][U2] = SOURCE_INTERP(data_in,wloc,iloc);  }

  myH5_read_basic( file_id, "v3", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][U3] = SOURCE_INTERP(data_in,wloc,iloc);  }

#if( !HYDRO_ONLY ) 
  myH5_read_basic( file_id, "B1", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][B1] = SOURCE_INTERP(data_in,wloc,iloc);  }

  myH5_read_basic( file_id, "B2", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][B2] = SOURCE_INTERP(data_in,wloc,iloc);  }

  myH5_read_basic( file_id, "B3", H5T_NATIVE_DOUBLE, data_in ) ;
  LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; p[i][j][k][B3] = SOURCE_INTERP(data_in,wloc,iloc);  }
#endif

  DEALLOC_ARRAY( data_in , n_data );

//no-emf-for-now   
//no-emf-for-now     /**************************************************************************************************
//no-emf-for-now        Interpolate EMFs to cell corners: 
//no-emf-for-now           -- essentially do the same but for the cell corners;  
//no-emf-for-now           -- read in special  x_emf[] coordinates; 
//no-emf-for-now           -- 1) read in x[1-3] for all use; 
//no-emf-for-now              2) loop over each grid function and interpolate one at a time; 
//no-emf-for-now              3) after interpolating from Cartesian to numerical coordinates, 
//no-emf-for-now                  transform components from Cartesian to numerical coordinate system;
//no-emf-for-now     **************************************************************************************************/
//no-emf-for-now     int nemf1  = nsrc1+1;
//no-emf-for-now     int nemf2  = nsrc2+1;
//no-emf-for-now     int nemf3  = nsrc3+1;
//no-emf-for-now     int n_data2 = nemf1*nemf2*nemf3;
//no-emf-for-now     int n1e_new,n2e_new,n3e_new;
//no-emf-for-now   
//no-emf-for-now   #define NEW_LOOP for(i=N1S;i<=n1e_new ;i++) for(j=N2S;j<=n2e_new;j++) for(k=N3S;k<=n3e_new;k++)
//no-emf-for-now   
//no-emf-for-now     ALLOC_ARRAY( data_in , n_data2 );
//no-emf-for-now   
//no-emf-for-now     set_interp_weights(interp_ind, interp_weights,n1e_new,n2e_new,n3e_new,nemf1,nemf2,nemf3);
//no-emf-for-now   
//no-emf-for-now     myH5_read_basic( file_id, "emf1", H5T_NATIVE_DOUBLE, data_in ) ;
//no-emf-for-now     NEW_LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; emf[0][i][j][k] = SOURCE_INTERP(data_in,wloc,iloc); }
//no-emf-for-now   
//no-emf-for-now     myH5_read_basic( file_id, "emf2", H5T_NATIVE_DOUBLE, data_in ) ;
//no-emf-for-now     NEW_LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; emf[1][i][j][k] = SOURCE_INTERP(data_in,wloc,iloc); }
//no-emf-for-now   
//no-emf-for-now     myH5_read_basic( file_id, "emf3", H5T_NATIVE_DOUBLE, data_in ) ;
//no-emf-for-now     NEW_LOOP { iloc=interp_ind[i][j][k];  wloc=interp_weights[i][j][k]; emf[2][i][j][k] = SOURCE_INTERP(data_in,wloc,iloc); }
//no-emf-for-now   
//no-emf-for-now     DEALLOC_ARRAY( data_in , n_data2 );

  /**************************************************************************************************
    Now transform everything, knowing that we know p[] has not been set over all of NEW_LOOP's elements: 
  **************************************************************************************************/
  double u1,u2,u3,b1,b2,b3,gdet;
  double *primloc;
  struct of_coord *coords; 
  double dxp_dxc[NDIM][NDIM];

  if( (top_type_choice_src != TOP_CARTESIAN) && (TOP_TYPE_CHOICE != TOP_SPHERICAL) ) { 
    fprintf(stdout,"%s(): problem with topologies here:  %d %d \n", __func__, top_type_choice_src,TOP_TYPE_CHOICE); 
    fflush(stdout); fail(FAIL_BASIC,0); 
  }

  //-no-emf-for-now  NEW_LOOP { 
  LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);    
    invert_matrix2(coords->dxc_dxp,dxp_dxc,&gdet);
    primloc = p[i][j][k];
    u1 = primloc[U1]; 
    u2 = primloc[U2]; 
    u3 = primloc[U3]; 

    primloc[U1] = 
      dxp_dxc[  1][1] * u1 + 
      dxp_dxc[  1][2] * u2 + 
      dxp_dxc[  1][3] * u3 ;
    
    primloc[U2] = 
      dxp_dxc[  2][1] * u1 + 
      dxp_dxc[  2][2] * u2 + 
      dxp_dxc[  2][3] * u3 ;
    
    primloc[U3] = 
      dxp_dxc[  3][1] * u1 + 
      dxp_dxc[  3][2] * u2 + 
      dxp_dxc[  3][3] * u3 ;
    
#if( !HYDRO_ONLY ) 
    b1 = primloc[B1]; 
    b2 = primloc[B2]; 
    b3 = primloc[B3]; 
    primloc[B1] = 
      dxp_dxc[  1][1] * b1 + 
      dxp_dxc[  1][2] * b2 + 
      dxp_dxc[  1][3] * b3 ;
    
    primloc[B2] = 
      dxp_dxc[  2][1] * b1 + 
      dxp_dxc[  2][2] * b2 + 
      dxp_dxc[  2][3] * b3 ;
    
    primloc[B3] = 
      dxp_dxc[  3][1] * b1 + 
      dxp_dxc[  3][2] * b2 + 
      dxp_dxc[  3][3] * b3 ;
#endif

#if(PRINT_DEBUG_INFO)
    if( (i==(N1S+3)) && (j==(N2S)) && (k==(N3S+3)) ) { 
      fprintf(stdout,"weights = "); 
      for(l=0;l<8;l++) { fprintf(stdout," %16.10e ",interp_weights[i][j][k][l]);  }
      fprintf(stdout,"\n"); 
      fprintf(stdout,"indices = "); 
      for(l=0;l<8;l++) { fprintf(stdout," %16d ",interp_ind[i][j][k][l]);  }
      fprintf(stdout,"\n"); 
      fprintf(stdout,"prims = "); 
      for(l=0;l<8;l++) { fprintf(stdout," %16.10e ",primloc[l]);  }
      fprintf(stdout,"\n"); 
      fflush(stdout);
    }
#endif

    /* note that EMF's do NOT transform like vectors:  

        emf[0-2] = emf_k = -(1/sqrt{-g}) \epsilon_{tijk} b^i u^j 

        -- so you have to multiply by source \sqrt{-g}, transform to local coords, then divide by local \sqrt{-g} 

     */ 
//-no-emf-for-now     u1 = emf[0][i][j][k];
//-no-emf-for-now     u2 = emf[1][i][j][k];
//-no-emf-for-now     u3 = emf[2][i][j][k];
//-no-emf-for-now 
//-no-emf-for-now     emf[0][i][j][k] = 
//-no-emf-for-now       dxp_dxc[  1][1] * u1 + 
//-no-emf-for-now       dxp_dxc[  1][2] * u2 + 
//-no-emf-for-now       dxp_dxc[  1][3] * u3 ;
//-no-emf-for-now     
//-no-emf-for-now     emf[1][i][j][k] = 
//-no-emf-for-now       dxp_dxc[  2][1] * u1 + 
//-no-emf-for-now       dxp_dxc[  2][2] * u2 + 
//-no-emf-for-now       dxp_dxc[  2][3] * u3 ;
//-no-emf-for-now     
//-no-emf-for-now     emf[2][i][j][k] = 
//-no-emf-for-now       dxp_dxc[  3][1] * u1 + 
//-no-emf-for-now       dxp_dxc[  3][2] * u2 + 
//-no-emf-for-now       dxp_dxc[  3][3] * u3 ;

  }

  DEALLOC_4D_ARRAY( interp_weights , N1TOT , N2TOT, N3TOT, 8 ); 
  DEALLOC_4D_ARRAY( interp_ind     , N1TOT , N2TOT, N3TOT, 8 ); 

#undef SOURCE_INTERP

  TRACE_END;
  return;
}

/****************************************************************************************

 set_interp_weights(): 
 ---------------------
   -- determines simulation cell indices used for the interpolation of source data
         at each point of the destination system;
   -- records the weights used for the interpolation also;
   -- routine borrowed and modified from bothros;
   -- assumes "weights" is pre-allocated; 

     Arguments: 
           X[]          = location of interpolation;
   xp1,xp2,xp3          = array of parent grid function's coordinates;
          sdata         = grid function to be interpolated, living on Xgrid[];
          dX[]          = grid discretizations in the various directions;
          
****************************************************************************************/
static void set_interp_weights( int ****src_indices, double ****weights, 
				int n1e_loc, int n2e_loc, int n3e_loc, 
				int n1, int n2, int n3 )
{
  int i,j,k,l;
  int kx1min,kx2min,kx3min,k1,k2,k3,k4,k5,k6,k7,k8;
  double dx1,dx2,dx3,ddx1,ddx2,ddx3,delx_1, delx_2, delx_3;
  double *xc;
  double xp[NDIM];
  int *ind_loc;
  double *weights_loc;
  struct of_coord *coords; 

  int n23 = n2*n3;

  for(i=N1S;i<=n1e_loc ;i++) for(j=N2S;j<=n2e_loc;j++) for(k=N3S;k<=n3e_loc;k++) { 
	
	delx_1 = (xc[1] - begx_src[1]) * inv_dx_src[1] ;
	delx_2 = (xc[2] - begx_src[2]) * inv_dx_src[2] ;
	delx_3 = (xc[3] - begx_src[3]) * inv_dx_src[3] ;

	/*  Find the indices of the nearest sim_data[COL_X*] coord's smaller than or equal to X[]: */
	/*  k1 -> (i,j),  k2 -> (i+1,j), k3 -> (i+1,j+1), k4 -> (i,j+1) */
	kx1min = (int) delx_1 ; 
	kx2min = (int) delx_2 ; 
	kx3min = (int) delx_3 ; 

	 /* Check for out bounds indices, and set the closest allowable value: */
	if( kx1min >= (n1-1) ) {  kx1min = n1 - 2 ; }   
	if( kx2min >= (n2-1) ) {  kx2min = n2 - 2 ; } 
	if( kx3min >= (n3-1) ) {  kx3min = n3 - 2 ; } 
	if( kx1min  <  0 )     {  kx1min = 0      ; }   
	if( kx2min  <  0 )     {  kx2min = 0      ; } 
	if( kx3min  <  0 )     {  kx3min = 0      ; } 

	/* Set the indices of the other cells used in the interpolation: */
	k1 = kx3min + n3*(kx2min + n2*kx1min);  // corner
	k2 = k1 + n23;                          // corner + x1
	k3 = k2 + n3 ;                          // corner + x1 + x2   
	k4 = k1 + n3 ;                          // corner      + x2   
	k5 = k1 + 1  ;                          // corner           + x3 
	k6 = k2 + 1  ;                          // corner + x1      + x3 
	k7 = k3 + 1  ;                          // corner + x1 + x2 + x3  
	k8 = k4 + 1  ;                          // corner      + x2 + x3  
	
	ind_loc     = src_indices[i][j][k]; 
	weights_loc =     weights[i][j][k]; 

	ind_loc[0] = k1;
	ind_loc[1] = k2;
	ind_loc[2] = k3;
	ind_loc[3] = k4;
	ind_loc[4] = k5;
	ind_loc[5] = k6;
	ind_loc[6] = k7;
	ind_loc[7] = k8;
  
	dx1 = delx_1 - kx1min; 
	dx2 = delx_2 - kx2min; 
	dx3 = delx_3 - kx3min; 

	ddx1 = 1.-dx1;
	ddx2 = 1.-dx2;
	ddx3 = 1.-dx3;

	weights_loc[0] = ddx1*ddx2*ddx3;
	weights_loc[1] =  dx1*ddx2*ddx3;
	weights_loc[2] =  dx1* dx2*ddx3;
	weights_loc[3] = ddx1* dx2*ddx3;
	weights_loc[4] = ddx1*ddx2* dx3;
	weights_loc[5] =  dx1*ddx2* dx3;
	weights_loc[6] =  dx1* dx2* dx3;
	weights_loc[7] = ddx1* dx2* dx3;

#if(PRINT_DEBUG_INFO)
	if( (i==(N1S+3)) && (j==(N2S)) && (k==(N3S+3)) ) { 
	  fprintf(stdout,"set_interp():  weights = "); 
	  for(l=0;l<8;l++) { fprintf(stdout," %16.10e ",weights_loc[l]);  }
	  fprintf(stdout,"\n"); 
	  fprintf(stdout,"set_interp():  indices = "); 
	  for(l=0;l<8;l++) { fprintf(stdout," %16d ",ind_loc[l]);  }
	  fprintf(stdout,"\n"); 
	  fprintf(stdout,"set_interp():  dx1-3 = %16.10e %16.10e %16.10e \n",dx1,dx2,dx3); 
	  fprintf(stdout,"set_interp(): ddx1-3 = %16.10e %16.10e %16.10e \n",ddx1,ddx2,ddx3); 
	  fprintf(stdout,"set_interp(): xc     = %16.10e %16.10e %16.10e \n",xc[1],xc[2],xc[3]); 
	  fprintf(stdout,"set_interp(): kxmin  = %16d %16d %16d \n",kx1min,kx2min,kx3min);
	  fflush(stdout);
	}
#endif

      }

  return;
}

/***********************************************************************************/
/***********************************************************************************
  init_data(): 
  ---------
   -- calculates the initial distribution of the MHD fields;
   -- driver routine for various initial data prescriptions;

***********************************************************************************/
void init_data(void) 
{

  TRACE_BEG;

  int i,j,k;

  read_and_interp_all_funcs();

  /*******************************************************************************
    Set the magnetic field :  
  *******************************************************************************/
  /* Correct bad points and setup boundary values since we will require them for B^i */
  fixup(p) ;
  bounds(p,0) ;

#if( HYDRO_ONLY ) 
  LOOP { 
    p[i][j][k][B1] = p[i][j][k][B2] = p[i][j][k][B3] =  0.; 
  }
#endif 

  fixup(p) ;    
  bounds(p,0) ;

  TRACE_END;
}


/***********************************************************************/
/***********************************************************************
  set_special_coord(): 
  --------------------
   -- set coordinate specific (global) quantities;
***********************************************************************/
void set_special_coord(void)
{
  int i;
  
  /* Spacetime parameters : */
  //  M = m_bh_tot = 1.;

  /* BBH parameters: */
#if( BBH_SPACETIME ) 
  //  m_bh1 = 0.5*M;                            /* Mass of left-most BH */
  //  m_bh2 = 0.5*M;                            /* Mass of right-most BH */
  //  initial_bbh_separation = 20.*m_bh_tot;  /* Initial separtion of BHs in units of total mass */

  r_horizon1 = m_bh1;             /* Radius of the horizon about BH1 in Harmonic coordinates */
  r_horizon2 = m_bh2;             /* Radius of the horizon about BH2 in Harmonic coordinates */

  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;

  //  a = 0.0;          /* Spin of the black hole in units of M */ 
  asq = a*a;
  r_isco    = 0.;
  r_horizon = 0.;

//  Rout     = 13.*initial_bbh_separation;        /* Radial extent of the grid           */
  Rout = 300. ;
  n_within_horizon = 0;            /* number of cells within horizon */
  //  Rin       = Rin_calc();

  Rin = 0.75 * (initial_bbh_separation) ; /*  put inner boundary at 2. times the binary orbit's radius */
  //  R0       = 0.9*Rin;
  R0 = 0.;

  /* Time at which to start shrinking the binary,  set to 0. if you want to shrink from the beginning and 
     set to an impossibly large number if you never want to shrink : */
  t_shrink_bbh = 0.; 
  phi_0_bbh    = 0.;

#else 
  /* Single BH Parameters: */
  //  a = 0.; 
  asq = a*a;
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  Rout     = 8.;         /* Radial extent of the grid           */
  n_within_horizon = 0;   /* number of cells within horizon */
  Rin      = 1.;
  
  r_horizon1 = r_horizon2 = r_horizon;
  m_bh1 = m_bh2 = M;
  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;
  initial_bbh_separation = 0.;

#endif

  h_slope   = 0.13;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

#if( EQUATORIAL_RUN ) 
  th_cutout = 0.;
#else 
  th_cutout = 2.*SMALL;     
  /* 0.02655 gives a cutout of ~0.045Pi for h_slope = 0.35 */
  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */ 
  //  th_cutout = 0.0* M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
#endif                 

  /* New way : */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */

  /* Old way: */
//  th_beg    = 0.;               /* beg/end set the limits of theta  and */
//  th_end    = M_PI;             /*  can override th_cutout              */
//  th_length = M_PI;             /* Angular size in theta of the grid    */
  
  /* COORD_DIAGONAL3 specific parameters (note h_slope is overloaded for dtheta_min/(pi/2) */
  diag3_exponent = 9;                            /* note that this should be an odd integer greater than 1  */
  diag3_factor   = 1. - 2*th_cutout/M_PI - h_slope;  /* Temporary variable */
  if( (diag3_exponent < 3) || (diag3_exponent % 2 == 0 )  ) { fprintf(stderr,"set_special_coord(): bad value for  diag3_exponent  = %d \n",diag3_exponent); fflush(stderr); exit(4112); }


  n_phi_avg = 100; 
  inv_n_phi_avg = 1./n_phi_avg; 
  d_phi_avg = 2.*M_PI * inv_n_phi_avg; 

#if( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  this-initial-data-is-not-setup-for-COORD_DIAGONAL2
#endif 

  SET_GEOM_COORD_TIME_FUNCS(t,1);

  if( myid == printer_pid ) { 
    fprintf(stdout,"a                = %28.18e \n", a);
    fprintf(stdout,"R0               = %28.18e \n", R0               );
    fprintf(stdout,"Rin              = %28.18e \n", Rin              );
    fprintf(stdout,"Rout             = %28.18e \n", Rout             );
    fprintf(stdout,"n_within_horizon = %6d     \n", n_within_horizon );
    fprintf(stdout,"h_slope          = %28.18e \n", h_slope          );
    fprintf(stdout,"X1_slope         = %28.18e \n", X1_slope         );
    fprintf(stdout,"X1_0             = %28.18e \n", X1_0             );
    fprintf(stdout,"r_isco           = %28.18e \n", r_isco           );
    fprintf(stdout,"r_horizon        = %28.18e \n", r_horizon        );
    fflush(stdout);
  }

  return;
}

#undef PRINT_DEBUG_INFO

#else

void init(void)
{
  fprintf(stdout,"init(): %s  relies on HDF5, so please include it!! \n",__FILE__); 
  fflush(stdout); 
  fail(FAIL_BASIC,0); 
}

#endif
