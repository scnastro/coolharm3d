/******************************************************************************************************************************/
/******************************************************************************************************************************
  interp_data.c:
 --------------

     -- contains all the routines for the interpolator. 

     -- there are two main driver routines:

          read_and_interp_all_funcs()          :   serial,   non-MPI version  for small problems
          read_and_interp_all_funcs_parallel() :   parallel,     MPI version  for   big problems


     -- both basically work the same way, but the parallel version breaks the process up into 
         small pieces so that each piece can fit into memory.  The number of the "batches" 
         is controlled by the compile-time variable, "nbatches", 
         set in read_and_interp_all_funcs_parallel(void).  

     -- The procedure works in the following way: 

           1) First you generate "destination" coordinates, which are arrays of x1, x2, x3 
              (or physical coorinates) that we will interpolate data ONTO (i.e. the DESTINATION). 
              These coordinates are assumed to be written to a HDF5 file just like harm3d 
              normally does, so it is easiest to just setup a code version with your new 
              coordinates (you would have to do this anyway), use some initialization routine 
              and run it so that you can generate a dump file with the x[1-3] coordinate arrays. 
              By default, the code looks a file "./rdump_start_other.h5"  for the destination coordinates, 
              so you can name your file this if you want.   Or, you can name it something else and set
              the environment variable DEST_INTERP_FILENAME to the path/name of the file.  Note that you 
              have to set the variable and run the code on the same line if doing a MPI run; see 
              the "testinterp" make option in "makefile."

           2) This program figures out the type of destination coordinates they are by reading the "TOP_TYPE_CHOICE"
              variable in the hdf5 file.   It then sets up the transformation matrices so that it can 
              transform from the source data (in memory) to the destination data, and transform 
              the velocity & magnetic field components into the new basis. 

           3) The "interpolation" run is initialized by an "init" routine or with a checkpoint file.  The code used 
              for the interpolation run should be that of the original run, from which you are interpolating data. 
              It should be setup to use the same coorindates and same metric used in that original run, 
              Use the same code and run it the same ways as you did before (with same command-line arguments). 
              a) If using a checkpoint file, it should be from this original run using the same code as you are using 
                 for the interpolation.  It should be data at the time you want the interpolation to occur.  If you 
                 do not have a regular checkpoint file (e.g., rdump_0.h5 or rdump_1.h5) then you can convert a 
                 regular 3-d dump file into one (see README file).  
              b) If you are using an "init" routine to start with, make sure you are starting at the right time. 
                 For instance, if the phase of the binary orbit is important, the phase is set by the initial time, 
                 t=0 usually implies phase=0. 

            4) In decs.h, set   
                     #define INTERPOLATE_DATA (1) 

            5) Recompile and run with as you did before to generate the source data.  You may kill the job once it 
               starts timestepping as by then it has already interpolated the data.  The interpolated data is stored 
               in the HDF5 file containing the destination coordinates.  

            6) Once you have destination data, using the code setup for the new run, in the destination coordinates, set in decs.h:
                     #define INTERPOLATE_DATA (0) 
                     #define RUN_FROM_INTERP_DATA (1)

               and recompile.  Run with the interpolated data file as your "rdump_start.h5" file.  You're done!  



     -- I have placed "USER:"  comments in this file near places where you may need to change 
         things. 

 ******************************************************************************************************************************/
/******************************************************************************************************************************/
#include "decs.h"
#include "metric.h"

#if( MAKE_HDF5 ) 
#include <hdf5.h>
#include <hdf5_hl.h>

#define PRINT_DEBUG_INFO (0)  
#define PRINT_NEWTON_RAPHSON_INFO (0)  

#define PERFORM_XP_OF_X_TEST (0)  

#define MISSING_VALUE_RHO (0.)

#define USE_EXTRAPOLATION (0)

/* Default name of the file that contains the destination coordinates, onto which we interpolate data: */
#define FILENAME       "rdump_start_other.h5"

/* The "TEMP" file is used to save data temporarily when doing a parallel interpolation: */
#define TEMP_FILENAME  "interp-tmp-file.h5"

/* Names of the grid functions we use in the "TEMP" file : */
#define PRIM_DEST_TMP_NAME  "prim_dest_tmp"
#define MASK_TMP_NAME "mask_tmp"

/* Rotate src data so that BHs are rotated clockwise by angle to the x-axis
   -- Intended so that BHs are placed back on the x-axis for SUPERIMPOSE_MINI_DISKS
   -- ONLY WORKS WHEN XP3 = X3 IN SRC COORDINATES */
#define PHASE_SHIFT (1)
int first_phase_calc = 1;

#define FTYPE double     /* your choice of floating-point data type */
#define MAX_NEWT_ITER 100     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-16    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER 2
#define GNR_SMALL (1.e-18)


static char      *filename;
static char *temp_filename;
static hid_t file_id, tmpfile_id;
static int top_type_choice_dest, coord_type_choice_dest;

static int general_newton_raphson_1d( double *x_phys, double *xp  );
static int general_newton_raphson_3d_orig( double *x_phys, double *xp  );
static int general_newton_raphson_3d_monotonic( double *x_phys, double *xp, short int *minimize, short int enforce_monotonicity);
static int general_newton_raphson_3d_select( double *x_phys, double *xp, short int *minimize);
static int general_newton_raphson_3d_general( double *x_phys, double *xp);

static void interp_prims( double **coord_list, double **prim_dest, usint *exists, int ndest);
static double **get_list_of_pts( char *filename, int *npts );
static double **get_list_of_pts2( char *filename, int *npts, int *npts_tot, int ibatch, int nbatches );
static void transform_dest_to_src_coords( double **coord_list, int npts ) ;
static void dump_interp_prims( char *filename, double **prim_dest, usint *exists, int npts, int npts_tot, int ibatch, int nbatches );
static void fix_interp_prims_dump( char *temp_filename, char *filename, int npts_tot );

extern void   myH5_read_scalar2( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
extern void     myH5_read_basic( hid_t loc_id, char *name, hid_t type_id, void *value ) ;
extern hsize_t *myH5_get_dims( hid_t loc_id, char *name,  int *ndims );
extern double  ranc(int iseed); 
extern int LU_decompose( double A[][NDIM], int permute[] );
extern void LU_substitution( double A[][NDIM], double B[], int permute[] );
extern hid_t    myH5_Fopen(char *name);

int ndims;
hsize_t *dims;

double theta_global;

int id_head;

void func_diagonal( double *x2, double *xp, double *dx, double *resid, double **jac, double *f, double *df, double *misc, int n, int n_misc );
void phys_val_diagonal(double *x2); 

#define phys_val  phys_val_diagonal 
#define func_resid  func_diagonal 




/***********************************************************************************/
/***********************************************************************************
  identity_phys_transform():
  -------------------------
   -- just copies the coordinate array; 
***********************************************************************************/
void  identity_phys_transform(double *x, double *xsrc )
{
  int i; 

  DLOOP1 { x[i] = xsrc[i]; }

  return;
}


/***********************************************************************************/
/***********************************************************************************
  identity_dxdest_dxsrc():
  -------------------------
   -- identify matrix for the transformation;
***********************************************************************************/
void  identity_dxdest_dxsrc(double *x, double *xphys, double *xsrc, double dxdest_dxsrc[][NDIM] )
{
  int i,j; 

  DLOOP1 { x[i] = xphys[i]; }
  dx_dxp_calc( xphys, xsrc, dxdest_dxsrc ); 

  return;
}


/***********************************************************************************/
/***********************************************************************************
  xp_of_x(): 
  ---------
   -- returns with xp (primed, numerical coordinates) from x (physical coordinates) 
       for any coordinate system; 
   -- developed from version in bothros; 
   -- returns with non-zero value if there was a problem finding the solution; 
   -- some coordinate systems use "periodic" functions that are not invertible, 
      meaning that there is not a unique solution.  The degenrate x() values lie 
      in the ghost zones.  So, we only assume that xp() solutions lie in the physical 
      region and not in the ghost zones.  
   -- this routine should uses the entering value of xp[] as a guess in the interative procedure; 
   -- this routine is meant to be used with a wrapper routine that supplies it with 
       good guesses (possibly using a lookup table); 
***********************************************************************************/
int  xp_of_x(double *xp, double *x )
{
  int i,j,k;
  double r,th,phi;
  double xp_orig[NDIM]; 
  double *misc_vars;
  int retval;
  struct of_coord *coords; 

  static int local_first_time = 1; 

#if( COORD_TYPE_CHOICE == COORD_DIAGONAL2         || \
     COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD || \
     COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN       )
  fprintf(stdout,"%s():COORD_TYPE_CHOICE == %d is not supported !! \n",__func__,COORD_TYPE_CHOICE); 
  fflush(stdout);  fail(FAIL_BASIC,0);
#endif

  retval = 0 ; 

  DLOOP1 { xp_orig[i] = xp[i]; } 

  /*******************************************************************************************************
     Time coordinate is always assumed to be the same, for all coordinate systems so far :
  *******************************************************************************************************/
  xp[0] = x[0]; 


  /*******************************************************************************************************/

#if( COORD_TYPE_CHOICE == COORD_IDENTITY ) 
  xp[1] = x[1];   xp[2] = x[2];  xp[3] = x[3];
  return(0);
#endif

  /*******************************************************************************************************/

#if( COORD_TYPE_CHOICE == COORD_DIAGONAL || \
     COORD_TYPE_CHOICE == COORD_MIXED    || \
     COORD_TYPE_CHOICE == COORD_DIAGONAL3    )

  retval = general_newton_raphson_1d( x, xp);
  //  retval = general_newton_raphson_3d( x, xp);

#endif

#if( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL )
  short int minimize[NDIM];
  minimize[0] = minimize[1] = minimize[2] = minimize[3] = 0;
  minimize[2] = 1; 
  retval = general_newton_raphson_3d_select( x, xp, minimize);
  minimize[2] = 0; 
  minimize[1] = 1; 
  retval = general_newton_raphson_3d_select( x, xp, minimize);
  minimize[1] = 0; 
  minimize[3] = 1; 
  retval = general_newton_raphson_3d_select( x, xp, minimize);

#endif


  if( retval ) { 
    fprintf(stderr,"%s(): Inversion failure x      :  %28.20e %28.20e %28.20e %28.20e \n", __func__,x[0],x[1],x[2],x[3]);
    fprintf(stderr,"%s(): Inversion failure xp_orig:  %28.20e %28.20e %28.20e %28.20e \n", __func__,xp_orig[0],xp_orig[1],xp_orig[2],xp_orig[3]);
    fprintf(stderr,"%s(): Inversion failure xp     :  %28.20e %28.20e %28.20e %28.20e \n", __func__,xp[0],xp[1],xp[2],xp[3]);
    fflush(stderr);
  }

  return( retval ) ;

}

/***********************************************************************************/
/***********************************************************************************
  xp_of_x_wrapper(): 
  ---------
   -- wrapper routine for xp_of_x()  to handle exceptions and automatically 
       supply it with good guesses;
***********************************************************************************/
int  xp_of_x_wrapper(double *xp, double *x, int use_guess )
{
  int i,j,k;
  double x_tmp[NDIM]; 

  DLOOP1 { x_tmp[i] = x[i]; } 

  /*******************************************************************************************************
     Create lookup table for theta coordinate:
  *******************************************************************************************************/
  //--HERE-HERE-HERE --   optimize and implement this later:
#if( COORD_TYPE_CHOICE == COORD_DIAGONAL          ||\
     COORD_TYPE_CHOICE == COORD_DIAGONAL2         ||\
     COORD_TYPE_CHOICE == COORD_DIAGONAL3         ||\
     COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD ||\
     COORD_TYPE_CHOICE == COORD_MIXED                )

//   if( local_first_time) { 
//     i = (int) ((N1S + N1E)/2); 
//     k = (int) ((N3S + N3E)/2); 
//     for(j=N2S; j<=(N2E+1) ; j++ ) { 
//       get_coord(i,j,k,CORN,ncurr,coords);
//       
//     }
//     
//     local_first_time = 0 ; 
//   }

  /* just use a decent guess for now: */ 

  xp[1] = log(x[1]-R0);

  if( use_guess ) { 
    double th = x[2];   
    x_tmp[2] = PERIODIC(th,0.,M_PI);
    xp[2] = x_tmp[2]/M_PI; 
  }

  xp[3] = x[3];

#elif(  COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL    )

  xp[3] = x[3]/(2.*M_PI);

  if( use_guess ) { 
    double th = x[2];   
    x_tmp[2] = PERIODIC(th,0.,M_PI);
    xp[2] = x_tmp[2]/M_PI; 
  }

#endif

  return(xp_of_x(xp,x_tmp));
}

/***********************************************************************************/
/***********************************************************************************
  test_xp_of_x():
  ---------
   -- routine to test performance of xp_of_x inversion;
   -- requires that x = x(xp) is known already; 
   -- we perturb the solution to use it as a new guess to find a "new" solution;
   -- returns with non-zero value if inversion of perturbed guess fails; 
   -- if print_out is set, then we print the comparison of the answer to the 
       solution with the perturbed guess; 
***********************************************************************************/
int  test_xp_of_x(struct of_coord *coords, double perturbation_amplitude, int print_out )
{
  int i,retval; 
  double xp_g[NDIM],xp_g_orig[NDIM],reldiff[NDIM];
  double *xp, *x;
  double dx_dxp[NDIM][NDIM];
  
  xp = coords->xp;
  x  = coords->x;
  dx_dxp_calc(x,xp,dx_dxp);

  /* Make the guess be a random perturbation of the given solution : */ 

  DLOOP1 { 
    xp_g_orig[i] = xp_g[i] = xp[i] * ( 1. + perturbation_amplitude * ( ranc(0) - 0.5 ) ) ; 
  }

  retval = xp_of_x_wrapper(xp_g,x,0);

  DLOOP1 { reldiff[i] = REL_DIFF_FUNC(xp_g[i],xp[i]); }

  if( print_out ) { 
    fprintf(stdout,"x         = %26.16e %26.16e %26.16e %26.16e \n",x[0],x[1],x[2],x[3]); 
    fprintf(stdout,"xp_guess  = %26.16e %26.16e %26.16e %26.16e \n",xp_g_orig[0],xp_g_orig[1],xp_g_orig[2],xp_g_orig[3]); 
    fprintf(stdout,"xp_exact  = %26.16e %26.16e %26.16e %26.16e \n",xp[0],xp[1],xp[2],xp[3]); 
    fprintf(stdout,"xp_numsol = %26.16e %26.16e %26.16e %26.16e \n",xp_g[0],xp_g[1],xp_g[2],xp_g[3]); 
    fprintf(stdout,"reldiff   = %26.16e %26.16e %26.16e %26.16e \n",reldiff[0],reldiff[1],reldiff[2],reldiff[3]); 
//     fprintf(stdout,"dxdxpold  = %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n",
// 	    dx_dxp[1][1],
// 	    dx_dxp[1][2],
// 	    dx_dxp[1][3],
// 	    dx_dxp[2][1],
// 	    dx_dxp[2][2],
// 	    dx_dxp[2][3],
// 	    dx_dxp[3][1],
// 	    dx_dxp[3][2],
// 	    dx_dxp[3][3]);
//     fprintf(stdout,"dxdxpnew  = %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n",
// 	    coords->dx_dxp[1][1],
// 	    coords->dx_dxp[1][2],
// 	    coords->dx_dxp[1][3],
// 	    coords->dx_dxp[2][1],
// 	    coords->dx_dxp[2][2],
// 	    coords->dx_dxp[2][3],
// 	    coords->dx_dxp[3][1],
// 	    coords->dx_dxp[3][2],
// 	    coords->dx_dxp[3][3]);
    fflush(stdout); 
  }

  return(retval);
}

/******************************************************************************
  func_diagonal(): 
  ----------------
     -- routine to calculate the residual and jacobian used by 
           the Newton-Raphson routine general_newton_raphson(), which is 
           used to find X2 from theta of the COORD_DIAGONAL coordinate system;
***********************************************************************************/
void func_diagonal(double *x2, double *xp, double *dx, double *resid, double **jac, 
		   double *f, double *df, double *misc_vars, int n, int n_misc )
{

  double x[NDIM]; 
  double dx_dxp[NDIM][NDIM]; 
  
  xp[2] = x2[0]; 
  x_of_xp(x,xp);
  dx_dxp_calc(x,xp,dx_dxp);

  resid[0]  = x[2] - misc_vars[0];
  jac[0][0] = dx_dxp[2][2]; 

  dx[0] = -resid[0] / jac[0][0];
  *df = -resid[0]*resid[0];
  *f = -0.5*(*df);

}

/******************************************************************************
  phys_val_diagonal():
  -------------------
     -- adjust the given argument (xp2) to lie in its valid range, closest to the original value; 
***********************************************************************************/
void phys_val_diagonal(double *xp2)
{
  double xp2_new;
  xp2_new = (*xp2 > 1.) ? (1. - GNR_SMALL)  :  fabs(*xp2);
  *xp2 = xp2_new; 
  return;
}


/***********************************************************************************/
/***********************************************************************************
  read_and_interp_all_funcs_parallel(): 
  ---------
   -- a distributed version of read_and_interp_all_funcs() that processes the 
       interpolation to the "destination points" in bite-sized portions ; 
   -- also adds communication support;

   -- responsible for reading all grid function data and interpolating it 
      to the local grid and coordinate system; 

***********************************************************************************/
void read_and_interp_all_funcs_parallel(void)
{
  TRACE_BEG;
  
  /**************************************************************************************************
     Interpolate rho, uu, v[1-3], B[1-3] to local cell centers:
        -- 1) read in x[1-3] for all use; 
           2) loop over each grid function and interpolate one at a time; 
           3) after interpolating from Cartesian to numerical coordinates, 
               transform components from Cartesian to numerical coordinate system;
  **************************************************************************************************/
  int i,j,k,l;
  int npts, npts_tot;
  usint *exists;
  double **prim_dest;
  double **coord_list; 
  herr_t status;

  /* Seed random number generator : */ 
  ranc(myid+1); 
  int id_head = out_pid[OUT_HDF5];

  /* USER:  You may want to change this number.  Set it to be large enough so that the 
   batches fit into memory.  
   Rule of thumb: set it to be the number of processors in the new run. */
  int nbatches = 2400;

  fprintf(stdout,"%s(): using  nbatches=[%d]  and  id_head=[%d] \n",__func__,nbatches,id_head); fflush(stdout); 

#if( PERFORM_XP_OF_X_TEST )
  struct of_coord *coords; 
  LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    test_xp_of_x( coords, 0.5, 1 );
    //    fail(FAIL_BASIC,0); 
  }
#endif

  /* Find the list of coordinates that we need to interpolate data to */
  filename = getenv("DEST_INTERP_FILENAME"); 
  if( filename == NULL ) { 
    filename = (char *) malloc(strlen(FILENAME)+1); 
    sprintf(filename,"%s",FILENAME); 
  }
  temp_filename = (char *) malloc(strlen(TEMP_FILENAME)+1); 
  sprintf(temp_filename,"%s",TEMP_FILENAME); 

  fprintf(stdout,"%s():  Using  file      = %s \n",__func__,filename); 
  fprintf(stdout,"%s():  Using  temp_file = %s \n",__func__,temp_filename); 
  fflush(stdout);

  for( i = 0 ; i < nbatches ; i++ ) { 

    coord_list = get_list_of_pts2(filename,&npts, &npts_tot, i, nbatches);  

    /* Transform given list of coordinates to source's numerical coordinates:  */
    transform_dest_to_src_coords(coord_list, npts);

    ALLOC_2D_ARRAY( prim_dest , npts, NP ); 
    ALLOC_ARRAY( exists , npts); 

    interp_prims(coord_list,prim_dest,exists,npts);

    dump_interp_prims( temp_filename, prim_dest, exists, npts, npts_tot, i, nbatches );

    DEALLOC_2D_ARRAY(coord_list,npts,NDIM*2); 
    DEALLOC_2D_ARRAY( prim_dest , npts, NP ); 
    DEALLOC_ARRAY( exists , npts); 
  }

  if( myid == id_head ) { 
    fix_interp_prims_dump(temp_filename,filename,npts_tot); 
  }

  DEALLOC_ARRAY( dims, ndims );

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
  
  /**************************************************************************************************
     Interpolate rho, uu, v[1-3], B[1-3] to local cell centers:
        -- 1) read in x[1-3] for all use; 
           2) loop over each grid function and interpolate one at a time; 
           3) after interpolating from Cartesian to numerical coordinates, 
               transform components from Cartesian to numerical coordinate system;
  **************************************************************************************************/
  int i,j,k,l;
  int npts;
  usint *exists;
  double *data_out; 
  double **prim_dest;
  double **coord_list; 
  herr_t status;

  /* Seed random number generator : */ 
  ranc(myid+1); 
  int id_head = out_pid[OUT_HDF5];

#if( PERFORM_XP_OF_X_TEST )
  struct of_coord *coords; 
  LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    test_xp_of_x( coords, 0.5, 1 );
    //    fail(FAIL_BASIC,0); 
  }
#endif

  
  /* Find the list of coordinates that we need to interpolate data to */
  filename = getenv("DEST_INTERP_FILENAME"); 
  if( filename == NULL ) { 
    filename = (char *) malloc(strlen(FILENAME)+1); 
    sprintf(filename,"%s",FILENAME); 
  }

  fprintf(stdout,"%s():  Using  file      = %s \n",__func__,filename); 
  fflush(stdout);

  coord_list = get_list_of_pts(filename,&npts);  

  /* Transform given list of coordinates to source's numerical coordinates:  */
  transform_dest_to_src_coords(coord_list, npts);

  ALLOC_ARRAY( data_out , npts );
  ALLOC_ARRAY( exists , npts );
  ALLOC_2D_ARRAY( prim_dest , npts, NP ); 

  interp_prims(coord_list,prim_dest,exists,npts);

  /* Now print the data: */ 
  file_id = myH5_Fopen(filename);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][RHO];  }
  status = H5LTmake_dataset_double(file_id,"rho_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][UU];  }
  status = H5LTmake_dataset_double(file_id,"uu_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][U1];  }
  status = H5LTmake_dataset_double(file_id,"v1_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][U2];  }
  status = H5LTmake_dataset_double(file_id,"v2_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][U3];  }
  status = H5LTmake_dataset_double(file_id,"v3_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][B1];  }
  status = H5LTmake_dataset_double(file_id,"B1_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][B2];  }
  status = H5LTmake_dataset_double(file_id,"B2_dest",ndims, dims, data_out);

  for(i=0;i<npts;i++) {  data_out[i] = prim_dest[i][B3];  }
  status = H5LTmake_dataset_double(file_id,"B3_dest",ndims, dims, data_out);

  DEALLOC_ARRAY( data_out , npts );
  DEALLOC_ARRAY( exists , npts );
  DEALLOC_2D_ARRAY( prim_dest , npts, NP ); 


  DEALLOC_2D_ARRAY( coord_list     , npts  , NDIM*2);

  H5Fclose(    file_id);
  
  DEALLOC_ARRAY( dims, ndims );

  TRACE_END;
  return;
}

/****************************************************************************************

 interp_prims(): 
 ---------------------
   -- determines simulation cell indices used for the interpolation of source data
         at each point of the destination system;
   -- records the weights used for the interpolation also;
   -- routine borrowed and modified from bothros;
   -- assumes "weights" is pre-allocated; 
   -- include the ghost zones as they hold relevant data for interpolating, e.g., 
        necessary when doing azimuthal symmetry boundary conditions;

     Arguments: 
           X[]          = location of interpolation;
   xp1,xp2,xp3          = array of parent grid function's coordinates;
          sdata         = grid function to be interpolated, living on Xgrid[];
          dX[]          = grid discretizations in the various directions;
          
****************************************************************************************/
static void interp_prims( double **coord_list, double **prim_dest, usint *exists, int ndest)
{
  int i,j,k,l,n,ind;
  short unsigned int extrapolate;
  double dx1,dx2,dx3,ddx1,ddx2,ddx3,delx_1, delx_2, delx_3;
  double *xsrc;
  double *xphys;
  double xdest[NDIM];
  double weights_loc[8];
  double *prim_src[8];
  double *prim_loc;
  double dxdest_dxsrc[NDIM][NDIM];

  void (*dxdest_dxsrc_calc)( double *xdest , double *xphys, double *xsrc, double dxdest_dxsrc[][NDIM] ); 

  int n1 = N1TOT-2;
  int n2 = N2TOT-2;
  int n3 = N3TOT-2;

  double begx[NDIM];

  TRACE_BEG; 

  /* Coordinates of first cell center: */ 
  begx[1] = startx[1] + (0.5+globalpos[1])*dx[1]; 
  begx[2] = startx[2] + (0.5+globalpos[2])*dx[2]; 
  begx[3] = startx[3] + (0.5+globalpos[3])*dx[3];

  /* Transform between the physical coordinate systems : */
  if( top_type_choice_dest != TOP_TYPE_CHOICE ) { 
    if( (top_type_choice_dest == TOP_SPHERICAL) && 
	(TOP_TYPE_CHOICE      == TOP_CARTESIAN) ) { 
      dxdest_dxsrc_calc = xspher_of_xcart_special3; 
    }
    else if( (top_type_choice_dest == TOP_CARTESIAN) && 
	     (TOP_TYPE_CHOICE     == TOP_SPHERICAL) ) { 
      dxdest_dxsrc_calc = xcart_of_xspher_special3;
    }
    else { 
      fprintf(stdout,"%s(): top_type combination not yet implemented!! :  src=[%d]  dest=[%d] \n", 
	      __func__,top_type_choice_dest, TOP_TYPE_CHOICE); 
      fflush(stdout);  fail(FAIL_BASIC,0); 
    }
  }
  else { 
    dxdest_dxsrc_calc = identity_dxdest_dxsrc;
  }

  /*  Loop over list of destination points, interpolating and transforming vectors along the way :   */
  for(ind = 0 ; ind < ndest; ind++ ) { 
    extrapolate  = 0;
    exists[ind] = 1; 

    prim_loc = prim_dest[ind]; 
    xphys = coord_list[ind]      ;
    xsrc  = coord_list[ind]+NDIM ;

#if( PHASE_SHIFT )
    extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
    struct of_bbh_traj bbh_traj;
    double xbh1, ybh1, BH_PHASE;

    get_bbh_traj_data(&bbh_traj);
    xbh1 = bbh_traj.xi1x;  ybh1 = bbh_traj.xi1y;
    BH_PHASE = atan2(ybh1,xbh1);
    if( BH_PHASE < 0 ) { BH_PHASE += 2. * M_PI; }
    /* Rotate the destination coordinates counter-clockwise by BH_PHASE */
    if( myid == printer_pid && first_phase_calc ) { first_phase_calc = 0; fprintf(stdout,"Checkpoint time is %e adjusting by %e degree\n",t,BH_PHASE*180./M_PI); fflush(stdout); }
    xsrc[3] += BH_PHASE;
    if( xsrc[3] < 0 ) { xsrc[3] += 2. * M_PI; }
    else if( xsrc[3] > 2. * M_PI ) { xsrc[3] -= 2. * M_PI; }
#endif

    /* Distance in each dimension of interpolated point from the first cell center of the local subdomain w.r.t. the grid spacing: */
    delx_1 = (xsrc[1] - begx[1]) * invdx[1] ;
    delx_2 = (xsrc[2] - begx[2]) * invdx[2] ;
    delx_3 = (xsrc[3] - begx[3]) * invdx[3] ;

    /*  Find the indices of the nearest sim_data[COL_X*] coord's smaller than or equal to X[]: */
    /*  k1 -> (i,j),  k2 -> (i+1,j), k3 -> (i+1,j+1), k4 -> (i,j+1) */
    i = (int) floor(delx_1) ; 
    j = (int) floor(delx_2) ; 
    k = (int) floor(delx_3) ; 

    dx1 = delx_1 - i; 
    dx2 = delx_2 - j; 
    dx3 = delx_3 - k; 

    /* Check for out bounds indices, and set the closest allowable value: Do not use ghost zones: */
    if( fabs(dx1) > SMALL ) { 
      if( i  >  n1  ) {  i = n1   ; extrapolate = 1 ; }   
      if( i  <  0   ) {  i = 0    ; extrapolate = 1 ; }   
    }
    if( fabs(dx2) > SMALL ) { 
      if( j  >  n2  ) {  j = n2   ; extrapolate = 1 ; } 
      if( j  <  0   ) {  j = 0    ; extrapolate = 1 ; } 
    }
    if( fabs(dx3) > SMALL ) { 
      if( k  >  n3  ) {  k = n3   ; extrapolate = 1 ; } 
      if( k  <  0   ) {  k = 0    ; extrapolate = 1 ; } 
    }


#if( !USE_EXTRAPOLATION ) 
    /*****************************************************************************
      Set the primitives to a floor state and assume that another subdomain will 
       have the necessary data to overwrite it,  or leave it to be set by the 
       fixup() prescription.  
    ******************************************************************************/ 
    if( extrapolate ) { 

      /* double xtmp[NDIM] ; */
      /* x_of_xp( xtmp, xsrc ) ; */

      /* fprintf(stdout,"extrapolate0: %d : %d : ijk = [%d,%d,%d]  %26.16e  %26.16e  %26.16e\n",myid,ind,i,j,k,xtmp[1],xtmp[2],xtmp[3]);  fflush(stdout); */
      /* fprintf(stdout,"extrapolate1: %d : %d : ijk = [%d,%d,%d]  %26.16e  %26.16e  %26.16e\n",myid,ind,i,j,k,dx1,dx2,dx3);  fflush(stdout); */
      /* fprintf(stdout,"extrapolate2: %d : %d : ijk = [%d,%d,%d]  %26.16e  %26.16e  %26.16e\n",myid,ind,i,j,k,begx[1],begx[2],begx[3]);  fflush(stdout); */
      /* fprintf(stdout,"extrapolate3: %d : %d : ijk = [%d,%d,%d]  %26.16e  %26.16e  %26.16e\n",myid,ind,i,j,k,xsrc[1],xsrc[2],xsrc[3]);  fflush(stdout); */
      /* fprintf(stdout,"extrapolate4: %d : %d : ijk = [%d,%d,%d]  %26.16e  %26.16e  %26.16e\n",myid,ind,i,j,k,begx[1]+n1*dx[1],begx[2]+n2*dx[2],begx[3]+n3*dx[3]);  fflush(stdout); */

      exists[ind] = 0; 
      prim_loc[RHO] = MISSING_VALUE_RHO;
      prim_loc[ UU] = MISSING_VALUE_RHO;
      prim_loc[ U1] = 0.;
      prim_loc[ U2] = 0.;
      prim_loc[ U3] = 0.;
      prim_loc[ B1] = 0.;
      prim_loc[ B2] = 0.;
      prim_loc[ B3] = 0.;
      continue;
    }
#endif

    /* Set the indices of the other cells used in the interpolation: */
    prim_src[0] = p[i  ][j  ][k  ];   /* corner		          */
    prim_src[1] = p[i+1][j  ][k  ];   /* corner + x1	      	  */
    prim_src[2] = p[i+1][j+1][k  ];   /* corner + x1 + x2         */
    prim_src[3] = p[i  ][j+1][k  ];   /* corner      + x2         */
    prim_src[4] = p[i  ][j  ][k+1];   /* corner           + x3    */
    prim_src[5] = p[i+1][j  ][k+1];   /* corner + x1      + x3    */
    prim_src[6] = p[i+1][j+1][k+1];   /* corner + x1 + x2 + x3    */
    prim_src[7] = p[i  ][j+1][k+1];   /* corner      + x2 + x3    */

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
  
    PLOOP {  
      dx1 = 0. ;   /* temporary  variable */
      for(n=0;n<8;n++) { dx1 += weights_loc[n] * prim_src[n][l]; }
      prim_loc[l] = dx1 ; 
    }

    /* Transform vector components into destination's physical : */ 

    dxdest_dxsrc_calc( xdest, xphys, xsrc, dxdest_dxsrc) ; 

    dx1 = 
      prim_loc[U1] * dxdest_dxsrc[1][1] + 
      prim_loc[U2] * dxdest_dxsrc[1][2] + 
      prim_loc[U3] * dxdest_dxsrc[1][3] ;
    
    dx2 = 
      prim_loc[U1] * dxdest_dxsrc[2][1] + 
      prim_loc[U2] * dxdest_dxsrc[2][2] + 
      prim_loc[U3] * dxdest_dxsrc[2][3] ;
    
    dx3 = 
      prim_loc[U1] * dxdest_dxsrc[3][1] + 
      prim_loc[U2] * dxdest_dxsrc[3][2] + 
      prim_loc[U3] * dxdest_dxsrc[3][3] ;
    

    ddx1 = 
      prim_loc[B1] * dxdest_dxsrc[1][1] + 
      prim_loc[B2] * dxdest_dxsrc[1][2] + 
      prim_loc[B3] * dxdest_dxsrc[1][3] ;
    
    ddx2 = 
      prim_loc[B1] * dxdest_dxsrc[2][1] + 
      prim_loc[B2] * dxdest_dxsrc[2][2] + 
      prim_loc[B3] * dxdest_dxsrc[2][3] ;
    
    ddx3 = 
      prim_loc[B1] * dxdest_dxsrc[3][1] + 
      prim_loc[B2] * dxdest_dxsrc[3][2] + 
      prim_loc[B3] * dxdest_dxsrc[3][3] ;
    
      prim_loc[U1] =  dx1;  
      prim_loc[U2] =  dx2;  
      prim_loc[U3] =  dx3; 
      prim_loc[B1] = ddx1;  
      prim_loc[B2] = ddx2;  
      prim_loc[B3] = ddx3; 

  }

  TRACE_END; 

  return;
}

/***********************************************************************************/
/***********************************************************************************
  get_list_of_pts():
  -----------------
  -- Reads a list of coordinates (destination) to which we will interpolate;
  -- this routine is basically for testing purposes as will be passing data via MPI 
     in actual runs; 
***********************************************************************************/
static double **get_list_of_pts( char *filename, int *npts )
{
  int i, j, k, nloc;
  double *data_in;
  double **coord_list;
  struct of_coord *coords; 
  
  TRACE_BEG; 

  /******************************************************************************
    Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  /* Open existing file : wrapper handles error checking */ 
  file_id = myH5_Fopen(filename);

  myH5_read_scalar2(file_id,"/Header/Macros/TOP_TYPE_CHOICE" 	 ,H5T_NATIVE_INT,&top_type_choice_dest           );
  myH5_read_scalar2(file_id,"/Header/Macros/COORD_TYPE_CHOICE"   ,H5T_NATIVE_INT,&coord_type_choice_dest         );

  fprintf(stdout,"%s(): using   top_type_choice_dest = %d \n",__func__,top_type_choice_dest); 
  fprintf(stdout,"%s(): using coord_type_choice_dest = %d \n",__func__,coord_type_choice_dest); 
  fflush(stdout);

  /* Find dimensions: */
  dims = myH5_get_dims( file_id , "x1", &ndims );

  nloc = 1; 
  for(i=0;i<ndims;i++) {  nloc *= dims[i];  }

  /* Allocate arrays: */
  ALLOC_ARRAY(data_in,nloc); 
  ALLOC_2D_ARRAY(coord_list,nloc,NDIM*2); 
  for(i=0;i<nloc;i++) { coord_list[i][0] = t; }

  /* Read in coordinates and assign coord_list one by one: */
  myH5_read_basic( file_id, "x1", H5T_NATIVE_DOUBLE, data_in ) ;
  for(i=0;i<nloc;i++) { coord_list[i][1] = data_in[i]; }

  myH5_read_basic( file_id, "x2", H5T_NATIVE_DOUBLE, data_in ) ;
  for(i=0;i<nloc;i++) { coord_list[i][2] = data_in[i]; }

  myH5_read_basic( file_id, "x3", H5T_NATIVE_DOUBLE, data_in ) ;
  for(i=0;i<nloc;i++) { coord_list[i][3] = data_in[i]; }

  *npts = nloc; 

  DEALLOC_ARRAY(data_in,nloc); 
  H5Fclose(    file_id);

  TRACE_END; 
  return(coord_list);
}

/***********************************************************************************/
/***********************************************************************************
  get_list_of_pts2():
  -----------------
  -- version that only reads in a part or batch of the points;
  -- Reads a list of coordinates (destination) to which we will interpolate;
  -- this routine is basically for testing purposes as will be passing data via MPI 
     in actual runs; 
***********************************************************************************/
static double **get_list_of_pts2( char *filename, int *npts, int *npts_tot, int ibatch, int nbatches )
{
  int i, j, k;
  int npts_tmp, ibeg, iend ;
  double *data_in;
  double **coord_list;

  struct of_coord *coords; 

  static int nloc;
  double *data;
  double **sendv;

  TRACE_BEG; 

  if( ibatch >= nbatches ) { 
    fprintf(stdout,"%s(): ibatch must be smaller than nbatches :  %d %d \n",__func__,ibatch,nbatches);
    fflush(stdout);  fail(FAIL_BASIC,0); 
  }

  /******************************************************************************
    Open hdf5 object (e.g. file, dataspace, dataset, ...)
  ******************************************************************************/ 
  if( myid == id_head ) { 
    file_id = myH5_Fopen(filename);
  }

  if( ibatch == 0 ) { 
    myH5_read_scalar2(file_id,"/Header/Macros/TOP_TYPE_CHOICE"  ,H5T_NATIVE_INT,&top_type_choice_dest  );
    myH5_read_scalar2(file_id,"/Header/Macros/COORD_TYPE_CHOICE",H5T_NATIVE_INT,&coord_type_choice_dest);
    /* If phase shifting do not read in the time from the destination (we assume t = 0) */
#if( !PHASE_SHIFT )
    myH5_read_scalar2(file_id,"/Header/Grid/t"                  ,H5T_NATIVE_DOUBLE,&t);
#endif
    /* Find dimensions: do this once so that we don't create a memory leak with dims[] */
    dims = myH5_get_dims( file_id , "x1", &ndims );

    nloc = 1; 
    for(i=0;i<ndims;i++) {  nloc *= dims[i];  }
  }

  /* npts_tmp = nloc / (nbatches-1);  */
  /* ibeg = ibatch * npts_tmp;  */
  /* iend = ibeg + npts_tmp - 1; */
  /* i = nloc-1; */
  /* iend = MIN(iend,i); */
  /* npts_tmp = iend - ibeg + 1;  */

  npts_tmp = nloc / nbatches ; 
  ibeg = ibatch * npts_tmp; 
  iend = ibatch != nbatches-1 ?
    ibeg + npts_tmp - 1 : nloc - 1 ;
  //    (ibatch + 1) * npts_tmp - 1: nloc - 1 ;

  npts_tmp = iend - ibeg + 1; 

  if( myid == id_head ) { 
    fprintf(stdout,"%s(): ibatch=[%6d/%6d]  using   top_type_choice_dest = %d \n",__func__,ibatch,nbatches,top_type_choice_dest); 
    fprintf(stdout,"%s(): ibatch=[%6d/%6d]  using coord_type_choice_dest = %d \n",__func__,ibatch,nbatches,coord_type_choice_dest); 
    fprintf(stdout,"%s(): ibatch=[%6d/%6d]  using indices [%10d,%10d]  of %d \n",__func__,ibatch,nbatches,ibeg,iend,nloc);
    fflush(stdout);
  }

  *npts     = npts_tmp;
  *npts_tot = nloc;

  ALLOC_2D_ARRAY(sendv,NDIM,npts_tmp); 

  /******************************************************************************
    Read in the coordinate data :
  ******************************************************************************/ 
  if( myid == id_head ) { 
    ALLOC_ARRAY(data_in,nloc); 

    j = 0; 
    for(i=ibeg;i<=iend;i++) { sendv[0][j++] = t;}

    /* Read in coordinates and assign coord_list one by one: */
    myH5_read_basic( file_id, "x1", H5T_NATIVE_DOUBLE, data_in ) ;
    j = 0; 
    for(i=ibeg;i<=iend;i++) { sendv[1][j++] = data_in[i]; }

    myH5_read_basic( file_id, "x2", H5T_NATIVE_DOUBLE, data_in ) ;
    j = 0; 
    for(i=ibeg;i<=iend;i++) { sendv[2][j++] = data_in[i]; }

    myH5_read_basic( file_id, "x3", H5T_NATIVE_DOUBLE, data_in ) ;
    j = 0; 
    for(i=ibeg;i<=iend;i++) { sendv[3][j++] = data_in[i]; }

    DEALLOC_ARRAY(data_in,nloc); 
    H5Fclose(    file_id);
  }

  /* Broadcast data to others processes : */
  DLOOP1 { sync_vect_from_rank(sendv[i],npts_tmp,id_head); }

  /* Transpose data to more usable form: */
  ALLOC_2D_ARRAY(coord_list,npts_tmp,NDIM*2); 
  for(i=0;i<npts_tmp;i++) { 
    coord_list[i][0] = sendv[0][i] ;  
    coord_list[i][1] = sendv[1][i] ;  
    coord_list[i][2] = sendv[2][i] ;  
    coord_list[i][3] = sendv[3][i] ;  
  }

  DEALLOC_2D_ARRAY(sendv,NDIM,npts_tmp); 

  TRACE_END; 
  return(coord_list);
}

/***********************************************************************************/
/***********************************************************************************
  get_list_of_pts_old():
  -----------------
  -- Reads a list of coordinates (destination) to which we will interpolate;
  -- this routine is basically for testing purposes as will be passing data via MPI 
     in actual runs; 
***********************************************************************************/
double **get_list_of_pts_old( int *npts )
{
  int i, j, k, ind, nloc;
  double *ctmp; 
  double **coord_list;
  struct of_coord *coords; 
  
  TRACE_BEG; 

  nloc = NTOT; 

  ALLOC_2D_ARRAY(coord_list,nloc,NDIM*2); 

  ind = 0 ; 
  ALL_LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    ctmp = coord_list[ind];
    ctmp[0] = t; 
    ctmp[1] = coords->x[1];
    ctmp[2] = coords->x[2];
    ctmp[3] = coords->x[3];
    ind++;
  }

  *npts = nloc; 

  TRACE_END; 
  return(coord_list);
}


/***********************************************************************************/
/***********************************************************************************
  dump_interp_prims():
  -----------------
  -- version that only reads in a part or batch of the points;
  -- Reads a list of coordinates (destination) to which we will interpolate;
  -- this routine is basically for testing purposes as will be passing data via MPI 
     in actual runs; 
***********************************************************************************/
static void dump_interp_prims( char *filename, double **prim_dest, usint *exists, int npts, int npts_tot, int ibatch, int nbatches )
{
  int i, j, k,l;
  int *mask;
  double *sendv;

  extern void collect_mask_values(double *val, int *mask, int n, int rank );
  extern hid_t    myH5_Fcreate_gather(char *name, int overwrite);

  static ulint offset=0;

  TRACE_BEG; 

  /******************************************************************************
    Check on parameters: 
  ******************************************************************************/ 
  if( ibatch >= nbatches ) { 
    fprintf(stdout,"%s(): ibatch must be smaller than nbatches :  %d %d \n",__func__,ibatch,nbatches);
    fflush(stdout);  fail(FAIL_BASIC,0); 
  }

  /* int npts_tmp = npts_tot / (nbatches-1);  */
  /* int ibeg = ibatch * npts_tmp;  */
  /* int iend = ibeg + npts_tmp - 1; */
  /* i = npts_tot-1; */
  /* iend = MIN(iend,i); */
  /* npts_tmp = iend - ibeg + 1;  */

  int npts_tmp = npts_tot / nbatches; 
  int ibeg = ibatch * npts_tmp; 
  int iend = ibatch != nbatches-1 ?
    ibeg + npts_tmp - 1 : npts_tot-1;
  //    (ibatch + 1) * npts_tmp - 1 : npts_tot-1;

  npts_tmp = iend - ibeg + 1; 

  if( npts != npts_tmp ) { 
    fprintf(stdout,"%s(): npts mismatch :  %d %d \n",__func__,npts,npts_tmp); 
    fflush(stdout);  fail(FAIL_BASIC,0); 
  }

  /******************************************************************************
   Format data to reduce, masking out missing data :
  ******************************************************************************/ 
  int nsend = npts*NP;
  ALLOC_ARRAY(sendv,nsend); 
  ALLOC_ARRAY(mask ,nsend); 

  j = 0; 
  for( i=0; i<npts; i++ ) { 
    k = exists[i];
    PLOOP { 
      sendv[j] = prim_dest[i][l]; 
      mask[j]  = k;
      j++;
    }
  }

  collect_mask_values(sendv,mask,nsend,id_head); 

  /*
  int ind ;
  for(ind = 0; ind < nsend; ind++) {
    fprintf(stdout, "id = %d, mask[%d] = %d \n", myid, ind, mask[ind]) ; fflush(stdout) ;
  }
  */
  
  /******************************************************************************
    Write unformatted data for later reformatting in the temporary location:
  ******************************************************************************/ 
  if( myid == id_head ) { 
    static char fnames[NP][200]; 
    static usint local_first_time = 1 ; 

    hsize_t filedims[1], memdims[1], hoffset[1]; 
    hid_t   dataset_id, filespace_id;

    fprintf(stdout,"%s(): dumping ibatch=[%6d/%6d] \n",__func__,(ibatch+1),nbatches);
    fflush(stdout);

    memdims[ 0] = npts;
    filedims[0] = npts_tot;
    hoffset[ 0] = offset;

    if( local_first_time ) { 
      file_id = myH5_Fcreate_gather(filename, 0);
      filespace_id = H5Screate_simple(1, filedims, NULL);
      PLOOP { 
	sprintf(fnames[l],"%s-%1d",PRIM_DEST_TMP_NAME,l);
	dataset_id = H5Dcreate(file_id,fnames[l], H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT);
	H5Dclose(dataset_id);
      }
      dataset_id = H5Dcreate(file_id,MASK_TMP_NAME , H5T_NATIVE_INT, filespace_id, H5P_DEFAULT);
      H5Dclose(dataset_id);
      H5Sclose(filespace_id);

      local_first_time = 0; 
    }
    else { 
      file_id = myH5_Fopen(filename);
    }

    hid_t  memspace_id = H5Screate_simple(1,  memdims, NULL);
    double *masktmp;
    double **primttmp;

    ALLOC_ARRAY(masktmp,npts);
    ALLOC_2D_ARRAY(primttmp,NP,npts);

    j=0;
    for( i=0; i<npts; i++ ) { 
      masktmp[i] = mask[j];
      PLOOP{ primttmp[l][i] = sendv[j++];  }
    }

    PLOOP { 
      dataset_id = H5Dopen(file_id, fnames[l]);
      filespace_id = H5Dget_space(dataset_id);  
      H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, hoffset, NULL, memdims, NULL); 
      H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id,  H5P_DEFAULT, primttmp[l]);
      H5Dclose(dataset_id);
    }

    dataset_id = H5Dopen(file_id, MASK_TMP_NAME ); 
    filespace_id = H5Dget_space(dataset_id);  
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, hoffset, NULL, memdims, NULL); 
    H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace_id, filespace_id,  H5P_DEFAULT, mask);
    H5Dclose(dataset_id);
    H5Sclose(filespace_id);
    H5Sclose(memspace_id);
    H5Fclose(file_id);

    DEALLOC_ARRAY(masktmp,npts);
    DEALLOC_2D_ARRAY(primttmp,NP,npts);

    offset += npts;
  }

  DEALLOC_ARRAY(sendv,nsend); 
  DEALLOC_ARRAY(mask ,nsend); 


  TRACE_END; 
  return;
}
/***********************************************************************************/
/***********************************************************************************
  fix_interp_prims_dump():
  -----------------
  -- after all the data has been interpolated and written out to the temporary file, 
        read it back in, format it, and write it to the original file; 
***********************************************************************************/
void fix_interp_prims_dump( char *temp_filename, char *filename, int npts_tot )
{
  int i, j, k,l;
  int *mask;
  double *data;
  char tmpfnames[NP][200]; 
  char fnames[NP][200]; 
  herr_t status;

  TRACE_BEG; 

  /*************************************************************************************
    Allocate memory and read in data : 
   *************************************************************************************/

  PLOOP { sprintf(tmpfnames[l],"%s-%1d",PRIM_DEST_TMP_NAME,l); }
  l = 0; 
  sprintf(fnames[l++],"rho_dest");
  sprintf(fnames[l++], "uu_dest");
  sprintf(fnames[l++], "v1_dest");
  sprintf(fnames[l++], "v2_dest");
  sprintf(fnames[l++], "v3_dest");
  sprintf(fnames[l++], "B1_dest");
  sprintf(fnames[l++], "B2_dest");
  sprintf(fnames[l++], "B3_dest");

  hid_t tmpfile_id = myH5_Fopen(temp_filename);
  hid_t    file_id = myH5_Fopen(filename);

  ALLOC_ARRAY(data,npts_tot); 
  ALLOC_ARRAY(mask,npts_tot); 

  PLOOP {
    status = H5LTread_dataset_double( tmpfile_id, tmpfnames[l] , data );
    if( status < 0 ) { 
      fprintf(stdout,"%s() %d %d %d : problem with reading dataset: %s  %d \n",__func__,myid,nstep,n_substep,tmpfnames[l],status); 
      fflush(stdout); fail(FAIL_HDF,0); 
    }
    H5LTmake_dataset_double(file_id, fnames[l], ndims, dims, data);
  }

  status = H5LTread_dataset_int( tmpfile_id, MASK_TMP_NAME , mask );
  if( status < 0 ) { 
    fprintf(stdout,"%s() %d %d %d : problem with reading dataset: %s  %d \n",__func__,myid,nstep,n_substep,MASK_TMP_NAME,status); 
    fflush(stdout); fail(FAIL_HDF,0); 
  }
  H5LTmake_dataset_int(file_id, MASK_TMP_NAME, ndims, dims, mask);

  H5Fclose( tmpfile_id);
  H5Fclose(    file_id);
  
  DEALLOC_ARRAY(data,npts_tot); 
  DEALLOC_ARRAY(mask,npts_tot); 

  TRACE_END; 
  return;
}

/***********************************************************************************/
/***********************************************************************************
  transform_dest_to_src_coords():
  ------------------------------
  -- Transform given list of coordinates to source's numerical coordinates: 
***********************************************************************************/
static void transform_dest_to_src_coords( double **coord_list, int npts ) 
{
  int i; 
  double *x, xp[NDIM]; 
  double xorig[NDIM], *xsrc;
  
  void (*phys_coord_transform)( double *xdest , double *xsrc ); 

  TRACE_BEG; 

  /* Transform between the physical coordinate systems : */
  if( top_type_choice_dest != TOP_TYPE_CHOICE ) { 
    if( (top_type_choice_dest == TOP_SPHERICAL) && 
	(TOP_TYPE_CHOICE      == TOP_CARTESIAN) ) { 
      phys_coord_transform = xcart_of_xspher_only; 
      fprintf(stdout,"Using xcart_of_xspher_only \n"); fflush(stdout);
    }
    else if( (top_type_choice_dest == TOP_CARTESIAN) && 
	     (TOP_TYPE_CHOICE     == TOP_SPHERICAL) ) { 
      phys_coord_transform = xspher_of_xcart_only; 
      fprintf(stdout,"Using xspher_of_xcart_only \n"); fflush(stdout);
    }
    else { 
      fprintf(stdout,"%s(): top_type combination not yet implemented!! :  src=[%d]  dest=[%d] \n", 
	      __func__,top_type_choice_dest, TOP_TYPE_CHOICE); 
      fflush(stdout);  fail(FAIL_BASIC,0); 
    }
  }
  else { 
    phys_coord_transform = identity_phys_transform;
  }

  /* Now transform from local (dest) physical coordinates to local (dest) numerical coordinates : */
  for( i = 0 ; i < npts; i++ ) { 
    x    = coord_list[i]     ;
    xsrc = coord_list[i]+NDIM;
    xorig[0] = x[0];
    xorig[1] = x[1];
    xorig[2] = x[2];
    xorig[3] = x[3];

#if( PRINT_DEBUG_INFO )
      fprintf(stdout,"%s(): xcart = %26.16e %26.16e %26.16e %26.16e \n",__func__,xsrc[0],xsrc[1],xsrc[2],xsrc[3]); fflush(stdout);
#endif
      
      phys_coord_transform(x,xorig); 
      xp_of_x_wrapper(xsrc,x,1);

#if( PRINT_DEBUG_INFO )
      fprintf(stdout,"%s(): xsphr = %26.16e %26.16e %26.16e %26.16e \n",__func__,xsrc[0],xsrc[1],xsrc[2],xsrc[3]); fflush(stdout);
#endif

  }
  
  TRACE_END; 
  return;
}

/********************************************************************** 
/************************************************************

  general_newton_raphson_1d(): 

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  X2 from theta 
        by ensuring that   0 <= X2 <= 1  (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson_1d( double *x_phys, double *xp  )
{
  double x, df, dx, x_old, resid, jac,errx, x_orig;
  double x_tmp[NDIM];
  double dx_dxp[NDIM][NDIM];
  int    n_iter, id, jd, i_extra, doing_extra;

  int   keep_iterating;


  // Initialize various parameters and variables:
  i_extra = doing_extra = 0;
  x = xp[2];

  if( (x < 0.) || (x > 1.) ) {   x = ranc(0);  }  /* make sure initial guess is physical */

  x_old = x_orig = x;

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    xp[2] = x;
    x_of_xp(x_tmp,xp);
    dx_dxp_calc(x_tmp,xp,dx_dxp);

    resid  = x_tmp[2] - x_phys[2];

    //    jac = dx_dxp[2][2]; 
    //    dx = -resid / jac;

    dx = -resid / dx_dxp[2][2];
    df = resid*resid;

    /* Save old values before calculating the new: */
    errx = 0.;
    x_old = x;

    /* don't use line search : */
    x += dx;

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx = (x==0.) ?  fabs(dx) : fabs(dx/x);

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    if( x < 0. ) { 
      x = ranc(0) * x_old;  /* random number between old value and 0. */
    }
    else if( x > 1. ) { 
      x = 1. + ranc(0)*(x_old-1.);   /* random number between old value and 1. */
    }
    
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (errx <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( (errx <= NEWT_TOL) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    fprintf(stdout,"n,xo,x,dx,df,res,errx = %10d %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e \n",n_iter,x_old,x,dx,df,resid,errx); fflush(stdout);

    n_iter++;

  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
  if( finite(df)==0    ) {  
    return(2);  
  }
  if( errx <= NEWT_TOL ) {  
    return(0);  
  }
  if( (errx <= MIN_NEWT_TOL) && (errx > NEWT_TOL) ){   
    return(0);  
  }
  if( errx > MIN_NEWT_TOL){
    return(1);
  }

  return(0);

}

/********************************************************************** 
/************************************************************

  general_newton_raphson_3d_orig(): 

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  X2 from theta and X3 from phi
        by ensuring that   0 <= X2 <= 1  (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson_3d_orig( double *x_phys, double *xp  )
{
  double df, dx[NDIM], x_old[NDIM], x_tmp[NDIM], x_orig[NDIM], resid[NDIM], errx, errx_try[NDIM], dx_dxp[NDIM][NDIM];
  double errx_old;
  int    i, n_iter, i_extra, doing_extra, keep_iterating, permute[NDIM];

  int n_survey = 200; 
  double x1_min=1e30;
  double errx_min=1e30;
  double x1_max;

  static int loc_iter = 0; 

#if( COORD_TYPE_CHOICE != COORD_WARPED_SPHERICAL )
  fprintf(stdout,"%s(): Probably should use another routine to perform coordinate inversion!  \n",__func__); 
  fflush(stdout); fail(FAIL_BASIC,0); 
#endif

  // Initialize various parameters and variables:
  errx = 1.e20 ; 
  i_extra = doing_extra = 0;

  if( (xp[2] < 0.) || (xp[2] > 1.) ) {   xp[2] = ranc(0);  }  /* make sure initial guess is physical */
  if( (xp[3] < 0.) || (xp[3] > 1.) ) {   xp[3] = ranc(0);  }  /* make sure initial guess is physical */
  

  DLOOP1 {  x_old[i] = x_orig[i] = xp[i] ; }

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    x_of_xp(x_tmp,xp);
    dx_dxp_calc(x_tmp,xp,dx_dxp);

    DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
    DLOOP1 {     dx[i] = -resid[i]            ; } 

    LU_decompose(dx_dxp,permute);
    LU_substitution( dx_dxp, dx, permute );

    df = 0. ; 
    DLOOP1 {  df += resid[i]*resid[i] ; } 
    
    /* Save old values before calculating the new: */
    errx = 0.;
    DLOOP1 {  x_old[i] = xp[i] ; }


    /* don't use line search : */
    DLOOP1 {  xp[i] += dx[i]  ; } 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx_old = errx; 
    errx = -1.e20;
    DLOOP1 { 
      errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]);
      if( errx_try[i] > errx ) { errx = errx_try[i]; } 
    }

    if( (errx > errx_old) && (n_iter > 5) ) {
      DLOOP1 {  xp[i] -= (ranc(0))*dx[i]; }
#if( PRINT_NEWTON_RAPHSON_INFO ) 
      fprintf(stdout,"not-using-full-step\n"); fflush(stdout);
#endif
    }

    if( errx_min > errx ) { 
      errx_min = errx;
      x1_min = xp[1]; 
    }

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    /*  I can't think of invalid values of xp,  Newton-Raphson should find the solution...  */
    if( xp[2] < 0. ) { 
      xp[2] = ranc(0) * x_old[2];  /* random number between old value and 0. */
    }
    else if( xp[2] > 1. ) { 
      xp[2] = 1. + ranc(0)*(x_old[2]-1.);   /* random number between old value and 1. */
    }

    if( xp[3] < 0. ) { 
      xp[3] = ranc(0) * x_old[3];  /* random number between old value and 0. */
    }
    else if( xp[3] > 1. ) { 
      xp[3] = 1. + ranc(0)*(x_old[3]-1.);   /* random number between old value and 1. */
    }

    
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (errx <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( (errx <= NEWT_TOL) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    fprintf(stdout,"xp1-resid-survey1-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],dx[1],x_tmp[1],resid[1],errx,1./dx_dxp[1][1]);
    fflush(stdout);

    fprintf(stdout,"xp      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,xp[0],xp[1],xp[2],xp[3]); 
    fprintf(stdout,"dx      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,dx[0],dx[1],dx[2],dx[3]); 
    fprintf(stdout,"resid   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,resid[0],resid[1],resid[2],resid[3]); 
    fprintf(stdout,"x_phys  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_phys[0],x_phys[1],x_phys[2],x_phys[3]); 
    fprintf(stdout,"x_tmp   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_tmp[0],x_tmp[1],x_tmp[2],x_tmp[3]); 
    fprintf(stdout,"errx_i  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,errx_try[0],errx_try[1],errx_try[2],errx_try[3]); 
    fprintf(stdout,"df,errx = %10d %26.16e  %26.16e \n",n_iter,df,errx); 
    fflush(stdout);
#endif

    n_iter++;

  }   // END of while(keep_iterating)

#if( PRINT_NEWTON_RAPHSON_INFO ) 
  x1_max = 1.01 * x1_min; 
  x1_min *= 0.99; 
  double dxp1 = (x1_max - x1_min)/(n_survey-1);
  for( n_iter = 0 ; n_iter < n_survey; n_iter++ ) { 
    xp[1] = x1_min + n_iter*dxp1; 
    x_of_xp(x_tmp,xp);    
    dx_dxp_calc(x_tmp,xp,dx_dxp);
    DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
    DLOOP1 {     dx[i] = -resid[i]            ; } 

    LU_decompose(dx_dxp,permute);
    LU_substitution( dx_dxp, dx, permute );

    errx_try[1]  = (xp[1]==0.) ?  fabs(dx[1]) : fabs(dx[1]/xp[1]);
    fprintf(stdout,"xp1-resid-survey2-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],dx[1],x_tmp[1],resid[1],errx_try[1],1./dx_dxp[1][1]);
    fflush(stdout);
  }
  loc_iter++;
#endif

  /*  Check for bad untrapped divergences : */
  if( finite(df)==0    ) {  
    return(2);  
  }
  if( errx <= NEWT_TOL ) {  
    return(0);  
  }
  if( (errx <= MIN_NEWT_TOL) && (errx > NEWT_TOL) ){   
    return(0);  
  }
  if( errx > MIN_NEWT_TOL){
    return(1);
  }

  return(0);

}
/********************************************************************** 
/************************************************************

  general_newton_raphson_3d_monotonic(): 

    -- special version of general_newton_raphson_3d_select()
            that checks for non-monotocity in the residual; 

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  X2 from theta and X3 from phi
        by ensuring that   0 <= X2 <= 1  (look near "METHOD specific:")
        by ensuring that   0 <= X3 <= 2*Pi  (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson_3d_monotonic( double *x_phys, double *xp, short int *minimize, short int enforce_monotonicity)
{
  double df, dx[NDIM], x_old[NDIM], x_tmp[NDIM], x_orig[NDIM], resid[NDIM], errx, errx_try[NDIM], dx_dxp[NDIM][NDIM];
  double errx_old;
  int    i, n_iter, i_extra, doing_extra, keep_iterating, permute[NDIM];

  double bbox[NDIM][2];

#if( PRINT_NEWTON_RAPHSON_INFO ) 
  double x1_min=1e30;
  double x3_min=1e30;
  double errx_min=1e30;
  double x1_max;
  double x3_max;
  static int loc_iter = 0; 
#endif


#if( COORD_TYPE_CHOICE != COORD_WARPED_SPHERICAL )
  fprintf(stdout,"%s(): Probably should use another routine to perform coordinate inversion!  \n",__func__); 
  fflush(stdout); fail(FAIL_BASIC,0); 
#endif

  // Initialize various parameters and variables:
  errx = 1.e20 ; 
  errx_old = 2*errx;
  i_extra = doing_extra = 0;

  /*   METHOD specific: */
  if( minimize[2] && ( (xp[2] < 0.) || (xp[2] > 1.) ) ) {   xp[2] = ranc(0);  }  /* make sure initial guess is physical */
  if( minimize[3] && ( (xp[3] < 0.) || (xp[3] > 1.) ) ) {   xp[3] = ranc(0);  }  /* make sure initial guess is physical */
  
  DLOOP1 {  x_old[i] = x_orig[i] = xp[i] ; }

  DLOOP1 {  bbox[i][0] = xp[i];  bbox[i][1] = 2.*fabs(bbox[i][0]); }

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    x_of_xp(x_tmp,xp);
    dx_dxp_calc(x_tmp,xp,dx_dxp);
    DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
    DLOOP1 {     dx[i] = -resid[i]            ; } 
    DLOOP1 { if( !minimize[i] ) {   dx[i] = 0. ; } }
    //   LU_decompose(dx_dxp,permute);     LU_substitution( dx_dxp, dx, permute );
    DLOOP1 {     dx[i] /= dx_dxp[i][i]            ; } 
    
    /* assume that x(xp) is monotonic for all xp[i]  */
    /*   METHOD specific: -- not all systems are necessarily monotonic */
    if( enforce_monotonicity ) { 
      DLOOP1 if( minimize[i] ) {  
	if( x_tmp[i] > x_phys[i]  ) { 
	  if( bbox[i][1]  > xp[i]  ) { bbox[i][1] = xp[i]; }
	}
	else { 
	  if( bbox[i][0]  < xp[i]  ) { bbox[i][0] = xp[i]; }
	}
	if( bbox[i][0] > bbox[i][1] ) { 
	  keep_iterating = 0; 
	  fprintf(stdout,"exitting from monotonicity criterion:  %26.16e  %26.16e \n", bbox[i][0], bbox[i][1]); fflush(stdout);
	}
      }
    }

    /* Save old values before calculating the new: */
    DLOOP1 {  x_old[i] = xp[i] ; }

    /* don't use line search : */
    DLOOP1 {  xp[i] += dx[i]  ; } 

    df = 0. ; 
    DLOOP1 {  df += resid[i]*resid[i] ; } 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx_old = errx; 
    errx = -1.e20;
    DLOOP1 { 
      errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]);
      if( errx_try[i] > errx)  { errx = errx_try[i]; } 
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    if( errx_min > errx ) { 
      errx_min = errx;
      x1_min = xp[1]; 
      x3_min = xp[3];
    }
#endif

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    /*   METHOD specific: */
    /*  I can't think of invalid values of xp,  Newton-Raphson should find the solution...  */
    if( xp[2] < 0. ) { 
      xp[2] = ranc(0) * x_old[2];  /* random number between old value and 0. */
    }
    else if( xp[2] > 1. ) { 
      xp[2] = 1. + ranc(0)*(x_old[2]-1.);   /* random number between old value and 1. */
    }

    if( xp[3] < 0. ) { 
      xp[3] = ranc(0) * x_old[3];  /* random number between old value and 0. */
    }
    else if( xp[3] > 1. ) { 
      xp[3] = 1. + ranc(0)*(x_old[3]-1.);   /* random number between old value and 1. */
    }

    
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (errx <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( (errx <= NEWT_TOL) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    fprintf(stdout,"xp1-resid-survey1-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],xp[3],errx_try[1],errx_try[3],resid[1],resid[3],dx[1],dx[3]); 
    fflush(stdout);

    fprintf(stdout,"xp      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,xp[0],xp[1],xp[2],xp[3]); 
    fprintf(stdout,"dx      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,dx[0],dx[1],dx[2],dx[3]); 
    fprintf(stdout,"resid   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,resid[0],resid[1],resid[2],resid[3]); 
    fprintf(stdout,"x_phys  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_phys[0],x_phys[1],x_phys[2],x_phys[3]); 
    fprintf(stdout,"x_tmp   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_tmp[0],x_tmp[1],x_tmp[2],x_tmp[3]); 
    fprintf(stdout,"errx_i  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,errx_try[0],errx_try[1],errx_try[2],errx_try[3]); 
    fprintf(stdout,"df,errx = %10d %26.16e  %26.16e \n",n_iter,df,errx); 
    fflush(stdout);
#endif

    n_iter++;

  }   // END of while(keep_iterating)

#if( PRINT_NEWTON_RAPHSON_INFO ) 
   double eps=5.e-8; 
   x1_max = (1.+eps) * x1_min; 
   x1_min *= (1.-eps);
   x3_max = (1.+eps) * x3_min; 
   x3_min *= (1.-eps);
   int n_survey = 400; 
   double dxp1 = (x1_max - x1_min)/(n_survey-1);
   double dxp3 = (x3_max - x3_min)/(n_survey-1);
   int m_iter; 
   double dxdxp11;
   m_iter = n_survey/2;
   for( n_iter = 0 ; n_iter < n_survey; n_iter++ ) { 
     //     for( m_iter = 0 ; m_iter < n_survey; m_iter++ ) { 
       xp[1] = x1_min + n_iter*dxp1; 
       xp[3] = x3_min + m_iter*dxp3; 
       x_of_xp(x_tmp,xp);    
       dx_dxp_calc(x_tmp,xp,dx_dxp);
       DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
       DLOOP1 {     dx[i] = -resid[i]            ; } 
       dxdxp11 = dx_dxp[1][1];
       LU_decompose(dx_dxp,permute);
       LU_substitution( dx_dxp, dx, permute );
 
       DLOOP1 { errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]); }
       fprintf(stdout,"xp1-resid-survey2-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],xp[3],dx[1],dx[3],x_tmp[1],x_tmp[3]);
       fprintf(stdout,"xp1-resid-survey3-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,errx_try[1],errx_try[3],resid[1],resid[3],dxdxp11);
       fflush(stdout);
       //     }
   }

  loc_iter++;
#endif


  /*  Check for bad untrapped divergences : */
  if( finite(df)==0    ) {  
    return(2);  
  }
  if( errx <= NEWT_TOL ) {  
    return(0);  
  }
  if( (errx <= MIN_NEWT_TOL) && (errx > NEWT_TOL) ){   
    return(0);  
  }
  if( errx > MIN_NEWT_TOL){
    return(1);
  }

  return(0);

}

/********************************************************************** 
/************************************************************

  general_newton_raphson_3d_select(): 

    -- special version of general_newton_raphson_3d_general() that allows you 
        to select one dimension to minimize at a time.  

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  X2 from theta and X3 from phi
        by ensuring that   0 <= X2 <= 1  (look near "METHOD specific:")
        by ensuring that   0 <= X3 <= 2*Pi  (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson_3d_select( double *x_phys, double *xp, short int *minimize)
{
  double df, dx[NDIM], x_old[NDIM], x_tmp[NDIM], x_orig[NDIM], resid[NDIM], errx, errx_try[NDIM], dx_dxp[NDIM][NDIM];
  double errx_old;
  int    i, n_iter, i_extra, doing_extra, keep_iterating;

#if( PRINT_NEWTON_RAPHSON_INFO ) 
  double x1_min=1e30;
  double x3_min=1e30;
  double errx_min=1e30;
  double x1_max;
  double x3_max;
  static int loc_iter = 0; 
#endif


#if( COORD_TYPE_CHOICE != COORD_WARPED_SPHERICAL )
  fprintf(stdout,"%s(): Probably should use another routine to perform coordinate inversion!  \n",__func__); 
  fflush(stdout); fail(FAIL_BASIC,0); 
#endif

  // Initialize various parameters and variables:
  errx = 1.e20 ; 
  errx_old = 2*errx;
  i_extra = doing_extra = 0;

  /*   METHOD specific: */
  if( minimize[2] && ( (xp[2] < 0.) || (xp[2] > 1.) ) ) {   xp[2] = ranc(0);  }  /* make sure initial guess is physical */
  if( minimize[3] && ( (xp[3] < 0.) || (xp[3] > 1.) ) ) {   xp[3] = ranc(0);  }  /* make sure initial guess is physical */
  
  DLOOP1 {  x_old[i] = x_orig[i] = xp[i] ; }

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    x_of_xp(x_tmp,xp);
    dx_dxp_calc(x_tmp,xp,dx_dxp);
    DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
    DLOOP1 {     dx[i] = -resid[i]            ; } 
    DLOOP1 { if( !minimize[i] ) {   dx[i] = 0. ; } }
    DLOOP1 {     dx[i] /= dx_dxp[i][i]            ; } 
    
    /* Save old values before calculating the new: */
    DLOOP1 {  x_old[i] = xp[i] ; }

    /* don't use line search : */
    DLOOP1 {  xp[i] += dx[i]  ; } 

    df = 0. ; 
    DLOOP1 {  df += resid[i]*resid[i] ; } 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx_old = errx; 
    errx = -1.e20;
    DLOOP1 { 
      errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]);
      if( errx_try[i] > errx)  { errx = errx_try[i]; } 
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    if( errx_min > errx ) { 
      errx_min = errx;
      x1_min = xp[1]; 
      x3_min = xp[3];
    }
#endif

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    /*   METHOD specific: */
    /*  I can't think of invalid values of xp,  Newton-Raphson should find the solution...  */
    if( xp[2] < 0. ) { 
      xp[2] = ranc(0) * x_old[2];  /* random number between old value and 0. */
    }
    else if( xp[2] > 1. ) { 
      xp[2] = 1. + ranc(0)*(x_old[2]-1.);   /* random number between old value and 1. */
    }

    if( xp[3] < 0. ) { 
      xp[3] = ranc(0) * x_old[3];  /* random number between old value and 0. */
    }
    else if( xp[3] > 1. ) { 
      xp[3] = 1. + ranc(0)*(x_old[3]-1.);   /* random number between old value and 1. */
    }

    
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (errx <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( (errx <= NEWT_TOL) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    fprintf(stdout,"xp1-resid-survey1-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],xp[3],errx_try[1],errx_try[3],resid[1],resid[3],dx[1],dx[3]); 
    fflush(stdout);

    fprintf(stdout,"xp      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,xp[0],xp[1],xp[2],xp[3]); 
    fprintf(stdout,"dx      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,dx[0],dx[1],dx[2],dx[3]); 
    fprintf(stdout,"resid   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,resid[0],resid[1],resid[2],resid[3]); 
    fprintf(stdout,"x_phys  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_phys[0],x_phys[1],x_phys[2],x_phys[3]); 
    fprintf(stdout,"x_tmp   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_tmp[0],x_tmp[1],x_tmp[2],x_tmp[3]); 
    fprintf(stdout,"errx_i  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,errx_try[0],errx_try[1],errx_try[2],errx_try[3]); 
    fprintf(stdout,"df,errx = %10d %26.16e  %26.16e \n",n_iter,df,errx); 
    fflush(stdout);
#endif

    n_iter++;

  }   // END of while(keep_iterating)

#if( PRINT_NEWTON_RAPHSON_INFO ) 
   double eps=5.e-8; 
   x1_max = (1.+eps) * x1_min; 
   x1_min *= (1.-eps);
   x3_max = (1.+eps) * x3_min; 
   x3_min *= (1.-eps);
   int n_survey = 400; 
   double dxp1 = (x1_max - x1_min)/(n_survey-1);
   double dxp3 = (x3_max - x3_min)/(n_survey-1);
   int m_iter; 
   double dxdxp11;
   m_iter = n_survey/2;
   for( n_iter = 0 ; n_iter < n_survey; n_iter++ ) { 
     //     for( m_iter = 0 ; m_iter < n_survey; m_iter++ ) { 
       xp[1] = x1_min + n_iter*dxp1; 
       xp[3] = x3_min + m_iter*dxp3; 
       x_of_xp(x_tmp,xp);    
       dx_dxp_calc(x_tmp,xp,dx_dxp);
       DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
       DLOOP1 {     dx[i] = -resid[i]            ; } 
       dxdxp11 = dx_dxp[1][1];
       LU_decompose(dx_dxp,permute);
       LU_substitution( dx_dxp, dx, permute );
 
       DLOOP1 { errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]); }
       fprintf(stdout,"xp1-resid-survey2-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],xp[3],dx[1],dx[3],x_tmp[1],x_tmp[3]);
       fprintf(stdout,"xp1-resid-survey3-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,errx_try[1],errx_try[3],resid[1],resid[3],dxdxp11);
       fflush(stdout);
       //     }
   }

  loc_iter++;
#endif


  /*  Check for bad untrapped divergences : */
  if( finite(df)==0    ) {  
    return(2);  
  }
  if( errx <= NEWT_TOL ) {  
    return(0);  
  }
  if( (errx <= MIN_NEWT_TOL) && (errx > NEWT_TOL) ){   
    return(0);  
  }
  if( errx > MIN_NEWT_TOL){
    return(1);
  }

  return(0);

}


/********************************************************************** 
/************************************************************

  general_newton_raphson_3d_general(): 

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  X2 from theta and X3 from phi
        by ensuring that   0 <= X2 <= 1  (look near "METHOD specific:")
        by ensuring that   0 <= X3 <= 2*Pi  (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson_3d_general( double *x_phys, double *xp)
{
  double df, dx[NDIM], x_old[NDIM], x_tmp[NDIM], x_orig[NDIM], resid[NDIM], errx, errx_try[NDIM], dx_dxp[NDIM][NDIM];
  double errx_old;
  int    i, n_iter, i_extra, doing_extra, keep_iterating, permute[NDIM];

#if( PRINT_NEWTON_RAPHSON_INFO ) 
  double x1_min=1e30;
  double x3_min=1e30;
  double errx_min=1e30;
  double x1_max;
  double x3_max;
  static int loc_iter = 0; 
#endif


#if( COORD_TYPE_CHOICE != COORD_WARPED_SPHERICAL )
  fprintf(stdout,"%s(): Probably should use another routine to perform coordinate inversion!  \n",__func__); 
  fflush(stdout); fail(FAIL_BASIC,0); 
#endif

  // Initialize various parameters and variables:
  errx = 1.e20 ; 
  errx_old = 2*errx;
  i_extra = doing_extra = 0;

  DLOOP1 {  x_old[i] = x_orig[i] = xp[i] ; }

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    x_of_xp(x_tmp,xp);
    dx_dxp_calc(x_tmp,xp,dx_dxp);
    DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
    DLOOP1 {     dx[i] = -resid[i]            ; } 
    LU_decompose(dx_dxp,permute);     LU_substitution( dx_dxp, dx, permute );
    
    /* Save old values before calculating the new: */
    DLOOP1 {  x_old[i] = xp[i] ; }

    /* don't use line search : */
    DLOOP1 {  xp[i] += dx[i]  ; } 

    df = 0. ; 
    DLOOP1 {  df += resid[i]*resid[i] ; } 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx_old = errx; 
    errx = -1.e20;
    DLOOP1 { 
      errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]);
      if( errx_try[i] > errx)  { errx = errx_try[i]; } 
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    if( errx_min > errx ) { 
      errx_min = errx;
      x1_min = xp[1]; 
      x3_min = xp[3];
    }
#endif

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    /*   METHOD specific: */
    /*  I can't think of invalid values of xp,  Newton-Raphson should find the solution...  */
    if( xp[2] < 0. ) { 
      xp[2] = ranc(0) * x_old[2];  /* random number between old value and 0. */
    }
    else if( xp[2] > 1. ) { 
      xp[2] = 1. + ranc(0)*(x_old[2]-1.);   /* random number between old value and 1. */
    }

    if( xp[3] < 0. ) { 
      xp[3] = ranc(0) * x_old[3];  /* random number between old value and 0. */
    }
    else if( xp[3] > 1. ) { 
      xp[3] = 1. + ranc(0)*(x_old[3]-1.);   /* random number between old value and 1. */
    }

    
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (errx <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( (errx <= NEWT_TOL) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

#if( PRINT_NEWTON_RAPHSON_INFO ) 
    fprintf(stdout,"xp1-resid-survey1-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],xp[3],errx_try[1],errx_try[3],resid[1],resid[3],dx[1],dx[3]); 
    fflush(stdout);

    fprintf(stdout,"xp      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,xp[0],xp[1],xp[2],xp[3]); 
    fprintf(stdout,"dx      = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,dx[0],dx[1],dx[2],dx[3]); 
    fprintf(stdout,"resid   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,resid[0],resid[1],resid[2],resid[3]); 
    fprintf(stdout,"x_phys  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_phys[0],x_phys[1],x_phys[2],x_phys[3]); 
    fprintf(stdout,"x_tmp   = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,x_tmp[0],x_tmp[1],x_tmp[2],x_tmp[3]); 
    fprintf(stdout,"errx_i  = %10d %26.16e %26.16e %26.16e %26.16e\n",n_iter,errx_try[0],errx_try[1],errx_try[2],errx_try[3]); 
    fprintf(stdout,"df,errx = %10d %26.16e  %26.16e \n",n_iter,df,errx); 
    fflush(stdout);
#endif

    n_iter++;

  }   // END of while(keep_iterating)

#if( PRINT_NEWTON_RAPHSON_INFO ) 
   double eps=5.e-8; 
   x1_max = (1.+eps) * x1_min; 
   x1_min *= (1.-eps);
   x3_max = (1.+eps) * x3_min; 
   x3_min *= (1.-eps);
   int n_survey = 400; 
   double dxp1 = (x1_max - x1_min)/(n_survey-1);
   double dxp3 = (x3_max - x3_min)/(n_survey-1);
   int m_iter; 
   double dxdxp11;
   m_iter = n_survey/2;
   for( n_iter = 0 ; n_iter < n_survey; n_iter++ ) { 
     //     for( m_iter = 0 ; m_iter < n_survey; m_iter++ ) { 
       xp[1] = x1_min + n_iter*dxp1; 
       xp[3] = x3_min + m_iter*dxp3; 
       x_of_xp(x_tmp,xp);    
       dx_dxp_calc(x_tmp,xp,dx_dxp);
       DLOOP1 {  resid[i] = x_tmp[i] - x_phys[i] ; } 
       DLOOP1 {     dx[i] = -resid[i]            ; } 
       dxdxp11 = dx_dxp[1][1];
       LU_decompose(dx_dxp,permute);
       LU_substitution( dx_dxp, dx, permute );
 
       DLOOP1 { errx_try[i]  = (xp[i]==0.) ?  fabs(dx[i]) : fabs(dx[i]/xp[i]); }
       fprintf(stdout,"xp1-resid-survey2-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,xp[1],xp[3],dx[1],dx[3],x_tmp[1],x_tmp[3]);
       fprintf(stdout,"xp1-resid-survey3-%010d %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", loc_iter,errx_try[1],errx_try[3],resid[1],resid[3],dxdxp11);
       fflush(stdout);
       //     }
   }

  loc_iter++;
#endif


  /*  Check for bad untrapped divergences : */
  if( finite(df)==0    ) {  
    return(2);  
  }
  if( errx <= NEWT_TOL ) {  
    return(0);  
  }
  if( (errx <= MIN_NEWT_TOL) && (errx > NEWT_TOL) ){   
    return(0);  
  }
  if( errx > MIN_NEWT_TOL){
    return(1);
  }

  return(0);

}


#undef FTYPE
#undef NEWT_DIM
#undef MAX_NEWT_ITER
#undef NEWT_TOL
#undef MIN_NEWT_TOL
#undef EXTRA_NEWT_ITER 
#undef GNR_SMALL 

#undef PRINT_DEBUG_INFO
#undef PRINT_NEWTON_RAPHSON_INFO
#undef PERFORM_XP_OF_X_TEST
#undef MISSING_VALUE_RHO
#undef USE_EXTRAPOLATION
#undef phys_val  
#undef func_resid 


#else


#endif
