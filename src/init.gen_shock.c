
#include "decs.h"
#include "metric.h"



/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

    General slab-symmetric hydrodynamic shock in generalized coordinates;
      -- shock is specified in Cartesian coordinates but then is transformed 
         back to xp coordinates;

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
  double X[NDIM];

  void calc_all_geom(void) ;
  void set_special_coord(void);

  /* Set the parameters that define the grid and other constants : */ 
  gam = (5./3.) ;
  cour = 0.45 ;

  t = 0.0 ;                   /* Initial time */ 



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
  GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
  GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
  GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
  //	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
  GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  GridLength[1] = log((Rout-R0)/(Rin-R0)) ;              /* Length of X1 dimension */ 
  GridLength[2] = xi_diag2[n_diag2_lines]-xi_diag2[0]  ; /* Length of X2 dimension */ 
  GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
  GridLength[1] = log((coord_params->upsilon_r-coord_params->f0_r)/(1.-coord_params->f0_r)) ;              /* Length of X1 dimension */ 
  GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
  //	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
  GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
# elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
  GridLength[1] = 1.  ;             /* Length of X1 dimension */
  GridLength[2] = 1.  ;             /* Length of X2 dimension */
  GridLength[3] = 1.  ;             /* Length of X3 dimension */ 
# endif

  /*************** CARTESIAN  *******************************************/
#elif( TOP_TYPE_CHOICE == TOP_CARTESIAN )
# if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  GridLength[1] = 40. ;   /* Length of X1 dimension */ 
  GridLength[2] = 40. ;   /* Length of X1 dimension */ 
  GridLength[3] = 40. ;   /* Length of X1 dimension */ 
	
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
  dx[2] = 1./160;
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
  dx[0] = 1.e-3;
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

  

  return;

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
  M = m_bh_tot = 1.;

  a = 0.; 
  asq = a*a;
  r_isco    = 0.;
  r_horizon = 0.;
  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  Rout     = 20.;         /* Radial extent of the grid           */
  n_within_horizon = 0;   /* number of cells within horizon */
  //  Rin       = Rin_calc();
  Rin = 1.e-13 ; 
  
  r_horizon1 = r_horizon2 = r_horizon;
  m_bh1 = m_bh2 = M;
  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;
  initial_bbh_separation = 0.;

  h_slope   = 0.13;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.0* M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
  th_cutout = 0.;     
                                  /* 0.02655 gives a cutout of ~0.045Pi for h_slope = 0.35 */

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


/***********************************************************************************/
/***********************************************************************************
  init_data(): 
  ---------
   -- calculates the initial distribution of the MHD fields;
   -- driver routine for various initial data prescriptions;

***********************************************************************************/
void init_data(void) 
{
  int i, j, k ; 
  struct of_coord *coords;
 

  /* Loop over the grid  */
  LOOP {
    
    get_coord(i,j,k,CENT,ncurr,coords) ; 

    if( coords->xcart[1] < 10. ) { 
      p[i][j][k][RHO] = 1.;
      p[i][j][k][UU]  = 1./(gam - 1.) ;
      p[i][j][k][U1]  = 0. ;
      p[i][j][k][U2]  = 0. ;
      p[i][j][k][U3]  = 0. ;
      p[i][j][k][B1]  = 0. ;
      p[i][j][k][B2]  = 0. ;
      p[i][j][k][B3]  = 0. ;
    }
    else { 
      p[i][j][k][RHO] = 2.;
      p[i][j][k][UU]  = 2./(gam - 1.) ;
      p[i][j][k][U1]  = 0. ;
      p[i][j][k][U2]  = 0. ;
      p[i][j][k][U3]  = 0. ;
      p[i][j][k][B1]  = 0. ;
      p[i][j][k][B2]  = 0. ;
      p[i][j][k][B3]  = 0. ;
    }
  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

}

