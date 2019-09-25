
#include "decs.h"

// set to 1 to see the Bondi parameters
#define LTRACE   1
#define LTRACE2  0
#define LTRACE3  0   /* for find_bondi_solution */

//Newton-Raphson parameters:
#define NEWT_DIM_B        (1      )
#define MAX_NEWT_ITER_B   (30     )    /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL_B        (1.0e-15)    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL_B    (1.0e-10)    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER_B (2      )
#define SMALL_BONDI (1.e-20)

double Mdot, rs, vs_sq, vs, cs_sq, cs, rhos, hs, K, Qdot, gamma_eos, r_sol;

static void bondi_resid(double x[], double dx[], double resid[], double jac[][NEWT_DIM_B], 
			double *f, double *df, int n );

static int gnr_bondi( double x[], int n, 
		      void (*funcd) (double [], double [], double [], 
				     double [][NEWT_DIM_B], double *, 
				     double *, int) );

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

    Bondi solution using Bondi solution routines written by Scott Noble

***********************************************************************************/

void init(void)
{
  double ftmp;
  void init_data(void) ;

  /* Do some checks on the pre-compiler definitions : */ 
  ftmp = (double) N1TOT;  ftmp *= N2TOT  ;  ftmp *= N3TOT  ;  ftmp *= NP; 
  if( ftmp > ( (double) UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Max. no. of elements too large !!! N = %g \n", ftmp);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  init_base();  /* Set global grid parameters and static work arrays */

  init_data();  /* Set MHD grid functions */

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

	/* Set the parameters that define the grid and other constants : */ 
	gam = (4./3.) ;
	cour = 0.4 ;

	t = 0. ;                     /* Initial time */ 

	/* Coordinate dependent quantities : */ 
	set_special_coord();

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 100. ;         /* Length of X0 dimension (evolution period) */ 
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	GridLength[1] = Rout-Rin  ;   /* Length of X1 dimension */ 
	GridLength[2] = th_length ;   /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	GridLength[1] = log((coord_params->upsilon_r-coord_params->f0_r)/(1.-coord_params->f0_r)) ; /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
#else
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                     /* Length of X2 dimension */ 
#endif
	GridLength[3] = 2.*M_PI/4. ;                 /* Length of X3 dimension */ 
#if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[1] = 1.  ;              /* Length of X1 dimension */
        GridLength[2] = 1.  ;                    /* Length of X2 dimension */
	GridLength[3] = 1. ;                 /* Length of X3 dimension */ 
#endif


	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];
#if(EQUATORIAL_RUN)
 	dx[2] = 1.e-20;
#endif

	/**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
	**************************************************************************/
	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	startx[1] = Rin;                /* Set the Physical Minimum boundary  */
	startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	startx[1] = 0.;        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
#else 
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
#endif
#if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )

	startx[1] = 1.e-10;        /* Set the Physical Minimum boundary  */

# if(EQUATORIAL_RUN)
	startx[2] = 0.5-0.5*dx[2] ; 
# else 
 	startx[2] = 0.            ; 
# endif 

#endif
	startx[3] =  0.;               /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	DT_out[OUT_ASCII]   = GridLength[0]/10. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/40. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = GridLength[0]/200. ;      /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	n_restart = 100;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	dx[0] = 1.e-5;
	calc_all_geom() ;  

	/**************************************************************************
	  Print out basic grid parameters: 
	**************************************************************************/
	fprintf(stdout,"Length: %10.4e %10.4e %10.4e %10.4e\n", 
		GridLength[0], GridLength[1], GridLength[2], GridLength[3]);
	fprintf(stdout,"dx    : %10.4e %10.4e %10.4e %10.4e\n", 
		dx[0], dx[1], dx[2], dx[3]);
	fprintf(stdout,"startx: %10.4e %10.4e %10.4e %10.4e\n", 
		startx[0], startx[1], startx[2], startx[3]);
	fflush(stdout);

	/* in case you want to test the metric and coord. transf. routines */
	//test_geom();   

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
  int i, j, k, l, g ; 
  double X[NDIM], prim[NP];
  
  int init_dsandvels(int i, int j, int k, double *xp, double *pr);


  /* Komissarov's explosion problem */
  N2ALL_LOOP N3ALL_LOOP N1ALL_LOOP {
    coord(i,j,k,CENT,X) ; 
    
    init_dsandvels( i, j, k, X, prim);

    PLOOP p[i][j][k][l] = prim[l] ; 

  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

  /* Only need to use reflective boundaries for U2 : 
     this is silly since U2 is 0 from the start. */
//   GLOOP   N1ALL_LOOP   N3ALL_LOOP  p[i][N2S-1-g][k][U2] *= -1.;
//   GLOOP   N1ALL_LOOP   N3ALL_LOOP  p[i][N2E+1+g][k][U2] *= -1.;

  /* if we are not setting boundary conditions ever, we need to set ph[] too */ 
  ALL_LOOP PLOOP ph[i][j][k][l] = p[i][j][k][l];

}


/***********************************************************************/
/***********************************************************************
  set_special_coord(): 
  --------------------
   -- set coordinate specific (global) quantities;
***********************************************************************/
void set_special_coord(void)
{
  
  a = 0.;          /* Spin of the black hole in units of M */ 

  R0       = 0.;          /* Offset in Radius from 1.            */
  Rout     = 100.;         /* Radial extent of the grid           */
  n_within_horizon = 5;   /* number of cells within horizon */

#if(COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD||COORD_TYPE_CHOICE==COORD_DIAGONAL3)
  h_slope  = 1.;      
#else
  h_slope  = 0.;      /* Severity of the focusing            */
#endif

  th_cutout = SMALL * M_PI;       /* Angular size to excise from the axis */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */

  X1_slope = 1.;       /* Severity of transition              */
  X1_0     = log(10.);     /* Location of transition in X1 units  */ 
  
  
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  Rin       = Rin_calc();

  diag3_exponent = 3;                            /* note that this should be an odd integer greater than 1  */
  diag3_factor   = 1. - 2*th_cutout/M_PI - h_slope;  /* Temporary variable */
  if( (diag3_exponent < 3) || (diag3_exponent % 2 == 0 )  ) { fprintf(stderr,"set_special_coord(): bad value for  diag3_exponent  = %d \n",diag3_exponent); fflush(stderr); exit(4112); }


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


/***********************************************************************/
/***********************************************************************
        SPECIAL BONDI FUNCTIONS 
***********************************************************************/

/***********************************************************************************
 ***********************************************************************************/
int init_dsandvels(int i, int j, int k, double *xp, double *pr)
{

  int iit; 

  int set_bondi_parameters( double M_in, double Mdot_in,  double rs_in, double gam );
  int find_bondi_solution( double r, double *rho, double *u, double *v ) ;
  
  int retval ;

  double Mbh, Mdot, r_sonicpt, rhotmp;
  static double *r_bondi, *rho_bondi, *u_bondi, *v_bondi;

  double x[NDIM], r,th, d_spt, d_spt_min;
  
  static int  first_one  = 1;
  static int  ir_spt ; 
  static int j_old = -1;
  static int k_old = -1;
  int calc_bondi;
  

  double gcov[NDIM][NDIM], ucon_ks[NDIM], ucon_bl[NDIM];
  int gp_bak, nall;

  
  /********************************************************
    Do the following only once: 
       -- set parameters, allocate temp. arrays 
       -- calc. solution vs. r
   ********************************************************/
  // set Parameters that determine the Bondi solution:

#if( COORD_RADIUS_DIM == 1 ) 
  calc_bondi = first_one ;
#else
  calc_bondi = (j_old != j) || (k_old != k) ;
#endif

  if( calc_bondi ) { 
    fprintf(stdout,"calculating bondi... \n");

    if(   r_bondi != NULL ) { FREE(   r_bondi); }
    if( rho_bondi != NULL ) { FREE( rho_bondi); }
    if(   u_bondi != NULL ) { FREE(   u_bondi); }
    if(   v_bondi != NULL ) { FREE(   v_bondi); }
    
    first_one = 0; 
    Mbh       = 1.;
    Mdot      = 4*M_PI;
    r_sonicpt = 8.;

    nall = totalsize[1]+2*NG;
    if( (retval=set_bondi_parameters( Mbh, Mdot, r_sonicpt, gam)) ) { 
      fprintf(stdout,"init_dsandvels: Problem with set_bondi(), ret = %d\n", retval);
      return( 1 ) ; 
    }

    // Allocate temp. arrays to save time since Bondi is spherically symmetric
    if( (r_bondi = (double *) calloc(nall, sizeof(double))) == NULL ) { 
      fprintf(stdout,"init_dsandvels: Cannot alloc r_bondi \n");
      return(1);
    }
    if( (rho_bondi = (double *) calloc(nall, sizeof(double))) == NULL ) { 
      fprintf(stdout,"init_dsandvels: Cannot alloc rho_bondi \n");
      free(r_bondi);
      return(1);
    }
    if( (u_bondi = (double *) calloc(nall, sizeof(double))) == NULL ) { 
      fprintf(stdout,"init_dsandvels: Cannot alloc u_bondi \n");
      free(r_bondi);      free(rho_bondi);
      return(1);
    }
    if( (v_bondi = (double *) calloc(nall, sizeof(double))) == NULL ) { 
      fprintf(stdout,"init_dsandvels: Cannot alloc v_bondi \n");
      free(r_bondi);      free(rho_bondi);      free(u_bondi);
      return(1);
    }

    d_spt_min = r_sonicpt;
    ir_spt = 0;

    gp_bak = globalpos[1];   
    globalpos[1] = 0;   /* So that we can loop over all global radii */
    for( iit = 0; iit < nall; iit++ ) { 
      coord(iit, j, k, CENT, xp);
      x_of_xp(x , xp);
      r_bondi[iit] = x[RR]; 
      d_spt = fabs(r_bondi[iit] - r_sonicpt);
      if( d_spt < d_spt_min ) {  d_spt_min = d_spt ; ir_spt = iit;  }
    }

    globalpos[1] = gp_bak;

    // get Bondi solution for all r:
    rhotmp = -1.;  // start with guess
    for( iit = (ir_spt+1); iit < nall; iit++ ) { 
      if( (retval=find_bondi_solution( r_bondi[iit], &rhotmp, &(u_bondi[iit]), &(v_bondi[iit]) )) ) { 
	fprintf(stdout,"init_dsandvels: Problem with find_bondi() at ir = %d,  r = %g \n",iit,r_bondi[iit]);
	return(1);
      }
      rho_bondi[iit] = rhotmp;
    }

    rhotmp = -1.;  // start with guess
    for( iit = ir_spt; iit >= 0; iit-- ) { 
      if( (retval=find_bondi_solution( r_bondi[iit], &rhotmp, &(u_bondi[iit]), &(v_bondi[iit]) )) ) { 
	fprintf(stdout,"init_dsandvels: Problem with find_bondi() at ir = %d,  r = %g \n",iit,r_bondi[iit]);
	return(1);
      }
      rho_bondi[iit] = rhotmp;
    }

//     if( myid == printer_pid ) { 
//       for( iit = 0 ; iit < nall ; iit++) { 
// 	fprintf(stdout,"bondsol: %28.18e  %28.18e  %28.18e  %28.18e \n", 
// 		r_bondi[iit],rho_bondi[iit],u_bondi[iit], v_bondi[iit]); 
// 	fflush(stdout);
//       }
//     }    
  }  // end if(first_one)


  j_old = j;
  k_old = k;

  /********************************************************
    Do the following for each call or point:
       -- set the primitive variables for the right r
   ********************************************************/

  pr[RHO] =  rho_bondi[i+globalpos[1]];
  pr[UU]  =  u_bondi[i+globalpos[1]];

  x_of_xp( x, xp ); 

  ucon_bl[TT] = ucon_bl[TH] = ucon_bl[PH] = 0.;
  ucon_bl[RR] = -v_bondi[i+globalpos[1]];

  /* Find time component of 4-velocity: */
  bl_gcov_func( x, gcov );      
  setutcon( ucon_bl, gcov );

  //  fprintf(stdout,"bondi2: %10.4e  %10.4e  %10.4e  \n",  r, ucon_bl[TT], ucon_bl[RR]) ;

  /* Transform into computation coordinates (need to first do BL->KS): */
  bl_to_ks_con( x, ucon_bl, ucon_ks ); 

  //  fprintf(stdout,"bondi3: %10.4e  %10.4e  %10.4e  \n",  r, ucon_ks[TT], ucon_ks[RR]) ;

  transform_rank1con( x, xp, ucon_ks ); 

  //  fprintf(stdout,"bondi4: %10.4e  %10.4e  %10.4e  \n",  r, ucon_ks[TT], ucon_ks[RR]) ;

  /* get metric in numerical coordinates */
  ks_gcon_func( x, gcov );      /* gcov is now the contravariant metric */
  transform_rank2con( x, xp, gcov );      /* gcov is now the contravariant metric */

  /* Calculate the primitive variables : */
  ucon2pr( pr, ucon_ks, gcov ); 
  pr[ B1] = 0. ;
  pr[ B2] = 0. ;
  pr[ B3] = 0. ;

  // If we are the last point then free up the arrays like a good little C program: 
  if( (i == (N1TOT-1)) && (j == (N2TOT-1)) && (k == (N3TOT-1)) ) { 
    free(r_bondi);      free(rho_bondi);      free(u_bondi);    free(v_bondi);
  }

  return(0); // need to transform to Kerr-Schild

}



/***************************************************************************

set_bondi_parameters():

   -- finds the values of the hydro. quantities  at the sonic point, which
         serves as a reference point for the conservation equations given 
         in Shapiro and Teukolsky equations (G.21,G.22).  

   -- The sonic point values are then used in find_bondi_solution() to determine
       the hydro. quantities at an arbitrary radius;

   -- the "boundary conditions" that uniquely determine the Bondi solution
       are the radius of the sonic point, "rs",  and the mass accretion rate, 
       Mdot;

   -- the Bondi solution here in is the isentropic, spherically symmetric, 
       perfect fluid solution to Einstein's equations.  That is, we only 
       assume an r-dependence, there's a in-going radial velocity only, 
       and the EOS are :  P = (G-1)*rho   and   P = K rho^G  
       where  K = const.  and  G is the adiabatic constant "gam".  

***************************************************************************/
int set_bondi_parameters( double M_in, double Mdot_in,  double rs_in, double gam )
{

  int retval=0;
  double my_pi, checkp; 
  double  gtemp; 

  /* Set the solution-determining parameters:    */
  M    = M_in;
  Mdot = Mdot_in;
  rs   = rs_in;
  gamma_eos = gam;

  
  /* Calculate the hydro. quantities: */
  vs_sq  =  M / ( 2. * rs ) ; 
  vs     =  sqrt(vs_sq); 

  cs_sq  =  M / ( 2.*rs - 3.*M ) ; 
  cs     =  sqrt(cs_sq);

  rhos   =  Mdot / ( 4. * M_PI * vs * rs * rs ) ; 

  gtemp  =  gam - 1.;
  hs     =  1. / ( 1. - cs_sq / (gam - 1.) );

  K      = hs * cs_sq * pow( rhos, (-gtemp) ) / gam ; 

  Qdot   = hs * hs * ( 1. - 3. * vs_sq ) ;

  gamma_eos = gam;

  if( cs_sq > (gam - 1.) ) { retval = -1 ; } 

#if( LTRACE ) 
	fprintf(stdout,"rs   = : %28.20e \n", rs) ;
	fprintf(stdout,"urs  = : %28.20e \n", vs) ;
	fprintf(stdout,"rhos = : %28.20e \n", rhos) ;
	fprintf(stdout,"K    = : %28.20e \n", K) ;
#endif

#if( LTRACE2 ) 
  fprintf(stdout, "\n#######################################################\n");
  fprintf(stdout, "Bondi Solution Parameters1: \n");
  fprintf(stdout, "------------------------- \n\n");
  fprintf(stdout, "M  = %28.20e     Mdot = %28.20e     rs   = %28.20e  \n",M,Mdot,rs);
  fprintf(stdout, "vs = %28.20e     cs   = %28.20e     rhos = %28.20e  \n",vs,cs,rhos);
  fprintf(stdout, "hs = %28.20e     K    = %28.20e     Qdot = %28.20e   \n",hs,K,Qdot);
  fprintf(stdout, "gam= %28.20e     r_sol= %28.20e       \n",gamma_eos, r_sol);
  fprintf(stdout, "#######################################################\n\n");
#endif

  return(retval);
  
}


/***************************************************************************/
/***************************************************************************

find_bondi_solution():

   -- essentially just calls gnr_bondi() to find the solution for
      the density at a given radius, given the parameters calculated from 
      set_bondi_parameters();

   -- after the density is found, the density (rho), internal energy densit (u), magnitude 
        of the radial component of the 4-velocity (v) are returned to the calling 
        routine; 

   -- requires r = radius at which we want solution

   -- note that v is a magnitude, so the user will have to set u^r = -v ;

   -- if there is an error in finding the solution, it returns the error 
       status from the root-finding routine.  See documentation of the 
       gnr_bondi() for further details;

***************************************************************************/
int find_bondi_solution( double r, double *rho, double *u, double *v )
{

  int  retval=0;
  int  ntries = 10000;
  int  itry;

  double rhotmp, rho_guess;
  double dr;

  
  /************************************************************************/
  /* Find the initial guess for the newton iterations:                    */
  /* Take the sonic point values if we have no better guess (when rho<0)  */
  /************************************************************************/

  rhotmp = (sqrt(Qdot) - 1.) * (gamma_eos - 1.) / ( gamma_eos * K );
  rho_guess = pow( rhotmp , (1./(gamma_eos - 1.)) );

  if( *rho < 0. ) {  
    if( r > 0.9*rs && r < 1.1*rs ) { 
      *rho = rhos;
    }
    else { 
      *rho = rho_guess; 
    }
  } 

  // set global variables needed by residual function:
  r_sol = r ; 


  // Use Newton's method to find rho:
  retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     


  // first try guess if failure 
  if( retval ) { 
    *rho = rho_guess;
    retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     
  }

  // If we were unsure about the guess and solver fails, then creep from known solution to desired point:
  if( retval ) { 

    dr = (r - rs)/(1.*(ntries-1));
    
    *rho = rhos;   // start with sonic point value and near sonic point

    // go gradually away from sonic point toward location where we want the solution:
    r_sol = rs ;
    for( itry = 1; itry < ntries; itry++ ) { 
      r_sol  += dr;

      retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     

#if(LTRACE3) 
      fprintf(stderr,"find_bondi: rsol, rho, rhos = %28.20e %28.20e %28.20e \n", r_sol, *rho, rhos);
      fflush(stderr);
#endif    

      if( retval ) { 
	fprintf(stderr,"find_bondi: Incr. guess failed, decrease dfactor, retval = %d, r = %g, itry = %d \n", retval, r, itry);
	fflush(stderr);
	return(10);
      }
    }
    
    // No try where we want the solution:
    r_sol = r ; 
    retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid );  

    if( retval ) { 
      fprintf(stderr,"find_bondi: Final Incr. guess failed, decrease dfactor??, retval = %d \n", retval);
      fflush(stderr);
      return(11);
    }

  }

  rhotmp = *rho;

  // Calculate other quantities:
  *u = K * pow( rhotmp, gamma_eos ) / (gamma_eos - 1.);

  *v = Mdot / ( 4. * M_PI * r * r * rhotmp );

  if( *u <= 0. ) { retval = -1; } 

#if( LTRACE2 ) 
  fprintf(stdout, "\n#######################################################\n");
  fprintf(stdout, "Bondi Solution Parameters2: \n");
  fprintf(stdout, "------------------------- \n\n");
  fprintf(stdout, "M  = %28.20e     Mdot = %28.20e     rs   = %28.20e  \n",M,Mdot,rs);
  fprintf(stdout, "vs = %28.20e     cs   = %28.20e     rhos = %28.20e  \n",vs,cs,rhos);
  fprintf(stdout, "hs = %28.20e     K    = %28.20e     Qdot = %28.20e   \n",hs,K,Qdot);
  fprintf(stdout, "gam= %28.20e     r_sol= %28.20e       \n",gamma_eos, r_sol);
  fprintf(stdout, "#######################################################\n\n");
#endif
  
  return( retval ) ;

}  


/******************************************************************************

bondi_resid():

     -- routine to calculate the residual and jacobian used by 
           the Newton-Raphson routine general_newton_raphson(), which is 
           used to find X2 from theta;

***********************************************************************************/
static void bondi_resid(double x[], double dx[], double resid[], double jac[][NEWT_DIM_B], 
               double *f, double *df, int n )
{
  double v, vp, h, hp, term;


  hp  =  K * gamma_eos * pow( x[0], (gamma_eos - 2.) );   //   dh/drho 
	 
  h   =  1. +  hp * x[0] / ( gamma_eos - 1. ); 
  	 
  v   =  Mdot / ( 4. * M_PI * r_sol * r_sol * x[0] );    
  	 
  vp  =  -v / x[0];   //  dv/drho


  term = 1. - 2.*M/r_sol + v*v;

  resid[0]  =  -Qdot  +  h * h * term;
  
  jac[0][0] =  2. * h *( hp*term + h*v*vp );

  dx[0] = -resid[0] / jac[0][0];
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);


#if( LTRACE2 ) 
  fprintf(stdout, "\n#######################################################\n");
  fprintf(stdout, "Bondi Solution Parameters3: \n");
  fprintf(stdout, "------------------------- \n\n");
  fprintf(stdout, "M  = %28.20e     Mdot = %28.20e     rs   = %28.20e  \n",M,Mdot,rs);
  fprintf(stdout, "vs = %28.20e     cs   = %28.20e     rhos = %28.20e  \n",vs,cs,rhos);
  fprintf(stdout, "hs = %28.20e     K    = %28.20e     Qdot = %28.20e   \n",hs,K,Qdot);
  fprintf(stdout, "gam= %28.20e     r_sol= %28.20e       \n",gamma_eos, r_sol);
  fprintf(stdout, "#######################################################\n\n");
#endif

#if( LTRACE3 ) 
  //fprintf(stderr,"hp,  h, v,   vp = %28.20e %28.20e %28.20e %28.20e\n", hp, h, v, vp);
  fprintf(stderr,"x, res, jac, dx = %28.20e %28.20e %28.20e %28.20e\n\n", x[0], resid[0],jac[0][0], dx[0]);
#endif 

  return;

}

/**********************************************************************/
/************************************************************

  gnr_bondi():

    -- should be just like the routine general_newton_raphson() in utoprim*.c 
       except the "physicality" condition is different;

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  rho  second Bondi Conservation eq. 
        by ensuring that  rho > 0    (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int gnr_bondi( double x[], int n, 
	       void (*funcd) (double [], double [], double [], 
			      double [][NEWT_DIM_B], double *, 
			      double *, int) )
{
  double f, df, dx[NEWT_DIM_B], x_old[NEWT_DIM_B], resid[NEWT_DIM_B], 
    jac[NEWT_DIM_B][NEWT_DIM_B];
  double errx, x_orig[NEWT_DIM_B];
  int    n_iter, id, jd, i_extra, doing_extra;

  int   keep_iterating, i_increase;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* don't use line search : */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at relative error in indep. variable: */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = (x[0] == 0.) ? (SMALL_BONDI)  :  fabs(x[0]);
    

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL_B) && (doing_extra == 0) && (EXTRA_NEWT_ITER_B > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL_B)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER_B) || (n_iter >= (MAX_NEWT_ITER_B-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0)  ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL_B){
    fprintf(stderr,"newt: errx = %28.20e \n", errx); fflush(stderr);
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL_B) && (fabs(errx) > NEWT_TOL_B) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL_B ){
    return(0);
  }

  return(0);

}


#undef LTRACE 
#undef LTRACE2
#undef LTRACE3
#undef NEWT_DIM_B        
#undef MAX_NEWT_ITER_B   
#undef NEWT_TOL_B        
#undef MIN_NEWT_TOL_B    
#undef EXTRA_NEWT_ITER_B 
#undef SMALL_BONDI 


