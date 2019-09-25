
#include "decs.h"


static double rho_start, beta_start, K_adiab, B_mag;


/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

    Initial data consists of a uniform density/pressure/magnetic-field embedded in a 
    post-Newtonian approximation for a binary black hole spacetime. 

    This is done in Cartesian coordinates.  

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

  if( n_cells_glob > ( UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Num. of global points too large !!! N = %ld \n", 
	    n_cells_glob);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

#if( (METRIC_TYPE_CHOICE != METRIC_GENERAL_DYNAMIC) || (TOP_TYPE_CHOICE != TOP_CARTESIAN) || (COORD_TYPE_CHOICE != COORD_IDENTITY) )  
  //  wont-compile10
#endif

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
	gam = (5./3.) ;
	cour = 0.8;

	/* Coordinate dependent quantities : */ 
	set_special_coord();

	t = 0. ;                     /* Initial time */ 

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 10000. ;          /* Length of X0 dimension (evolution period) */ 
//--both	GridLength[1] = 1.5*initial_bbh_separation  ;   /* Length of X1 dimension */ 
//--both	GridLength[2] = 1.5*initial_bbh_separation  ;   /* Length of X2 dimension */ 
//--both	GridLength[3] = 1.5*initial_bbh_separation  ;   /* Length of X3 dimension */ 
	GridLength[1] = 0.25*initial_bbh_separation  ;   /* Length of X1 dimension */ 
	GridLength[2] = 0.25*initial_bbh_separation  ;   /* Length of X2 dimension */ 
	GridLength[3] = 0.25*initial_bbh_separation  ;   /* Length of X3 dimension */ 

	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];

	/**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
	**************************************************************************/
	/* Set the Physical Minimum boundary  */
	startx[0] =  0.;                                              
//--both 	SDLOOP1 startx[i] = -0.5*GridLength[i];
	SDLOOP1 startx[i] = -0.5*GridLength[i];
	startx[1] += 0.5*initial_bbh_separation;

	/* Adjust to numerical minimum boundary */
	SDLOOP1 startx[i] -= NG*dx[i];  

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = 1.e-5 ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = 1.;         /* radiative flux dumps       */
	n_restart = 1000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	/* Need to set dx[0] before geom is calculated as we calculate the substep and 
	   future step spacetimes now */
	dx[0] = 1.e-3;    
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

	fprintf(stdout,"rhomin,minlimt = %28.18e %28.18e \n", RHOMIN, RHOMINLIMIT);
	fprintf(stdout,"uumin ,minlimt = %28.18e %28.18e \n", UUMIN,  UUMINLIMIT);
	
	fflush(stdout);


	/* in case you want to test the metric and coord. transf. routines */
	//	test_geom();   

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
  double X[NDIM];
  double r;
  
  void init_prim( int i, int j, int k, double *xp, double *pr);

  /*******************************************************************************
     Set the hydrodynamic quantities (rho,uu,v^i) : 
  *******************************************************************************/
  //  ranc(myid+1);

  LOOP {
    coord(i,j,k,CENT,X) ; 
    init_prim( i, j, k, X, p[i][j][k]);
  }

  /*******************************************************************************
    Set the magnetic field :  
  *******************************************************************************/
  /* Correct bad points and setup boundary values since we will require them for B^i */
  fixup(p) ;
  bounds(p,0) ;

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
  initial_bbh_separation = 20.;   /* Initial separtion of BHs in units of total mass */
  m_bh1 = 0.5;                    /* Mass of left-most BH */
  m_bh2 = 0.5;                    /* Mass of right-most BH */
  m_bh_tot = m_bh1 + m_bh2;       /* Total mass of BHs     */

  r_horizon1 = m_bh1;             /* Radius of the horizon about BH1 */
  r_horizon2 = m_bh2;             /* Radius of the horizon about BH2 */
  
  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;

  /* MHD parameters : */
  rho_start       = 1.;
  K_adiab         = 1.e-2;
  beta_start      = 1.e-2;
  B_mag           = sqrt( 2. * K_adiab * pow(rho_start,gam) / beta_start );


  fprintf(stdout,"init():   initial_bbh_separation = %28.18e \n", initial_bbh_separation);
  fprintf(stdout,"init():   m_bh1                  = %28.18e \n", m_bh1                 );
  fprintf(stdout,"init():   m_bh2                  = %28.18e \n", m_bh2                 );
  fprintf(stdout,"init():   m_bh_tot               = %28.18e \n", m_bh_tot              );
  fprintf(stdout,"init():   rho_start              = %28.18e \n", rho_start             );
  fprintf(stdout,"init():   K_adiab                = %28.18e \n", K_adiab               );
  fprintf(stdout,"init():   beta_start             = %28.18e \n", beta_start            );
  fprintf(stdout,"init():   B_mag                  = %28.18e \n", B_mag                 );
  fflush(stdout);

  return;
}


/***********************************************************************************/
/***********************************************************************************
  init_prim(): 
  ------------------
   -- This is supposed to follow the routine in inits.f of the Hawley et al. code; 
       -- inconsistencies are due to not knowing what various grid functions represent 
          in inits.f  (there is very little documentation in that file); 
   -- when in doubt, I used   
             De Villiers, Hawley, Krolik  ApJ 599 (2003)

   -- The disk parameters mentioned in this paper are used below; 
   -- Uses BL coordinates until the very end and then transforms the velocities to KS; 
***********************************************************************************/
void init_prim(int i, int j, int k, double *xp, double *pr)
{

  pr[RHO] = rho_start  ;	
  pr[UU ] = K_adiab * pow(pr[RHO],gam) / (gam-1.)  ;	
  pr[U1 ] = 0.  ;	
  pr[U2 ] = 0.  ;	
  pr[U3 ] = 0.  ;	
  pr[B1 ] = 0.  ;	
  pr[B2 ] = 0.  ;	
  pr[B3 ] = 0.  ;	

#if( !HYDRO_ONLY )
  pr[B3 ] = B_mag;      /* Uniform in z */
#endif

  return;

}

