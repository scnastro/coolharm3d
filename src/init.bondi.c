
#include "decs.h"




/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

    Bondi solution using Bondi solution routines written by Gammie, McKinney & Toth

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

#if( (COORD_TYPE_CHOICE != COORD_MIXED)  && (COORD_TYPE_CHOICE != COORD_DIAGONAL ) )
  fprintf(stderr,"init_base(): Wrong COORD_TYPE_CHOICE =  %d \n", COORD_TYPE_CHOICE );
  fflush(stderr);
  fail(FAIL_BASIC,0);
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
	gam = (4./3.) ;
	cour = 0.8;

	/* Coordinate dependent quantities : */ 
	set_special_coord();

	t = 0. ;                     /* Initial time */ 

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 10000. ;         /* Length of X0 dimension (evolution period) */ 
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	GridLength[1] = Rout-Rin  ;   /* Length of X1 dimension */ 
	GridLength[2] = th_length ;   /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	//	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ;              /* Length of X1 dimension */ 
	GridLength[2] = xi_diag2[n_diag2_lines]-xi_diag2[0]  ; /* Length of X2 dimension */ 
#endif
	GridLength[3] = M_PI/2. ;                 /* Length of X3 dimension */ 

	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];

	/**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
	**************************************************************************/
	startx[0] =  0.;                /* Set the Physical Minimum boundary  */
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	startx[1] = Rin;                /* Set the Physical Minimum boundary  */
	startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = 1. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = 1.;         /* radiative flux dumps       */
	n_restart = 1000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	dx[0] = 1.e-5;
	calc_all_geom() ;  
	dt_global_min = find_min_dt();

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
  int i, j, k, l ; 
  double X[NDIM], prim[NP];
  
  void get_bondi_data( double X[NDIM], double prim[NP] );

//  if( METRIC_TYPE_CHOICE != METRIC_KS_SPHERICAL ) { 
//    fprintf( stderr, "This initial data type was written for KS spherical coordinates !\n");
//    fflush(stderr); 
//    fail( FAIL_BASIC,0 ) ; 
//  }

  ALL_LOOP {
    coord(i,j,k,CENT,X) ; 
    
    get_bondi_data( X, prim );

    PLOOP p[i][j][k][l] = prim[l] ; 

  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */
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
  Rout     = 20.;         /* Radial extent of the grid           */
  n_within_horizon = 5;   /* number of cells within horizon */
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  Rin       = Rin_calc();

  h_slope  = 0.0;      /* Severity of the focusing            */

  th_cutout = M_PI*SMALL;       /* Angular size to excise from the axis */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_end-th_beg;    /* Angular size in theta of the grid    */

  X1_slope = 1.;       /* Severity of transition              */
  X1_0     = 1.e3;     /* Location of transition in X1 units  */ 
  
  
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
  fprintf(stdout,"th_cutout        = %28.18e \n", th_cutout        );
  fprintf(stdout,"th_beg           = %28.18e \n", th_beg           );
  fprintf(stdout,"th_end           = %28.18e \n", th_end           );
  fprintf(stdout,"th_length        = %28.18e \n", th_length        );
  fflush(stdout);


  /* Check for unreasonable values : */
  if( n_within_horizon <= 1 ) { 
    fprintf(stderr,"set_special_coord(): n_within_horizon must be greater than 1 , %d \n", 
	    n_within_horizon ); 
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

  return;
}


/***********************************************************************/
/***********************************************************************
        SPECIAL BONDI FUNCTIONS 
***********************************************************************/

#define MDOT (-1.)

/***********************************************************************/
/***********************************************************************
  get_bondi_data(): 
  -----------------
      -- returns the bondi solution at a point specified by xp in 
         numerical coordinates; 
      -- this routine also specifies parameters that define the solution;
      -- exact same method as used in HARM2D's tests;
      -- assumes that we are using KS_SPHERICAL coordinates;
***********************************************************************/
void get_bondi_data( double xp[NDIM], double prim[NP] )
{
  int i,j,k,l;

  static int first_call = 1;
  static double rs, urs, rhos, K, ps, edot;
  double r, th, rho, u;
  double ucon_bl[NDIM], ucon_ks[NDIM], gcov[NDIM][NDIM], x[NDIM];
  double hawley_bondi(double K,double gam,double edot,double mdot, double rs, double r) ;

  /* solve bondi eigenvalue problem */
  /* assume GM = c = Mdot = 1 */
  if( first_call ) { 

    rs = 8. ;		/* pick sonic radius */
    urs = -1./sqrt(2.*rs) ;	/* find ur at sonic point */
				/* find entropy constant */
    rhos = MDOT/(rs*rs*urs) ;	
    K = pow(2.,(1.-gam)/2.)*(gam - 1.)*
      pow(-MDOT/(rs*sqrt(rs)),1.-gam)/
      (gam*(2. - 2.*rs + gam*(-3. + 2.*rs))) ;
    ps = K*pow(rhos,gam) ;
    /* find Edot */
    edot = rs*rs*(rhos + ps*gam/(gam - 1.))*sqrt(1. - 2./rs
						 + urs*urs)*urs ;

    fprintf(stdout,"sonic radius: %g\n", rs) ;
    fprintf(stdout,"ur at sonic radius: %g\n", urs) ;
    fprintf(stdout,"edot: %g\n", edot) ;
    fflush(stdout);
    first_call = 0 ; 
  }

  x_of_xp( x, xp );
  r  = x[RR];
  th = x[TH];

  //fprintf(stdout,"coords: %10.4e %10.4e \n", r, th);  fflush(stdout);

  /* Contravariant 4-velocity components in BL coordinates: */
//  if( r <= r_horizon ) { 
//    PLOOP prim[l] = 0.;
//    prim[RHO]   = -MDOT/(r_horizon*r_horizon);
//    prim[ UU]   = K*pow(rho,gam)/(gam - 1.) ;
//    return;
//  }
    
  ucon_bl[TT] = 0.;
  ucon_bl[RR] = hawley_bondi(K,gam,edot,MDOT,rs,r) ;
  ucon_bl[TH] = 0. ;
  ucon_bl[PH] = 2.*a/(r*r*r) ;

  rho = MDOT/(r*r*ucon_bl[RR]) ;
  u   = K*pow(rho,gam)/(gam - 1.) ;

  fprintf(stdout,"bondi: %10.4e %10.4e %10.4e %10.4e %10.4e \n", 
	  r, rho, u, ucon_bl[RR], ucon_bl[PH]);  fflush(stdout);

  /* Find time component of 4-velocity: */
  bl_gcov_func( x, gcov );      
  setutcon( ucon_bl, gcov );

  fprintf(stdout,"bondi2: %10.4e  %10.4e  %10.4e  \n",  r, ucon_bl[TT], ucon_bl[RR]) ;

  /* Transform into computation coordinates (need to first do BL->KS): */
  bl_to_ks_con( x, ucon_bl, ucon_ks ); 

  fprintf(stdout,"bondi3: %10.4e  %10.4e  %10.4e  \n",  r, ucon_ks[TT], ucon_ks[RR]) ;

  transform_rank1con( x, xp, ucon_ks ); 

  fprintf(stdout,"bondi4: %10.4e  %10.4e  %10.4e  \n",  r, ucon_ks[TT], ucon_ks[RR]) ;

  /* get metric in numerical coordinates */
  ks_gcon_func( x, gcov );      /* gcov is now the contravariant metric */
  transform_rank2con( x, xp, gcov );      /* gcov is now the contravariant metric */

  /* Calculate the primitive variables : */
  prim[RHO] = rho ;
  prim[ UU] = u ;
  prim[ U1] = ucon_ks[RR] - gcov[TT][RR] * ucon_ks[TT] / gcov[TT][TT] ; 
  prim[ U2] = ucon_ks[TH] - gcov[TT][TH] * ucon_ks[TT] / gcov[TT][TT] ; 
  prim[ U3] = ucon_ks[PH] - gcov[TT][PH] * ucon_ks[TT] / gcov[TT][TT] ; 
  prim[ B1] = 0. ;
  prim[ B2] = 0. ;
  prim[ B3] = 0. ;

  return;
}

#define DEL 	2.e-5 
#define TOL	1.e-12

double hawley_bondi(double K,double gam,double edot,double mdot, double rs, double r) 
{
	double ur,x1,x2,y1,y2 ;
	double bondi_func(double K, double gam, double edot, double mdot, 
			double r, double ur) ;

	if(r < rs) 
		ur = -pow(r,-0.5) ;	/* guess */
	else
		ur = -0.5*pow(r,-1.5) ;	/* guess */

	x1 = ur ;
	x2 = ur*(1. + DEL) ;

	y1 = bondi_func(K,gam,edot,MDOT,r,x1) ;	
	y2 = bondi_func(K,gam,edot,MDOT,r,x2) ;	
	/*
	fprintf(stderr,"%g %g %g %g %g\n",r,x1,x2,y1,y2) ;
	*/

	while(fabs(y1) > TOL) {
		x1 = (x2*y1 - x1*y2)/(y1 - y2) ;
		x2 = x1*(1. + DEL) ;

		y1 = bondi_func(K,gam,edot,MDOT,r,x1) ;	
		y2 = bondi_func(K,gam,edot,MDOT,r,x2) ;	
		/*
		fprintf(stderr,"%g %g %g %g %g\n",r,x1,x2,y1,y2) ;
		*/
	}

	return(x1) ;
}
#undef TOL
#undef DEL

double bondi_func(double K, double gam, double edot, double mdot, double r, 
		double ur) 
{
	return(-edot + sqrt(1. - 2./r + ur*ur)*(MDOT +
		gam*K*r*r*pow(MDOT/(r*r*ur),gam)*ur/(gam - 1.))) ; 
}


#undef MDOT




