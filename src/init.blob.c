
#include "decs.h"

double rin, rdiskout, beta, rhomax=0, umax=0, bsq_max=0;

double  ranc(int iseed); 

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

Localized sphere at a given phi orbitting at 4-velocity of  equatorial keplerian 
  rate of its center. 


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
	cour = 0.9;


	/* Coordinate dependent quantities : */ 
	set_special_coord();

	t = 0. ;                     /* Initial time */ 

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 500. ;         /* Length of X0 dimension (evolution period) */ 
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	GridLength[1] = Rout-Rin  ;   /* Length of X1 dimension */ 
	GridLength[2] = th_length ;   /* Length of X2 dimension */ 
#else
	GridLength[1] = log((Rout-R0)/(Rin-R0)) ; /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                     /* Length of X2 dimension */ 
#endif
	GridLength[3] = 0.5*M_PI ;                 /* Length of X3 dimension */ 

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
#else 
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
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

	fprintf(stdout,"rhomin,minlimt = %28.18e %28.18e \n", RHOMIN, RHOMINLIMIT);
	fprintf(stdout,"uumin ,minlimt = %28.18e %28.18e \n", UUMIN,  UUMINLIMIT);
	
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
  double X[NDIM];
  
  void init_prim( int i, int j, int k, double *xp, double *pr);
  void set_mag_field( void ) ; 

  if( METRIC_TYPE_CHOICE != METRIC_KS_SPHERICAL ) { 
    fprintf( stderr, "This initial data type was written for KS spherical coordinates !\n");
    fflush(stderr); 
    fail( FAIL_BASIC,0 ) ; 
  }

  /*******************************************************************************
     Set the hydrodynamic quantities (rho,uu,v^i) : 
  *******************************************************************************/
  ranc(0);
  LOOP {
    coord(i,j,k,CENT,X) ; 
    
    init_prim( i, j, k, X, p[i][j][k]);

  }

  /*******************************************************************************
    Normalize the densities: 
  *******************************************************************************/
  mpi_global_max(&rhomax); 
  mpi_global_max(&umax); 
  fprintf(stdout, "init_data(): orig rhomax = %28.18e \n", rhomax); 
  fprintf(stdout, "init_data(): orig   umax = %28.18e \n",   umax);   fflush(stdout); 

  LOOP {
    p[i][j][k][RHO]  /=  rhomax; 
    p[i][j][k][ UU]  /=  rhomax; 
  }
  umax /= rhomax; 
  rhomax = 1. ; 

  fprintf(stdout, "init_data(): new  rhomax = %28.18e \n", rhomax); 
  fprintf(stdout, "init_data(): new    umax = %28.18e \n", umax);   fflush(stdout); 


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
#else 
  set_mag_field(); 
#endif 

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
  
  a = 0.0;          /* Spin of the black hole in units of M */ 

  R0       = 0.;          /* Offset in Radius from 1.            */
  Rout     = 120.;         /* Radial extent of the grid           */
  n_within_horizon = 10;   /* number of cells within horizon */

  h_slope  = 0.13;      /* Severity of the focusing            */

  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
  th_beg    = th_cutout;        /* beg/end set the limits of theta  and */
  th_end    = M_PI-th_cutout;   /*  can override th_cutout              */
  th_length = th_beg-th_end;    /* Angular size in theta of the grid    */

  X1_slope = 1.;       /* Severity of transition              */
  X1_0     = log(1.e6);     /* Location of transition in X1 units  */ 
  
  
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  //  Rin       = Rin_calc();
  Rin       = 1.9;

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

  return;
}


/***********************************************************************/
/***********************************************************************
        SPECIAL TORUS FUNCTIONS 
***********************************************************************/



/***********************************************************************************/
/***********************************************************************************
  init_prim(): 
  ------------------
   This is pretty much the same routine as in HARM2D to set Fishbone & Moncrief 
    initial data. 
 ***********************************************************************************/
void init_prim(int i, int j, int k, double *xp, double *pr)
{
  int  prograde;
  double x[NDIM], exp_f;
  double del_blob;
  double ucon[NDIM], ucon2[NDIM];

  /* for disk interior */
  double A,B,C,r_minus_3halfs, kappa, r,th,ph,xo,yo,zo,rhoflr,uuflr,rhoscal,uuscal;

  static double pr_max[NP], x_max[NDIM], rmax, xm,ym,zm, xpmax[NDIM];

  static int first_call = 1;


  /* Set the parameters that define the torus : */
  del_blob = 12.;
  kappa = 1.e-3 ;       /* Constant entropy factor in polytropic EOS used at t=0 */
  beta = 1.e2 ;         /* Plasma beta parameter */


  if( first_call ) { 
    xpmax[0] = 0.;  
    xpmax[1] = startx[1] + 0.8*GridLength[1];  
    xpmax[2] = startx[2] + 0.5*GridLength[2] + NG*dx[2];
    xpmax[3] = 0.75*GridLength[3];
    x_of_xp(x_max,xpmax);
    
    rmax = x_max[RR];       /* Radius of pressure maximum    */
    r_minus_3halfs = 1./(rmax*sqrt(rmax));
    B = 1.0  +  a * r_minus_3halfs ;
    C = 2.*B - 1.  -  3./rmax;
    C = 1./sqrt(C);
    ucon[TT] = B * C; 
    ucon[RR] = 0.;
    ucon[TH] = 0.;
    ucon[PH] = C * r_minus_3halfs;
    ucon2[TT] = ucon[TT]; ucon2[RR] = ucon[RR]; ucon2[TH] = ucon[TH]; ucon2[PH] = ucon[PH]; 
    transform_rank1con( x_max, xpmax, ucon ); 
    ucon2pr( pr_max, ucon, gcon[i][j][CENT] );

    r = rmax ; th = x_max[TH]; ph = x_max[PH];
    xm = r*sin(th)*cos(ph);
    ym = r*sin(th)*sin(ph);
    zm = r*cos(th);

    fprintf(stdout,"rmax = %28.18e \n", rmax); fflush(stdout);

    first_call = 0; 
  }

  x_of_xp( x, xp ); 
  r = x[RR];  th = x[TH]; ph = x[PH];
  xo = r*sin(th)*cos(ph);
  yo = r*sin(th)*sin(ph);
  zo = r*cos(th);


  exp_f = (xm-xo)*(xm-xo) + (ym-yo)*(ym-yo) + (zm-zo)*(zm-zo);
  if( sqrt(exp_f) < 5*del_blob ) { 
    pr[RHO] = exp(-exp_f/del_blob);
    pr[UU] =  kappa * pow(pr[RHO], gam) / (gam - 1.);
    ucon2[TT] = ucon[TT]; ucon2[RR] = ucon[RR]; ucon2[TH] = ucon[TH]; ucon2[PH] = ucon[PH]; 
    transform_rank1con( x_max, xpmax, ucon ); 
    ucon2pr( pr, ucon, gcon[i][j][CENT] );
    if( pr[RHO] > rhomax ) {  rhomax = pr[RHO]; }
    if( pr[UU]  > umax   ) {  umax   = pr[UU];  }
  }
  else{
    pr[RHO] = pr[UU] = pr[U1] = pr[U2] =  pr[U3] = 0.;
  }    
    
    


//  fprintf(stdout,"xm,ym,zm = %g %g %g \n",xm,ym,zm);
//  fprintf(stdout,"xo,yo,zo = %g %g %g \n",xo,yo,zo);
//  fprintf(stdout,"exp_f,del=  %g %g \n",exp_f,del_blob);
//  fprintf(stdout,"pr =  %g %g %g %g %g \n",pr[RHO],pr[UU],pr[U1],pr[U2],pr[U3]);


  return;

}

/***********************************************************************************/
/***********************************************************************************
  set_mag_field():
  ------------------
   -- Sets the magnetic field; 
   -- Currently only supports poloidal fields derived from the azimuthal component 
      of the vector potential; 
   -- Before the calculation, does a fixup() and sets the boundaries since 
      the derivatives of A_\phi require ghost cells (field follows density contours); 
   -- We assume below that the magnetic field is independent of X3 ; 
***********************************************************************************/
void set_mag_field(void )
{
  int i,j,k,l;

  double rho_avg, q, bsq_ij, beta_act, norm ; 
  double xp[NDIM],x[NDIM], fval,rmin;
  struct of_geom *geom;  
  
  /*************************************************************************************
    Calculate the azimuthal component of the vector potential, A_\phi: 
       --  Use F[] as a work array to hold A_\phi . 
       --  We assume here that rho is initially independent of X3 
       --  Need to calc. A_\phi on first ghost cells for finite diff. calculation for B^i 
   *************************************************************************************/
  rmin = 3.;

  k = N3S; 
  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) { 

    coord(i,j,k,CENT,xp);
    x_of_xp(x,xp);

    fval = x[RR]*cos(x[TH]); 

    F[i][j][0][0][0] = (fval > rmin)  ?  fval  :  rmin ;

  }
    
  /*************************************************************************************
    Calculate the poloidal field components given A_\phi
       --  Note that we can take derivatives w.r.t. numerical coordinates here because 
           the connection is symmetric in its lower indices, we just have to use 
           the numerical metric determinant; 
  *************************************************************************************/
  N1_LOOP  N2_LOOP { 
    get_geometry(i,j,0,CENT,ncurr,geom); 
    /* flux-ct */
    p[i][j][0][B1] = -( F[i  ][j  ][0][0][0] 
		      - F[i  ][j+1][0][0][0] 
		      + F[i+1][j  ][0][0][0] 
                      - F[i+1][j+1][0][0][0] 
		       )/(2.*dx[2]*geom->g) ;

    p[i][j][0][B2] = (  F[i  ][j  ][0][0][0] 
		      + F[i  ][j+1][0][0][0] 
		      - F[i+1][j  ][0][0][0] 
		      - F[i+1][j+1][0][0][0] 
			)/(2.*dx[1]*geom->g) ;
    
    p[i][j][0][B3] = 0. ;
    
    bsq_ij = bsq_calc(p[i][j][0],geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }

  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) {  F[i][j][0][0][0] = 0.; }


  mpi_global_max( &bsq_max );
  fprintf(stdout,"initial bsq_max = %28.18e \n", bsq_max); fflush(stdout); 


  /************************************************************************************
     Normalize magnetic field to match specified value of beta (using umax);
   ************************************************************************************/
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
  fprintf(stdout,"initial beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);

  norm = sqrt(beta_act/beta) ;
  bsq_max = 0. ;
  
  N1_LOOP  N2_LOOP { 
    get_geometry(i,j,0,CENT,ncurr,geom) ;

    p[i][j][0][B1] *= norm ;
    p[i][j][0][B2] *= norm ;

    bsq_ij = bsq_calc(p[i][j][0],geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }

  mpi_global_max( &bsq_max );
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
  fprintf(stdout,"final beta: %28.18e (should be %28.18e)\n",beta_act,beta); fflush(stdout);

  /************************************************************************************
     Copy B^i  from X3=0 slice to other slices X3 slices : 
  ************************************************************************************/
  ALL_LOOP  for( l = B1 ; l <= B3 ; l++ ) { 
    p[i][j][k][l] = p[i][j][0][l] ; 
  }
    

  return;

}

