
#include "decs.h"

static double rhomax, umax;

static void mm08(double r, double r_s, double h_o_r, double M, double Sigma0, double delta, double xi, double a_sep, double a, double alpha, double *c_s_sq, double *Sigma, double *P, double *v_r, double *Omega, double *Omega_k);

/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 
 
   -- Post-Newtonian version of MacFadyen and Milosavljevic (2008) simulation; 
      We will refer to this paper as MM08

   -- Initial data for a cicumbinary disk in an evolving spacetime determined by 
       post-Newtonian approximations.  The initial data consists of an isothermal
       disk surrounding a binary.  With the assumption that the disk is isothermal and
       heat from dissipation is immediately radiated, we can assume the disk is independent 
       of theta.  MM08 solved the height-integrated equations, or the z-integrated equations. 

   --  The disk starts at approximately twice the initial 
       separation of the binary.  The numerical grid starts outside the binary's orbit, 
       so the black holes live outside the grid.  The inner boundary condition is imposed
       to being inflow conditions.  The coordinates are spherical coordinates. 

   -- MM08 parameters:
          -- rmin = a
          -- rout = 100 a 
          -- no flow onto the grid from the radial boundaries is allowed 
          -- Sigma = Sigma_0 * (r/r_s)^(-delta) * exp( -(r/r_s)^(-xi) )
                -- Sigma_0 = arbitrary (set it to 1)
                -- r_s = 10 * a
                -- delta = 3
                -- xi = 2 
          -- the initial angular velocity is corrected for radial pressure gradient and the initial 
              radial velocity is corrected for the viscous drift 
          -- v_phi is corrected from the torque 
               Omega = Omega_K^2 (1 + (3*a^2/(16*r^2)))^2 +  (dP/dr) / (r Sigma)
              where P = vertically integrated pressure = 2*h*p  
                    p = normal pressure

          -- v_r = 2*(r^2*Omega*Sigma)^(-1) d(r^2 sigma_{r \phi})/dr
               where sigma_{r \phi} is the only non-vanishing component of the viscous stree tensor:
               -- sigma_{r \phi} = Sigma(r) * \nu(r)  * ( r d(v_phi/r)/dr + (dv_r/dphi)/r ) 
               -- \nu(r) = alpha c_s(r)^2 / Omega_K(r)
               -- alpha = 0.01

          -- Resolution:
              dr   (r=rmin) =  0.039 a 
              dphi (r=rmin) =  0.0078 * (2*pi)

              dr   (r=rmax) =  0.31 a 
              dphi (r=rmax) =  0.062  * (2*pi)


            


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

// #if( (METRIC_TYPE_CHOICE != METRIC_GENERAL_DYNAMIC) || (TOP_TYPE_CHOICE != TOP_SPHERICAL) )  
//   initial-data-is-inconsistent-with-metric-and-coordinates
// #endif

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

	/* Set the parameters that define the grid and other constants : */ 
	gam = (2.) ;
	cour = 0.8;


	/* Coordinate dependent quantities : */ 
	set_special_coord();

	t = 0. ;                     /* Initial time */ 

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 100000. ;         /* Length of X0 dimension (evolution period) */ 
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
	GridLength[3] = M_PI ;                 /* Length of X3 dimension */ 

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
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	startx[1] = Rin;                /* Set the Physical Minimum boundary  */
	startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5 - 0.5*dx[2] ;
# else 
 	startx[2] = 0.            ; 
# endif
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	startx[1] = 1.e-10;        /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5-0.5*dx[2] ; 
# else 
 	startx[2] = 0.            ; 
# endif 
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 50. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = 50. ;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = 1.;         /* radiative flux dumps       */
	DT_out[OUT_MIN_DT]  = 100.;
	n_restart = 1000;                /* number of time steps between restart dumps */
	

	/**************************************************************************
	  Calculate the static metric grid functions 
	**************************************************************************/
	dx[0] = 1.e-5;
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


	/* in case you want to test the metric and coord. transf. routines */
	//	test_geom();   

  TRACE_END;

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
  struct of_coord *coords;
  
  void init_prim( int i, int j, int k, struct of_coord *coords, double *pr);
  void set_mag_field( void ) ; 

  TRACE_BEG;

  /*******************************************************************************
     Set the hydrodynamic quantities (rho,uu,v^i) : 
  *******************************************************************************/
  ALL_LOOP {
    get_coord(i,j,k,CENT,ncurr,coords);	
    init_prim( i, j, k, coords, p[i][j][k]);
  }

  /*******************************************************************************
    Normalize the densities: 
  *******************************************************************************/
  mpi_global_max(&rhomax); 
  mpi_global_max(&umax); 
  fprintf(stdout, "init_data(): orig rhomax   = %28.18e \n", rhomax); 
  fprintf(stdout, "init_data(): orig   umax   = %28.18e \n",   umax);   fflush(stdout); 
  
//  LOOP {
//    p[i][j][k][RHO]  /=  rhomax; 
//    p[i][j][k][ UU]  /=  rhomax; 
//  }
//  umax /= rhomax; 
//  rhomax = 1. ; 
//
//  fprintf(stdout, "init_data(): new  rhomax = %28.18e \n", rhomax); 
//  fprintf(stdout, "init_data(): new    umax = %28.18e \n", umax);   fflush(stdout); 


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
  M = m_bh_tot = 1.;
  m_bh1 = M/2.;                            /* Mass of left-most BH */
  m_bh2 = M/2.;                            /* Mass of right-most BH */
  initial_bbh_separation = 100.*m_bh_tot;  /* Initial separtion of BHs in units of total mass */

  r_horizon1 = m_bh1;             /* Radius of the horizon about BH1 in Harmonic coordinates */
  r_horizon2 = m_bh2;             /* Radius of the horizon about BH2 in Harmonic coordinates */
  
  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;

  a = 0.0;          /* Spin of the black hole in units of M */ 
  r_isco    = 0.;
  r_horizon = 0.;

  n_within_horizon = 0;            /* number of cells within horizon */

  Rin      = 1. * (initial_bbh_separation) ; 
  Rout     = 100.*(initial_bbh_separation);        /* Radial extent of the grid           */
  R0       = 0.9*Rin;

  h_slope   = 0.0729419909586;           /* Severity of the focusing             */
  X1_slope  = 1.;       /* Severity of transition              */
  X1_0      = log(1.e6);     /* Location of transition in X1 units  */ 

  //  th_cutout = 0.02655 * M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.0* M_PI;     /* Angular size to excise from the axis */
  //  th_cutout = 0.045 * M_PI;     /* Angular size to excise from the axis */
  th_cutout = M_PI*SMALL;     /* Angular size to excise from the axis */
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


#if( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
  this-initial-data-is-not-setup-for-COORD_DIAGONAL2
#endif 

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
        SPECIAL DISK FUNCTIONS 
***********************************************************************/

/***********************************************************************************/
/***********************************************************************************
  init_prim(): 
  ------------------
   -- This is supposed to follow the configuration described in MM08
   -- The disk parameters mentioned in MM08 are used below; 
***********************************************************************************/
void init_prim(int i, int j, int k, struct of_coord *coords, double *pr)
{
  double r,th, rho, uu, r_rat,ftmp1,ftmp2;
  double rhoflr, uuflr;
  struct of_geom geom_bl;
  struct of_geom *geom;

  static double Sigma_0, hor_tmp,h_o_r,r_s,delta,xi; 
  double Omega, Omega_k, Sigma, P, scale_height, c_s_sq, v_r, h; 
  double ucon[NDIM], ucon_bl[NDIM];

  static int first_call = 1;


  /***********************************************************************************************
    Set parameters used for all points : 
  **********************************************************************************************/  
  
  if( first_call ) { 
    alpha_visc = 10.;                          /* alpha used in the viscosity function */ 
    h_o_r      = 0.1;                           /* Disk scale height */
    r_s        = 10.*initial_bbh_separation;    /* Transition radius */
    Sigma_0    = 1.;                            /* Arbitrary density scale */
    delta      = 3.;                            /* exponent parameter */ 
    xi         = 2.;                            /* exponent parameter */
    /* Extra factor is to compensate for the fact that MM08 use H = Gaussian Std. Dev. instead of scaleheight;
       See Noble, Krolik & Hawley 2010 */
    hor_tmp = h_o_r / sqrt(0.5*M_PI);
    
    fprintf(stdout,"\n##################################################\n");
    fprintf(stdout,"  MM DISK PARAMETERS \n------------------------------------\n");
    fprintf(stdout,"\t  alpha_visc =  %28.18e \n",alpha_visc );
    fprintf(stdout,"\t  h_o_r      =  %28.18e \n",h_o_r      ); 
    fprintf(stdout,"\t  r_s        =  %28.18e \n",r_s        ); 
    fprintf(stdout,"\t  Sigma_0    =  %28.18e \n",Sigma_0    ); 
    fprintf(stdout,"\t  delta      =  %28.18e \n",delta      ); 
    fprintf(stdout,"\t  xi         =  %28.18e \n",xi        ); 
    fprintf(stdout,"\t  gam        =  %28.18e \n",gam        ); 

    first_call = 0;
  }
  
  /***********************************************************************************************
    Set the values in the disk :
  **********************************************************************************************/  

  /*  note we may have to use a different definition of radius if we are using PN data :*/
#if( METRIC_DYNAMIC_TYPE_CHOICE != METRIC_DYNAMIC_KS_SPHERICAL )
  coords->x[RR] += m_bh_tot;  /* Transform to BL coordinates */
#endif

  r            = coords->x[RR];
  h            = h_o_r * r; 
  scale_height = 2.*h;

#if( METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_KS_SPHERICAL )
  /* if we are using kerr, then we should assume the potential is that of a single black hole: */
  mm08(r, r_s, hor_tmp, M, Sigma_0, delta, xi, 0., a, alpha_visc, &c_s_sq, &Sigma, &P, &v_r, &Omega, &Omega_k);
#else
  mm08(r, r_s, hor_tmp, M, Sigma_0, delta, xi, initial_bbh_separation, a, alpha_visc, &c_s_sq, &Sigma, &P, &v_r, &Omega, &Omega_k);
#endif

  Katm[i]      = c_s_sq;
  nu_visc[i]   = alpha_visc * c_s_sq / Omega_k; 

  rho      = Sigma / scale_height;
  uu       = P / ( scale_height * (gam - 1.) );

  fprintf(stdout,"mmdata %d  %26.16e %26.16e %26.16e %26.16e %26.16e \n",i,r,rho,uu,c_s_sq,
	  target_temperature( r, hor_tmp )
	  ); fflush(stdout);


  get_special_geometry( coords, &geom_bl, METRIC_BL_SPHERICAL );
  ucon_bl[TT] = 1./sqrt( -(geom_bl.gcov[TT][TT] + v_r*v_r*geom_bl.gcov[RR][RR] + Omega * ( 2.*geom_bl.gcov[TT][PH] + Omega*geom_bl.gcov[PH][PH] ) ) );
  ucon_bl[RR] = v_r * ucon_bl[TT];
  ucon_bl[TH] = 0.;
  ucon_bl[PH] = Omega * ucon_bl[TT];

  bl_to_ks_con(coords->x,ucon_bl,ucon);  // transform from BL to KS coordinates 

  //  fprintf(stdout,"v_r(r=%26.16e) = %26.16e \n",r,v_r); fflush(stdout);

#if( METRIC_DYNAMIC_TYPE_CHOICE != METRIC_DYNAMIC_KS_SPHERICAL )
  coords->x[RR] -= m_bh_tot;  /* Transform back from BL coordinates */
#endif
  transform_rank1con2(coords->dxp_dx,ucon); // transform from (r,th,ph) to (x1,x2,x3)

  get_geometry(i,j,k,CENT,ncurr,geom); // get (x1,x2,x3) metric, from here on it's only (x1-3)

  rhoflr = RHOMIN*pow(r,RHOPOWER);
  uuflr  = UUMIN*pow(r, UUPOWER);

  if( (rho < rhoflr) || (uu < uuflr) ) {
    rho = rhoflr;
    uu  = uuflr;
    pr[U1] = 0.;
    pr[U2] = 0.;
    pr[U3] = 0.;
  }
  else { 
    ucon2pr( pr, ucon, geom->gcon );  // get prim. velocities;
  }
  
  pr[RHO] = rho;
  pr[UU]  = uu;

  if( rho > rhomax )  { rhomax = pr[RHO] ; } 
  if( uu  >   umax )  {   umax = pr[ UU] ; } 

  return;

}

/*************************************************************************************************************************/
/*************************************************************************************************************************
  mm08(): 
  ----------

    -- assumes the definition of constant scale height temperature defined in  target_temperature();
    -- assumes relativistic form of  Omega_K , i.e. with spin term; 

*************************************************************************************************************************/
void mm08(double r, double r_s, double h_o_r, double M, double Sigma0, double delta, double xi, double a_sep, double a, double alpha, double *c_s_sq, double *Sigma, double *P, double *v_r, double *Omega, double *Omega_k)
{ 

    double t1                 ;
    double t6 		       ;
    double t7 		       ;
    double t9 		       ;
    double t10 	       ;
    double t15 	       ;
    double t16 	       ;
    double t20 	       ;
    double cs_n 	       ;
    double t21 	       ;
    double t22 	       ;
    double t34 	       ;
    double dcs_dr_n 	       ;
    double t40 	       ;
    double t41 	       ;
    double t46 	       ;
    double t48 	       ;
    double t52 	       ;
    double t54 	       ;
    double t65 	       ;
    double t71 	       ;
    double t76 	       ;
    double t103 	       ;
    double t107 	       ;
    double t109 	       ;
    double d2cs_dr2_n 	       ;
    double Omega_k_n 	       ;
    double t116 	       ;
    double t117 	       ;
    double t119 	       ;
    double t120 	       ;
    double Sigma_n 	       ;
    double t121 	       ;
    double dSigma_dr_n        ;
    double t125 	       ;
    double d2Sigma_dr2_n      ;
    double t136 	       ;
    double P_n 	       ;
    double dP_dr_n 	       ;
    double t144 	       ;
    double d2P_dr2_n 	       ;
    double t151 	       ;
    double t152 	       ;
    double t156 	       ;
    double t159 	       ;
    double Omega_sq_n 	       ;
    double Omega_n 	       ;
    double t161 	       ;
    double dOmega_dr_n        ;
    double t166 	       ;
    double t172 	       ;
    double t173 	       ;
    double dOmega_sq_dr_n     ;
    double t194 	       ;
    double t206 	       ;
    double d2Omega_sq_dr2_n   ;
    double t223 	       ;
    double d2Omega_dr2_n      ;
    double t231 	       ;
    double nu_n 	       ;
    double dnu_dr_n 	       ;
    double v_phi_n 	       ;
    double t235 	       ;
    double t236 	       ;
    double sigma_rp_n 	       ;
    double dsigma_rp_dr_n     ;
    double v_r_n              ;   

    t1 = sqrt(r);
    t6 = r*r;
    t7 = (4.0*t1-3.0*a)*a-t6;
    t9 = h_o_r*h_o_r;
    t10 = t1*r;
    t15 = t10-3.0*t1+2.0*a;
    t16 = 1/t15;
    t20 = sqrt(-2.0*M_PI*t7*t9/t10*t16);
    cs_n = t20/2.0;
    t21 = t6*r;
    t22 = t1*t21;
    t34 = 1/r;
    dcs_dr_n = (t22+((18.0-11.0*r)*r+(-26.0*t1+9.0*t10+9.0*a)*a)*a)*t34*t16/
      t7*cs_n/2.0;
    t40 = t6*t6;
    t41 = t40*t40;
    t46 = t1*t40;
    t48 = t40*t6;
    t52 = t1*t6;
    t54 = t40*r;
    t65 = a*a;
    t71 = t65*t65;
    t76 = t65*a;
    t103 = t15*t15;
    t107 = t7*t7;
    t109 = 1/t6;
    d2cs_dr2_n = 3.0/4.0*(t1*t41-3.0*t1*t40*t21+(261.0*t46+22.0*t1*t48+144.0*
	 t22-756.0*t52-161.0*t1*t54+(45.0*t46-3768.0*t10-597.0*t22+2426.0*t52+(225.0*t10
      -911.0*t1+126.0*a)*t65)*t65)*t65+(2627.0*t71*a+((2682.0-1240.0*t65)*t76+((
    -1824.0+171.0*t65)*t76+((270.0+476.0*t65)*a+((-315.0-28.0*t65)*a+(150.0*a-23.0*
    a*r)*r)*r)*r)*r)*r)*r)/t103/t15/t107*t109*cs_n;
    Omega_k_n = sqrt(M/t21);
    t116 = r/r_s;
    t117 = pow(t116,-delta);
    t119 = pow(t116,-xi);
    t120 = exp(-t119);
    Sigma_n = Sigma0*t117*t120;
    t121 = t119*xi;
    dSigma_dr_n = (-delta+t121)*t34*Sigma_n;
    //    t125 = pow(t116,-2.0*xi);
    t125 = t119 * t119;
    d2Sigma_dr2_n = -((t119+(-t125+t119)*xi)*xi+(2.0*t121-1.0-delta)*delta)*
      t109*Sigma_n;
    t136 = cs_n*cs_n;
    P_n = Sigma_n*t136;
    dP_dr_n = cs_n*(dSigma_dr_n*cs_n+2.0*Sigma_n*dcs_dr_n);
    t144 = dcs_dr_n*dcs_dr_n;
    d2P_dr2_n = 4.0*dcs_dr_n*dSigma_dr_n*cs_n+2.0*Sigma_n*t144+d2Sigma_dr2_n*
      t136+2.0*cs_n*Sigma_n*d2cs_dr2_n;
    t151 = Omega_k_n*Omega_k_n;
    t152 = a_sep*a_sep;
    t156 = pow(1.0+3.0/16.0*t152*t109,2.0);
    t159 = 1/Sigma_n;
    Omega_sq_n = t151*t156+dP_dr_n*t34*t159;
    Omega_n = sqrt(Omega_sq_n);
    t161 = 1/Omega_n;
    t166 = dP_dr_n*dSigma_dr_n;
    t172 = Sigma_n*Sigma_n;
    t173 = t151*t172;
    dOmega_sq_dr_n = ((-128.0*dP_dr_n*Sigma_n+(128.0*d2P_dr2_n*Sigma_n-128.0*
       t166)*r)*t21+(-96.0*t173*t6-18.0*t173*t152)*t152)/t54/t172/128.0;
    dOmega_dr_n = t161*dOmega_sq_dr_n/2.0;
    t194 = dSigma_dr_n*dSigma_dr_n;
    t206 = t172*Sigma_n;
    d2Omega_sq_dr2_n = ((64.0*dP_dr_n*t172+(-64.0*d2P_dr2_n*t172+64.0*t166*
    Sigma_n+(-64.0*d2P_dr2_n*dSigma_dr_n*Sigma_n+64.0*dP_dr_n*t194-32.0*dP_dr_n*
     d2Sigma_dr2_n*Sigma_n)*r)*r)*t48+(384.0*t40*t206+(360.0*t6*t206+63.0*t206*t152)
		       *t152)*M)/t41/r/t206/32.0;
    t223 = dOmega_sq_dr_n*dOmega_sq_dr_n;
    d2Omega_dr2_n = (-t223+2.0*d2Omega_sq_dr2_n*Omega_sq_n)/Omega_n/
      Omega_sq_n/4.0;
    t231 = 1/Omega_k_n;
    nu_n = alpha*t136*t231;
    dnu_dr_n = 2.0*alpha*cs_n*t231*dcs_dr_n;
    v_phi_n = r*Omega_n;
    t235 = Sigma_n*nu_n;
    t236 = r*dOmega_dr_n;
    sigma_rp_n = t235*t236;
    dsigma_rp_dr_n = dSigma_dr_n*nu_n*t236+Sigma_n*dnu_dr_n*t236+t235*
      dOmega_dr_n+t235*r*d2Omega_dr2_n;
    v_r_n = 2.0*(2.0*r*sigma_rp_n+t6*dsigma_rp_dr_n)*t109*t161*t159;

    *c_s_sq = t136;
    *Sigma = Sigma_n;
    *P = P_n;
    *v_r = v_r_n;
    *Omega = Omega_n;
    *Omega_k = Omega_k_n; 

    return;
}
