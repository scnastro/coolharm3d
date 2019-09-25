
#include "decs.h"
#include "metric.h"


static double rhomax=0.,umax=0.,bsq_max=0.;
double  ranc(int iseed); 


/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

   -- This is the monopole field initial data that is used in 
             Komissarov MNRAS 350 (2004). 

   -- generalized to generalized coordinates and spacetimes; 
         -- using physics from init.monopole.c  but the modern coordinate system construction 
            from init.gen_disk.c    circa  Wed May 21 13:48:34 EDT 2014

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

#if( COORD_RADIUS_DIM == 1 )
//--HERE  fprintf(stdout,"init.gen_disk.c:  we should not be here.... use init.gen_axi_disk.c instead \n"); fflush(stdout); fail(FAIL_BASIC,0);
#endif


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
	gam = (5./3.) ;
	cour = 0.45;

	t = 0. ;                     /* Initial time */ 

	/* Coordinate dependent quantities : */ 
	set_special_coord();

	/**************************************************************************
          Length of each dimension : 
	**************************************************************************/
	GridLength[0] = 5.e4 ;         /* Length of X0 dimension (evolution period) */ 

	/* SPHERICAL CASES: */
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL ) 
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
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	GridLength[1] = log((coord_params->upsilon_r-coord_params->f0_r)/(1.-coord_params->f0_r)) ;              /* Length of X1 dimension */ 
	GridLength[2] = 1.  ;                    /* Length of X2 dimension */ 
	//	GridLength[2] = 1. - 2*th_cutout/M_PI ;          /* Length of X2 dimension */ 
#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[1] = 1.  ;              /* Length of X1 dimension */
        GridLength[2] = 1.  ;                    /* Length of X2 dimension */

#endif
	GridLength[3] = 2.*M_PI ;                 /* Length of X3 dimension */ 
#if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	GridLength[3] = 1. ;                 /* Length of X3 dimension */ 
#endif

	/**************************************************************************
          Grid discretization scales in each dimension : 
	**************************************************************************/
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];
#if(EQUATORIAL_RUN)
 	dx[2] = 1./160.;
#endif

	/**************************************************************************
          Starting coordinates in each dimension (global numerical boundary) :
	**************************************************************************/
	startx[0] =  0.;                /* Set the Physical Minimum boundary  */
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	startx[1] = Rin;                /* Set the Physical Minimum boundary  */
	startx[2] = th_beg;             /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5*GridLength[2]-0.5*dx[2] ; 
# else 
 	startx[2] = th_beg      ;
# endif 

#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
        //	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5 - 0.5*dx[2] ;
# else 
 	startx[2] = 0.            ; 
# endif
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL2 )
	startx[1] = log(Rin-R0);        /* Set the Physical Minimum boundary  */
	startx[2] = xi_diag2[0];        /* Set the Physical Minimum boundary  */
#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
	startx[1] = 0.;        /* Set the Physical Minimum boundary  */
	startx[2] = 0.;                 /* Set the Physical Minimum boundary  */
	//	startx[2] = th_cutout/M_PI ;       /* Set the Physical Minimum boundary  */

#elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
	startx[1] = 1.e-10;        /* Set the Physical Minimum boundary  */
# if(EQUATORIAL_RUN)
	startx[2] = 0.5-0.5*dx[2] ; 
# else 
 	startx[2] = 0.            ; 
# endif 
#endif
	startx[3] =  0.;                /* Set the Physical Minimum boundary  */
#endif  /* TOP_SPHERICAL */



	/* CARTESIAN CASES: */
#if( TOP_TYPE_CHOICE == TOP_CARTESIAN ) 
#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	GridLength[1] = 1. ;   /* Length of X1 dimension */ 
	GridLength[2] = 1. ;   /* Length of X1 dimension */ 
	GridLength[3] = 1. ;   /* Length of X1 dimension */ 

	SDLOOP1 startx[i] = -0.5*GridLength[i];
	startx[1] +=  initial_bbh_separation * m_bh1 / m_bh_tot;
	SDLOOP1 dx[i] =  GridLength[i] / totalsize[i];
#  if(EQUATORIAL_RUN)
	dx[2] = 1.e-10;
	startx[2] = -0.5*dx[2] ; 
#  endif 
#else 
	--not supported yet
#endif
#endif /* TOP_CARTESIAN */


	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	/**************************************************************************
          Output frequencies : 
	**************************************************************************/
	DT_out[OUT_ASCII]   = GridLength[0]/2. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/4. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 10. ;                	/* history file frequency */
	DT_out[OUT_SURFACE] = 30. ;                	/* surface file frequency */
	DT_out[OUT_HDF5]    = 4.e-1/20;                     /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 4.*DT_out[OUT_HDF5];      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_STAT];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	DT_out[OUT_MIN_DT]  = DT_out[OUT_HISTORY];
	n_restart = 10000;                /* number of time steps between restart dumps */
	

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
	if( myid == printer_pid ) {  
	  fprintf(stdout,"Length: %10.4e %10.4e %10.4e %10.4e\n", 
		  GridLength[0], GridLength[1], GridLength[2], GridLength[3]);
	  fprintf(stdout,"dx    : %10.4e %10.4e %10.4e %10.4e\n", 
		  dx[0], dx[1], dx[2], dx[3]);
	  fprintf(stdout,"startx: %10.4e %10.4e %10.4e %10.4e\n", 
		  startx[0], startx[1], startx[2], startx[3]);

	  fprintf(stdout,"rhomin,minlimt = %28.18e %28.18e \n", RHOMIN, RHOMINLIMIT);
	  fprintf(stdout,"uumin ,minlimt = %28.18e %28.18e \n", UUMIN,  UUMINLIMIT);
	
	  fflush(stdout);
	}

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
  double *ploc; 
  struct of_geom *geom;
  struct of_coord *coords;

  double Acov[NDIM];
  double A0 = 1.;
  double factor = 1.e3;

  double xtmp[NDIM];
  double *xspher;
  double dxsrc_dxdest[NDIM][NDIM];
  extern void curl_of_A( double ****A , double ****prim ) ;

  /*******************************************************************************
     Set the primitives, calculating the mag-field without assumptions : 
  *******************************************************************************/
  ranc(myid+1);
  for(i=N1S; i<=(N1E+1); i++)  for(j=N2S; j<=(N2E+1); j++) for(k=N3S; k<=(N3E+1); k++) { 

	get_coord(i,j,k,CORN,ncurr,coords);

#if	(TOP_TYPE_CHOICE == TOP_CARTESIAN) 
	xspher = xtmp; 
	xspher_of_xcart(xspher, coords->x, dxsrc_dxdest );
#elif	(TOP_TYPE_CHOICE == TOP_SPHERICAL) 
	xspher = coords->x;
#else 
	--option--not--supported-yet
#endif 

	/* Specify vector potential in covariant spherical coordinates: */
	Acov[0] = 0.;
	Acov[1] = 0.;
	Acov[2] = 0.;
	Acov[3] = -A0 * cos(xspher[TH]) ;


	/* transform to the appropriate physical coordinates if needed : */
#if	(TOP_TYPE_CHOICE == TOP_CARTESIAN) 
        transform_rank1cov2(dxsrc_dxdest,Acov);
#endif

	/* Transform to our numerical coordinates: */
        transform_rank1cov2(coords->dx_dxp,Acov);

	ploc = ph[i][j][k];
	ploc[B1] = Acov[1];
	ploc[B2] = Acov[2];
	ploc[B3] = Acov[3];
      }

  curl_of_A(ph,p);

  LOOP { 
	ploc = p[i][j][k]; 
	ploc[U1] = ploc[U2] = ploc[U3] = 0.;
	get_geometry(i,j,k,CENT,ncurr,geom);
	double Bsq = bsq_calc(ploc,geom) / ( geom->alpha * geom->alpha );
	ploc[RHO] = factor * Bsq;
	ploc[UU]  = ploc[RHO] / (gam-1.);
  }

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
  M = m_bh_tot = 1.;

  /* BBH parameters: */
#if( BBH_SPACETIME ) 
  m_bh1 = 0.5*M;                            /* Mass of left-most BH */
  m_bh2 = 0.5*M;                            /* Mass of right-most BH */
  initial_bbh_separation = 20.*m_bh_tot;  /* Initial separtion of BHs in units of total mass */

  r_horizon1 = m_bh1;             /* Radius of the horizon about BH1 in Harmonic coordinates */
  r_horizon2 = m_bh2;             /* Radius of the horizon about BH2 in Harmonic coordinates */

  rsq_horizon1 = r_horizon1*r_horizon1;
  rsq_horizon2 = r_horizon2*r_horizon2;

  a = 0.0;          /* Spin of the black hole in units of M */ 
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

#if( SET_TSHRINK )
  /* Time at which to start shrinking the binary,  set to 0. if you want to shrink from the beginning and 
     set to an impossibly large number if you never want to shrink : */
  t_shrink_bbh = 0.; 
  phi_0_bbh    = 0.;
#endif

#else 
  /* Single BH Parameters: */
  a = 0.; 
  asq = a*a;
  r_isco    = risco_calc(1);
  r_horizon = rhorizon_calc(1);
  R0       = 0.;             /* Offset in Radius from 1.   (HARM-like)         */
  Rout     = 300.;         /* Radial extent of the grid           */
  n_within_horizon = 5;   /* number of cells within horizon */
  Rin       = Rin_calc();
  
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

