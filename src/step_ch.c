
#include "decs.h"
#include "metric.h"

/* Global predicate arrays used to do reconstruction and find fluxes : */ 
static const int idir[] = {1,0,0};
static const int jdir[] = {0,1,0};
static const int kdir[] = {0,0,1};
static double max_char_speed[NDIM];
static double dt_char_min;


/* External routines : */
extern void reconstruct_fast( double ****p_h );
extern void flux_ct(            double *****F );
extern void flux_ct_para(       double *****F );
extern void flux_ct_para_fast(  double *****F );
extern void flux_ct_para_fast2( double *****F );
extern int  Utoprim_2d_fast(double *U, struct of_geom *geom, double *prim, double *gamma, double *bsq);
//extern int  Utoprim_2d(double U[NP],double gcov[NDIM][NDIM],double gcon[NDIM][NDIM],double gdet,double prim[NP]);
extern int  Utoprim_1d(double *U,double gcov[NDIM][NDIM],double gcon[NDIM][NDIM],double gdet,double *prim, double *gamma, double *bsq);
//extern int  Utoprim_2d_new(double U[NP],double gcov[NDIM][NDIM],double gcon[NDIM][NDIM],double gdet,double prim[NP]);
extern int Utoprim_1dfix1(double *U, struct of_geom *geom, double *prim, double *gamma, double *bsq, double K );
extern int Utoprim_1dvsq2fix1(double *U, struct of_geom *geom, double *prim, double *gamma, double *bsq, double K );

extern void source(double *ph,struct of_geom *geom,int ii,int jj,int kk,double *dU, double Dt);
extern void fixup_interp_prim(double ****pv);
extern void fixup_interp_v(double ****pv);
extern void fixup_global(double ****pv) ;
extern int fixup_entropy_eq(int i, int j, int k, double *prim_old, double *prim, struct of_geom *geom, struct of_geom *geom_old, double *gamma_out);
extern void dump_stat(void) ;
extern void advance_geometry( void );

static double timestep(void);
static void numerical_flux(double *****p_L, double *****p_R, double *****F);
static void recover_primitives(int i, int j, int k, double *pf, double *U, struct of_geom *geom, struct of_geom *geom_old, struct of_coord *coords);
static void recover_primitives_simple(int i, int j, int k, double *pf, double *U, struct of_geom *geom, struct of_coord *coords);

static void advance_viscous(double ****prim, double Dt);
static void advance(double ****p_i, double ****p_h, double Dt, double ****p_f);

/******************************************************************************************/
/******************************************************************************************
  step_ch(): 
 -----------
     -- top-level time-stepping driver routine; 
     -- responsible for handling the order in which the timestep, boundary conditions 
        and fixups are done;
******************************************************************************************/
void step_ch(void)
{
  int i,j,k,l;


  TRACE_BEG;
  
	/* evaluate timestep */
	/* half-step */
        
	/*******************************************************************************/
	/***  PREDICTOR STEP *******/
	/*******************************************************************************/
        
	/* Zero out failure array */
	if( failure_exists ) { 
	  failure_exists = 0 ; 
	  for(i=0;i<N_PFLAG_TYPES;i++) { failure_type[i] = 0; }
	  ALL_LOOP {  pflag[i][j][k] = 0; }
	}

        n_substep = 0;  /* NOT the last substep in full dt integration */
	n_U = (n_substep + 1)  % N0;

	ncurr = n_beg;  /* set the current time position to be the beginning  */

	ALL_LOOP PLOOP {  
	  p_old[i][j][k][l] = ph[i][j][k][l] = p[i][j][k][l];  /* Needed for Utoprim() guess */
	}

#if(DYNAMIC_SPACETIME) 
	BEG_TIMING(TIMER_TYPE_METRIC);
	set_general_conn(t);  /* Need to update the connection here */
	END_TIMING(TIMER_TYPE_METRIC);
#endif

#if( USE_KINEMATIC_VISCOSITY )
	bound_visc_stress_functions(n_substep,p);
#endif


	BEG_TIMING(TIMER_TYPE_ADVANCE);
	advance(p, p, 0.5*dx[0], ph) ;
	END_TIMING(TIMER_TYPE_ADVANCE);

	ncurr = n_mid;  /* set the current time position to be "halfstep"  */


#if( FIX_CUTOUT_GAMMA || FIX_CUTOUT_PRESSURE )
	fixup_global(ph);
#endif 

	BEG_TIMING(TIMER_TYPE_BOUNDS);
	bounds(ph,0);         // set boundaries so fixup's interpolation works

	if( failure_exists ) {   // fix/interpolate over bad points 
	  if( failure_type[PFLAG_INTERP_V] )    { fixup_interp_v(ph);    }
	  if( failure_type[PFLAG_INTERP_PRIM] ) { fixup_interp_prim(ph); }
	}

	if( boundary_phys_pflag || boundary_mpi_pflag ) { 
	  bounds(ph,1);         // set boundaries using new fixed points
	  boundary_phys_pflag = boundary_mpi_pflag = 0;
	}
	END_TIMING(TIMER_TYPE_BOUNDS);

	/* Zero out failure array */
	if( failure_exists ) { 
	  failure_exists = 0 ; 
	  for(i=0;i<N_PFLAG_TYPES;i++) { failure_type[i] = 0; } 
	  ALL_LOOP {  pflag[i][j][k] = 0; }
	}


#if( CALC_CURRENT )
	if( current_time() ) {
	  current_calc();  /* Uses pL and pR to calculate the 4-volume-centered current density */
	}
#endif

#if(DYNAMIC_SPACETIME) 
	BEG_TIMING(TIMER_TYPE_METRIC);
	set_general_conn(t+0.5*dx[0]);  /* Need to update the connection here */
	END_TIMING(TIMER_TYPE_METRIC);
#endif


	/*******************************************************************************/
	/***  CORRECTOR STEP *******/
	/*******************************************************************************/

	n_substep = 1; /* THE LAST substep in full dt integration */
	n_U = (n_substep + 1)  % N0;

#if( MAKE_RADFLUX ) 
	/* Needed to reset radiative flux dump array so we do not redump old data
	    since we only calculate radiative flux for bound cells:  */
	for(i=0;i<NCELLS;i++)  { coolflux[0][i] = 0.;  }
	for(i=0;i<NCELLS;i++)  { coolflux[1][i] = 0.;  }
	for(i=0;i<NCELLS;i++)  { coolflux[2][i] = 0.;  }
	for(i=0;i<NCELLS;i++)  { coolflux[3][i] = 0.;  }
#endif	

#if( USE_KINEMATIC_VISCOSITY )
	bound_visc_stress_functions(n_substep,ph);
#endif

	/* full-step */
	BEG_TIMING(TIMER_TYPE_ADVANCE);
	advance(p, ph, dx[0], p) ;
	END_TIMING(TIMER_TYPE_ADVANCE);

	ncurr = n_end;  /* set the current time position to be the end time  */

#if( FIX_CUTOUT_GAMMA || FIX_CUTOUT_PRESSURE )
	fixup_global(p);
#endif

	BEG_TIMING(TIMER_TYPE_BOUNDS);
	bounds(p,0);         // set boundaries so fixup's interpolation works

	if( failure_exists ) {   // fix/interpolate over bad points 
	  if( failure_type[PFLAG_INTERP_V] )    { fixup_interp_v(p);    }
	  if( failure_type[PFLAG_INTERP_PRIM] ) { fixup_interp_prim(p); }
	}
	
	if( boundary_phys_pflag || boundary_mpi_pflag ) { 
	  bounds(p,1);         // set boundaries using new fixed points
	  boundary_phys_pflag = boundary_mpi_pflag = 0;
	}
	END_TIMING(TIMER_TYPE_BOUNDS);

	/*******************************************************************************/
	/***  VISCOUS SOURCE STEP *******/
	/*******************************************************************************/
#if( USE_KINEMATIC_VISCOSITY )

	/* Zero out failure array */
	if( failure_exists ) { 
	  failure_exists = 0 ; 
	  for(i=0;i<N_PFLAG_TYPES;i++) { failure_type[i] = 0; } 
	  ALL_LOOP {  pflag[i][j][k] = 0; }
	}

	set_all_visc_stress_functions(2,p);

	/* full-step */
	advance_viscous(p,dx[0]) ;


#if( FIX_CUTOUT_GAMMA || FIX_CUTOUT_PRESSURE )
	fixup_global(p);
#endif

	bounds(p,0);         // set boundaries so fixup's interpolation works

	if( failure_exists ) {   // fix/interpolate over bad points 
	  if( failure_type[PFLAG_INTERP_V] )    { fixup_interp_v(p);    }
	  if( failure_type[PFLAG_INTERP_PRIM] ) { fixup_interp_prim(p); }
	}
	
	if( boundary_phys_pflag || boundary_mpi_pflag ) { 
	  bounds(p,1);         // set boundaries using new fixed points
	  boundary_phys_pflag = boundary_mpi_pflag = 0;
	}
#endif

	/* done! */
	t_old = t; 
	t += dx[0] ;
	
	/* Calculate what the next time interval should be : */ 
	dt_old = dx[0]; 
	dx[0]  = timestep();

	BEG_TIMING(TIMER_TYPE_METRIC);
	advance_geometry();  /* Read in the geometry for the next temporal substep (the present one is already set) */
	END_TIMING(TIMER_TYPE_METRIC);

	ncurr = n_beg;   /* set the current time position to be the beginning  */

	TRACE_END;

	return;
}

/******************************************************************************************/
/******************************************************************************************
  advance(): 
 -----------
     -- handles details pertaining to how a time step is performed, e.g. reconstruction, 
         flux calculation, and EOS evaluation; 
     -- timestep here means integrating the primitives from t to (t+Dt);
     -- arguments: 
           p_i = prim. var's at t
           p_h = prim. var's to be used for flux calculation 
           p_f = prim. var's at t+Dt  (on exit)
******************************************************************************************/
void advance(
	double ****p_i,
	double ****p_h,
	double Dt,
	double ****p_f
	)
{
        int i,j,k,l;
        long int ind_U; 
	double dU[NP],U[NP];

	struct of_state  q;
	struct of_geom   *geom_i,  *geom_h,  *geom_f;
	struct of_coord *coords;

  TRACE_BEG;

#if( NO_EVOLUTION )
	return;
#endif

	/****************************************************************
           Interpolate/reconstruct for left/right states at each face 
	*****************************************************************/
	/* find interpolation coefficients for primitive 
	   variable reconstruction.  Assume BCs are set */
	reconstruct_fast( p_h );

	/****************************************************************
          Find the fluxes in all directions : 
	*****************************************************************/
	numerical_flux( p_L, p_R, F );
	fix_flux();

	/****************************************************************
	  Constrain the fluxes to preserve  Div B = 0
	*****************************************************************/
#if(!HYDRO_ONLY)
# if(USE_FLUX_CT_PARA)
	//flux_ct_para( F ) ;
	//	flux_ct_para_fast( F ) ;
	flux_ct_para_fast2( F ) ;
# else 
	flux_ct( F ) ;
# endif
#endif

	/****************************************************************
	  Calculate special cooling functions or source term: 
	*****************************************************************/
#if( USE_COOLING_FUNCTION == 4 )
	setup_corona_cooling(p_h);
#endif
#if( USE_COOLING_FUNCTION == 5 )
	setup_corona_cooling2(p_h, Dt);
#endif

	/****************************************************************
	  Perform the time integration: 
	*****************************************************************/
	M_tot = 0.;
	ind_U = 0; 

	LOOP { 
#if( USE_MASK )
	  if( evol_mask[ncurr][i][j][k] == MASK_NORMAL )  { 
#else 
	  { 
#endif

	  get_geometry(i,j,k,CENT,n_beg,geom_i); 
#if( KEEP_CONSERVED_VARS ) 
	  PLOOP  { U[l] = U_gf[0][i][j][k][l]; }
	  /* calculate metric/coordinate quantities */
	  if( n_substep == 0 ) {
	    geom_h = geom_i; 
	    get_geometry(i,j,k,CENT,n_mid,geom_f);
	    get_coord(   i,j,k,CENT,n_mid,coords);	    
	  } 
	  else { 
	    get_geometry(i,j,k,CENT,n_mid,geom_h);
	    get_geometry(i,j,k,CENT,n_end,geom_f);
	    get_coord(   i,j,k,CENT,n_end,coords);	    
	  }
#else 	  
	  /* calculate conserved quantities */
	  if( n_substep == 0 ) {
	    geom_h = geom_i; 
	    get_geometry(i,j,k,CENT,n_mid,geom_f);
	    get_coord(   i,j,k,CENT,n_mid,coords);	    
	    get_state(  p_i[i][j][k], geom_i, &q );
	    primtoflux( p_i[i][j][k], &q, 0, geom_i, U) ;
	    PLOOP  { U_gf[0][i][j][k][l] = U[l] ; }
	  } 
	  else { 
	    get_geometry(i,j,k,CENT,n_mid,geom_h);
	    get_geometry(i,j,k,CENT,n_end,geom_f);
	    get_coord(   i,j,k,CENT,n_end,coords);	    
	    PLOOP  { U[l] = U_gf[0][i][j][k][l]; }
	  }
#endif
	    
	  source(p_h[i][j][k], geom_h, i, j, k, dU, Dt);
	    
#if( MAKE_STAT && DUMP_ALL_STAT ) 
	  PLOOP U_pre[l][ind_U] = U[l];
#endif
	    
	  /* evolve conserved quantities */
	  PLOOP {
	    U[l] += Dt*(
			- ( F[i+1][j  ][k  ][0][l] - F[i][j][k][0][l] ) * invdx[1] 
			- ( F[i  ][j+1][k  ][1][l] - F[i][j][k][1][l] ) * invdx[2] 
			- ( F[i  ][j  ][k+1][2][l] - F[i][j][k][2][l] ) * invdx[3] 
			+ dU[l] 
			);
	  }

#if( USE_PRESSURE_FLUX_FIX )
	  double ftmp = -Dt * geom_h->g; 
	  U[U1] +=  ftmp * (F2[i+1][j  ][k  ][0] - F2[i][j][k][0]) * invdx[1];
	  U[U2] +=  ftmp * (F2[i  ][j+1][k  ][1] - F2[i][j][k][1]) * invdx[2];
	  U[U3] +=  ftmp * (F2[i  ][j  ][k+1][2] - F2[i][j][k][2]) * invdx[3];
#endif
	    

	  PLOOP  { U_gf[n_U][i][j][k][l] = U[l]; } 


	  /**************************************************************************
               Find primitive var's from conserved var's, taking into special 
               conditions, most fixups (e.g. floor and GAMMAMAX), inversion checks...
	  **************************************************************************/
#if( USE_SIMPLE_EOS )
	  recover_primitives_simple(i,j,k, p_f[i][j][k], U, geom_f, coords);
#else 
	  recover_primitives(i,j,k, p_f[i][j][k], U, geom_f, geom_i, coords );
#endif

	  /* old-way, not conserved mass:  	  M_tot +=  geom_f->g * p_f[i][j][k][RHO];    */
	  M_tot +=  U[0] ;


#if( MAKE_STAT && DUMP_ALL_STAT ) 
	  if( n_substep == (N0-1) )  {  PLOOP U_out[l][ind_U] = U[l]; }
#endif
#if( MAKE_STAT2 ) 
	  S_stat2[0][ind_U] = dU[U1]; S_stat2[1][ind_U] = dU[U2]; 
#endif

	  } /* end either extra bracket or evol_mask condition */

	  /****************************************************************************************
            EVOLVE MAGNETIC FIELD IN EXCISION REGION: 
	    -- If evol_mask[] has a non-normal state, then we only update the induction equation : 
            -- still don't keep STAT data in excised region;
            -- don't need to calculate the source as the induction equation is always zero as far 
               as we are concerned for ow; 
	  ****************************************************************************************/
#if( USE_MASK && (!HYDRO_ONLY) )
	  else { 

# if( KEEP_CONSERVED_VARS ) 
	    BLOOP  { U[l] = U_gf[0][i][j][k][l]; }
# else 
	    /* calculate conserved quantities */
	    get_geometry(i,j,k,CENT,n_beg,geom_i); 
	    if( n_substep == 0 ) {
	      get_geometry(i,j,k,CENT,n_mid,geom_f);
	      BLOOP  { U_gf[0][i][j][k][l] = U[l] = p_i[i][j][k][l] * geom_i->g;  }
	    } 
	    else { 
	      get_geometry(i,j,k,CENT,n_end,geom_f);
	      BLOOP  { U[l] = U_gf[0][i][j][k][l]; }
	    }
# endif
	    /* evolve conserved quantities */
	    BLOOP {
	      U[l] += Dt*(
			  - ( F[i+1][j  ][k  ][0][l] - F[i][j][k][0][l] ) * invdx[1] 
			  - ( F[i  ][j+1][k  ][1][l] - F[i][j][k][1][l] ) * invdx[2] 
			  - ( F[i  ][j  ][k+1][2][l] - F[i][j][k][2][l] ) * invdx[3] 
			  );

	      U_gf[n_U][i][j][k][l] = U[l]; 
	      p_f[i][j][k][l] = U[l] * geom_f->g_inv;
	    }
	  }
#endif /* ( USE_MASK && (!HYDRO_ONLY) ) */

#if( (MAKE_STAT && DUMP_ALL_STAT) || MAKE_STAT2 ) 
	  ind_U++;   /* need to advance this even when we are excised */
#endif
	    
      }


  TRACE_END;

  return;

	/* done! */
}

/******************************************************************************************/
/******************************************************************************************
  advance_viscous(): 
 -----------
     -- adjusts the conserved variables to take into account the kinematic viscosity
     -- viscosity terms enter the EOM through the source 
     -- this step is only performed during the corrector step;
******************************************************************************************/
void advance_viscous(
		     double ****prim,
		     double Dt
		     )
{
#if( USE_KINEMATIC_VISCOSITY )
        int i,j,k,l;
	double U[NP];

	struct of_state  q;
	struct of_geom   *geom;
	struct of_coord *coords;

  TRACE_BEG;


	/****************************************************************
	  Calculate any viscous stresses for dissipative source terms: 
	*****************************************************************/
        calc_visc_stress_source();


	/****************************************************************
	  Perform the time integration: 
	*****************************************************************/
#if( USE_MASK )
	LOOP if( evol_mask[ncurr][i][j][k] == MASK_NORMAL )  { 
#else 
	LOOP { 
#endif

	  get_geometry(i,j,k,CENT,ncurr,geom); 
	  get_coord(   i,j,k,CENT,ncurr,coords);	    

#if( KEEP_CONSERVED_VARS ) 
	  PLOOP  { U[l] = U_gf[n_U][i][j][k][l]; }
#else 

	  /* calculate conserved quantities */
	  get_state(  prim[i][j][k], geom, &q );
	  primtoflux( prim[i][j][k], &q, 0, geom, U) ;
#endif
	    
	    
	  /* evolve conserved quantities */
	  for(l=0; l<NDIM; l++) {
	    U[l+UU] += Dt*( visc_source[i][j][k][l] );
	  }

	  PLOOP  { U_gf[n_U][i][j][k][l] = U[l]; } 

	    
	  /**************************************************************************
               Find primitive var's from conserved var's, taking into special 
               conditions, most fixups (e.g. floor and GAMMAMAX), inversion checks...
	  **************************************************************************/
#if( USE_SIMPLE_EOS )
	  recover_primitives_simple(i,j,k, prim[i][j][k], U, geom, coords);
#else 
	  recover_primitives(i,j,k, prim[i][j][k], U, geom,geom, coords);
#endif


	}


  TRACE_END;
#endif 

  return;

	/* done! */
}

/******************************************************************************************/
/******************************************************************************************
  timestep(): 
 -----------
     -- calculates the next time interval to integrate over;
     -- it is set in order to resolve fastest wave speeds present at current time step;
     -- prevents integration past final time;
     -- uses the global maximum wave speeds instead of the minimum time step since 
        the way we calculate the timestep here is not independent on how we decompose 
        the domains;
******************************************************************************************/
static double timestep(void)
{
  double dtsqij,dtt ;
  double cmax,cmin ;
  double c1, c2, c3;
  double dt_tmp;

  TRACE_BEG;

#if( NO_EVOLUTION ) 
  return(dx[0]);
#endif

#if( TIMESTEP_METHOD == TIMESTEP_ORIGINAL ) 
  /* Get the fastest wave speeds over all domains : */ 
  mpi_global_vmax(NDIM, max_char_speed);

  c1 = max_char_speed[1]; 
  c2 = max_char_speed[2]; 
  c3 = max_char_speed[3]; 

  dtsqij = sqrt( (1./3.)/(
		    (c1*c1)/(dx[1]*dx[1]) + 
		    (c2*c2)/(dx[2]*dx[2]) +
		    (c3*c3)/(dx[3]*dx[3]) 
		    ) 
		 );
#endif


#if( TIMESTEP_METHOD == TIMESTEP_MIN_CROSSING_TIME )
  dtsqij = dt_char_min;  /* note that the synchronization occurs at the bottom of this routine automatically */
#endif


#if( USE_LIGHT_SPEED ) 
  dt_tmp = MIN(dtmin_glob,dtsqij);
#else 
  dt_tmp = dtsqij;
#endif


#if( EQUATORIAL_RUN )
  static int local_first_time = 1; 
  if( local_first_time ) { 
    fprintf(stdout,"timestep():  Note that you need to set the timestep here to an appropriate value as equatorial this routine does not find correct dt in this case: line %d  of  %s  \n", __LINE__,__FILE__);
    fflush(stdout);
    local_first_time = 0;
  }

  static double dtt0 =  1.e-2; 
#if( BBH_SPACETIME && LARGE_CUTOUT )
  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  struct of_bbh_traj bbh_traj; 
  get_bbh_traj_data(&bbh_traj) ;
  dtt = dtt0 * pow( bbh_traj.r12 / initial_bbh_separation, 1.5) ;
#else
  dtt = dtt0; 
#endif // LARGE_CUTOUT

#else
  dtt = cour*dt_tmp ;
#endif


  /* Make the time step change smoothly : */
  if( dtt > SAFE*dx[0] ) { dtt = SAFE * dx[0] ; } 

  /* don't step beyond end of run */
  if(t + dtt > (GridLength[0]+startx[0]) )  dtt = GridLength[0] + startx[0] - t ;

  if (min_Dt_IC_cool < dtt) {
      fprintf(stdout, "min_Dt_IC_cool = %e < dtt\n", min_Dt_IC_cool); fflush(stdout);
//    dtt = min_Dt_IC_cool;
  }

  /* Get the smallest dt over all domains (to eliminate the possibility of differences at round-off level): */
  mpi_global_min( &dtt );

  TRACE_END;
  
  return(dtt) ;
}

/******************************************************************************************/
/******************************************************************************************
  numerical_flux(): 
 ------------------
     -- Using the reconstructed prim. variables, calculates the numerical flux used 
        when integrating the equations of motion;
     -- Responsible for calculating wave speeds and such needed for the 
        flux formula that is used here;
     -- Currently uses HLLE or LF-like flux formula 
     -- Assumes that dir = 0,1,2 = x1,x2,x3
******************************************************************************************/
static void numerical_flux( 
		    double *****p_L,
		    double *****p_R,
		    double *****F
		    )
{
  unsigned int i,j,k,l,i_s,j_s,k_s,i_e,j_e,k_e,d,dim,pos,ind;
  double Ul[NP],Ur[NP],Fl[NP],Fr[NP],btmp ;
  double *pl, *pr;
  double cminl, cminr, cmaxl, cmaxr, cmax, cmin, ctop;
  struct of_state  ql, qr;
  struct of_geom   *geom;

  TRACE_BEG;

  /* Reset max speeds for timestep calculation: */
  DLOOP1 max_char_speed[i] = 0. ; 
  dt_char_min = 1.e200;


  FACE_LOOP {

#if( USE_FLUX_CT_PARA ) 
    i_s = N1S + NG*(idir[d]-1);  j_s = N2S + NG*(jdir[d]-1);  k_s = N3S + NG*(kdir[d]-1);
    i_e = N1E + 1 + (NG-1)*(1-idir[d]);  j_e = N2E + 1 + (NG-1)*(1-jdir[d]);  k_e = N3E + 1 + (NG-1)*(1-kdir[d]);
#else 
    i_s = N1S + idir[d] - 1;  j_s = N2S + jdir[d] - 1;  k_s = N3S + kdir[d] - 1;
    i_e = N1E+1;  j_e = N2E+1;  k_e = N3E+1;
#endif


    dim = d + 1;  // Position in a 4-vector to which this direction corresponds
    pos = d;      // Position of the flux in the cell, assumes FACE1-3 = 0-2

    ind = 0;

    for(i=i_s;i<=i_e;i++)
      for(j=j_s;j<=j_e;j++) {
	for(k=k_s;k<=k_e;k++) {
//#if(USE_MASK) 
#if(0)
	if( evol_mask[ncurr][i][j][k] != MASK_EXCISED ) {
#else
	{
#endif
	    get_geometry(i,j,k,pos,ncurr,geom);  

	    /* Make local copy for efficiency's sake:  */ 
	    //-fast	  PLOOP pl[l] = p_L[i][j][k][d][l];
	    //-fast	  PLOOP pr[l] = p_R[i][j][k][d][l];
	    pl    = p_L[i][j][k][d];
	    pr    = p_R[i][j][k][d];

	    //	    /* Ensure that the B-field component perp. to the face is continuous: */
	    //	    btmp = 0.5*(pl[B1+d] + pr[B1+d]);
	    //	    pl[B1+d] = pr[B1+d] = btmp;

	    /* Calculate the fluxes and conserved variables : */
	    get_state( pl, geom, &ql );
	    get_state( pr, geom, &qr );
	    primtoflux( pl, &ql,   0, geom, Ul) ;
	    primtoflux( pr, &qr,   0, geom, Ur) ;
	    primtoflux( pl, &ql, dim, geom, Fl) ;
	    primtoflux( pr, &qr, dim, geom, Fr) ;

#if( USE_PRESSURE_FLUX_FIX )
	    Fl[U1+d] -= ql.ptot * geom->g; 
	    Fr[U1+d] -= qr.ptot * geom->g; 
	    F2[i][j][k][d] = 0.5*(ql.ptot + qr.ptot);
#endif

	    /* calculate maximum and minimum speeds */
	    //-fast	  vchar(pl, &ql, geom, dim, &cmaxl, &cminl);
	    //-fast	  vchar(pr, &qr, geom, dim, &cmaxr, &cminr);
	    vchar_fast(&ql, geom, dim, &cmaxl, &cminl);
	    vchar_fast(&qr, geom, dim, &cmaxr, &cminr);
	    cmax = fabs(MAX(MAX(0.,cmaxl),cmaxr)) ;
	    cmin = fabs(MAX(MAX(0.,-cminl),-cminr)) ;
	    ctop = MAX( cmin, cmax ) ; 

#if( TIMESTEP_METHOD == TIMESTEP_ORIGINAL ) 
	    if( ctop > max_char_speed[dim] ) 
#if(USE_MASK) 
	      if( evol_mask[ncurr][i][j][k] == MASK_NORMAL )
#endif
	      if( i >= N1S  && i <= N1E ) 
		if( j >= N2S  && j <= N2E ) 
		  if( k >= N3S  && k <= N3E ) { max_char_speed[dim] = ctop; } 
#endif

#if( TIMESTEP_METHOD == TIMESTEP_MIN_CROSSING_TIME )
	    if( ctop ) 
#if(USE_MASK) 
	      if( evol_mask[ncurr][i][j][k] == MASK_NORMAL )
#endif
	      if( i >= N1S  && i <= N1E ) 
		if( j >= N2S  && j <= N2E ) 
		  if( k >= N3S  && k <= N3E ) { 
		    btmp = dx[dim]/(ctop);
		    dt_char_min = MIN( dt_char_min ,  btmp  ) ;  
		  }
#endif	  

#if( OUTPUT_MAX_VCHAR )
	    max_vchar[i][j][k][d] = ctop; 
#endif

#if( USE_ENTROPY_EQ ) 
	    c_gf[i][j][k][d][0] = cmin;
	    c_gf[i][j][k][d][1] = cmax;
	    ucon_L[i][j][k][d][0] = geom->g * ql.ucon[0];   // LEFT  u^t     
	    ucon_L[i][j][k][d][1] = geom->g * ql.ucon[dim]; // LEFT  u^{face}
	    ucon_R[i][j][k][d][0] = geom->g * qr.ucon[0]  ; // RIGHT u^t     
	    ucon_R[i][j][k][d][1] = geom->g * qr.ucon[dim]; // RIGHT u^{face}
#endif

	    /* HLLE flux: */
#if( FLUX_TYPE_CHOICE == FLUX_HLLE ) 	
	    PLOOP {
	      F[i][j][k][d][l] = (
				  cmax*Fl[l] + cmin*Fr[l]
				  - cmax*cmin*(Ur[l] - Ul[l])
				  )/(cmax + cmin + SMALL) ;
	    }
#elif( FLUX_TYPE_CHOICE == FLUX_LF ) 	
	    /* Lax-Friedrichs flux:  */
	    PLOOP { 
	      F[i][j][k][d][l] = 0.5*(Fl[l] + Fr[l] - ctop*(Ur[l] - Ul[l])) ;
	    }
#endif 

#if( MAKE_STAT2 ) 
	    if( (i>=N1S) && (i<=N1E) && (j>=N2S) && (j<=N2E) && (k>=N3S) && (k<=N3E) ) { 
	      ind = ((k-NG) + (N3)*((j-NG) + (N2)*((i-NG)))); 
	      PLOOP { 
		FL_stat2[d][l][ind] = Fl[l];
		FR_stat2[d][l][ind] = Fr[l];
		UL_stat2[d][l][ind] = Ul[l];
		UR_stat2[d][l][ind] = Ur[l];
		ctop_stat2[d][ind]  = ctop;
	      }
	    }
#endif	  

	}  /* extra bracket or evol_mask conditional */
       } /* k  loop */
     }   /* ij loop */

  }   /* direction loop */ 


  TRACE_END;

  return;
}

/******************************************************************************************/
/******************************************************************************************
  recover_primitives():
 ------------------
     -- responsible for returning with a reasonable set of primitive variables; 
     -- first performs the inversion, then sees if they are valid (e.g. good
        exit status or that gamma is not too large)
     -- using gamma > GAMMAMAX2  for unphysicality test;
******************************************************************************************/
static void recover_primitives(int i, int j, int k, double *pf, double *U, 
			       struct of_geom *geom , struct of_geom *geom_old, struct of_coord *coords  )
{
  int l ;
  int ret,inv_err, ret_ee, use_ee=0;
  double pg[NP], gamma,bsq,pf2[NP];
  double Unew[NP],rd_err[NP];
  struct of_state qf;
  
  static int inv_count=0, prim_count=0;

  //  TRACE_BEG;

  /* Do not reset B-fields in prim vars (since they evolved fine) if FIXUP_TREE=3  */
#if( FIXUP_TREE == 3 )
# define REC_LOOP FLOOP
#else 
# define REC_LOOP PLOOP
#endif



  PLOOP { pg[l] = pf[l]; }

  
  /***************************************************************************************
    Invert  U  for  P   and make sure that the solution is valid: 
       -- try out alternative inversion schemes if we first fail 
       -- make sure that the inversion schemes return with a reasonable value of gamma 
   ***************************************************************************************/
  //-fast	    ret = Utoprim_2d(U, geom.gcov, geom.gcon, geom.g, p_f[i][j][k]);
  ret = Utoprim_2d_fast(U, geom, pf, &gamma, &bsq);

  /* Sometimes this inverter can return with what it thinks is a solution w/ gamma>1e4 but isn't */
  if( (!ret) && (check_gamma(gamma)==2) ) { /* check here so we can try for another solution: */
    REC_LOOP { pf[l] = pg[l]; }   
    ret = 1; 
    gamma = 1.;
  }   

  if( ret ) { 
    ret = Utoprim_1d(U, geom->gcov, geom->gcon, geom->g, pf, &gamma,&bsq);
    if( (!ret) && (check_gamma(gamma)==2) ) { 
      REC_LOOP { pf[l] = pg[l]; }   
      ret = 1; 
    } 
    if( ret ) { 
      gamma = 1.; /* Don't use bad gamma */
      //      PLOOP { fprintf(stderr,"BAD INVERSION: U[%1d] = %26.16e \n",l,U[l]); }
    }
  }

  p_gamma[i][j][k] = gamma;

#if( FIXUP_TREE == 1 ) 
  /* Check to see if everything is fine: */
  if( !ret ) { 
    if( !check_entropy_eq(coords->x[TH],pf[UU],bsq,pf[RHO]) ) { 
      if( !check_floor(pf,coords) ) { 
	if( !check_Tmax(pf) ) { 
	  if( !check_gamma(gamma) ) { 
	    fail(FAIL_USE_FULL_INV,0);
	    //	    TRACE_END;
	    return; 
	  }
	}
      }
    }
  }
# if( USE_ENTROPY_EQ ) 
  use_ee = 1;
# endif 
#endif 


#if( (FIXUP_TREE == 2) || (FIXUP_TREE == 3) ) 
# if( USE_ENTROPY_EQ ) 
  if( ret  || check_entropy_eq(coords->x[TH],pf[UU],bsq,pf[RHO]) ) {  use_ee = 1 ; } 
# endif
#endif
    

  /* Entropy fix : */

  if( !use_ee ) { 
    if( !ret ) { 
      fail(FAIL_USE_FULL_INV,0);
    }
    else { 
       fail( FAIL_UTOPRIM,1 ) ;
    }
    //    fprintf(stdout,"recover1 [%3d,%3d,%3d]   %6d %6d \n", i,j,k,ret, use_ee); fflush(stdout);
  }
  else { 
//    if( !ret ) { 
//      fprintf(stdout,"eetest(%d,%d,%d,%d) Pf : ",
//	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
//      PLOOP { fprintf(stdout," %26.16e",pf[l]); }
//      fprintf(stdout,"\n");
//    }
    if( !ret ) {  PLOOP pf2[l] = pf[l]; }
    else       {  fail( FAIL_UTOPRIM,0 ) ;}

    REC_LOOP { pf[l]  = pg[l]; }   
    ret_ee = fixup_entropy_eq(i,j,k,p_old[i][j][k],pf,geom,geom_old,&gamma);
    if( ret_ee ) { 
      if( ret ) { 
	fail( FAIL_UTOPRIM_EE, 1 ) ;  
	//	TRACE_END;
	return;
      }
      else { 
	fail( FAIL_UTOPRIM_EE, 0 ) ;  
	PLOOP {  pf[l] = pf2[l] ; }
	gamma = p_gamma[i][j][k];
	fail(FAIL_USE_FULL_INV,0);
      }
    }
  }

  //  fprintf(stdout,"recover2 [%3d,%3d,%3d]   %6d %6d \n", i,j,k,ret, use_ee); fflush(stdout);
  p_gamma[i][j][k] = gamma;

  /***************************************************************************************
    Output cases in which inversion yields large errors in U : 
  ****************************************************************************************/
#if( CHECK_INVERSION )
  if( !ret && (inv_count < 1000) ) { 
    inv_err = 0;
    get_state(  pf, geom, &qf );
    primtoflux( pf, &qf, 0, geom, Unew) ;
    HLOOP { 
      if( fabs(U[l]) > 1.e-14 ) { 
	rd_err[l] = fabs((Unew[l] - U[l])/U[l]); 
	if( rd_err[l] > INVERSION_THRESH ) { inv_err = 1; }
      }
      else { 
	rd_err[l] = 0.;
      }
    }
    if( inv_err ) { 
      inv_count++;
      fprintf(stdout,"INVERR(%d,%d,%d,%d) Pg   : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pg[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"INVERR(%d,%d,%d,%d) Pout : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pf[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"INVERR(%d,%d,%d,%d) U    : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",U[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"INVERR(%d,%d,%d,%d) Unew : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",Unew[l]); }
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
#endif


  /***************************************************************************************
    Output cases when any one of the primitive variables changes by a large amount : 
  ***************************************************************************************/
#if( CHECK_PRIM_CHANGE )
  if( !ret && (prim_count < 1000) ) { 
    inv_err = 0;
    HLOOP { 
      if( fabs(pg[l]) > 1.e-14 ) { 
	rd_err[l] = fabs((pf[l] - pg[l])/pg[l]); 
	if( rd_err[l] > PRIM_CHANGE_THRESH ) { inv_err++; }
      }
      else { 
	rd_err[l] = 0.;
      }
    }
    if( inv_err > 1 ) { 
      prim_count++;
      fprintf(stdout,"PRIMERR(%d,%d,%d,%d) U    : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",U[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"PRIMERR(%d,%d,%d,%d) Pg   : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pg[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"PRIMERR(%d,%d,%d,%d) Pout : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pf[l]); }
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
#endif


  /***************************************************************************************
    Now set the floor and GAMMAMAX 
  ***************************************************************************************/
  fixup_floor(pf,coords); 
  if( check_Tmax(pf)     )  {  fixup_Tmax(pf) ; }
  if( check_gamma(gamma) )  {  
    fail(FAIL_GAMMA_MAX,2);
    //    TRACE_END;
    return;
  }

  /* Store gamma for later use : */
  p_gamma[i][j][k] = gamma; 


#undef REC_LOOP

  //  TRACE_END;
  return;
}

/******************************************************************************************/
/******************************************************************************************
  recover_primitives_simple():
 ------------------
     -- like recover_primitives(), but uses the "fix" prim. var. solvers that assume a 
        simple EOS like isothermal or isentropic;
     -- note that you have to set the macros in "u2p_defs.h" to specify which simple EOS to use;
     -- responsible for returning with a reasonable set of primitive variables; 
     -- first performs the inversion, then sees if they are valid (e.g. good
        exit status or that gamma is not too large)
     -- using gamma > GAMMAMAX2  for unphysicality test;
******************************************************************************************/
static void recover_primitives_simple(int i, int j, int k, double *pf, double *U, 
				      struct of_geom *geom, struct of_coord *coords)
{
  int l ;
  int ret,inv_err, ret_ee, use_ee=0;
  double pg[NP], gamma,bsq,pf2[NP];
  double Unew[NP],rd_err[NP];
  struct of_state qf;
  
  static int inv_count=0, prim_count=0;

  //  TRACE_BEG;

  /* Do not reset B-fields in prim vars (since they evolved fine) if FIXUP_TREE=3  */
#if( FIXUP_TREE == 3 )
# define REC_LOOP FLOOP
#else 
# define REC_LOOP PLOOP
#endif



  PLOOP { pg[l] = pf[l]; }

  
  /***************************************************************************************
    Invert  U  for  P   and make sure that the solution is valid: 
       -- try out alternative inversion schemes if we first fail 
       -- make sure that the inversion schemes return with a reasonable value of gamma 
   ***************************************************************************************/
  //-fast	    ret = Utoprim_2d(U, geom.gcov, geom.gcon, geom.g, p_f[i][j][k]);
  ret = Utoprim_1dfix1(U, geom, pf, &gamma, &bsq,Katm[i]);

  /* Sometimes this inverter can return with what it thinks is a solution w/ gamma>1e4 but isn't */
  if( (!ret) && (check_gamma(gamma)==2) ) { /* check here so we can try for another solution: */
    REC_LOOP { pf[l] = pg[l]; }   
    ret = 1; 
    gamma = 1.;
  }   

  if( ret ) { 
    ret = Utoprim_1dvsq2fix1(U, geom, pf, &gamma,&bsq,Katm[i]);
    if( (!ret) && (check_gamma(gamma)==2) ) { 
      REC_LOOP { pf[l] = pg[l]; }   
      ret = 1; 
    } 
    if( ret ) { 
      gamma = 1.; /* Don't use bad gamma */
      //      PLOOP { fprintf(stderr,"BAD INVERSION: U[%1d] = %26.16e \n",l,U[l]); }
    }
  }

  p_gamma[i][j][k] = gamma;

#if( FIXUP_TREE == 1 ) 
  /* Check to see if everything is fine: */
  if( !ret ) { 
    if( !check_floor(pf,coords) ) { 
      if( !check_Tmax(pf) ) { 
	if( !check_gamma(gamma) ) { 
	  fail(FAIL_USE_FULL_INV,0);
	  //	  TRACE_END;
	  return; 
	}
      }
    }
  }
#endif 


  if( !ret ) { 
    fail(FAIL_USE_FULL_INV,0);
  }
  else { 
    fail( FAIL_UTOPRIM,1 ) ;
  }

  p_gamma[i][j][k] = gamma;

  /***************************************************************************************
    Output cases in which inversion yields large errors in U : 
  ****************************************************************************************/
#if( CHECK_INVERSION )
  if( !ret && (inv_count < 1000) ) { 
    inv_err = 0;
    get_state(  pf, geom, &qf );
    primtoflux( pf, &qf, 0, geom, Unew) ;
    HLOOP { 
      if( fabs(U[l]) > 1.e-14 ) { 
	rd_err[l] = fabs((Unew[l] - U[l])/U[l]); 
	if( rd_err[l] > INVERSION_THRESH ) { inv_err = 1; }
      }
      else { 
	rd_err[l] = 0.;
      }
    }
    if( inv_err ) { 
      inv_count++;
      fprintf(stdout,"INVERR(%d,%d,%d,%d) Pg   : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pg[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"INVERR(%d,%d,%d,%d) Pout : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pf[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"INVERR(%d,%d,%d,%d) U    : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",U[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"INVERR(%d,%d,%d,%d) Unew : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",Unew[l]); }
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
#endif


  /***************************************************************************************
    Output cases when any one of the primitive variables changes by a large amount : 
  ***************************************************************************************/
#if( CHECK_PRIM_CHANGE )
  if( !ret && (prim_count < 1000) ) { 
    inv_err = 0;
    HLOOP { 
      if( fabs(pg[l]) > 1.e-14 ) { 
	rd_err[l] = fabs((pf[l] - pg[l])/pg[l]); 
	if( rd_err[l] > PRIM_CHANGE_THRESH ) { inv_err++; }
      }
      else { 
	rd_err[l] = 0.;
      }
    }
    if( inv_err > 1 ) { 
      prim_count++;
      fprintf(stdout,"PRIMERR(%d,%d,%d,%d) U    : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",U[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"PRIMERR(%d,%d,%d,%d) Pg   : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pg[l]); }
      fprintf(stdout,"\n");

      fprintf(stdout,"PRIMERR(%d,%d,%d,%d) Pout : ",
	      nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
      PLOOP { fprintf(stdout," %26.16e",pf[l]); }
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
#endif


  /***************************************************************************************
    Now set the floor and GAMMAMAX 
  ***************************************************************************************/
  fixup_floor(pf,coords); 
  if( check_Tmax(pf)     )  {  fixup_Tmax(pf) ; }
  if( check_gamma(gamma) )  {  
    fail(FAIL_GAMMA_MAX,2);
    //    TRACE_END;
    return;
  }

  /* Store gamma for later use : */
  p_gamma[i][j][k] = gamma; 


#undef REC_LOOP

  //  TRACE_END;
  return;
}

/******************************************************************************************/
/******************************************************************************************
  check_for_crash():
 ------------------
     -- Probes the state of the simulation (gridfunctions, global variables) to see if 
        the simulation has effectively crashed; 

     -- Since the fixup routines are written to continue on when bad states are encountered,
        the code does not have an explicit condition on when to give up;  this routine 
        is meant to serve this purpose. 
 
     -- The types of things that indicate a "crash" are the following: 
 
          1) a sudden change in E_tot (amount controlled by "reldiff_E_bad") 

          2) when the number of the subdomain's cells that report inversion failures 
             since the last timestep exceeds some fraction ("nfail_fraction") of 
             the subdomain's total cell count;

        when integrating the equations of motion;
     -- Responsible for calculating wave speeds and such needed for the 
        flux formula that is used here;
     -- Currently uses HLLE or LF-like flux formula 
     -- Assumes that dir = 0,1,2 = x1,x2,x3
******************************************************************************************/
void check_for_crash( void )
{
  int i,j,k;
  unsigned long int n_bad;

  double fraction_bad;

  unsigned long int ncells = (N0*NCELLS); /* Number of possible instances of failure */
  static unsigned long int ntot1,ntot2, ntot1_old, ntot2_old,dntot1_old,dntot2_old,dntot1,dntot2;

  const double reldiff_E_bad = 1.e-3;  /* Rel. amount that E_tot must change to show crash */
  const double nfail_fraction = 0.5;   /* Fraction of domain's cells that must "fail" */ 

  static double E_tot_old = -1.;
  static int    first_time = 1; 
  
  TRACE_BEG;

  /*********************************************************************************
    Check to see if E_tot is suddenly increasing : 
  ***********************************************************************************/
//--bbh  if( myid == master_pid ) { 
//--bbh    if( first_time ) { 
//--bbh      E_tot_old = E_tot; 
//--bbh    }
//--bbh    else { 
//--bbh      if( REL_DIFF_FUNC(E_tot_old,E_tot) > reldiff_E_bad ) { 
//--bbh	fprintf(stderr,"\n\n#############################################################\n");
//--bbh	fprintf(stderr,"\t <<<<---- check_for_crash():  ---->>>> \n"); 
//--bbh	fprintf(stderr,"  Detecting sudden increase in E_tot :\n"); 
//--bbh	fprintf(stderr,"       nstep         = %d      \n", nstep        ); 
//--bbh	fprintf(stderr,"       t             = %g      \n", t            ); 
//--bbh	fprintf(stderr,"       E_old         = %28.18e \n", E_tot_old    ); 
//--bbh	fprintf(stderr,"       E             = %28.18e \n", E_tot        ); 
//--bbh	fprintf(stderr,"       reldiff_E_bad = %28.18e \n", reldiff_E_bad); 
//--bbh	fprintf(stderr,"\n\n Triggering exit.... \n");
//--bbh	fprintf(stderr,"#############################################################\n\n");
//--bbh	fflush(stderr);
//--bbh	fail(FAIL_BASIC,0); 
//--bbh      }
//--bbh      E_tot_old = E_tot;
//--bbh    }
//--bbh  }
  
  /*********************************************************************************
    See if a large fraction of this domains cells have inversion or fixup failures : 
  ***********************************************************************************/
  if( first_time ) { 
    ntot1 = ntot2 = dntot1 = dntot2 = 0; 
  }

  dntot1_old = dntot1; 
  dntot2_old = dntot2; 
  ntot1_old  = ntot1; 
  ntot2_old  = ntot2; 
  //  ntot1 = nfailtot[NFAIL_UTOPRIM]+nfailtot[NFAIL_UTOPRIM_EE];
  ntot1 = nfailtot[NFAIL_UTOPRIM_EE];
  ntot2 = nfailtot[NFAIL_INTERP_PRIM]+nfailtot[NFAIL_INTERP_V];

  /* Need to consider when nfailtot[] is reset */
  if( ntot1 < ntot1_old ) {  dntot1 = 0; }  
  else                    {  dntot1 = ntot1 - ntot1_old; }

  if( ntot2 < ntot2_old ) {  dntot2 = 0; } 
  else                    {  dntot2 = ntot2 - ntot2_old; }

  n_bad = MAX( dntot1 , dntot2 ); 
  fraction_bad = ((double) n_bad) / ((double) ncells); 

  if( fraction_bad  > nfail_fraction ) { 
    fprintf(stderr,"\n\n#############################################################\n");
    fprintf(stderr,"\t <<<<---- check_for_crash():   ---->>>> \n"); 
    fprintf(stderr,"  Reached critical fraction of cells that have failed this time step :\n"); 
    fprintf(stderr,"       nstep          = %d  \n", nstep        ); 
    fprintf(stderr,"       t              = %g  \n", t            ); 
    fprintf(stderr,"       n_bad          = %ld \n", n_bad        ); 
    fprintf(stderr,"       ncells         = %ld \n", ncells       );
    fprintf(stderr,"       ntot1          = %ld \n", ntot1        ); 
    fprintf(stderr,"       ntot1_old      = %ld \n", ntot1_old    ); 
    fprintf(stderr,"       ntot2          = %ld \n", ntot2        ); 
    fprintf(stderr,"       ntot2_old      = %ld \n", ntot2_old    ); 
    fprintf(stderr,"       nfail_fraction = %g \n", nfail_fraction); 
    fprintf(stderr,"\n\n Triggering exit.... \n");
    fprintf(stderr,"#############################################################\n\n");
    fprintf(stderr,"fails = %ld %ld \n", nfailtot[NFAIL_INTERP_PRIM],nfailtot[NFAIL_INTERP_V]);
    fflush(stderr);
    fflush(stderr);
    fail(FAIL_BASIC,0); 
  }

  first_time = 0; 

  TRACE_END;

  return;
}
