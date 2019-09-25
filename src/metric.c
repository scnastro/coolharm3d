/***********************************************************************************
***********************************************************************************/

#include "decs.h"
#include "metric.h"
#include "conn.h"

/* Sets the fraction of the numerical discretization to be used in the numerical derivatives : */
#define DEL_X_CONN  (1.e-2)

/***************************************************************************/
/***************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     G E N E R A L   M E T R I C    R O U T I N E S

 *****************************************************************************/

/* Handles error-checking for cos/sin evaluation: */ 
void my_sincos( double th, double *sth, double *cth)
{
  *sth = sin(th);
  *cth = cos(th);
  return;
}

/* Handles error-checking for cos/sin evaluation: */ 
static void trigvalues( double th, double *cth, double *sth)
{

  sincos(th, sth, cth); 

#if(COORDSINGFIX)
  if (fabs(*sth) < SINGSMALL) {
    if((*sth)>=0) *sth =  SINGSMALL;
    if((*sth)<0)  *sth = -SINGSMALL;
  }
#endif

}

/***************************************************************************
   gdet_func(): 
  ---------------
   -- returns the sqrt( -det(gcov) )  where gcov[][] is already defined; 
***************************************************************************/
static double gdet_func(double gcov[][NDIM]) 
{
  int i;
  int permute[NDIM]; 
  double gcovtmp[NDIM][NDIM];
  double detg;
  int LU_decompose( double A[][NDIM], int permute[] );

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j; 
  DLOOP2 {  gcovtmp[i][j] = gcov[i][j]; }
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
#endif

  if( LU_decompose( gcovtmp,  permute ) != 0  ) { 
    fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
    fail(FAIL_METRIC,0);
  }
  detg = 1.;
  DLOOP1 detg *= gcovtmp[i][i];
  return( sqrt(fabs(detg)) );

}

/****************************************************************************
  setutcon():
     -- find the contravariant time-component of a time-like vector 
        pointing forward in time;
****************************************************************************/
void setutcon(double *vcon, double gcov[][NDIM])
{
  int i,j;
  double d,b,c;
  
  d=gcov[TT][TT];

  b = gcov[TT][1]*vcon[1] + gcov[TT][2]*vcon[2] + gcov[TT][3]*vcon[3] ;

  c = gcov[1][1] * vcon[1] * vcon[1] 
    + gcov[2][2] * vcon[2] * vcon[2] 
    + gcov[3][3] * vcon[3] * vcon[3] 
    + 2.*(   gcov[1][2] * vcon[1] * vcon[2] 
	   + gcov[1][3] * vcon[1] * vcon[3] 
           + gcov[2][3] * vcon[2] * vcon[3] );

  c += 1. ;  /* vector is timelike */

  vcon[0]=(-b-sqrt(b*b-d*c))/(d);   /* sign for pointing forward in time */

  //  fprintf(stdout,"bcd = %10.4e %10.4e %10.4e \n", b,c,d); fflush(stdout);

  return;
}


/***************************************************************************/
/***************************************************************************
  conn_func():
  -----------

   -- calculates the connection coefficient by finite differencing 
      gcov (ala gcov_func() ) : 
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   --  where i = {1,2,3,4} corresponds to {t,r,theta,phi}

   !! could be optimized further by skipping the Killing directions (when 
          taking the derivative of the metric);
            -- and probably other things;

***************************************************************************/
static void conn_func(struct of_coord *coords, struct of_geom *geom, double ***conn)
{
        int i,j,k,l,n,ii,jj,kk ;
	double tmp[NDIM][NDIM][NDIM] ;
	double xp_h[NDIM],xp_l[NDIM] ;
	double gh[NDIM][NDIM] ;
	double gl[NDIM][NDIM] ;
	double delx, half_inv_delx;
	void gcov_func(struct of_coord *coords, double gcov[][NDIM]) ;
	struct of_coord coords_h, coords_l;

	ii = coords->i;
	jj = coords->j;
	kk = coords->k;

	/* Time component: */
	CONN_2ND_ORDER_TIME_DERIVATIVE ; 

	/* Space components: */
	for(k=1;k<NDIM;k++) {
		for(l=0;l<NDIM;l++) xp_h[l] = coords->xp[l] ;
		for(l=0;l<NDIM;l++) xp_l[l] = coords->xp[l] ;
		delx = DEL_X_CONN * dx[k] ;
		xp_h[k] += delx ;
		xp_l[k] -= delx ;
		coord_of_xp( xp_h , &coords_h );
		coord_of_xp( xp_l , &coords_l );
		half_inv_delx = 0.5/delx;
		gcov_func(&coords_h,gh) ;
		gcov_func(&coords_l,gl) ;

		for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {
		    conn[i][j][k] = (gh[i][j] - gl[i][j])*half_inv_delx;
		  }
	}

	DGCOV_TO_CONN ; 

	/* done! */
}

/***************************************************************************/
/***************************************************************************
  conn_func_fd2():
  -----------
   -- calculates the connection using existing metric data to 2nd order;
***************************************************************************/
static void conn_func_fd2(struct of_coord *coords, struct of_geom *geom, double ***conn)
{
        int i,j,k,l,n ;
	int ii,jj,kk;
	double delx, half_inv_delx;
	double tmp[NDIM][NDIM][NDIM] ;
	double gh[NDIM][NDIM],gl[NDIM][NDIM];
	void gcov_func(struct of_coord *coords, double gcov[][NDIM]) ;
	struct of_geom *geom_h, *geom_l;

	ii = coords->i;
	jj = coords->j;
	kk = coords->k;

	/* Time component: */
	CONN_2ND_ORDER_TIME_DERIVATIVE ; 

	/* Space components: */
	k = 1; 
	half_inv_delx = invdx[k];
	get_geometry(ii  ,jj  ,kk  ,FACE1,ncurr,geom_l);
	get_geometry(ii+1,jj  ,kk  ,FACE1,ncurr,geom_h);
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {
	    conn[i][j][k] = (geom_h->gcov[i][j] - geom_l->gcov[i][j]) * half_inv_delx;
	  }

	k = 2; 
	half_inv_delx = invdx[k];
	get_geometry(ii  ,jj  ,kk  ,FACE2,ncurr,geom_l);
	get_geometry(ii  ,jj+1,kk  ,FACE2,ncurr,geom_h);
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {
	    conn[i][j][k] = (geom_h->gcov[i][j] - geom_l->gcov[i][j]) * half_inv_delx;
	  }

	k = 3; 
	half_inv_delx = invdx[k];
	get_geometry(ii  ,jj  ,kk  ,FACE3,ncurr,geom_l);
	get_geometry(ii  ,jj  ,kk+1,FACE3,ncurr,geom_h);
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {
	    conn[i][j][k] = (geom_h->gcov[i][j] - geom_l->gcov[i][j]) * half_inv_delx;
	  }

	DGCOV_TO_CONN ; 
	/* done! */
}

/***************************************************************************/
/***************************************************************************
  conn_func_fd4():
  -----------
   -- calculates the connection using existing metric data to 4th order in 
       space, 2nd order in time (order determined by DEL_X_CONN
***************************************************************************/
static void conn_func_fd4(struct of_coord *coords, struct of_geom *geom, double ***conn)
{
        int i,j,k,l,n ;
	int ii,jj,kk;
	double delx, half_inv_delx;
	double tmp[NDIM][NDIM][NDIM] ;
	double gh[NDIM][NDIM],gl[NDIM][NDIM];
	void gcov_func(struct of_coord *coords, double gcov[][NDIM]) ;
	struct of_geom *geom_hh, *geom_ll, *geom_h, *geom_l;
	double w1, w2; 

	ii = coords->i;
	jj = coords->j;
	kk = coords->k;

	/* Time component: */
	CONN_2ND_ORDER_TIME_DERIVATIVE ; 

	/* Space components: */
	FOURTH_ORDER_METRIC_DERIVATIVE1 ; 
	FOURTH_ORDER_METRIC_DERIVATIVE2 ; 
	FOURTH_ORDER_METRIC_DERIVATIVE3 ; 

	DGCOV_TO_CONN ; 

	/* done! */
}

/***************************************************************************/
/***************************************************************************
  conn_func_fd4_lo2():
  -----------
   -- like conn_func_fd4() but does 4th order backward differencing in 
      x2 direction on the lower surface;
***************************************************************************/
static void conn_func_fd4_lo2(struct of_coord *coords, struct of_geom *geom, double ***conn)
{
        int i,j,k,l,n ;
	int ii,jj,kk;
	double delx, half_inv_delx;
	double tmp[NDIM][NDIM][NDIM] ;
	double gh[NDIM][NDIM],gl[NDIM][NDIM];
	void gcov_func(struct of_coord *coords, double gcov[][NDIM]) ;
	struct of_geom *geom_hh, *geom_ll, *geom_h, *geom_l;
	double w1, w2; 

	ii = coords->i;
	jj = coords->j;
	kk = coords->k;

	/* Time component: */
	CONN_2ND_ORDER_TIME_DERIVATIVE ; 

	/* Space components: */
	FOURTH_ORDER_METRIC_DERIVATIVE1 ; 

	/* Special backward differencing, we must be on the physical x2 boundary */
	k = 2;  
	get_geometry(ii  ,jj  ,kk  ,FACE2,ncurr,geom_ll );
	get_geometry(ii  ,jj+1,kk  ,FACE2,ncurr,geom_l  );
	get_geometry(ii  ,jj+1,kk  ,CENT ,ncurr,geom_h  );
	get_geometry(ii  ,jj+2,kk  ,FACE2,ncurr,geom_hh );
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {
	    conn[i][j][k] = w_2lo[0]*geom_ll->gcov[i][j] + w_2lo[1]*geom->gcov[i][j] + w_2lo[2]*geom_l->gcov[i][j] + w_2lo[3]*geom_h->gcov[i][j] + w_2lo[4]*geom_hh->gcov[i][j]; 
	}

	FOURTH_ORDER_METRIC_DERIVATIVE3 ; 

	DGCOV_TO_CONN ; 
	/* done! */
}


/***************************************************************************/
/***************************************************************************
  conn_func_fd4_hi2():
  -----------
   -- like conn_func_fd4() but does 4th order backward differencing in 
      x2 direction on the upper surface;
***************************************************************************/
static void conn_func_fd4_hi2(struct of_coord *coords, struct of_geom *geom, double ***conn)
{
        int i,j,k,l,n ;
	int ii,jj,kk;
	double delx, half_inv_delx;
	double tmp[NDIM][NDIM][NDIM] ;
	double gh[NDIM][NDIM],gl[NDIM][NDIM];
	void gcov_func(struct of_coord *coords, double gcov[][NDIM]) ;
	struct of_geom *geom_hh, *geom_ll, *geom_h, *geom_l;
	double w1, w2; 

	ii = coords->i;
	jj = coords->j;
	kk = coords->k;

	/* Time component: */
	CONN_2ND_ORDER_TIME_DERIVATIVE ; 

	/* Space components: */
	FOURTH_ORDER_METRIC_DERIVATIVE1 ; 

	/* Special backward differencing, we must be on the physical x2 boundary */
	k = 2;  
	get_geometry(ii  ,jj-1,kk  ,FACE2,ncurr,geom_ll );
	get_geometry(ii  ,jj-1,kk  ,CENT ,ncurr,geom_l  );
	get_geometry(ii  ,jj  ,kk  ,FACE2,ncurr,geom_h  );
	get_geometry(ii  ,jj+1,kk  ,FACE2,ncurr,geom_hh );
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {
	    conn[i][j][k] = w_2hi[0]*geom_ll->gcov[i][j] + w_2hi[1]*geom_l->gcov[i][j] + w_2hi[2]*geom_h->gcov[i][j] + w_2hi[3]*geom->gcov[i][j] + w_2hi[4]*geom_hh->gcov[i][j]; 
	}

	FOURTH_ORDER_METRIC_DERIVATIVE3 ; 

	DGCOV_TO_CONN ; 
	/* done! */
}


/* Lowers a contravariant rank-1 tensor to a covariant one */
void lower(double *ucon, struct of_geom *geom, double *ucov)
{

	ucov[0] = geom->gcov[0][0]*ucon[0] 
		+ geom->gcov[0][1]*ucon[1] 
		+ geom->gcov[0][2]*ucon[2] 
		+ geom->gcov[0][3]*ucon[3] ;
	ucov[1] = geom->gcov[1][0]*ucon[0] 
		+ geom->gcov[1][1]*ucon[1] 
		+ geom->gcov[1][2]*ucon[2] 
		+ geom->gcov[1][3]*ucon[3] ;
	ucov[2] = geom->gcov[2][0]*ucon[0] 
		+ geom->gcov[2][1]*ucon[1] 
		+ geom->gcov[2][2]*ucon[2] 
		+ geom->gcov[2][3]*ucon[3] ;
	ucov[3] = geom->gcov[3][0]*ucon[0] 
		+ geom->gcov[3][1]*ucon[1] 
		+ geom->gcov[3][2]*ucon[2] 
		+ geom->gcov[3][3]*ucon[3] ;

        return ;
}

/* Raises a covariant rank-1 tensor to a contravariant one */
void raise(double *ucov, struct of_geom *geom, double *ucon)
{

	ucon[0] = geom->gcon[0][0]*ucov[0] 
		+ geom->gcon[0][1]*ucov[1] 
		+ geom->gcon[0][2]*ucov[2] 
		+ geom->gcon[0][3]*ucov[3] ;
	ucon[1] = geom->gcon[1][0]*ucov[0] 
		+ geom->gcon[1][1]*ucov[1] 
		+ geom->gcon[1][2]*ucov[2] 
		+ geom->gcon[1][3]*ucov[3] ;
	ucon[2] = geom->gcon[2][0]*ucov[0] 
		+ geom->gcon[2][1]*ucov[1] 
		+ geom->gcon[2][2]*ucov[2] 
		+ geom->gcon[2][3]*ucov[3] ;
	ucon[3] = geom->gcon[3][0]*ucov[0] 
		+ geom->gcon[3][1]*ucov[1] 
		+ geom->gcon[3][2]*ucov[2] 
		+ geom->gcon[3][3]*ucov[3] ;

        return ;
}

/********************************************************************************************
  advance_geometry():
 ----------------
   -- reads in the geometry functions for the next time step and next halfstep;
   -- does something only if DYNAMIC_SPACETIME  is set;
   -- t_beg is the starting time time, t_mid the half-step time, and t_end the final time to which 
      we are integrating;  we need to set t_mid and t_end times; 
   -- As we advance in time, the "new timestep" will become the present time step.  In order 
      to not copy the new data into the old data's space, we can just keep track of which 
      index houses the new and old data and flip this index when we advance in time.  This 
      routine handles that bookkeeping.

********************************************************************************************/
void advance_geometry( void )
{

  TRACE_BEG;

#if(DYNAMIC_SPACETIME) 
  double t_beg, t_mid, t_end; 
  int     n_tmp;
  void set_general_geometry( int n, double t_input ) ;

  t_beg = t;
  t_mid = t_beg + 0.5*dx[0];
  t_end = t_beg + dx[0];

  /* swap t_end -> t_beg,  keep n_mid the same  : */
  n_tmp = n_beg;   n_beg = n_end;   n_end = n_tmp; 

  /* Read in the new n_mid and n_end levels : */
  set_general_geometry(n_mid,t_mid);
  set_general_geometry(n_end,t_end);

#endif

  TRACE_END;
  return;
}

/***************************************************************************
/***********************************************************************
  null_time_funcs_setup():
  -------------
   -- filler routine for those "dynamic" spacetimes that don't use a 
      "*_setup()" routine to set the quantities only dependent on time; 
***********************************************************************/
void null_time_funcs_setup(double t)
{
#if( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW || \
     METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_DROP      || \
     METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND    || \
     METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND || \
     METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_SIMPLE_NEWTONIAN      )
  nz_params[37] = t;  /* at least set the time coordinate */
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )
  bbh_params->t = t;
#else
  nz_params[36] = t;  /* at least set the time coordinate */
#endif
  return;
}

/***********************************************************************/
/***********************************************************************
  risco_calc():
  -------------
   -- returns the radius of the (equatorial) 
      innermost stable circular orbit ; 
   -- if "is_prograde" is non-zero, then return risco of a prograde
      orbit, otherwise return that of a retrograde orbit; 
   -- "a" must already be set;
***********************************************************************/
double risco_calc( int do_prograde )
{
  double Z1,Z2,sign,term1,term2 ;

  sign = (do_prograde) ? 1. : -1. ; 

  term1 = pow(1. + a,1./3.);
  term2 = pow(1. - a,1./3.);
  
  Z1 = 1. + term1*term2*(term1 + term2);

  Z2 = sqrt(3.*asq + Z1*Z1) ;

  return( 3. + Z2-sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2))  );

}
//Dennis add function below
double risco_calc_general( int do_prograde, double spin, double mass )
{
  double Z1,Z2,sign,term1,term2 ;
  double spinsq;

  spinsq = spin*spin;
  sign = (do_prograde) ? 1. : -1. ; 

  term1 = pow(1. + spin,1./3.);
  term2 = pow(1. - spin,1./3.);
  
  Z1 = 1. + term1*term2*(term1 + term2);

  Z2 = sqrt(3.*spinsq + Z1*Z1) ;

  return( mass * (3. + Z2-sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)))  );

}
//Dennis added function above
/***********************************************************************/
/***********************************************************************
  rhorizon_calc(): 
  --------------------
   -- returns the value of the horizon's radius ;
   -- "a" must already be set;
   -- if "pos_sign" is non-zero, then the radius of the outermost
      horizon ("r_+") is returned ;  otherwise the innermost "r_-" one 
      is returned;
***********************************************************************/
double rhorizon_calc(int pos_sign)
{
  double sign;
  
  sign = (pos_sign) ? 1. : -1.; 

  return(1. + sign*sqrt((1.-a)*(1.+a)) );
}


/***************************************************************************/
/***************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     B O Y E R - L I N D Q U I S T   M E T R I C    R O U T I N E S

 *****************************************************************************/
/* sqrt( - det(gcov) ) */ 
void bl_gdet_func( double *x, double *gdet )
{
  double a2,r2 ;
  double cth,sth, cos_2th; 
  double r,th;

  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cth, &sth);
  cos_2th = (cth-sth)*(cth+sth);
  a2 = asq ;
  r2 = r*r ;

  *gdet = r2*fabs(sth)*(1. + 0.5*(a2/r2)*(1. + cos_2th)) ;

  return;
}

/* covariant metric  */
void bl_gcov_func( double *x, double gcov[][NDIM])
{
  int i,j,k ;
  double sth,cth,s2,a2,r2,DD,mu ;
  double r,th;
  double Mass, Spin;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  DLOOP2 {  gcov[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcov[0][i] = 0. ; } 
#endif

  /* Generic Mass and Spin */
  double t1,t2,t3,t4,t5,t7,t8,t13,t14,t18,t24;
  Spin = a;
  if     ( MASS_TYPE == MBH1 ) { Mass = m_bh1; }
  else if( MASS_TYPE == MBH2 ) { Mass = m_bh2; }
  else                         { Mass = M;    }
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cth, &sth);

  t1 = Mass*r;
  t2 = r*r;
  t3 = Spin*Spin;
  t4 = cth;
  t5 = t4*t4;
  t7 = t2+t3*t5;
  t8 = 1/t7;
  t13 = sth;
  t14 = t13*t13;
  t18 = 2.0*Mass*Spin*r*t14*t8;
  t24 = t2*t2;
  gcov[0][0] = -1.0+2.0*t1*t8;
  gcov[0][1] = 0.0;
  gcov[0][2] = 0.0;
  gcov[0][3] = -t18;
  gcov[1][0] = 0.0;
  gcov[1][1] = t7/(t2-2.0*t1+t3);
  gcov[1][2] = 0.0;
  gcov[1][3] = 0.0;
  gcov[2][0] = 0.0;
  gcov[2][1] = 0.0;
  gcov[2][2] = t7;
  gcov[2][3] = 0.0;
  gcov[3][0] = -t18;
  gcov[3][1] = 0.0;
  gcov[3][2] = 0.0;
  gcov[3][3] = t14*t8*(t24+((2.0-t14)*t2+(1.0-t14)*t3)*t3+2.0*t3*t14*t1);
  return;
}

/* covariant metric */
void bl_gcov_func_old( double *x, double gcov[][NDIM])
{
  int i,j,k ;
  double sth,cth,s2,a2,r2,DD,mu ;
  double r,th;
  double Mass, Spin;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  DLOOP2 {  gcov[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcov[0][i] = 0. ; } 
#endif
  /* Assumes M = 1 */
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cth, &sth);

  s2 = sth*sth ;
  a2 = asq ;
  r2 = r*r ;
  DD = 1. - 2./r + a2/r2 ;
  mu = 1. + a2*cth*cth/r2 ;
	
  gcov[TT][TT] = -(1. - 2./(r*mu)) ;
  gcov[TT][3] = -2.*a*s2/(r*mu) ;
  gcov[3][TT] = gcov[TT][3] ;
  gcov[1][1] = mu/DD ;
  gcov[2][2] = r2*mu ;
  gcov[3][3] = r2*s2*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu)) ;

  return;
}

/* contravariant metric  */
void bl_gcon_func( double *x, double gcon[][NDIM])
{
  int i, j,k ;
  double sth,cth,a2,r2,r3,DD,mu ;
  double r,th;
  double Mass,Spin;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  DLOOP2 {  gcon[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcon[0][i] = 0. ; } 
#endif

  /* Generic Mass and Spin */
  Spin = a;
  if     ( MASS_TYPE == MBH1 ) { Mass = m_bh1; }
  else if( MASS_TYPE == MBH2 ) { Mass = m_bh2; }
  else                         { Mass = M;    }
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cth, &sth);

  double t1,t2,t3,t4,t7,t8,t17,t20,t21,t22,t28,t32,t33;
  t1 = r*r;
  t2 = t1*t1;
  t3 = cth;
  t4 = t3*t3;
  t7 = Spin*Spin;
  t8 = t7*t4;
  t17 = 1/(t1+t8);
  t20 = 2.0*Mass*r;
  t21 = -t1+t20-t7;
  t22 = 1/t21;
  t28 = 2.0*Spin*Mass*r*t17*t22;
  t32 = sth;
  t33 = t32*t32;
  gcon[0][0] = -(-t2+((-1.0-t4)*t1-t8)*t7+2.0*(t4-1.0)*r*t7*Mass)*t17*t22;
  gcon[0][1] = 0.0;
  gcon[0][2] = 0.0;
  gcon[0][3] = t28;
  gcon[1][0] = 0.0;
  gcon[1][1] = -t21*t17;
  gcon[1][2] = 0.0;
  gcon[1][3] = 0.0;
  gcon[2][0] = 0.0;
  gcon[2][1] = 0.0;
  gcon[2][2] = t17;
  gcon[2][3] = 0.0;
  gcon[3][0] = t28;
  gcon[3][1] = 0.0;
  gcon[3][2] = 0.0;
  gcon[3][3] = (-t1-t8+t20)*t22*t17/t33;
  return;
}


/* contravariant metric  */
void bl_gcon_func_old( double *x, double gcon[][NDIM])
{
  int i, j,k ;
  double sth,cth,a2,r2,r3,DD,mu ;
  double r,th;
  double Mass,Spin;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  DLOOP2 {  gcon[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcon[0][i] = 0. ; } 
#endif

  /* Assumes M = 1 */
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cth, &sth);

  a2 = asq ;
  r2 = r*r ;
  r3 = r2*r ;
  DD = 1. - 2./r + a2/r2 ;
  mu = 1. + a2*cth*cth/r2 ;

  gcon[TT][TT] = -1. - 2.*(1. + a2/r2)/(r*DD*mu) ;
  gcon[TT][3] = -2.*a/(r3*DD*mu) ;
  gcon[3][TT] = gcon[TT][3] ;
  gcon[1][1] = DD/mu ;
  gcon[2][2] = 1./(r2*mu) ;
  gcon[3][3] = (1. - 2./(r*mu))/(r2*sth*sth*DD) ;
  return;
}

/* connection   */
void bl_conn_func( struct of_coord *coords, double ***connp)
{
  double r  = coords->r;
  double th = coords->x[TH];

  /* General Mass and Spin */
  double Mass, Spin;
  double cth, sth;

  Spin = a;
  if     ( MASS_TYPE == MBH1 ) { Mass = m_bh1; }
  else if( MASS_TYPE == MBH2 ) { Mass = m_bh2; }
  else                         { Mass = M;     }

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  trigvalues( th, &cth, &sth);
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
  double t17,t18,t22,t23,t25,t26,t29;
  double t36,t38,t39,t40,t41,t42,t45,t48;
  double t51,t53,t59,t62,t70,t73,t74;
  double t86,t92,t93,t95,t97;
  double t102,t103,t109,t112,t128,t131;
  double t144,t148,t151,t159,t166;
  double conn[NDIM][NDIM][NDIM];
  t1 = r*r;
  t2 = Spin*Spin;
  t3 = cth;
  t4 = t3*t3;
  t5 = t2*t4;
  t6 = -t1+t5;
  t7 = Mass*t6;
  t8 = t1+t2;
  t9 = t1*t1;
  t10 = t9*t1;
  t17 = t4*t4;
  t18 = t2*t17;
  t22 = ((-2.0*t4-1.0)*t9+((-2.0-t4)*t4*t1-t18)*t2)*t2;
  t23 = t9*r;
  t25 = t1*r;
  t26 = t4*t25;
  t29 = t17*r*t2;
  t36 = 1/(-t10+t22+(2.0*t23+(4.0*t26+2.0*t29)*t2)*Mass);
  t38 = t7*t8*t36;
  t39 = t2*Mass;
  t40 = t39*r;
  t41 = sth;
  t42 = t3*t41;
  t45 = 2.0*t1*t4+t18;
  t48 = 1/(t9+t45*t2);
  t51 = 2.0*t40*t42*t48;
  t53 = -1.0-t4;
  t59 = t41*t41;
  t62 = (-3.0*t9+(t53*t1+t5)*t2)*Mass*Spin*t59*t36;
  t70 = 2.0*Mass*t2*Spin*r*t59*t41*t3*t48;
  t73 = -t1+2.0*Mass*r-t2;
  t74 = t73*Mass;
  t86 = 1/(t10+(3.0*t9*t4+(3.0*t17*t1+t17*t4*t2)*t2)*t2);
  t92 = t74*Spin*t59*t6*t86;
  t93 = 1/t73;
  t95 = 1/(t1+t5);
  t97 = 1.0-t4;
  t102 = t95*t2;
  t103 = t102*t42;
  t109 = (-2.0*t26-t29)*t2;
  t112 = -t97*t4*t2;
  t128 = 2.0*t8*Spin*Mass*r*t3*t41*t86;
  t131 = t95*r;
  t144 = Spin*Mass*t6*t36;
  t148 = 1/t41;
  t151 = 2.0*t3*Spin*Mass*r*t48*t148;
  t159 = (-t23+t109+(2.0*t9+(-t53*t1+t112)*t2)*Mass)*t36;
  t166 = t3*(-t9-t45*t2-2.0*t97*r*t39)*t48*t148;
  conn[0][0][0] = 0.0;
  conn[0][0][1] = t38;
  conn[0][0][2] = -t51;
  conn[0][0][3] = 0.0;
  conn[0][1][0] = t38;
  conn[0][1][1] = 0.0;
  conn[0][1][2] = 0.0;
  conn[0][1][3] = -t62;
  conn[0][2][0] = -t51;
  conn[0][2][1] = 0.0;
  conn[0][2][2] = 0.0;
  conn[0][2][3] = t70;
  conn[0][3][0] = 0.0;
  conn[0][3][1] = -t62;
  conn[0][3][2] = t70;
  conn[0][3][3] = 0.0;
  conn[1][0][0] = t74*t6*t86;
  conn[1][0][1] = 0.0;
  conn[1][0][2] = 0.0;
  conn[1][0][3] = -t92;
  conn[1][1][0] = 0.0;
  conn[1][1][1] = -t93*t95*(t97*r*t2+t7);
  conn[1][1][2] = -t103;
  conn[1][1][3] = 0.0;
  conn[1][2][0] = 0.0;
  conn[1][2][1] = -t103;
  conn[1][2][2] = t73*t95*r;
  conn[1][2][3] = 0.0;
  conn[1][3][0] = -t92;
  conn[1][3][1] = 0.0;
  conn[1][3][2] = 0.0;
  conn[1][3][3] = -t73*t59*(-t23+t109+(t97*t1+t112)*t2*Mass)*t86;
  conn[2][0][0] = -2.0*t40*t42*t86;
  conn[2][0][1] = 0.0;
  conn[2][0][2] = 0.0;
  conn[2][0][3] = t128;
  conn[2][1][0] = 0.0;
  conn[2][1][1] = -t102*t42*t93;
  conn[2][1][2] = t131;
  conn[2][1][3] = 0.0;
  conn[2][2][0] = 0.0;
  conn[2][2][1] = t131;
  conn[2][2][2] = -t103;
  conn[2][2][3] = 0.0;
  conn[2][3][0] = t128;
  conn[2][3][1] = 0.0;
  conn[2][3][2] = 0.0;
  conn[2][3][3] = t42*(-t10+t22+(-4.0*t97*t25+2.0*(-1.0+t17)*r*t2)*t2*Mass)*t86;
  conn[3][0][0] = 0.0;
  conn[3][0][1] = t144;
  conn[3][0][2] = -t151;
  conn[3][0][3] = 0.0;
  conn[3][1][0] = t144;
  conn[3][1][1] = 0.0;
  conn[3][1][2] = 0.0;
  conn[3][1][3] = t159;
  conn[3][2][0] = -t151;
  conn[3][2][1] = 0.0;
  conn[3][2][2] = 0.0;
  conn[3][2][3] = -t166;
  conn[3][3][0] = 0.0;
  conn[3][3][1] = t159;
  conn[3][3][2] = -t166;
  conn[3][3][3] = 0.0; 

  transform_connection2(coords,conn,connp);

  return;
}

/* connection   */
void bl_conn_func_old( struct of_coord *coords, double ***connp)
{
  double r  = coords->r;
  double th = coords->x[TH];

  /* Assumes M = 1 */
  double t1;
  double t10;
  double t101;
  double t106;
  double t11;
  double t111;
  double t12;
  double t127;
  double t13;
  double t14;
  double t16;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t24;
  double t28;
  double t31;
  double t37;
  double t4;
  double t41;
  double t47;
  double t48;
  double t5;
  double t50;
  double t57;
  double t58;
  double t6;
  double t67;
  double t68;
  double t7;
  double t71;
  double t72;
  double t75;
  double t79;
  double t82;
  double t9;
  double conn[NDIM][NDIM][NDIM];

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  trigvalues( th, &t5, &t21 );


  conn[0][0][0] = 0.0;
  t1 = asq;
  t2 = r*r;
  t4 = M*(t1+t2);
  //  t5 = cos(th);
  t6 = a*t5;
  t7 = -r+t6;
  t9 = r+t6;
  t10 = t5*t5;
  t11 = t1*t10;
  t12 = t2+t11;
  t13 = t12*t12;
  t14 = 1/t13;
  t16 = M*r;
  t18 = -t2+2.0*t16-t1;
  t19 = 1/t18;
  t20 = t9*t14*t19;
  conn[0][0][1] = t4*t7*t20;
  //  t21 = sin(th);
  t24 = M*t1;
  conn[0][0][2] = -2.0*t21*r*t5*t24*t14;
  conn[0][0][3] = 0.0;
  conn[0][1][0] = conn[0][0][1];
  conn[0][1][1] = 0.0;
  conn[0][1][2] = 0.0;
  t28 = t2*t2;
  t37 = t21*t21;
  t31 = (-1.0-t10)*t2;
  conn[0][1][3] = -(-3.0*t28+(t31+t11)*t1)*M*a*t37*t14*t19;
  conn[0][2][0] = conn[0][0][2];
  conn[0][2][1] = 0.0;
  conn[0][2][2] = 0.0;
  t41 = t16*t1;
  conn[0][2][3] = 2.0*t41*a*t5*t37*t21*t14;
  conn[0][3][0] = 0.0;
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  conn[0][3][3] = 0.0;
  t47 = t18*M;
  t48 = 1/t12;
  t50 = t14*t48;
  conn[1][0][0] = t47*t7*t9*t50;
  conn[1][0][1] = 0.0;
  conn[1][0][2] = 0.0;
  conn[1][0][3] = -t47*a*t37*t7*t9*t50;
  conn[1][1][0] = 0.0;
  t57 = t48;
  t58 = 1.0-t10;
  conn[1][1][1] = -t57*(t58*r*t1+(-t2+t11)*M)*t19;
  t67 = t5*t21;
  t68 = t57*t1*t67;
  conn[1][1][2] = -t68;
  conn[1][1][3] = 0.0;
  conn[1][2][0] = 0.0;
  conn[1][2][1] = conn[1][1][2];
  conn[1][2][2] = t18*r*t57;
  conn[1][2][3] = 0.0;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = 0.0;
  conn[1][3][2] = 0.0;
  t71 = t28*r;
  t72 = t2*r;
  t75 = t10*t10;
  t79 = (-2.0*t72*t10-r*t75*t1)*t1;
  t82 = (t75-t10)*t1;
  conn[1][3][3] = -t18*t37*(-t71+t79+(t58*t2+t82)*t1*M)*t50;
  conn[2][0][0] = -2.0*t41*t67*t50;
  conn[2][0][1] = 0.0;
  conn[2][0][2] = 0.0;
  conn[2][0][3] = 2.0*t4*a*r*t5*t21*t50;
  conn[2][1][0] = 0.0;
  conn[2][1][1] = -t68*t19;
  conn[2][1][2] = t57*r;
  conn[2][1][3] = 0.0;
  conn[2][2][0] = 0.0;
  conn[2][2][1] = conn[2][1][2];
  conn[2][2][2] = conn[1][2][1];
  conn[2][2][3] = 0.0;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = 0.0;
  conn[2][3][2] = 0.0;
  t101 = 2.0*t10;
  t106 = t1*t75;
  t111 = -t58;
  conn[2][3][3] = t67*(-t28*t2+((-1.0-t101)*t28+((-t101-t75)*t2-t106)*t1)*t1+
		       (4.0*t111*t72+(-2.0+2.0*t75)*r*t1)*t1*M)*t50;
  conn[3][0][0] = 0.0;
  conn[3][0][1] = M*a*t7*t20;
  t127 = 1/t21;
  conn[3][0][2] = -2.0*t16*a*t5*t14*t127;
  conn[3][0][3] = 0.0;
  conn[3][1][0] = conn[3][0][1];
  conn[3][1][1] = 0.0;
  conn[3][1][2] = 0.0;
  conn[3][1][3] = (-t71+t79+(2.0*t28+(-t31+t82)*t1)*M)*t14*t19;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = 0.0;
  conn[3][2][2] = 0.0;
  conn[3][2][3] = -t5*(-t28+(-2.0*t10*t2-t106)*t1+2.0*t111*r*t24)*t14*t127;
  conn[3][3][0] = 0.0;
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = 0.0;

  transform_connection2(coords,conn,connp);

  return;
}


/******************************************************************************
  bl_dxc_dxs_calc():
  ----------------------
       -- calculates the transformation matrix  Lambda^\hat{a}_a  defined:

      x^\hat{a}[Cartesian]  = \Lambda^\hat{a}_a  x^a[Spherical] 

            where dxc_dxs[i][j] = \Lambda^i_j 

        for  Boyer-Lindquist coordinates.

 ******************************************************************************/
void bl_dxc_dxs_calc(double *x_cart, double *x_spher, double dxc_dxs[][NDIM] )
{
  int i; 
  double r, th, ph, rterm, dr; 
  double sth,cth,sph,cph;
  

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  dxc_dxs[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  dxc_dxs[0][i] = 0. ; } 
#endif
  
  r  = x_spher[RR]; 
  th = x_spher[TH];
  ph = x_spher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  rterm = sqrt( r*r + asq );
  dr = r / rterm; 
  
  dxc_dxs[TT][TT] = 1.                    ;  // dt/dt
  dxc_dxs[XX][RR] = cph * sth * dr        ;  // dx/dr
  dxc_dxs[XX][TH] = rterm * cph * cth     ;  // dx/dtheta
  dxc_dxs[XX][PH] = -(x_cart[YY] - y0_bh) ;  // dx/dphi
  dxc_dxs[YY][RR] = sph * sth * dr        ;  // dy/dr
  dxc_dxs[YY][TH] = rterm * sph * cth     ;  // dy/dtheta
  dxc_dxs[YY][PH] = x_cart[XX] - x0_bh    ;  // dy/dphi
  dxc_dxs[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dxs[ZZ][TH] = -r*sth                ;  // dz/dtheta
  //  dxc_dxs[ZZ][PH] = 0.                    ;  // dz/dphi
  
  return;

}

/******************************************************************************
  bl_dxs_dxc_calc():
  ----------------------
       -- calculates the transformation matrix  Lambda^a_\hat{a}  defined:

      x^a[Spherical]  = \Lambda^a_\hat{a}  x^\hat{a}[Cartesian] 

            where dxs_dxc[i][j] = \Lambda^i_j 

        for  Boyer-Lindquist coordinates.

      -- note that the "dr/d?" and "dth/d?" derivatives are the same 
         as KS, the "dphi/d?" are different, however. 

******************************************************************************/
void bl_dxs_dxc_calc(double *x_spher, double *x_cart, double dxs_dxc[][NDIM] )
{
  int i; 
  double r, xx,yy,zz,rsq,rhosq,term_dr,term_dth,term_dph,dr_dx,dr_dy,dr_dz;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  dxs_dxc[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  dxs_dxc[0][i] = 0. ; } 
#endif
  
  r  = x_spher[RR]; 
  xx = x_cart[XX] - x0_bh; 
  yy = x_cart[YY] - y0_bh; 
  zz = x_cart[ZZ] - z0_bh; 

  rsq = r*r; 
  rhosq = xx*xx + yy*yy + zz*zz ;
  term_dr = 1./( r * (2.*rsq -  rhosq + asq  ) );
  dr_dx  = xx * rsq       * term_dr; 
  dr_dy  = yy * rsq       * term_dr; 
  dr_dz  = zz * (rsq+asq) * term_dr; 

  term_dth =  1./(r*sqrt(rsq-zz*zz));
  term_dph = 1./(xx*xx + yy*yy);

  dxs_dxc[TT][TT] = 1.                                 ;  // dt/dt
  dxs_dxc[RR][XX] = dr_dx                              ;  // dr/dx
  dxs_dxc[RR][YY] = dr_dy                              ;  // dr/dy
  dxs_dxc[RR][ZZ] = dr_dz                              ;  // dr/dz
  dxs_dxc[TH][XX] = term_dth * zz * dr_dx              ;  // dth/dx
  dxs_dxc[TH][YY] = term_dth * zz * dr_dy              ;  // dth/dy
  dxs_dxc[TH][ZZ] = term_dth * ( zz * dr_dz  -  r )    ;  // dth/dz
  dxs_dxc[PH][XX] = -yy*term_dph                       ;  // dph/dx
  dxs_dxc[PH][YY] =  xx*term_dph                       ;  // dph/dy
  //  dxs_dxc[PH][ZZ] = 0.                                 ;  // dph/dz

  return;

}


/***************************************************************************/
/***************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     K E R R - S C H I L D   M E T R I C    R O U T I N E S (spherical coordinates)

 *****************************************************************************/
/******************************************************************************
  ks_gcov_func():
  ----------
       -- covariant form of the metric in Kerr-Schild coordinates
 ******************************************************************************/
void ks_gcov_func(double *x, double gcov[][NDIM])
{
  double r,th;
  double t10,  t12,  t13,  t2 ,  t3 ,  t4 ,  t5 ,  t7 ,  t9 ;

  /* Set up the coordinates : */ 
  r  = x[RR];
  th = x[TH];
  trigvalues(th, &t4, &t12);

  /* Calculate the KS metric using derivation with MAPLE: */
  t2 = r*r;
  t3 = asq;
  //  t4 = cos(th);
  t5 = t4*t4;
  t7 = t2+t3*t5;
  t9 = M*r/t7;
  t10 = 2.0*t9;
  gcov[0][0] = t10-1.0;
  gcov[0][1] = t10;
  gcov[0][2] = 0.0;
  //  t12 = sin(th);
  t13 = t12*t12;
  gcov[0][3] = -2.0*t9*a*t13;
  gcov[1][0] = t10;
  gcov[1][1] = 1.0+t10;
  gcov[1][2] = 0.0;
  gcov[1][3] = -a*gcov[1][1]*t13;
  gcov[2][0] = 0.0;
  gcov[2][1] = 0.0;
  gcov[2][2] = t7;
  gcov[2][3] = 0.0;
  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][2] = 0.0;
  gcov[3][3] = t13*(t2+t3*(1.0+2.0*t9*t13));
  
  return;
}


/******************************************************************************
  ks_gcon_func():
  ----------
       -- contravariant form of the metric in Kerr-Schild coordinates
 ******************************************************************************/
void ks_gcon_func(double *x, double gcon[][NDIM])
{
  double r,th;
  double t1 ,   t10,   t13,   t14,   t2 ,   t3 ,   t4 ,   t5 ,   t8 ;

  /* Set up the coordinates : */ 
  r  = x[RR];
  th = x[TH];
  trigvalues(th, &t4, &t13);

  /* Calculate the KS metric using derivation with MAPLE: */
  t1 = M*r;
  t2 = r*r;
  t3 = asq;
  //  t4 = cos(th);
  t5 = t4*t4;
  t8 = 1/(t2+t3*t5);
  t10 = 2.0*t1*t8;
  gcon[0][0] = -1.0-t10;
  gcon[0][1] = t10;
  gcon[0][2] = 0.0;
  gcon[0][3] = 0.0;
  gcon[1][0] = t10;
  gcon[1][1] = (t2-2.0*t1+t3)*t8;
  gcon[1][2] = 0.0;
  gcon[1][3] = a*t8;
  gcon[2][0] = 0.0;
  gcon[2][1] = 0.0;
  gcon[2][2] = t8;
  gcon[2][3] = 0.0;
  gcon[3][0] = 0.0;
  gcon[3][1] = gcon[1][3];
  gcon[3][2] = 0.0;
  //  t13 = sin(th);
  t14 = t13*t13;
  gcon[3][3] = t8/t14;

  return;
}


/*********************************************************************************************
   Scott's MKS connection that can be used for any transformation between r,th <-> X1,X2 :
     -- assumes that only X1,X2 are different from (r,th,ph)
*********************************************************************************************/
void ks_conn_func_1(struct of_coord *coords, double ***conn )
{
  int i, j, k, l;
  double r,th,sigma,dx_dxp_dxp[NDIM][NDIM][NDIM];

  double t1,   t10,   t102,   t1024,   t1035,   t1037,   t104,   t11,   t114,   t116,   t119;
  double  t12,   t121,   t123,   t126,   t129,   t130,   t132,   t14,   t148,   t149,   t15,   t152;
  double t154,   t156,   t157,   t158,   t159,   t161,   t169,   t17,   t171,   t172,   t175,   t177;
  double t185,   t2,   t203,   t204,   t208,   t209,   t21,   t210,   t212,   t214,   t22,   t221;
  double   t222,   t224,   t227,   t23,   t236,   t24,   t240,   t241,   t242,   t243,   t245,   t246;
  double    t247,   t248,   t25,   t250,   t251,   t258,   t26,   t260,   t261,   t263,   t264,   t271;
  double   t273,   t275,   t276,   t278,   t28,   t280,   t281,   t283,   t284,   t285,   t286,   t288;
  double    t289,   t29,   t297,   t299,   t3,   t30,   t300,   t303,   t305,   t306,   t308,   t309;
  double    t31,   t310,   t313,   t314,   t320,   t325,   t327,   t328,   t329,   t330,   t333,   t336;
  double    t338,   t34,   t340,   t342,   t344,   t346,   t35,   t358,   t361,   t363,   t366,   t367;
  double    t368,   t370,   t372,   t375,   t38,   t380,   t381,   t384,   t385,   t387,   t39,   t399;
  double    t4,   t40,   t400,   t402,   t404,   t405,   t406,   t408,   t409,   t41,   t411,   t412;
  double    t418,   t42,   t421,   t425,   t428,   t431,   t432,   t433,   t434,   t437,   t440,   t442;
  double    t448,   t451,   t453,   t454,   t459,   t462,   t467,   t469,   t480,   t481,   t486,   t487;
  double    t488,   t491,   t492,   t498,   t501,   t504,   t507,   t508,   t510,   t512,   t52,   t521;
  double    t528,   t53,   t530,   t553,   t556,   t56,   t57,   t588,   t60,   t607,   t627,   t628;
  double    t63,   t630,   t631,   t632,   t634,   t636,   t637,   t64,   t651,   t652,   t654,   t656;
  double    t657,   t659,   t661,   t662,   t670,   t673,   t675,   t677,   t686,   t689,   t7,   t712;
  double    t74,   t748,   t75,   t78,   t793,   t794,   t795,   t799,   t8,   t800,   t801,   t803;
  double    t806,  t807,   t813,   t816,   t822,   t83,   t831,   t84,   t845,   t86,   t89,   t891; 
  double    t90,   t91,   t916,   t917,   t920,   t924,   t928,   t940,   t946,   t968,   t97, t970,t991;


  // set coordinate transformation matrices: 
  r  = coords->x[RR];
  th = coords->x[TH];
  dx_dxp_dxp_calc( coords->x, coords->xp, dx_dxp_dxp) ; 


  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  trigvalues( th, &t3, &t29 );

  // the connection coefficients: 

  //  t1 = rf(X1,X2);
  t1 = r;
  //  t2 = thf(X1,X2);
  t2 = th;
  //  t3 = cos(t2);
  t4 = a*t3;
  t7 = (-t1+t4)*(t1+t4);
  t8 = M*M;
  t10 = t1*t1;
  t11 = asq;
  t12 = t3*t3;
  t14 = t10+t11*t12;
  t15 = t14*t14;
  t17 = 1/t15/t14;
  conn[0][0][0] = -2.0*t7*t8*t1*t17;
  t21 = t10*t10;
  //  t22 = diff(rf(X1,X2),X1);
  t22 = coords->dx_dxp[RR][1];
  t23 = t21*t22;
  t24 = t10*t1;
  t25 = M*t24;
  t26 = t25*t22;
  t28 = t24*t3;
  //  t29 = sin(t2);
  //  t30 = diff(thf(X1,X2),X1);
  t30 = coords->dx_dxp[TH][1];
  t31 = t29*t30;
  t34 = M*t1;
  t35 = t22*t12;
  t38 = t12*t12;
  t39 = t38*t22;
  t40 = t29*t1;
  t41 = t12*t3;
  t42 = t41*t30;
  conn[0][0][1] = -M*(-t23-2.0*t26+(2.0*t28*t31+2.0*t34*t35+(t39+2.0*t40*t42)*t11)*t11)*t17;
  //  t52 = diff(rf(X1,X2),X2);
  t52 = coords->dx_dxp[RR][2];
  t53 = t21*t52;
  //  t56 = diff(thf(X1,X2),X2);
  t56 = coords->dx_dxp[TH][2];
  t57 = t29*t56;
  t60 = t52*t12;
  t63 = t38*t52;
  t64 = t41*t1;
  conn[0][0][2] = -M*(-t53-2.0*t25*t52+(2.0*t28*t57+2.0*t34*t60+(t63+2.0*t64*t57)*t11)*t11)*t17;
  t74 = -1.0+t3;
  t75 = 1.0+t3;
  t78 = t7*t74*t75;
  conn[0][0][3] = -2.0*t78*a*t1*t8*t17;
  conn[0][1][0] = conn[0][0][1];
  t83 = t30*t30;
  t84 = t21*t10;
  t86 = t22*t22;
  t89 = t24*t22;
  t90 = t3*t29;
  t91 = t90*t30;
  t97 = t86*t12;
  t102 = t22*t29;
  t104 = t64*t102;
  conn[0][1][1] = -2.0*(t83*t84-t21*t86-t25*t86+(2.0*t89*t91+2.0*t83*t21*t12+
						 t34*t97+(t83*t10*t38+t38*t86+2.0*t104*t30)*t11)*t11)*M*t17;
  t114 = t22*t52;
  t116 = t30*t56;
  t119 = t24*t52;
  t121 = t114*t12;
  t123 = t21*t12;
  t126 = t90*t56;
  t129 = t52*t29;
  t130 = t129*t30;
  t132 = t10*t38;
  conn[0][1][2] = -2.0*(-t25*t114+t116*t84-t23*t52+(t119*t91+t34*t121+2.0*
						    t116*t123+t89*t126
						    +(t39*t52+t64*t130+t116*t132+t104*t56)*t11)*t11)*M*t17;
  t148 = 2.0*t28*t30;
  t149 = t102*t12;
  t152 = t41*t24;
  t154 = 2.0*t152*t30;
  t156 = 2.0*t64*t30;
  t157 = t39*t29;
  t158 = t30*t1;
  t159 = t38*t3;
  t161 = 2.0*t158*t159;
  t169 = a*t29*t17;
  conn[0][1][3] = -(2.0*t25*t102+t23*t29+(-t148-2.0*t34*t149+t154+(-t156-t157+t161)*t11)*t11)*M*t169;
  conn[0][2][0] = conn[0][0][2];
  conn[0][2][1] = conn[0][1][2];
  t171 = t52*t52;
  t172 = t171*M;
  t175 = t56*t56;
  t177 = t1*t12;
  t185 = t52*t41;
  conn[0][2][2] = -2.0*(-t172*t24-t171*t21+t175*t84+(t172*t177+2.0*t119*t126+
						     2.0*t175*t21*t12
						     +(t171*t38+2.0*t185*t40*t56+t175*t10*t38)*t11)*t11)*M*t17;
  t203 = 2.0*t152*t56;
  t204 = t129*t12;
  t208 = 2.0*t28*t56;
  t209 = t63*t29;
  t210 = t56*t1;
  t212 = 2.0*t210*t159;
  t214 = 2.0*t64*t56;
  conn[0][2][3] = (-2.0*t25*t129-t53*t29+(-t203+2.0*t34*t204+t208+(t209-t212+t214)*t11)*t11)*M*t169;
  conn[0][3][0] = conn[0][0][3];
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  t221 = t21*t1;
  t222 = t24*t12;
  t224 = t10*t12;
  t227 = t1*t38;
  t236 = (-t221+(-2.0*t222+(-t224+t10)*M+(-t227+(t38-t12)*M)*t11)*t11)*t74*t75;
  conn[0][3][3] = -2.0*t236*t34*t17;
  t240 = t56*M;
  t241 = t240*t24;
  t242 = 2.0*t241;
  t243 = t56*t21;
  t245 = 2.0*t240*t177;
  t246 = t56*t10;
  t247 = t246*t12;
  t248 = t1*t3;
  t250 = 2.0*t248*t129;
  t251 = t56*t12;
  t258 = t30*t52;
  t260 = 1/(-t22*t56+t258);
  t261 = t260*t17;
  conn[1][0][0] = -(-t242+t243+(t245-t247+t250+t246-t251*t11)*t11)*M*t261;
  t263 = M*t22;
  t264 = t56*t38;
  t271 = (-t242+(t245-t247+t250+t246+(-t251+t264)*t11)*t11)*t260*t17;
  conn[1][0][1] = -t263*t271;
  t273 = M*t52;
  conn[1][0][2] = -t273*t271;
  t275 = t24*t29;
  t276 = t240*t275;
  t278 = t119*t3;
  t280 = t243*t29;
  t281 = t40*t12;
  t283 = 2.0*t240*t281;
  t284 = t29*t12;
  t285 = t246*t284;
  t286 = t52*t3;
  t288 = 2.0*t286*t1;
  t289 = t57*t10;
  t297 = a*t29*t260*t17;
  conn[1][0][3] = (-2.0*t276+2.0*t278+t280+(t283-t285+t288+t289-t56*t11*t284)*t11)*M*t297;
  conn[1][1][0] = conn[1][0][1];
  //  t299 = diff(diff(thf(X1,X2),X1),X1);
  t299 = dx_dxp_dxp[TH][1][1];
  t300 = t52*t299;
  t303 = t258*t221*t22;
  t305 = t56*t84;
  //  t306 = diff(diff(rf(X1,X2),X1),X1);
  t306 = dx_dxp_dxp[RR][1][1] ;
  t308 = t83*t56;
  t309 = t21*t24;
  t310 = t308*t309;
  t313 = 2.0*t308*t84;
  t314 = t56*t24;
  t320 = t30*t22;
  t325 = t258*t89*t12;
  t327 = t308*t221;
  t328 = t52*t83;
  t329 = t90*t21;
  t330 = t328*t329;
  t333 = t306*t12;
  t336 = t221*t12;
  t338 = 2.0*t308*t336;
  t340 = 4.0*t308*t123;
  t342 = t248*t29;
  t344 = 2.0*t52*t86*t342;
  t346 = t56*t86;
  t358 = t306*t38;
  t361 = t1*t22;
  t363 = t258*t361*t38;
  t366 = 2.0*t308*t222;
  t367 = t24*t38;
  t368 = t308*t367;
  t370 = t41*t29*t10;
  t372 = 2.0*t328*t370;
  t375 = 2.0*t308*t132;
  t380 = t159*t29;
  t381 = t380*t56;
  t384 = t308*t227;
  t385 = t38*t12;
  t387 = t328*t380;
  conn[1][1][1] = -(-t300*t84-2.0*t303+t305*t306-t310+(-t243*t86+t313-2.0*t314*t86*M)*M
		    +(-2.0*t320*t3*t280-4.0*t325-t327+t330-3.0*t300*t123+3.0*t243*t333
		      -t338+(t340+t344-t246*t97+t346*t10+2.0*t210*t97*M)*M
		      +(-4.0*t320*t41*t289-3.0*t300*t132+3.0*t246*t358-2.0*t363-t366-t368+t372
			+(-t346*t12+t375+2.0*t346*t38)*M
			+(-2.0*t320*t381-t384-t300*t385+t387+t56*t306*t385)*t11)*t11)*t11)*t260*t17;
  //  t399 = diff(diff(thf(X1,X2),X1),X2);
  t399 = dx_dxp_dxp[TH][1][2]; 
  t400 = t52*t399;
  //  t402 = diff(diff(rf(X1,X2),X1),X2);
  t402 = dx_dxp_dxp[RR][1][2]; 
  t404 = t30*t175;
  t405 = t404*t309;
  t406 = t171*t30;
  t408 = t56*t221;
  t409 = t114*t408;
  t411 = 2.0*t404*t84;
  t412 = t52*t56;
  t418 = t402*t12;
  t421 = t404*t221;
  t425 = t114*t314*t12;
  t428 = 2.0*t404*t336;
  t431 = t22*t175;
  t432 = t431*t329;
  t433 = t10*t22;
  t434 = t433*t12;
  t437 = 4.0*t404*t123;
  t440 = 2.0*t171*t22*t342;
  t442 = t210*M;
  t448 = 2.0*t370*t431;
  t451 = t114*t210*t38;
  t453 = 2.0*t404*t222;
  t454 = t402*t38;
  t459 = t404*t367;
  t462 = 2.0*t404*t132;
  t467 = t404*t227;
  t469 = t431*t380;
  conn[1][1][2] = (t400*t84-t305*t402+t405+t406*t221+t409+(-t411+t412*t23+2.0*t114*t241)*M
		   +(-3.0*t243*t418+t421+3.0*t400*t123+2.0*t425+t428+2.0*t406*t222+
		     t432+(t412*t434-t437-t440-t412*t433-2.0*t121*t442)*M
		     +(t448+t406*t227+t451+t453-3.0*t246*t454+3.0*t400*t132+t459
		       +(t114*t251-t462-2.0*t114*t264)*M+(t467+t400*t385+t469
							  -t56*t402*t385)*t11)*t11)*t11)*t260*t17;
  t480 = t286*t21;
  t481 = t57*t221;
  t486 = 2.0*t185*t10;
  t487 = t57*t222;
  t488 = 2.0*t487;
  t491 = t57*t227;
  t492 = t52*t159;
  t498 = (-t491+t492+(-t57*t12+t57*t38)*M)*t11;
  t501 = t480-t481+(2.0*t278-2.0*t276)*M+(t486-t488+(-t285+t289+t288+t283)*M+t498)*t11;
  conn[1][1][3] = t501*t22*t297;
  conn[1][2][0] = conn[1][0][2];
  conn[1][2][1] = conn[1][1][2];
  t504 = t171*t56;
  //  t507 = diff(diff(thf(X1,X2),X2),X2);
  t507 = dx_dxp_dxp[TH][2][2];
  t508 = t52*t507;
  //  t510 = diff(diff(rf(X1,X2),X2),X2);
  t510 = dx_dxp_dxp[RR][2][2];
  t512 = t175*t56;
  t521 = t512*t221;
  t528 = t52*t175;
  t530 = t510*t12;
  t553 = t510*t38;
  t556 = t512*t24;
  conn[1][2][2] = -(-2.0*t504*t221-t508*t84+t305*t510-t512*t309
		    +(2.0*t512*t84-t504*t21-2.0*t504*t25)*M
		    +(-2.0*t521*t12-3.0*t508*t123-4.0*t504*t222-t521-t528*t329+3.0*t243*t530
		      +(-t504*t224+t504*t10+2.0*t342*t171*t52+4.0*t512*t21*t12+2.0*t12*t171*t442)*M
		      +(-3.0*t508*t132-2.0*t504*t227-2.0*t528*t370+3.0*t246*t553-2.0*t556*t12-t556*t38
			+(2.0*t512*t10*t38-t504*t12+2.0*t504*t38)*M
			+(-t528*t380-t512*t1*t38-t508*t385+t56*t510*t385)*t11)*t11)*t11)*t260*t17;
  conn[1][2][3] = t501*t52*t297;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = conn[1][1][3];
  conn[1][3][2] = conn[1][2][3];
  t588 = t84*t29;
  t607 = t29*t38;
  conn[1][3][3] = -(-t56*t309*t29+t286*t84+2.0*t240*t588
		    +(-t481-2.0*t408*t284+t480+2.0*t185*t21+(-4.0*t185*t24+t280+3.0*t243*t284+4.0*t278
							     +(-2.0*t314*t29+2.0*t487)*M)*M
		      +(t492*t10+t486-t314*t607-t488+(-2.0*t492*t1+t289+t288-2.0*t285+3.0*t246*t607
						      +(-2.0*t491+2.0*t210*t284)*M)*M+t498)*t11)*t11)*t29*t261;
  t627 = t30*t21;
  t628 = t30*M;
  t630 = 2.0*t628*t24;
  t631 = t30*t10;
  t632 = t631*t12;
  t634 = 2.0*t628*t177;
  t636 = 2.0*t361*t90;
  t637 = t30*t12;
  conn[2][0][0] = -(-t627+t630+(t632-t634-t631-t636+t637*t11)*t11)*M*t261;
  t651 = (-t630+(t634+t636+t631-t632+(-t637+t30*t38)*t11)*t11)*t260*t17;
  conn[2][0][1] = t263*t651;
  conn[2][0][2] = t273*t651;
  t652 = t628*t275;
  t654 = t89*t3;
  t656 = t627*t29;
  t657 = t631*t284;
  t659 = 2.0*t628*t281;
  t661 = 2.0*t361*t3;
  t662 = t31*t10;
  conn[2][0][3] = (2.0*t652-2.0*t654-t656+(t657-t659-t661-t662+t30*t11*t284)*t11)*M*t297;
  conn[2][1][0] = conn[2][0][1];
  t670 = t30*t86;
  t673 = t30*t84;
  t675 = t22*t299;
  t677 = t83*t30;
  t686 = t677*t221;
  t689 = t83*t22;
  t712 = t677*t24;
  conn[2][1][1] = -(2.0*t670*t221-t673*t306+t675*t84+t677*t309+(-2.0*t677*t84+t627*t86+2.0*t670*t25)*M
		    +(2.0*t686*t12+t689*t329-3.0*t627*t333+t686+4.0*t670*t222+3.0*t675*t123
		      +(-2.0*t86*t22*t1*t90+t631*t97-t670*t10-4.0*t677*t21*t12-2.0*t637*t86*t34)*M
		      +(2.0*t712*t12+t712*t38+2.0*t670*t227+3.0*t675*t132-3.0*t631*t358+2.0*t689*t370
			+(t670*t12-2.0*t677*t10*t38-2.0*t670*t38)*M
			+(t675*t385-t30*t306*t385+t689*t380+t677*t1*t38)*t11)*t11)*t11)*t260*t17;
  t748 = t22*t399;
  conn[2][1][2] = -(t303+t310-t673*t402+t748*t84+t346*t221
		    +(-t313+t258*t23+2.0*t258*t26)*M
		    +(t327+t330+2.0*t325-3.0*t627*t418+2.0*t346*t222+3.0*t748*t123+t338
		      +(-t340-t258*t433-t344+t258*t434-2.0*t60*t30*t361*M)*M
		      +(t346*t227+t372+t363+t366+t368-3.0*t631*t454+3.0*t748*t132
			+(t258*t35-t375-2.0*t258*t39)*M+(t387+t748*t385+t384
							 -t30*t402*t385)*t11)*t11)*t11)*t260*t17;
  t793 = t31*t221;
  t794 = t3*t22;
  t795 = t794*t21;
  t799 = t31*t222;
  t800 = 2.0*t799;
  t801 = t22*t41;
  t803 = 2.0*t801*t10;
  t806 = t31*t227;
  t807 = t22*t159;
  t813 = (-t806+t807+(t31*t38-t31*t12)*M)*t11;
  t816 = -t793+t795+(2.0*t654-2.0*t652)*M+(-t800+t803+(t661-t657+t662+t659)*M+t813)*t11;
  conn[2][1][3] = -t816*t22*t297;
  conn[2][2][0] = conn[2][0][2];
  conn[2][2][1] = conn[2][1][2];
  t822 = t22*t507;
  t831 = t3*t56;
  t845 = t41*t56;
  conn[2][2][2] = -(-t673*t510+2.0*t409+t405+t822*t84+(t406*t21-t411+2.0*t406*t25)*M
		    +(-3.0*t627*t530+t421-t432+2.0*t130*t831*t21+t428+3.0*t822*t123+4.0*t425
		      +(t406*t224-t406*t10-t440-t437-2.0*t406*t177*M)*M
		      +(4.0*t130*t845*t10+2.0*t451+3.0*t822*t132+t453+t459-t448-3.0*t631*t553
			+(t406*t12-t462-2.0*t406*t38)*M+(2.0*t258*t381+t822*t385-t30*t510*t385
							 +t467-t469)*t11)*t11)*t11)*t260*t17;
  conn[2][2][3] = -t816*t52*t297;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = conn[2][1][3];
  conn[2][3][2] = conn[2][2][3];
  t891 = t30*t24;
  conn[2][3][3] = (t794*t84+2.0*t628*t588-t30*t309*t29
		   +(t795-t793-2.0*t30*t221*t284+2.0*t801*t21
		     +(3.0*t627*t284+t656-4.0*t801*t24+4.0*t654+(-2.0*t891*t29+2.0*t799)*M)*M
		     +(t807*t10+t803-t891*t607-t800+(3.0*t631*t607-2.0*t807*t1+t661-2.0*t657+t662
						     +(2.0*t158*t284-2.0*t806)*M)*M+t813)*t11)*t11)*t29*t261;
  t917 = a*M;
  conn[3][0][0] = -t7*t917*t17;
  t920 = t102*t10;
  t924 = 1/t29;
  t916 = t924*t17;
  conn[3][0][1] = -t917*(-t920+t148+(t156+t149)*t11)*t916;
  t928 = t129*t10;
  conn[3][0][2] = -t917*(t208-t928+(t214+t204)*t11)*t916;
  conn[3][0][3] = -t78*M*t11*t17;
  conn[3][1][0] = conn[3][0][1];
  t940 = t83*t29;
  t946 = t86*t29;
  t968 = t916;
  conn[3][1][1] = -(t940*t221+2.0*t794*t627+(4.0*t794*t891-t946*t10)*M
		    +(4.0*t801*t631+2.0*t940*t222+(4.0*t361*t42+t946*t12)*M+(t940*t227+2.0*t807*t30)*t11)
		    *t11)*a*t968;
  t970 = t29*t221;
  t991 = t52*t1;
  conn[3][1][2] = -(t116*t970+t286*t627+t794*t243
		    +(2.0*t286*t891-t102*t52*t10+2.0*t794*t314)*M
		    +(2.0*t185*t631+2.0*t801*t246+2.0*t116*t275*t12+(2.0*t361*t845+2.0*t991*t42+t102*t60)*M
		      +(t116*t40*t38+t492*t30+t807*t56)*t11)*t11)*a*t968;
  t1024 = t38*t41;
  conn[3][1][3] = (t3*t30*t84+t970*t22+(3.0*t42*t21+2.0*t275*t35
					+(t102*t224+t148-t154-t920)*M
					+(t40*t39+3.0*t159*t30*t10+(t156-t161+t149-t157)*M
					  +t1024*t30*t11)*t11)*t11)*t924*t17;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = conn[3][1][2];
  t1035 = t175*t29;
  t1037 = t171*t29;
  conn[3][2][2] = -(2.0*t286*t243+t1035*t221+(-t1037*t10+4.0*t286*t314)*M
		    +(4.0*t185*t246+2.0*t1035*t222+(t1037*t12+4.0*t991*t845)*M
		      +(t1035*t227+2.0*t492*t56)*t11)*t11)*a*t968;
  conn[3][2][3] = -(-t970*t52-t831*t84+(-2.0*t275*t60-3.0*t845*t21
					+(t203-t208-t129*t224+t928)*M
					+(-t40*t63-3.0*t159*t56*t10+(t209+t212-t204-t214)*M
					  -t1024*t56*t11)*t11)*t11)*t924*t17;
  conn[3][3][0] = conn[3][0][3];
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = -t236*a*t17;

  return;

}

/*********************************************************************************************
   ks_conn_func_2():
  -----------------
      -- this version of the KS-spherical connection does not assume anything about coordinate 
           transformation to numerical coordinates, but uses the transform_connection2() 
           routine; 
*********************************************************************************************/
void ks_conn_func_2(struct of_coord *coords, double ***connp )
{
  double r  = coords->r;
  double th = coords->x[TH];
  double t1;
  double t10;
  double t11;
  double t110;
  double t112;
  double t114;
  double t115;
  double t12;
  double t127;
  double t129;
  double t13;
  double t14;
  double t15;
  double t157;
  double t16;
  double t164;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t26;
  double t29;
  double t3;
  double t31;
  double t34;
  double t37;
  double t4;
  double t44;
  double t49;
  double t5;
  double t59;
  double t6;
  double t60;
  double t61;
  double t64;
  double t68;
  double t69;
  double t7;
  double t74;
  double t79;
  double t8;
  double t80;
  double t81;
  double t82;
  double t87;
  double t9;
  double t90;
  double t97;
  double t98;
  double conn[NDIM][NDIM][NDIM];

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  trigvalues( th, &t3, &t29 );

  t1 = M*M;
  t2 = t1*r;
  //  t3 = cos(th);
  t4 = a*t3;
  t5 = -r+t4;
  t6 = r+t4;
  t8 = r*r;
  t9 = asq;
  t10 = t3*t3;
  t11 = t9*t10;
  t12 = t8+t11;
  t13 = t12*t12;
  t7 = 1/t13;
  t14 = 1/t12;
  t15 = t7*t14;
  t16 = t5*t6*t15;
  conn[0][0][0] = -2.0*t2*t16;
  t19 = M*t5;
  t20 = M*r;
  t21 = 2.0*t20;
  t22 = t8+t11+t21;
  t23 = t6*t22;
  conn[0][0][1] = -t19*t23*t15;
  t26 = t7;
  //  t29 = sin(th);
  t31 = t9*t3*t29;
  conn[0][0][2] = -2.0*t20*t26*t31;
  t34 = t29*t29;
  t37 = t34*t5*t6*t15;
  conn[0][0][3] = 2.0*t2*a*t37;
  conn[0][1][0] = conn[0][0][1];
  conn[0][1][1] = -2.0*t19*t6*(t8+t11+t20)*t15;
  conn[0][1][2] = conn[0][0][2];
  t44 = M*a;
  conn[0][1][3] = t44*t5*t23*t34*t15;
  conn[0][2][0] = conn[0][1][2];
  conn[0][2][1] = conn[0][2][0];
  t49 = t14;
  conn[0][2][2] = -2.0*M*t8*t49;
  conn[0][2][3] = 2.0*t34*t29*t3*M*t9*a*r*t26;
  conn[0][3][0] = conn[0][0][3];
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  t59 = t8*t8;
  t60 = t59*r;
  t61 = t8*r;
  t64 = t10*t10;
  t68 = (-2.0*t61*t10-r*t64*t9)*t9;
  t69 = 1.0-t10;
  t74 = (t69*t8+(t64-t10)*t9)*t9;
  t79 = (-t60+t68+t74*M)*t15;
  conn[0][3][3] = 2.0*t20*t34*t79;
  t80 = -t8+t21-t9;
  t81 = t80*M;
  conn[1][0][0] = t81*t16;
  t82 = -t69;
  t87 = t6*M*t15;
  conn[1][0][1] = (t82*t9+t21)*t5*t87;
  conn[1][0][2] = 0.0;
  conn[1][0][3] = -t81*a*t37;
  conn[1][1][0] = conn[1][0][1];
  t90 = 2.0*t10;
  conn[1][1][1] = (t8+(-1.0+t90)*t9+t21)*t5*t87;
  conn[1][1][2] = -t31*t49;
  t97 = a*t34;
  t98 = r*t9;
  conn[1][1][3] = -t97*(-t60+t68+(t74+(-2.0*t61+2.0*t98*t10)*M)*M)*t15;
  conn[1][2][0] = 0.0;
  conn[1][2][1] = conn[1][1][2];
  conn[1][2][2] = t80*r*t49;
  conn[1][2][3] = 0.0;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = conn[1][1][3];
  conn[1][3][2] = 0.0;
  t110 = t79;
  conn[1][3][3] = -t80*t34*t110;
  t112 = t9*M;
  t114 = t3*t29;
  t115 = t114*t15;
  conn[2][0][0] = -2.0*t112*r*t115;
  conn[2][0][1] = conn[2][0][0];
  conn[2][0][2] = 0.0;
  conn[2][0][3] = 2.0*(t9+t8)*M*a*r*t3*t29*t15;
  conn[2][1][0] = conn[2][0][1];
  conn[2][1][1] = conn[2][1][0];
  conn[2][1][2] = t49*r;
  t127 = t64*t9;
  t129 = (2.0*t8*t10+t127)*t9;
  conn[2][1][3] = a*(t59+t129+(2.0*t61+2.0*t98)*M)*t115;
  conn[2][2][0] = 0.0;
  conn[2][2][1] = conn[2][1][2];
  conn[2][2][2] = conn[1][2][1];
  conn[2][2][3] = 0.0;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = conn[2][1][3];
  conn[2][3][2] = 0.0;
  conn[2][3][3] = t114*(-t59*t8+((-1.0-t90)*t59+((-t90-t64)*t8-t127)*t9)*t9+(4.0*t82*t61+(-2.0+2.0*t64)*r*t9)*t9*M)*t15;
  conn[3][0][0] = -t44*t16;
  conn[3][0][1] = conn[3][0][0];
  t157 = 1/t29;
  conn[3][0][2] = -2.0*t44*r*t3*t157*t26;
  conn[3][0][3] = t112*t34*t16;
  conn[3][1][0] = conn[3][0][1];
  conn[3][1][1] = conn[3][1][0];
  t164 = t157*t26;
  conn[3][1][2] = -t4*t22*t164;
  conn[3][1][3] = -t110;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = conn[3][1][2];
  conn[3][2][2] = -a*t49*r;
  conn[3][2][3] = -t3*(-t59-t129+2.0*t82*r*t112)*t164;
  conn[3][3][0] = conn[3][0][3];
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = t97*t110;

  transform_connection2(coords,conn,connp);
  return;

}


/******************************************************************************
  ks_gdet_func():
  ----------
       -- returns the sqrt( - det(g_uv) )  for the KS metric;
 ******************************************************************************/
void ks_gdet_func(double *x, double *gdet )
{
  bl_gdet_func( x, gdet );
  
  return;
}


/******************************************************************************
  bl_to_ks_cov():
  ----------
       -- transforms a covariant vector in BL coordinates to KS coordinates;
 ******************************************************************************/
void bl_to_ks_cov(double *x, double blcov[], double kscov[] )
{
  double delta_m1, dtBL_drKS, dphiBL_drKS;
  
  delta_m1 = 1./( x[RR] * ( x[RR] - 2*M ) + asq );
  dtBL_drKS   = -2*M*x[RR] * delta_m1; 
  dphiBL_drKS = -a * delta_m1;

  kscov[TT] = blcov[TT];
  kscov[RR] = blcov[RR]  +  blcov[TT] * dtBL_drKS  +  blcov[PH] * dphiBL_drKS;
  kscov[TH] = blcov[TH];
  kscov[PH] = blcov[PH];
  
  return;
}


/******************************************************************************
  bl_to_ks_con():
  ----------
       -- transforms a contravariant vector in BL coordinates to KS coordinates;
 ******************************************************************************/
void bl_to_ks_con(double *x, double blcon[], double kscon[] )
{
  double delta_m1, dtKS_drBL, dphiKS_drBL;
  
  delta_m1 = 1./( x[RR] * ( x[RR] - 2*M ) + asq );
  dtKS_drBL   = 2*M*x[RR] * delta_m1; 
  dphiKS_drBL = a * delta_m1;

  kscon[TT] = blcon[TT]  +  blcon[RR] * dtKS_drBL;
  kscon[RR] = blcon[RR];
  kscon[TH] = blcon[TH];
  kscon[PH] = blcon[PH]  +  blcon[RR] * dphiKS_drBL;
  
  return;
}


/******************************************************************************
  ks_to_bl_cov():
  ----------
       -- transforms a covariant vector in KS coordinates to BL coordinates;
 ******************************************************************************/
void ks_to_bl_cov(double *x, double kscov[], double blcov[] )
{
  double delta_m1, dtKS_drBL, dphiKS_drBL;
  
  delta_m1 = 1./( x[RR] * ( x[RR] - 2*M ) + asq );
  dtKS_drBL   = 2*M*x[RR] * delta_m1; 
  dphiKS_drBL = a * delta_m1;

  blcov[TT] = kscov[TT];
  blcov[RR] = kscov[RR]  +  kscov[TT] * dtKS_drBL  +  kscov[PH] * dphiKS_drBL;
  blcov[TH] = kscov[TH];
  blcov[PH] = kscov[PH];
  
  return;
}


/******************************************************************************
  ks_to_bl_con():
  ----------
       -- transforms a contravariant vector in KS coordinates to BL coordinates;
 ******************************************************************************/
void ks_to_bl_con(double *x, double kscon[], double blcon[] )
{
  double delta_m1, dtBL_drKS, dphiBL_drKS;
  
  delta_m1 = 1./( x[RR] * ( x[RR] - 2*M ) + asq );
  dtBL_drKS   = -2*M*x[RR] * delta_m1; 
  dphiBL_drKS = -a * delta_m1;

  blcon[TT] = kscon[TT]  +  kscon[RR] * dtBL_drKS;
  blcon[RR] = kscon[RR];
  blcon[TH] = kscon[TH];
  blcon[PH] = kscon[PH]  +  kscon[RR] * dphiBL_drKS;
  
  return;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     M I N K O W S K I  (spherical)  M E T R I C    R O U T I N E S

 *****************************************************************************/

/******************************************************************************
  mink_gcov_func():
  ----------
       -- Minkowski form of the covariant metric:
 ******************************************************************************/
void mink_spherical_gcov_func(double *x, double gcov[][NDIM])
{
  int i;
  double r,th;
  double t1, cthtmp;
  
  /* Set up the coordinates : */ 
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cthtmp, &t1);
  //t1 = sin(th);

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  gcov[i][j] = 0. ; } 
#else
  for(i = 0; i < NDIM*NDIM; i++) gcov[0][i] = 0. ; 
#endif

  gcov[0][0] = -1.;
  gcov[1][1] =  1.;
  gcov[2][2] = r*r;
  gcov[3][3] = gcov[2][2]*t1*t1;

  return;
}



/******************************************************************************
  mink_gcon_func():
  ---------------
       -- general routine to calculate the contravariant form of the metric
 ******************************************************************************/
void mink_spherical_gcon_func(double *x, double gcon[][NDIM])
{
  int i;
  double r,th;
  double t1, t2, t3, cthtmp;
  
  /* Set up the coordinates : */ 
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cthtmp, &t2);

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  gcon[i][j] = 0. ; } 
#else
  for(i = 0; i < NDIM*NDIM; i++) gcon[0][i] = 0. ; 
#endif

  //t2 = sin(th);
  t1 = r*r;
  t3 = t2*t2;
  gcon[0][0] = -1.;
  gcon[1][1] =  1.;
  gcon[2][2] = 1./t1;
  gcon[3][3] = gcon[2][2]/t3;

  return;

}

/******************************************************************************
  mink_gdet_func():
  ----------
       -- Determinant of Minkowski covariant metric:
 ******************************************************************************/
void mink_spherical_gdet_func(double *x, double *gdet)
{
  int i;
  double r,th;
  double t1, cthtmp, gdet_tmp;
  
  /* Set up the coordinates : */ 
  r  = x[RR];
  th = x[TH];
  trigvalues( th, &cthtmp, &t1);
  //t1 = sin(th);

  *gdet = r*r*fabs(t1); 
  return;
}



/******************************************************************************
  mink_oonn_func():
  -----------------
       -- general routine to calculate the connection coefficients
       -- (need to implement xp,x  coordinate transformation)
 ******************************************************************************/
void mink_spherical_conn_func(struct of_coord *coords, double ***connp)
{
  int i,j,k;
  double t2, t3,  t5;
  double r, th;
  double conn[NDIM][NDIM][NDIM];
  double *x = coords->x; 
  
  DLOOP3 {  conn[i][j][k] = 0. ; } 

  r  = x[RR];
  th = x[TH];
  trigvalues( th, &t2, &t5);
  //  t2 = cos(th);
  //  t5 = sin(th);

  t3 = t2*t2;
  conn[1][2][2] = -r;
  conn[1][3][3] = r*(-1.0+t3);
  conn[2][1][2] = 1./r;
  conn[2][2][1] = conn[2][1][2];
  conn[2][3][3] = -t5*t2;
  conn[3][1][3] = conn[2][2][1];
  conn[3][2][3] = 1./t5*t2;
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];

//  fprintf(stdout,"CONN "); 
//  DLOOP3 {  fprintf(stdout,"%10.4e ",conn[i][j][k]); }
//  fprintf(stdout,"\n "); 

  transform_connection2(coords,conn,connp);

  return;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     M I N K O W S K I  (cartesian)  M E T R I C    R O U T I N E S

 *****************************************************************************/

/******************************************************************************
  mink_cartesian_gcov_func():
  ----------
       -- Minkowski form of the covariant metric:
 ******************************************************************************/
void mink_cartesian_gcov_func(double *x, double gcov[][NDIM])
{
  int i;
  
#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  gcov[i][j] = mink[i][j] ; } 
#else
  for(i = 0; i < NDIM*NDIM; i++) gcov[0][i] = mink[0][i] ; 
#endif


  return;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     K E R R - S C H I L D   (cartesian)  M E T R I C    R O U T I N E S

 *****************************************************************************/

/******************************************************************************
  ks_cart_gcov_gcon_func():
  ----------
  -- calculates both the contravariant and covariant forms of the 
            metric in Cartesian Kerr-Schild coordinates
 ******************************************************************************/
void ks_cart_gcov_gcon_func(double *x, double gcov[][NDIM], double gcon[][NDIM])
{
  int i, j;
  double xx,yy,zz,r,rsq;
  double cth,  t3;
  double two_Hf, lcov[NDIM];

  /* Set up the coordinates : */ 
  xx = x[XX] - x0_bh; 
  yy = x[YY] - y0_bh; 
  zz = x[ZZ] - z0_bh; 
  t3 = 0.5*(xx*xx + yy*yy + zz*zz - asq); 
  rsq = t3 + sqrt( t3*t3 + asq*zz*zz );
  r   = sqrt(rsq);
  cth = zz / r; 

  t3 = 1./(rsq + asq); 
  
  lcov[TT] = 1.;
  lcov[XX] = t3 * ( r * xx + a*yy );
  lcov[YY] = t3 * ( r * yy - a*xx );
  lcov[ZZ] = cth; 

  two_Hf = 2.*M*r / ( rsq + asq*cth*cth); 

  DLOOP2 {  gcov[i][j] = mink[i][j] + two_Hf * lcov[i] * lcov[j] ;  }

  lcov[TT] = -1.; 
  DLOOP2 {  gcon[i][j] = mink[i][j] - two_Hf * lcov[i] * lcov[j] ;  }

  return;
}


/******************************************************************************
  ks_cart_gcov_func():
  ----------
       -- covariant form of the metric in Cartesian Kerr-Schild coordinates
 ******************************************************************************/
void ks_cart_gcov_func(double *x, double gcov[][NDIM])
{
  int i, j;
  double xx,yy,zz,r,rsq;
  double cth,  t3;
  double two_Hf, lcov[NDIM];

  /* Set up the coordinates : */ 
  xx = x[XX] - x0_bh; 
  yy = x[YY] - y0_bh; 
  zz = x[ZZ] - z0_bh; 
  t3 = 0.5*(xx*xx + yy*yy + zz*zz - asq); 
  rsq = t3 + sqrt( t3*t3 + asq*zz*zz );
  r   = sqrt(rsq);
  cth = zz / r; 

  t3 = 1./(rsq + asq); 
  
  lcov[TT] = 1.;
  lcov[XX] = t3 * ( r * xx + a*yy );
  lcov[YY] = t3 * ( r * yy - a*xx );
  lcov[ZZ] = cth; 

  two_Hf = 2.*M*r / ( rsq + asq*cth*cth); 

  DLOOP2 {  gcov[i][j] = mink[i][j] + two_Hf * lcov[i] * lcov[j] ;  }

  return;
}


/******************************************************************************
  ks_cart_gcon_func():
  ----------
       -- contravariant form of the metric in Cartesian Kerr-Schild coordinates
 ******************************************************************************/
void ks_cart_gcon_func(double *x, double gcon[][NDIM])
{
  int i, j;
  double xx,yy,zz,r,rsq;
  double cth,  t3;
  double two_Hf, lcon[NDIM];

  /* Set up the coordinates : */ 
  xx = x[XX] - x0_bh; 
  yy = x[YY] - y0_bh; 
  zz = x[ZZ] - z0_bh; 
  t3 = 0.5*(xx*xx + yy*yy + zz*zz - asq); 
  rsq = t3 + sqrt( t3*t3 + asq*zz*zz );
  r   = sqrt(rsq);
  cth = zz / r; 

  t3 = 1./(rsq + asq); 
  
  lcon[TT] = -1.;
  lcon[XX] = t3 * ( r * xx + a*yy );
  lcon[YY] = t3 * ( r * yy - a*xx );
  lcon[ZZ] = cth; 

  two_Hf = 2.*M*r / ( rsq + asq*cth*cth); 

  DLOOP2 {  gcon[i][j] = mink[i][j] - two_Hf * lcon[i] * lcon[j] ;  }

  return;
}



/******************************************************************************
  ks_cart_gdet_func():
  ----------
   -- for Cartesian KS coordinates;
 ******************************************************************************/
void ks_cart_gdet_func(double *x, double *gdet )
{
  *gdet = 1.; 
  return; 
}

/******************************************************************************
  ks_cart_to_ks_spher_pos():
  ----------
   -- transforms the position in Cartesian KS coordinates to Spherical KS coordinates;
 ******************************************************************************/
void ks_cart_to_ks_spher_pos(double *x_cart, double *x_spher)
{
  double xx,yy,zz,r,t3;

  xx = x_cart[XX] - x0_bh; 
  yy = x_cart[YY] - y0_bh; 
  zz = x_cart[ZZ] - z0_bh; 

  t3 = 0.5*(xx*xx + yy*yy + zz*zz - asq); 
  r = sqrt(   t3 + sqrt( t3*t3 + asq*zz*zz )  );

  x_spher[TT] = x_cart[TT];
  x_spher[RR] = r; 
  x_spher[TH] = acos(zz / r); 
  t3 = atan2( (yy*r - xx*a) , (xx*r + yy*a) ); 
  if( t3 < 0. ) {  t3 += 2.*M_PI; } 
  x_spher[PH] = t3;

  return;
}


/******************************************************************************
  ks_dxc_dxs_calc():
  ----------------------
       -- calculates the transformation matrix  Lambda^\hat{a}_a  defined:

      x^\hat{a}[Cartesian]  = \Lambda^\hat{a}_a  x^a[Spherical] 

            where dxc_dxs[i][j] = \Lambda^i_j 

           for Kerr-Schild coordinates 

 ******************************************************************************/
void ks_dxc_dxs_calc(double *x_cart, double *x_spher, double dxc_dxs[][NDIM] )
{
  int i; 
  double r, th, ph; 
  double sth,cth,sph,cph;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  dxc_dxs[i][j] = 0.; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  dxc_dxs[0][i] = 0. ; } 
#endif
  
  r  = x_spher[RR]; 
  th = x_spher[TH];
  ph = x_spher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

// #if(COORDSINGFIX)
//   if (fabs(sth) < SINGSMALL) {
//     if((sth)>=0) sth =  SINGSMALL;
//     if((sth)<0)  sth = -SINGSMALL;
//   }
// #endif
  
  dxc_dxs[TT][TT] = 1.                    ;  // dt/dt
  dxc_dxs[XX][RR] = cph * sth             ;  // dx/dr
  dxc_dxs[XX][TH] = (r*cph - a*sph) * cth ;  // dx/dtheta
  dxc_dxs[XX][PH] = -(x_cart[YY] - y0_bh) ;  // dx/dphi
  dxc_dxs[YY][RR] = sph * sth             ;  // dy/dr
  dxc_dxs[YY][TH] = (r*sph + a*cph) * cth ;  // dy/dtheta
  dxc_dxs[YY][PH] = x_cart[XX] - x0_bh    ;  // dy/dphi
  dxc_dxs[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dxs[ZZ][TH] = -r*sth                ;  // dz/dtheta
  //  dxc_dxs[ZZ][PH] = 0.                    ;  // dz/dphi
  
  return;

}

/******************************************************************************
  ks_dxs_dxc_calc();
  -----------------------
       -- calculates the transformation matrix  Lambda^a_\hat{a}  defined:

      x^a[Spherical]  = \Lambda^a_\hat{a}  x^\hat{a}[Cartesian] 

            where dxs_dxc[i][j] = \Lambda^i_j 

           for Kerr-Schild coordinates 

******************************************************************************/
void ks_dxs_dxc_calc(double *x_spher, double *x_cart, double dxs_dxc[][NDIM] )
{
  int i; 
  double r, th, ph, xx, yy, zz; 
  double rsq, rhosq; 
  double dr_dx, dr_dy, dr_dz;
  double term_dr, term_dth, term_dph1, term_dph2;

#if( USE_STRICT_ARRAY_BOUNDS ) 
  int j;
  DLOOP2 {  dxs_dxc[i][j] = 0.; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  dxs_dxc[0][i] = 0. ; } 
#endif
  
  r  = x_spher[RR]; 
  th = x_spher[TH];
  ph = x_spher[PH]; 
  xx = x_cart[XX] - x0_bh; 
  yy = x_cart[YY] - y0_bh; 
  zz = x_cart[ZZ] - z0_bh; 

  rsq = r*r; 
  rhosq = xx*xx + yy*yy + zz*zz ;
  term_dr = 1./( r * (2.*rsq -  rhosq + asq  ) );
  term_dph1 = (rsq+asq);
  dr_dx  = xx * rsq       * term_dr; 
  dr_dy  = yy * rsq       * term_dr; 
  dr_dz  = zz * term_dph1 * term_dr; 

  term_dth =  1./(r*sqrt(rsq-zz*zz));
  term_dph1 = a/term_dph1;
  term_dph2 = 1./(xx*xx + yy*yy);

  dxs_dxc[TT][TT] = 1.                                 ;  // dt/dt
  dxs_dxc[RR][XX] = dr_dx                              ;  // dr/dx
  dxs_dxc[RR][YY] = dr_dy                              ;  // dr/dy
  dxs_dxc[RR][ZZ] = dr_dz                              ;  // dr/dz
  dxs_dxc[TH][XX] = term_dth * zz * dr_dx                   ;  // dth/dx
  dxs_dxc[TH][YY] = term_dth * zz * dr_dy                   ;  // dth/dy
  dxs_dxc[TH][ZZ] = term_dth * ( zz * dr_dz  -  r )    ;  // dth/dz
  dxs_dxc[PH][XX] = term_dph1 * dr_dx  - yy * term_dph2;  // dph/dx
  dxs_dxc[PH][YY] = term_dph1 * dr_dy  + xx * term_dph2;  // dph/dy
  dxs_dxc[PH][ZZ] = term_dph1 * dr_dz                  ;  // dph/dz
  
  return;

}

/******************************************************************************
  ks_cart_to_ks_spher_cov():
  ----------
       -- transforms a covariant vector in Cartesian KS coordinates to Spherical KS coordinates;
 ******************************************************************************/
void ks_cart_to_ks_spher_cov(double *x_cart, double *x_spher, double ks_cart_cov[], double ks_spher_cov[] )
{
  int i, j ; 
  double dxc_dxs[NDIM][NDIM];

  ks_dxc_dxs_calc(x_cart, x_spher, dxc_dxs);
  DLOOP1 { ks_spher_cov[i] = 0.; } 
  DLOOP2 { ks_spher_cov[i] += dxc_dxs[j][i] * ks_cart_cov[j] ; }
  
  return;
}


/******************************************************************************
  ks_cart_to_ks_spher_con():
  ----------
       -- transforms a contravariant vector in Cartesian KS coordinates to Spherical KS coordinates;
 ******************************************************************************/
void ks_cart_to_ks_spher_con(double *x_cart, double *x_spher, double ks_cart_con[], double ks_spher_con[] )
{
  int i, j ; 
  double dxs_dxc[NDIM][NDIM];

  ks_dxs_dxc_calc(x_spher, x_cart, dxs_dxc);
  DLOOP1 { ks_spher_con[i] = 0.; } 
  DLOOP2 { ks_spher_con[i] += dxs_dxc[i][j] * ks_cart_con[j] ; }
  
  return;
}

/******************************************************************************
  ks_spher_to_ks_cart_pos():
  ----------
       -- transforms a covariant vector in Spherical KS coordinates to Cartesian KS coordinates;
 ******************************************************************************/
void ks_spher_to_ks_cart_pos(double *x_spher, double *x_cart)
{
  double r,th,ph,sth,cth,sph,cph;

  r  = x_spher[RR]; 
  th = x_spher[TH];
  ph = x_spher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

//  #if(COORDSINGFIX)
//    if (fabs(sth) < SINGSMALL) {
//      if((sth)>=0) sth =  SINGSMALL;
//      if((sth)<0)  sth = -SINGSMALL;
//    }
//  #endif

  x_cart[TT] = x_spher[TT];
  x_cart[XX] = (r * cph - a * sph) * sth   + x0_bh ; 
  x_cart[YY] = (r * sph + a * cph) * sth   + y0_bh ; 
  x_cart[ZZ] =          r          * cth   + z0_bh ;
  
  return;
}

/******************************************************************************
  ks_spher_to_ks_cart_cov():
  ----------
       -- transforms a covariant vector in Spherical KS coordinates to Cartesian KS coordinates;
 ******************************************************************************/
void ks_spher_to_ks_cart_cov(double *x_spher, double *x_cart, double ks_spher_cov[], double ks_cart_cov[] )
{
  int i, j ; 
  double dxs_dxc[NDIM][NDIM];

  ks_dxs_dxc_calc(x_spher, x_cart, dxs_dxc);
  DLOOP1 { ks_cart_cov[i] = 0.; } 
  DLOOP2 { ks_cart_cov[i] += dxs_dxc[j][i] * ks_spher_cov[j] ; }
  
  return;
}


/******************************************************************************
  ks_spher_to_ks_cart_con():
  ----------
       -- transforms a contravariant vector in Spherical KS coordinates to Cartesian KS coordinates;
 ******************************************************************************/
void ks_spher_to_ks_cart_con(double *x_spher, double *x_cart, double ks_spher_con[], double ks_cart_con[] )
{
  int i, j ; 
  double dxc_dxs[NDIM][NDIM];

  ks_dxc_dxs_calc(x_cart, x_spher, dxc_dxs);
  DLOOP1 { ks_cart_con[i] = 0.; } 
  DLOOP2 { ks_cart_con[i] += dxc_dxs[i][j] * ks_spher_con[j] ; }
  
  return;
}

#if( USE_STRICT_ARRAY_BOUNDS )
/******************************************************************************
  ks_cart_conn_func():
  ----------
    -- calculates the connection coefficients for Cartesian KS coordinates;
    -- responsible for transforming it finally to numerical coordinates; 
******************************************************************************/
void ks_cart_conn_func(struct of_coord *coords, double ***connp )
{
  int i, j, k;
  double conn[NDIM][NDIM][NDIM];

  double t1,   t105,   t11,   t116,   t123,   t126,   t127,   t13,   t130,   t132,   t135,   t136,   t137; 
  double   t14,   t140,   t141,   t142,   t144,   t148,   t151,   t153,   t157,   t16,   t164,   t170, t177;    
  double   t178,   t179,   t180,   t181,   t182,   t183,   t185,   t188,   t189,   t19,   t190,   t191,   t192;   
  double   t194,   t195,   t196,   t197,   t2,   t200,   t202,   t203,   t204,   t205,   t207,   t208,   t209; 
  double   t213,   t215,   t216,   t217,   t218,   t22,   t222,   t223,   t224,   t225,   t226,   t227,   t231; 
  double   t233,   t237,   t238,   t239,   t240,   t244,   t245,   t247,   t248,   t249,   t250,   t254,   t256; 
  double   t26,   t260,   t262,   t264,   t266,   t267,   t269,   t272,   t274,   t276,   t277,   t28,   t281; 
  double   t282,   t283,  t285,   t288,   t289,   t290,   t292,   t295,   t296,   t297,   t298,   t299,   t3; 
  double   t300,   t301,   t305,   t307,   t308,   t311,   t312,   t314,   t316,   t319,   t32,   t320,   t321; 
  double   t323,   t326,   t327,   t328,   t329,   t33,   t330,   t331,   t332,   t335,   t336,   t338,   t339;   
  double   t340,   t341,   t344,   t346,   t352,   t354,   t356,   t357,   t359,   t362,   t364,   t366,   t367; 
  double   t370,   t372,   t374,   t375,   t378,   t379,   t38,   t380,   t382,   t385,   t386,   t387,   t388; 
  double   t389,   t390,   t393,   t394,   t396,   t397,   t398,   t399,   t4,   t403,   t405,   t41,   t410; 
  double   t412,   t414,   t415,   t417,   t420,   t421,   t423,   t427,   t428,   t43,   t439,   t440,   t446; 
  double   t451,   t452,   t458,   t473,   t486,   t487,   t488,   t5,   t505,   t506,   t507,   t525,   t528; 
  double   t531,   t535,   t538,   t542,   t545,   t549,   t56,   t57,   t575,   t580,   t584,   t587,   t59;
  double   t591,   t594,   t598,   t61,   t613,   t617,   t624,   t63,   t631,   t639,   t644,   t65,   t650; 
  double   t654,   t657,   t668,   t67,   t69,   t7,   t70,   t71,   t72,   t722,   t724,   t725,   t74;
  double   t746,   t76,   t775,   t83,   t84,   t85,   t86,   t9,   t92,   t96,   t99;

  double x,y,z;
  double rhosq,rf,rsq,Hf,lx,ly,lz;
  double s1,s2,s3,s4;
  {
    x = coords->x[XX] - x0_bh ; 
    y = coords->x[YY] - y0_bh ; 
    z = coords->x[ZZ] - z0_bh ; 

    t1 = x*x;
    t2 = y*y;
    t3 = z*z;
    t4 = t1+t2+t3;  // rhosq
    t5 = a*a;
    s1 = t4-t5;
    //    t7 = pow(t4-t5,2.0);
    t7 = s1*s1;
    t9 = 4.0*t5*t3;
    t11 = sqrt(t7+t9);
    t41 = t11+s1;
    t13 = 0.5*t41; // rsq
    t14 = sqrt(t13);             // r 
    //    t16 = z/t14;
    t56 = 1/t14;
    t16 = z*t56;
    t26 = acos(t16);   // theta 
    t19 = y*t14-x*a;
    t22 = x*t14+y*a;
    //    eqs2[4] = atan(t19/t22);  // phi
    //    t27 = cos(t26);
    t28 = t16*t16;
    t188 = M*t14/(t13+t5*t28);   // H
    t32 = t13+t5;
    t33 = 1/t32;
    t178 = t22*t33;              // lx
    t181 = t19*t33;              // ly 
    //    t184 = t16;                  // lz
    //    t35 = t4*t4;
    t38 = t5*t5;
    //    t40 = t11;
    //    t41 = t11+t4-t5;
    //    t34 = 2.*x;  // drho_dx
    //    t45 = 2.*y;  // drho_dy
    //    t48 = 2.*z;  // drho_dz
    t43 = 1./t11;
    t57 = x*t41*t43;  // drsq_dx
    t59 = y*t41*t43;  // drsq_dy
    t61 = (t41+2.0*t5)*z*t43;  // drsq_dz
    t63 = t56*t57*0.5;  // dr_dx
    t72 = t56*t59*0.5;  // dr_dy
    t74 = t56*t61*0.5;  // dr_dz
    t65 = t14*t14;
    t67 = 1/t65;
    t69 = sqrt((t65-t3)*t67);
    t70 = 1/t69;
    t71 = t70*t67;
    t86  = z*t63*t71;  // dth_dx
    t105 = z*t72*t71;  // dth_dy
    t76  = t14-z*t74;
    t116 = -t76*t70*t67;  // dth_dz
    t83 = t14*t5;
    t84 = sin(t26);
    t85 = t16*t84;
    t92 = t13*t13;
    t96 = t28*t28;
    t99 = -M/(t92+2.0*t13*t5*t28+t38*t96);
    s2 = t13 + t5*t28;
    s3 = 2.*t83*t85;
    s4 = t13 + t5;
    t177 = (-t63*s2+t14*t57-s3*t86 )*t99;   // dH_dx
    t180 = (-t72*s2+t14*t59-s3*t105)*t99;  // dH_dy
    t183 = (-t74*s2+t14*t61-s3*t116)*t99;  // dH_dz
    t123 = z*t67;
    t194 = -t123*t63;  // dlz_dx 
    t215 = -t123*t72;  // dlz_dy
    t393 = t76*t67;    // dlz_dz
    t126 = t14*t13;
    t127 = x*t63;
    t130 = t57*x;
    t132 = t57*y;
    t135 = t32*t32;
    t136 = 1/t135;
    t244 = (t126+t83+t127*(t13+t5)-t130*t14-t132*a)*t136;  // dlx_dx
    t137 = x*t72;
    t140 = a*t13;
    t141 = t5*a;
    t142 = t59*x;
    t144 = t59*y;
    t189 = -(-t137*s4-t140-t141+t142*t14+t144*a)*t136;  // dlx_dy
    t148 = x*t74;
    t151 = t61*x;
    t153 = t61*y;
    t195 = -(-t148*s4+t151*t14+t153*a)*t136;            // dlx_dz
    t157 = y*t63;
    t190 = -(-t157*s4+t140+t141+t132*t14-t130*a)*t136;  // dly_dx
    t164 = y*t72;
    t335 = (t126+t83+t164*s4-t144*t14+t142*a)*t136;     // dly_dy
    t170 = y*t74;
    t216 = -(-t170*s4+t153*t14-t151*a)*t136;            // dly_dz
    t179 = t177*t178;
    t182 = t180*t181;
    t185 = t183*t16;
    conn[0][0][ 0] = -(t179+t182+t185)*2.*t188;
    t191 = -t189+t190;
    t192 = 2.0*t191;
    t196 = t194-t195;
    t197 = 2.0*t196;
    t200 = t188*t188;
    t202 = t181*t181;
    t203 = t177*t202;
    t204 = t180*t178;
    t205 = t204*t181;
    t207 = t177*t28;
    t208 = t183*t178;
    t209 = t208*t16;
    conn[0][0][ 1] = (t192*t181+t197*t16)*t200+2.*(-t177+t203-t205+t207-t209)*t188-t177;
    t213 = -t192;
    t217 = t215-t216;
    t218 = 2.0*t217;
    t222 = t178*t178;
    t223 = t222*t180;
    t224 = t179*t181;
    t225 = t28*t180;
    t226 = t183*t181;
    t227 = t226*t16;
    conn[0][0][ 2] = (t213*t178+t218*t16)*t200+2.*(-t180+t223-t224+t225-t227)*t188-t180;
    t231 = -t197;
    t233 = -t218;
    t237 = t222*t183;
    t238 = t179*t16;
    t239 = t202*t183;
    t240 = t182*t16;
    conn[0][0][ 3] = (t231*t178+t233*t181)*t200+2.*(-t183+t237-t238+t239-t240)*t188-t183;
    conn[0][1][ 0] = conn[0][0][ 1];
    t245 = t222*t244;
    t247 = 4.0*t191;
    t248 = t247*t181;
    t249 = 4.0*t196;
    t250 = t249*t16;
    t254 = t202*t244;
    t256 = t28*t244;
    t260 = t222*t178;
    t262 = 2.0*t260*t177;
    t264 = -2.0*t182-2.0*t185;
    t266 = -t177+t207+t203;
    t267 = 4.0*t266;
    t269 = 2.0*t244;
    t272 = 2.0*t179;
    conn[0][1][ 1] = (t245+(t248+t250)*0.25*t178-t244+t254+t256)*4.*t200
      +(t262+t264*t222+t267*t178-t269)*t188-t272;
    t274 = 2.0*t189+2.0*t190;
    t276 = t218*t178;
    t277 = t197*t181;
    t281 = 4.0*t222*t189;
    t282 = 2.0*t190;
    t283 = 2.0*t189;
    t285 = 4.0*t202*t190;
    t288 = t177*t181;
    t289 = t204+t288;
    t290 = 2.0*t289;
    t292 = t178*t181;
    t295 = 2.0*t204;
    t296 = t202*t181;
    t297 = t296*t177;
    t298 = 2.0*t297;
    t299 = 2.0*t288;
    t300 = t260*t180;
    t301 = 2.0*t300;
    conn[0][1][ 2] = (t274*t28+(t276+t277)*t16+t281-t282-t283+t285)*t200
      +(t290*t28-2.0*t185*t292-t189-t190-t295+t298-t299+t301)*t188-t204-t288;
    t305 = 2.0*t195+2.0*t194;
    t307 = t233*t178;
    t308 = t192*t16;
    t311 = 2.0*t195;
    t312 = 2.0*t194;
    t314 = 4.0*t222*t195;
    t316 = 4.0*t28*t194;
    t319 = t177*t16;
    t320 = t208+t319;
    t321 = 2.0*t320;
    t323 = t178*t16;
    t326 = 2.0*t208;
    t327 = 2.0*t319;
    t328 = t28*t16;
    t329 = t328*t177;
    t330 = 2.0*t329;
    t331 = t260*t183;
    t332 = 2.0*t331;
    conn[0][1][ 3] = (t305*t202+(t307+t308)*t181-t311-t312+t314+t316)*t200
      +(t321*t202-2.0*t182*t323-t195-t194-t326-t327+t330+t332)*t188-t208-t319;
    conn[0][2][ 0] = conn[0][0][ 2];
    conn[0][2][ 1] = conn[0][1][ 2];
    t336 = t202*t335;
    t338 = -t247;
    t339 = t338*t178;
    t340 = 4.0*t217;
    t341 = t340*t16;
    t344 = t28*t335;
    t346 = t222*t335;
    t352 = 2.0*t296*t180;
    t354 = -2.0*t185-2.0*t179;
    t356 = -t180+t225+t223;
    t357 = 4.0*t356;
    t359 = 2.0*t335;
    t362 = 2.0*t182;
    conn[0][2][2] = (t336+(t339+t341)*0.25*t181+t344+t346-t335)*4.*t200
      +(t352+t354*t202+t357*t181-t359)*t188-t362;
    t364 = 2.0*t215+2.0*t216;
    t366 = t231*t181;
    t367 = t213*t16;
    t370 = 2.0*t216;
    t372 = 4.0*t28*t215;
    t374 = 4.0*t202*t216;
    t375 = 2.0*t215;
    t378 = t180*t16;
    t379 = t226+t378;
    t380 = 2.0*t379;
    t382 = t181*t16;
    t385 = 2.0*t226;
    t386 = t296*t183;
    t387 = 2.0*t386;
    t388 = 2.0*t378;
    t389 = t328*t180;
    t390 = 2.0*t389;
    conn[0][2][3] = (t364*t222+(t366+t367)*t178-t370+t372+t374-t375)*t200
      +(t380*t222-2.0*t179*t382-t216-t215-t385+t387-t388+t390)*t188-t226-t378;
    conn[0][3][0] = conn[0][0][ 3];
    conn[0][3][1] = conn[0][1][ 3];
    conn[0][3][2] = conn[0][2][ 3];
    t394 = t28*t393;
    t396 = -t249;
    t397 = t396*t178;
    t398 = -t340;
    t399 = t398*t181;
    t403 = t202*t393;
    t405 = t222*t393;
    t410 = 2.0*t328*t183;
    t412 = -2.0*t179-2.0*t182;
    t414 = -t183+t237+t239;
    t415 = 4.0*t414;
    t417 = 2.0*t393;
    t420 = 2.0*t185;
    conn[0][3][3] = (t394+(t397+t399)*0.25*t16-t393+t403+t405)*4.*t200
      +(t410+t412*t28+t415*t16-t417)*t188-t420;
    t421 = t177*t222;
    t423 = -t264;
    conn[1][0][0] = (2.0*t421+t423*t178)*t188-t177;
    t427 = t213*t181;
    t428 = t231*t16;
    conn[1][0][1] = (t427+t428)*t178*t200+(t423*t222-2.0*t266*t178)*t188;
    t439 = t233*t16;
    t440 = t439*t178;
    t446 = 2.0*t180-2.0*t225+2.0*t227;
    conn[1][0][2] = (t192*t222+t440)*t200+(-t301+2.0*t421*t181+t446*t178-t190+t189)*t188+t204-t288;
    t451 = t218*t181;
    t452 = t451*t178;
    t458 = 2.0*(t183-t239+t240);
    conn[1][0][3] = (t197*t222+t452)*t200+(-t332+2.0*t421*t16+t458*t178+t195-t194)*t188+t208-t319;
    conn[1][1][0] = conn[1][0][1];
    t473 = t222*t222;
    conn[1][1][1] = (-4.0*t244*t260+(t338*t181+t396*t16)*t222
		+4.*(t244-t256-t254)*t178)*t200
      +(-2.0*t177*t473+t423*t260-t267*t222+2.0*t244*t178)*t188+t421;
    t486 = t428*t181;
    t487 = -t274;
    t488 = t487*t28;
    conn[1][1][2] = (-4.0*t189*t260+t439*t222+(-t285+t486+t488+t283+t282)*t178)*t200
      +(-2.0*t180*t473+t446*t222+2.*(t189+t288-t297-t207*t181)*t178)*t188+t223;
    t505 = -t305;
    t506 = t505*t202;
    t507 = t367*t181;
    conn[1][1][3] = (-4.0*t195*t260+t451*t222+(t506+t507+t312+t311-t316)*t178)*t200+
      (-2.0*t183*t473+t458*t222+2.*(t195+t319-t203*t16-t329)*t178)*t188+t237;
    conn[1][2][0] = conn[1][0][2];
    conn[1][2][1] = conn[1][1][2];
    t525 = -4.0*(t188*t181*t180+t200*t335);
    t528 = t188*t177;
    t531 = t248*t200+2.0*t528*t202;
    t535 = t188*t183;
    t538 = t399*t200+2.0*t535*t202;
    t542 = 4.*(-t336+t335)*t200;
    t545 = (t359+4.0*t182-t352)*t188;
    conn[1][2][2] = t525*t260+t531*t222+(t525*t28+t538*t16+t542+t545+t362)*t178+t427*t188-t203;
    t549 = -t364;
    conn[1][2][3] = (t549*t260+(t277+t308)*t222+(t375+t370-t374-t372)*t178)*t200
      +(-t380*t260+2.0*t421*t382+(t216+t385+t215-t387-t390+t388)*t178-t196*t181-t191*t16)*t188
      +t379*t178-t288*t16;
    conn[1][3][0] = conn[1][0][3];
    conn[1][3][1] = conn[1][1][3];
    conn[1][3][2] = conn[1][2][3];
    t575 = -4.0*(t188*t16*t183+t200*t393);
    t580 = t250*t200+2.0*t528*t28;
    t584 = t188*t180;
    t587 = t341*t200+2.0*t584*t28;
    t591 = 4.*(t393-t394)*t200;
    t594 = (4.0*t185+t417-t410)*t188;
    conn[1][3][3] = t575*t260+t580*t222+(t575*t202+t587*t181+t591+t594+t420)*t178+t428*t188-t207;
    t598 = t180*t202;
    conn[2][0][0] = 2.*(t224+t598+t227)*t188-t180;
    conn[2][0][1] = (t213*t202+t486)*t200+(2.*(t598+t227)*t178-t298
				      +2.*(t177-t207)*t181+t190-t189)*t188+t288-t204;
    t613 = t192*t178;
    t617 = -t354;
    conn[2][0][2] = (t613+t439)*t181*t200+(t617*t202-2.0*t356*t181)*t188;
    t624 = t277*t178;
    t631 = 2.0*t183-2.0*t237+2.0*t238;
    conn[2][0][3] = (t624+t218*t202)*t200+(-t387+2.0*t598*t16+t631*t181+t216-t215)*t188+t226-t378;
    conn[2][1][0] = conn[2][0][1];
    t639 = -4.0*(t188*t178*t177+t200*t244);
    t644 = t339*t200+2.0*t584*t222;
    t650 = t397*t200+2.0*t535*t222;
    t654 = 4.*(-t245+t244)*t200;
    t657 = (4.0*t179+t269-t262)*t188;
    conn[2][1][1] = t639*t296+t644*t202+(t639*t28+t650*t16+t654+t657+t272)*t181+t613*t188-t223;
    t668 = t202*t202;
    conn[2][1][2] = (-4.0*t190*t296+t428*t202+(-t281+t440+t488+t283+t282)*t181)*t200
      +(-2.0*t177*t668+2.*(t209+t177-t207)*t202+2.*(t204+t190-t300-t225*t178)*t181)*t188+t203;
    conn[2][1][3] = (t505*t296+(t276+t367)*t202+(t312-t314-t316+t311)*t181)*t200
      +(-t321*t296+2.0*t598*t323+(-t330-t332+t194+t326+t327+t195)*t181-t217*t178+t191*
	t16)*t188+t320*t181-t204*t16;
    conn[2][2][0] = conn[2][0][2];
    conn[2][2][1] = conn[2][1][2];
    conn[2][2][2] = (-4.0*t335*t296+(t247*t178+t398*t16)*t202
		+4.*(t335-t344-t346)*t181)*t200+(-2.0*t180*t668+t617*t296-t357*t202+2.0*t181*t335)*t188+t598;
    t722 = t197*t178;
    t724 = t549*t222;
    t725 = t308*t178;
    conn[2][2][3] = (-4.0*t216*t296+t722*t202+(t724+t725+t370+t375-t372)*t181)*t200+
      (-2.0*t183*t668+t631*t202+2.*(t216+t378-t223*t16-t389)*t181)*t188+t239;
    conn[2][3][0] = conn[2][0][3];
    conn[2][3][1] = conn[2][1][3];
    conn[2][3][2] = conn[2][2][3];
    conn[2][3][3] = t575*t296+t587*t202+(t575*t222+t580*t178+t591+t594+t420)*t181+t439*t188-t225;
    t746 = t183*t28;
    conn[3][0][0] = 2.*(t238+t240+t746)*t188-t183;
    conn[3][0][1] = (t507+t231*t28)*t200+(2.*(t240+t746)*t178-t330+2.*(t177-t203)*t16+t194-t195)*t188+t319-t208;
    conn[3][0][2] = (t725+t233*t28)*t200+(2.*(t238+t746)*t181-t390+2.*(t180-t223)*t16+t215-t216)*t188+t378-t226;
    t775 = -t412;
    conn[3][0][3] = (t722+t451)*t16*t200+(t775*t28-2.0*t414*t16)*t188;
    conn[3][1][0] = conn[3][0][1];
    conn[3][1][1] = t639*t328+t650*t28+(t639*t202+t644*t181+t654+t657+t272)*t16+t722*t188-t237;
    conn[3][1][2] = (t487*t328+(t307+t366)*t28+(t283-t285+t282-t281)*t16)*t200
      +(-t290*t328+2.0*t746*t292+(t299+t295+t190-t301+t189-t298)*t16+t217*t178+t196*
	t181)*t188+t289*t16-t208*t181;
    t96 = t28*t28;
    conn[3][1][3] = (-4.0*t194*t328+t213*t28*t181+(-t314+t452+t506+t311+t312)*t16)
      *t200+(-2.0*t177*t96+2.*(-t203+t177+t205)*t28+2.*(t208-t239*t178+t194-t331)*t16)*t188+t207;
    conn[3][2][0] = conn[3][0][2];
    conn[3][2][1] = conn[3][1][2];
    conn[3][2][2] = t525*t328+t538*t28+(t525*t222+t531*t178+t542+t545+t362)*t16+t451*t188-t239;
    conn[3][2][3] = (-4.0*t215*t328+t192*t28*t178+(t724+t624+t370+t375-t374)*t16)*
      t200+(-2.0*t180*t96+2.*(t224-t223+t180)*t28+2.*(t226-t386+t215-t237*t181)*t16)*t188+t225;
    conn[3][3][0] = conn[3][0][3];
    conn[3][3][1] = conn[3][1][3];
    conn[3][3][2] = conn[3][2][3];
    conn[3][3][3] = (-4.0*t393*t328+(t249*t178+t340*t181)*t28+4.*(-t405+t393-t403)*t16)*t200
      +(-2.0*t183*t96+t775*t328-t415*t28+2.0*t16*t393)*t188+t746;

    /*********************************************************************************************
       Up until this poin the connection is in non-primed coordinates.  
       Now we transform it to  primed/numerical coordinates. 
    *********************************************************************************************/

    transform_connection2(coords,conn,connp);

    return;
  }
}

#else
/******************************************************************************
  ks_cart_conn_func():
  ----------
    -- calculates the connection coefficients for Cartesian KS coordinates;
    -- responsible for transforming it finally to numerical coordinates; 
******************************************************************************/
void ks_cart_conn_func(struct of_coord *coords, double ***connp )
{
  int i, j, k;
  double conn[NDIM][NDIM][NDIM];

  double t1,   t105,   t11,   t116,   t123,   t126,   t127,   t13,   t130,   t132,   t135,   t136,   t137; 
  double   t14,   t140,   t141,   t142,   t144,   t148,   t151,   t153,   t157,   t16,   t164,   t170, t177;    
  double   t178,   t179,   t180,   t181,   t182,   t183,   t185,   t188,   t189,   t19,   t190,   t191,   t192;   
  double   t194,   t195,   t196,   t197,   t2,   t200,   t202,   t203,   t204,   t205,   t207,   t208,   t209; 
  double   t213,   t215,   t216,   t217,   t218,   t22,   t222,   t223,   t224,   t225,   t226,   t227,   t231; 
  double   t233,   t237,   t238,   t239,   t240,   t244,   t245,   t247,   t248,   t249,   t250,   t254,   t256; 
  double   t26,   t260,   t262,   t264,   t266,   t267,   t269,   t272,   t274,   t276,   t277,   t28,   t281; 
  double   t282,   t283,  t285,   t288,   t289,   t290,   t292,   t295,   t296,   t297,   t298,   t299,   t3; 
  double   t300,   t301,   t305,   t307,   t308,   t311,   t312,   t314,   t316,   t319,   t32,   t320,   t321; 
  double   t323,   t326,   t327,   t328,   t329,   t33,   t330,   t331,   t332,   t335,   t336,   t338,   t339;   
  double   t340,   t341,   t344,   t346,   t352,   t354,   t356,   t357,   t359,   t362,   t364,   t366,   t367; 
  double   t370,   t372,   t374,   t375,   t378,   t379,   t38,   t380,   t382,   t385,   t386,   t387,   t388; 
  double   t389,   t390,   t393,   t394,   t396,   t397,   t398,   t399,   t4,   t403,   t405,   t41,   t410; 
  double   t412,   t414,   t415,   t417,   t420,   t421,   t423,   t427,   t428,   t43,   t439,   t440,   t446; 
  double   t451,   t452,   t458,   t473,   t486,   t487,   t488,   t5,   t505,   t506,   t507,   t525,   t528; 
  double   t531,   t535,   t538,   t542,   t545,   t549,   t56,   t57,   t575,   t580,   t584,   t587,   t59;
  double   t591,   t594,   t598,   t61,   t613,   t617,   t624,   t63,   t631,   t639,   t644,   t65,   t650; 
  double   t654,   t657,   t668,   t67,   t69,   t7,   t70,   t71,   t72,   t722,   t724,   t725,   t74;
  double   t746,   t76,   t775,   t83,   t84,   t85,   t86,   t9,   t92,   t96,   t99;

  double x,y,z;
  double rhosq,rf,rsq,Hf,lx,ly,lz;
  double s1,s2,s3,s4;
  {
    x = coords->x[XX] - x0_bh ; 
    y = coords->x[YY] - y0_bh ; 
    z = coords->x[ZZ] - z0_bh ; 

    t1 = x*x;
    t2 = y*y;
    t3 = z*z;
    t4 = t1+t2+t3;  // rhosq
    t5 = a*a;
    s1 = t4-t5;
    //    t7 = pow(t4-t5,2.0);
    t7 = s1*s1;
    t9 = 4.0*t5*t3;
    t11 = sqrt(t7+t9);
    t41 = t11+s1;
    t13 = 0.5*t41; // rsq
    t14 = sqrt(t13);             // r 
    //    t16 = z/t14;
    t56 = 1/t14;
    t16 = z*t56;
    t26 = acos(t16);   // theta 
    t19 = y*t14-x*a;
    t22 = x*t14+y*a;
    //    eqs2[4] = atan(t19/t22);  // phi
    //    t27 = cos(t26);
    t28 = t16*t16;
    t188 = M*t14/(t13+t5*t28);   // H
    t32 = t13+t5;
    t33 = 1/t32;
    t178 = t22*t33;              // lx
    t181 = t19*t33;              // ly 
    //    t184 = t16;                  // lz
    //    t35 = t4*t4;
    t38 = t5*t5;
    //    t40 = t11;
    //    t41 = t11+t4-t5;
    //    t34 = 2.*x;  // drho_dx
    //    t45 = 2.*y;  // drho_dy
    //    t48 = 2.*z;  // drho_dz
    t43 = 1./t11;
    t57 = x*t41*t43;  // drsq_dx
    t59 = y*t41*t43;  // drsq_dy
    t61 = (t41+2.0*t5)*z*t43;  // drsq_dz
    t63 = t56*t57*0.5;  // dr_dx
    t72 = t56*t59*0.5;  // dr_dy
    t74 = t56*t61*0.5;  // dr_dz
    t65 = t14*t14;
    t67 = 1/t65;
    t69 = sqrt((t65-t3)*t67);
    t70 = 1/t69;
    t71 = t70*t67;
    t86  = z*t63*t71;  // dth_dx
    t105 = z*t72*t71;  // dth_dy
    t76  = t14-z*t74;
    t116 = -t76*t70*t67;  // dth_dz
    t83 = t14*t5;
    t84 = sin(t26);
    t85 = t16*t84;
    t92 = t13*t13;
    t96 = t28*t28;
    t99 = -M/(t92+2.0*t13*t5*t28+t38*t96);
    s2 = t13 + t5*t28;
    s3 = 2.*t83*t85;
    s4 = t13 + t5;
    t177 = (-t63*s2+t14*t57-s3*t86 )*t99;   // dH_dx
    t180 = (-t72*s2+t14*t59-s3*t105)*t99;  // dH_dy
    t183 = (-t74*s2+t14*t61-s3*t116)*t99;  // dH_dz
    t123 = z*t67;
    t194 = -t123*t63;  // dlz_dx 
    t215 = -t123*t72;  // dlz_dy
    t393 = t76*t67;    // dlz_dz
    t126 = t14*t13;
    t127 = x*t63;
    t130 = t57*x;
    t132 = t57*y;
    t135 = t32*t32;
    t136 = 1/t135;
    t244 = (t126+t83+t127*(t13+t5)-t130*t14-t132*a)*t136;  // dlx_dx
    t137 = x*t72;
    t140 = a*t13;
    t141 = t5*a;
    t142 = t59*x;
    t144 = t59*y;
    t189 = -(-t137*s4-t140-t141+t142*t14+t144*a)*t136;  // dlx_dy
    t148 = x*t74;
    t151 = t61*x;
    t153 = t61*y;
    t195 = -(-t148*s4+t151*t14+t153*a)*t136;            // dlx_dz
    t157 = y*t63;
    t190 = -(-t157*s4+t140+t141+t132*t14-t130*a)*t136;  // dly_dx
    t164 = y*t72;
    t335 = (t126+t83+t164*s4-t144*t14+t142*a)*t136;     // dly_dy
    t170 = y*t74;
    t216 = -(-t170*s4+t153*t14-t151*a)*t136;            // dly_dz
    t179 = t177*t178;
    t182 = t180*t181;
    t185 = t183*t16;
    conn[0][0][ 0] = -(t179+t182+t185)*2.*t188;
    t191 = -t189+t190;
    t192 = 2.0*t191;
    t196 = t194-t195;
    t197 = 2.0*t196;
    t200 = t188*t188;
    t202 = t181*t181;
    t203 = t177*t202;
    t204 = t180*t178;
    t205 = t204*t181;
    t207 = t177*t28;
    t208 = t183*t178;
    t209 = t208*t16;
    conn[0][0][ 1] = (t192*t181+t197*t16)*t200+2.*(-t177+t203-t205+t207-t209)*t188-t177;
    t213 = -t192;
    t217 = t215-t216;
    t218 = 2.0*t217;
    t222 = t178*t178;
    t223 = t222*t180;
    t224 = t179*t181;
    t225 = t28*t180;
    t226 = t183*t181;
    t227 = t226*t16;
    conn[0][0][ 2] = (t213*t178+t218*t16)*t200+2.*(-t180+t223-t224+t225-t227)*t188-t180;
    t231 = -t197;
    t233 = -t218;
    t237 = t222*t183;
    t238 = t179*t16;
    t239 = t202*t183;
    t240 = t182*t16;
    conn[0][0][ 3] = (t231*t178+t233*t181)*t200+2.*(-t183+t237-t238+t239-t240)*t188-t183;
    conn[0][0][ 4] = conn[0][0][ 1];
    t245 = t222*t244;
    t247 = 4.0*t191;
    t248 = t247*t181;
    t249 = 4.0*t196;
    t250 = t249*t16;
    t254 = t202*t244;
    t256 = t28*t244;
    t260 = t222*t178;
    t262 = 2.0*t260*t177;
    t264 = -2.0*t182-2.0*t185;
    t266 = -t177+t207+t203;
    t267 = 4.0*t266;
    t269 = 2.0*t244;
    t272 = 2.0*t179;
    conn[0][0][ 5] = (t245+(t248+t250)*0.25*t178-t244+t254+t256)*4.*t200
      +(t262+t264*t222+t267*t178-t269)*t188-t272;
    t274 = 2.0*t189+2.0*t190;
    t276 = t218*t178;
    t277 = t197*t181;
    t281 = 4.0*t222*t189;
    t282 = 2.0*t190;
    t283 = 2.0*t189;
    t285 = 4.0*t202*t190;
    t288 = t177*t181;
    t289 = t204+t288;
    t290 = 2.0*t289;
    t292 = t178*t181;
    t295 = 2.0*t204;
    t296 = t202*t181;
    t297 = t296*t177;
    t298 = 2.0*t297;
    t299 = 2.0*t288;
    t300 = t260*t180;
    t301 = 2.0*t300;
    conn[0][0][ 6] = (t274*t28+(t276+t277)*t16+t281-t282-t283+t285)*t200
      +(t290*t28-2.0*t185*t292-t189-t190-t295+t298-t299+t301)*t188-t204-t288;
    t305 = 2.0*t195+2.0*t194;
    t307 = t233*t178;
    t308 = t192*t16;
    t311 = 2.0*t195;
    t312 = 2.0*t194;
    t314 = 4.0*t222*t195;
    t316 = 4.0*t28*t194;
    t319 = t177*t16;
    t320 = t208+t319;
    t321 = 2.0*t320;
    t323 = t178*t16;
    t326 = 2.0*t208;
    t327 = 2.0*t319;
    t328 = t28*t16;
    t329 = t328*t177;
    t330 = 2.0*t329;
    t331 = t260*t183;
    t332 = 2.0*t331;
    conn[0][0][ 7] = (t305*t202+(t307+t308)*t181-t311-t312+t314+t316)*t200
      +(t321*t202-2.0*t182*t323-t195-t194-t326-t327+t330+t332)*t188-t208-t319;
    conn[0][0][ 8] = conn[0][0][ 2];
    conn[0][0][ 9] = conn[0][0][ 6];
    t336 = t202*t335;
    t338 = -t247;
    t339 = t338*t178;
    t340 = 4.0*t217;
    t341 = t340*t16;
    t344 = t28*t335;
    t346 = t222*t335;
    t352 = 2.0*t296*t180;
    t354 = -2.0*t185-2.0*t179;
    t356 = -t180+t225+t223;
    t357 = 4.0*t356;
    t359 = 2.0*t335;
    t362 = 2.0*t182;
    conn[0][0][10] = (t336+(t339+t341)*0.25*t181+t344+t346-t335)*4.*t200
      +(t352+t354*t202+t357*t181-t359)*t188-t362;
    t364 = 2.0*t215+2.0*t216;
    t366 = t231*t181;
    t367 = t213*t16;
    t370 = 2.0*t216;
    t372 = 4.0*t28*t215;
    t374 = 4.0*t202*t216;
    t375 = 2.0*t215;
    t378 = t180*t16;
    t379 = t226+t378;
    t380 = 2.0*t379;
    t382 = t181*t16;
    t385 = 2.0*t226;
    t386 = t296*t183;
    t387 = 2.0*t386;
    t388 = 2.0*t378;
    t389 = t328*t180;
    t390 = 2.0*t389;
    conn[0][0][11] = (t364*t222+(t366+t367)*t178-t370+t372+t374-t375)*t200
      +(t380*t222-2.0*t179*t382-t216-t215-t385+t387-t388+t390)*t188-t226-t378;
    conn[0][0][12] = conn[0][0][ 3];
    conn[0][0][13] = conn[0][0][ 7];
    conn[0][0][14] = conn[0][0][11];
    t394 = t28*t393;
    t396 = -t249;
    t397 = t396*t178;
    t398 = -t340;
    t399 = t398*t181;
    t403 = t202*t393;
    t405 = t222*t393;
    t410 = 2.0*t328*t183;
    t412 = -2.0*t179-2.0*t182;
    t414 = -t183+t237+t239;
    t415 = 4.0*t414;
    t417 = 2.0*t393;
    t420 = 2.0*t185;
    conn[0][0][15] = (t394+(t397+t399)*0.25*t16-t393+t403+t405)*4.*t200
      +(t410+t412*t28+t415*t16-t417)*t188-t420;
    t421 = t177*t222;
    t423 = -t264;
    conn[0][0][16] = (2.0*t421+t423*t178)*t188-t177;
    t427 = t213*t181;
    t428 = t231*t16;
    conn[0][0][17] = (t427+t428)*t178*t200+(t423*t222-2.0*t266*t178)*t188;
    t439 = t233*t16;
    t440 = t439*t178;
    t446 = 2.0*t180-2.0*t225+2.0*t227;
    conn[0][0][18] = (t192*t222+t440)*t200+(-t301+2.0*t421*t181+t446*t178-t190+t189)*t188+t204-t288;
    t451 = t218*t181;
    t452 = t451*t178;
    t458 = 2.0*(t183-t239+t240);
    conn[0][0][19] = (t197*t222+t452)*t200+(-t332+2.0*t421*t16+t458*t178+t195-t194)*t188+t208-t319;
    conn[0][0][20] = conn[0][0][17];
    t473 = t222*t222;
    conn[0][0][21] = (-4.0*t244*t260+(t338*t181+t396*t16)*t222
		+4.*(t244-t256-t254)*t178)*t200
      +(-2.0*t177*t473+t423*t260-t267*t222+2.0*t244*t178)*t188+t421;
    t486 = t428*t181;
    t487 = -t274;
    t488 = t487*t28;
    conn[0][0][22] = (-4.0*t189*t260+t439*t222+(-t285+t486+t488+t283+t282)*t178)*t200
      +(-2.0*t180*t473+t446*t222+2.*(t189+t288-t297-t207*t181)*t178)*t188+t223;
    t505 = -t305;
    t506 = t505*t202;
    t507 = t367*t181;
    conn[0][0][23] = (-4.0*t195*t260+t451*t222+(t506+t507+t312+t311-t316)*t178)*t200+
      (-2.0*t183*t473+t458*t222+2.*(t195+t319-t203*t16-t329)*t178)*t188+t237;
    conn[0][0][24] = conn[0][0][18];
    conn[0][0][25] = conn[0][0][22];
    t525 = -4.0*(t188*t181*t180+t200*t335);
    t528 = t188*t177;
    t531 = t248*t200+2.0*t528*t202;
    t535 = t188*t183;
    t538 = t399*t200+2.0*t535*t202;
    t542 = 4.*(-t336+t335)*t200;
    t545 = (t359+4.0*t182-t352)*t188;
    conn[0][0][26] = t525*t260+t531*t222+(t525*t28+t538*t16+t542+t545+t362)*t178+t427*t188-t203;
    t549 = -t364;
    conn[0][0][27] = (t549*t260+(t277+t308)*t222+(t375+t370-t374-t372)*t178)*t200
      +(-t380*t260+2.0*t421*t382+(t216+t385+t215-t387-t390+t388)*t178-t196*t181-t191*t16)*t188
      +t379*t178-t288*t16;
    conn[0][0][28] = conn[0][0][19];
    conn[0][0][29] = conn[0][0][23];
    conn[0][0][30] = conn[0][0][27];
    t575 = -4.0*(t188*t16*t183+t200*t393);
    t580 = t250*t200+2.0*t528*t28;
    t584 = t188*t180;
    t587 = t341*t200+2.0*t584*t28;
    t591 = 4.*(t393-t394)*t200;
    t594 = (4.0*t185+t417-t410)*t188;
    conn[0][0][31] = t575*t260+t580*t222+(t575*t202+t587*t181+t591+t594+t420)*t178+t428*t188-t207;
    t598 = t180*t202;
    conn[0][0][32] = 2.*(t224+t598+t227)*t188-t180;
    conn[0][0][33] = (t213*t202+t486)*t200+(2.*(t598+t227)*t178-t298
				      +2.*(t177-t207)*t181+t190-t189)*t188+t288-t204;
    t613 = t192*t178;
    t617 = -t354;
    conn[0][0][34] = (t613+t439)*t181*t200+(t617*t202-2.0*t356*t181)*t188;
    t624 = t277*t178;
    t631 = 2.0*t183-2.0*t237+2.0*t238;
    conn[0][0][35] = (t624+t218*t202)*t200+(-t387+2.0*t598*t16+t631*t181+t216-t215)*t188+t226-t378;
    conn[0][0][36] = conn[0][0][33];
    t639 = -4.0*(t188*t178*t177+t200*t244);
    t644 = t339*t200+2.0*t584*t222;
    t650 = t397*t200+2.0*t535*t222;
    t654 = 4.*(-t245+t244)*t200;
    t657 = (4.0*t179+t269-t262)*t188;
    conn[0][0][37] = t639*t296+t644*t202+(t639*t28+t650*t16+t654+t657+t272)*t181+t613*t188-t223;
    t668 = t202*t202;
    conn[0][0][38] = (-4.0*t190*t296+t428*t202+(-t281+t440+t488+t283+t282)*t181)*t200
      +(-2.0*t177*t668+2.*(t209+t177-t207)*t202+2.*(t204+t190-t300-t225*t178)*t181)*t188+t203;
    conn[0][0][39] = (t505*t296+(t276+t367)*t202+(t312-t314-t316+t311)*t181)*t200
      +(-t321*t296+2.0*t598*t323+(-t330-t332+t194+t326+t327+t195)*t181-t217*t178+t191*
	t16)*t188+t320*t181-t204*t16;
    conn[0][0][40] = conn[0][0][34];
    conn[0][0][41] = conn[0][0][38];
    conn[0][0][42] = (-4.0*t335*t296+(t247*t178+t398*t16)*t202
		+4.*(t335-t344-t346)*t181)*t200+(-2.0*t180*t668+t617*t296-t357*t202+2.0*t181*t335)*t188+t598;
    t722 = t197*t178;
    t724 = t549*t222;
    t725 = t308*t178;
    conn[0][0][43] = (-4.0*t216*t296+t722*t202+(t724+t725+t370+t375-t372)*t181)*t200+
      (-2.0*t183*t668+t631*t202+2.*(t216+t378-t223*t16-t389)*t181)*t188+t239;
    conn[0][0][44] = conn[0][0][35];
    conn[0][0][45] = conn[0][0][39];
    conn[0][0][46] = conn[0][0][43];
    conn[0][0][47] = t575*t296+t587*t202+(t575*t222+t580*t178+t591+t594+t420)*t181+t439*t188-t225;
    t746 = t183*t28;
    conn[0][0][48] = 2.*(t238+t240+t746)*t188-t183;
    conn[0][0][49] = (t507+t231*t28)*t200+(2.*(t240+t746)*t178-t330+2.*(t177-t203)*t16+t194-t195)*t188+t319-t208;
    conn[0][0][50] = (t725+t233*t28)*t200+(2.*(t238+t746)*t181-t390+2.*(t180-t223)*t16+t215-t216)*t188+t378-t226;
    t775 = -t412;
    conn[0][0][51] = (t722+t451)*t16*t200+(t775*t28-2.0*t414*t16)*t188;
    conn[0][0][52] = conn[0][0][49];
    conn[0][0][53] = t639*t328+t650*t28+(t639*t202+t644*t181+t654+t657+t272)*t16+t722*t188-t237;
    conn[0][0][54] = (t487*t328+(t307+t366)*t28+(t283-t285+t282-t281)*t16)*t200
      +(-t290*t328+2.0*t746*t292+(t299+t295+t190-t301+t189-t298)*t16+t217*t178+t196*
	t181)*t188+t289*t16-t208*t181;
    t96 = t28*t28;
    conn[0][0][55] = (-4.0*t194*t328+t213*t28*t181+(-t314+t452+t506+t311+t312)*t16)
      *t200+(-2.0*t177*t96+2.*(-t203+t177+t205)*t28+2.*(t208-t239*t178+t194-t331)*t16)*t188+t207;
    conn[0][0][56] = conn[0][0][50];
    conn[0][0][57] = conn[0][0][54];
    conn[0][0][58] = t525*t328+t538*t28+(t525*t222+t531*t178+t542+t545+t362)*t16+t451*t188-t239;
    conn[0][0][59] = (-4.0*t215*t328+t192*t28*t178+(t724+t624+t370+t375-t374)*t16)*
      t200+(-2.0*t180*t96+2.*(t224-t223+t180)*t28+2.*(t226-t386+t215-t237*t181)*t16)*t188+t225;
    conn[0][0][60] = conn[0][0][51];
    conn[0][0][61] = conn[0][0][55];
    conn[0][0][62] = conn[0][0][59];
    conn[0][0][63] = (-4.0*t393*t328+(t249*t178+t340*t181)*t28+4.*(-t405+t393-t403)*t16)*t200
      +(-2.0*t183*t96+t775*t328-t415*t28+2.0*t16*t393)*t188+t746;

    /*********************************************************************************************
       Up until this poin the connection is in non-primed coordinates.  
       Now we transform it to  primed/numerical coordinates. 
    *********************************************************************************************/

    transform_connection2(coords,conn,connp);

    return;
  }
}

#endif  /*  #if( USE_STRICT_ARRAY_BOUNDS )  */


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     G E N E R A L  (A R B I T R A R Y)   M E T R I C    R O U T I N E S

 *****************************************************************************/

/******************************************************************************
  general_gcon_func():
  ----------
       -- place holder for new gcon routine
       -- used on in special situations, see gcov_func() and set_general_geometry for actual usage
 ******************************************************************************/
void general_gcon_func( double *x, double gcon[][NDIM] ) 
{
  double gcov[NDIM][NDIM];
  double detg;

  /* Some dynamical metrics may also calculate inverse metric too, but so far they do not: */
  general_gcov_func(x,gcov);
  gcon_func(gcov,gcon,&detg);

  return;
}

/******************************************************************************
  set_general_geometry():
  ----------
    -- responsible for retrieving/calculating gcov, gcon and conn at all 
       necessary places for the given time;
    -- used by both METRIC_GENERAL_STATIC and METRIC_GENERAL_DYNAMIC
    -- also responsible for setting ncon,beta,alpha
    -- saves the geometry into the  "n"-th spot in the geometry arrays;
******************************************************************************/
void set_general_geometry( int n, double t_input ) 
{

#define GCON_EPS (-1.e-3)
  
  int i,j,k,pos;
  double *x;
  struct of_geom *geom; 
  struct of_coord *coords;

  TRACE_BEG;

  ncurr = n; 

  SET_GEOM_COORD_TIME_FUNCS(t_input,1) ;

#if(DYNAMIC_COORDINATES)
  extern void  calc_all_coord(int n, double t_now ) ;
  calc_all_coord(n,t_input);
#endif

# if( BBH_SPACETIME && USE_MASK ) 
  extern void set_excision_mask( short unsigned int ***evol_mask_loc );
  set_excision_mask(evol_mask[ncurr]); 
  double max_gdet = -1.e200;
# endif

  GEOM_LOOP { 

    //    fprintf(stdout,"set_general_geometry-all %6d   %6d  %6d  %6d %6d %6d \n",myid,n,pos,i,j,k); fflush(stdout);
    
    get_geometry(i,j,k,pos,n,geom);
    get_coord(   i,j,k,pos,n,coords);

#if(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_35PN_NZ || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_FAST || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NR_FZ || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_DROP || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_SIMPLE_NEWTONIAN || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP)

    gcov_func(coords, geom->gcov);
    gcon_func(geom->gcov,geom->gcon,&(geom->g));

# if( USE_MASK ) 
    //    if( evol_mask[ncurr][i][j][k] == MASK_EXCISED ) { geom->g = 1.; }
    if( evol_mask[ncurr][i][j][k] == MASK_BUFFER ) { 
      if( geom->g > max_gdet ) { max_gdet = geom->g ; } 
    }
# endif

 
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_KS_SPHERICAL  )
    x  = coords->x;
    ks_gcov_func( x, geom->gcov );
    ks_gcon_func( x, geom->gcon );
    ks_gdet_func( x, &(geom->g));
    transform_all(coords,geom);

#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_MINK_CARTESIAN  )
    x  = coords->x;
    mink_cartesian_gcov_func(x,geom->gcov); 
    mink_cartesian_gcon_func(x,geom->gcon); 
    geom->g = 1.;
    transform_all(coords,geom);

#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_MINK_SPHERICAL  )
    x  = coords->x;
    mink_spherical_gcov_func( x, geom->gcov );
    mink_spherical_gcon_func( x, geom->gcon );
    mink_spherical_gdet_func( x, &(geom->g));
    transform_all(coords,geom);

#else
    fprintf(stderr,"general_gcon_func(): option not implemented yet:  %d \n",METRIC_DYNAMIC_TYPE_CHOICE); 
    fflush(stderr);
    fail( FAIL_METRIC,0 );
#endif


    /* Set auxiliary, derived metric functions : */ 
    geom->g_inv  = 1./geom->g;
    if( geom->gcon[0][0] >= GCON_EPS )  { 
      geom->gcon[0][0] = GCON_EPS;
    }
    geom->ncon[0] = sqrt(-geom->gcon[0][0]) ; // 1/alpha
    geom->alpha   =  1./geom->ncon[0];
    geom->ncon[1] = -geom->alpha * geom->gcon[0][1];
    geom->ncon[2] = -geom->alpha * geom->gcon[0][2];
    geom->ncon[3] = -geom->alpha * geom->gcon[0][3];
    geom->beta[0] =  -geom->gcon[0][1]/geom->gcon[0][0];
    geom->beta[1] =  -geom->gcon[0][2]/geom->gcon[0][0];
    geom->beta[2] =  -geom->gcon[0][3]/geom->gcon[0][0];

  }

# if( USE_MASK ) 
  mpi_global_max(&max_gdet); 

  GEOM_LOOP { 
    if( evol_mask[ncurr][i][j][k] == MASK_EXCISED ) { 
      geom->g = max_gdet; 
      geom->g_inv  = 1./geom->g;
    }
  }
# endif

#undef GCON_EPS 

  TRACE_END;
  return;
}

/******************************************************************************
  set_general_conn():
  ----------
    -- calculates the connection or Christoffel symbols for the "general" metric;
    -- used by both METRIC_GENERAL_STATIC and METRIC_GENERAL_DYNAMIC
 ******************************************************************************/
void set_general_conn( double t_input )
{
  int i,j,k,l,jbeg,jend;
  struct of_geom *geom;
  struct of_coord *coords;

  TRACE_BEG;

#if( USE_KINEMATIC_VISCOSITY ) 
	fprintf(stdout,"set_general_conn():  USE_KINEMATIC_VISCOSITY may not be compatible with finite differencing of the metric due to the fact that x2 is reflected when setting the theta coordinates \n"); 
	fflush(stdout);
	fail(FAIL_BASIC,0); 
#endif

#if( NDC != 0 )
	fprintf(stdout,"set_general_conn():  NDC must be zero \n"); 
	fflush(stdout);
	fail(FAIL_BASIC,0); 
#endif

 DLOOP1 { w1_v[i] =    invdx[i] / 6. ; }
 DLOOP1 { w2_v[i] = 4.*invdx[i] / 3. ; }


 k = 2; 
 i = 0; 
 w_2lo[i++] = -0.5 * invdx[k]      ;     /* i-1/2 */
 w_2lo[i++] =   -5.* invdx[k] / 3. ;     /* i     */ 
 w_2lo[i++] =    3.* invdx[k]      ;     /* i+1/2 */
 w_2lo[i++] =       -invdx[k]      ;     /* i+1   */
 w_2lo[i++] =        invdx[k] / 6. ;     /* i+3/2 */

 i = 0; 
 w_2hi[i++] =       -invdx[k] / 6. ;    /* i-3/2 */
 w_2hi[i++] =        invdx[k]      ;    /* i-1   */
 w_2hi[i++] =   -3.* invdx[k]      ;    /* i-1/2 */
 w_2hi[i++] =    5.* invdx[k] / 3. ;    /* i     */ 
 w_2hi[i++] =  0.5 * invdx[k]      ;    /* i+1/2 */


  if( dt_conn > DEL_X_CONN*dx[0] ) { 
    fprintf(stderr,"set_general_conn():  dt_conn is too big!      %26.16e %26.16e \n",dt_conn,dx[0]); fflush(stderr);
  }

#if(DYNAMIC_SPACETIME)
  SET_GEOM_COORD_TIME_FUNCS((t_input-dt_conn),0);
  SET_GEOM_COORD_TIME_FUNCS((t_input        ),1);
  SET_GEOM_COORD_TIME_FUNCS((t_input+dt_conn),2);
#endif //dynamic spacetime


#if( (!DYNAMIC_COORDINATES) && (METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_KS_SPHERICAL) )
  GDUMP_LOOP  { 
        get_coord(i,j,k,CENT,ncurr,coords);   ks_conn_func( coords, conn[CONN_ID(i,j,k)] ); 
  }
    

#else 

# if( (CONN_METHOD_CHOICE == CONN_METHOD_4TH_ORDER) )

  if( N2S == N2E ) { 
    GDUMP_LOOP  { 
      get_coord(i,j,k,CENT,ncurr,coords);   get_geometry(i,j,k,CENT,ncurr,geom);  conn_func_fd4(   coords,geom, conn[CONN_ID(i,j,k)]);
    }
  } 
  else { 
    if( bc_pid[2][BCDN] == BC_PHYS ) { 
      if( bc_pid[2][BCUP] == BC_PHYS ) { 
	/* lo & hi  */
	j = N2S; 
	N1_LOOP  N3_LOOP { 
	  get_coord(i,j,k,CENT,ncurr,coords);   get_geometry(i,j,k,CENT,ncurr,geom);  conn_func_fd4_lo2(coords,geom, conn[CONN_ID(i,j,k)]);
	}
	N1_LOOP  for(j=N2S+1;j<N2E ;j++) N3_LOOP { 
	    get_coord(i,j,k,CENT,ncurr,coords); get_geometry(i,j,k,CENT,ncurr,geom);  conn_func_fd4(    coords,geom, conn[CONN_ID(i,j,k)]);
	  }
	j = N2E; 
	N1_LOOP N3_LOOP { 
	  get_coord(i,j,k,CENT,ncurr,coords);   get_geometry(i,j,k,CENT,ncurr,geom);  conn_func_fd4_hi2(coords,geom, conn[CONN_ID(i,j,k)]);
	}
      }
      else { 
	/* lo */
	j = N2S; 
	N1_LOOP  N3_LOOP { 
	  get_coord(i,j,k,CENT,ncurr,coords);   get_geometry(i,j,k,CENT,ncurr,geom); conn_func_fd4_lo2(coords,geom, conn[CONN_ID(i,j,k)]);
	}
	N1_LOOP  for(j=N2S+1;j<=N2E ;j++) N3_LOOP { 
	    get_coord(i,j,k,CENT,ncurr,coords); get_geometry(i,j,k,CENT,ncurr,geom); conn_func_fd4(    coords,geom, conn[CONN_ID(i,j,k)]);
	  }
      }
    }
    else { 
      if( bc_pid[2][BCUP] == BC_PHYS ) { 
	/* hi  */
	N1_LOOP for(j=N2S;j<N2E ;j++) N3_LOOP { 
	    get_coord(i,j,k,CENT,ncurr,coords); get_geometry(i,j,k,CENT,ncurr,geom); conn_func_fd4(    coords,geom, conn[CONN_ID(i,j,k)]);
	  }
	j = N2E; 
	N1_LOOP  N3_LOOP { 
	  get_coord(i,j,k,CENT,ncurr,coords);   get_geometry(i,j,k,CENT,ncurr,geom); conn_func_fd4_hi2(coords,geom, conn[CONN_ID(i,j,k)]);
	}
      }
      else { 
	/* neither */
	GDUMP_LOOP  { 
	  get_coord(i,j,k,CENT,ncurr,coords);   get_geometry(i,j,k,CENT,ncurr,geom); conn_func_fd4(    coords,geom, conn[CONN_ID(i,j,k)]);
	}
      }
    }
  } /* N2S = N2E */

# else
  l = 0;
  CONN2_LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    get_geometry(i,j,k,CENT,ncurr,geom);

#  if( CONN_METHOD_CHOICE == CONN_METHOD_ORIG )
    conn_func(coords,geom, conn[l]);

#  elif( CONN_METHOD_CHOICE == CONN_METHOD_2ND_ORDER )
    conn_func_fd2(coords,geom, conn[l]);

#  else
    conn_func(coords,geom, conn[l]);

#  endif

    l++;
  }

# endif
#endif


#if(DYNAMIC_SPACETIME)
  SET_GEOM_COORD_TIME_POINTERS(1);
#endif

  TRACE_END;
  return;
}

/******************************************************************************
  init_general_metric():
  ----------
    -- initializes the arrays needed for a general spacetime
    -- used by both METRIC_GENERAL_STATIC and METRIC_GENERAL_DYNAMIC
    -- keep in mind that this routine is used to restart from checkpoint files; 
 ******************************************************************************/
void init_general_metric( void ) 
{
  double t_beg, t_mid, t_end; 

  TRACE_BEG;

#if(DYNAMIC_SPACETIME)
# if( N0_GEOM < 3 ) 
    fprintf(stderr,"init_general_metric():  We need more than 2 time storage spots !   %d \n",N0_GEOM); fflush(stderr);
    fail( FAIL_METRIC, 0);
# endif
    
  /* We need a guess at least for the initial time step : */
  if( fabs(dx[0]) < SMALL ) { 
    fprintf(stderr,"init_general_metric():  dx[0] needs to be set first! \n"); fflush(stderr);
    fail( FAIL_METRIC, 0);
  }

  n_beg = 0;
  n_mid = n_beg+1;
  n_end = n_mid+1;
  ncurr = n_beg;

  t_beg = t;
  t_mid = t_beg + 0.5*dx[0];
  t_end = t_beg + dx[0];

  set_general_geometry(n_beg,t_beg);
  set_general_geometry(n_mid,t_mid);
  set_general_geometry(n_end,t_end);

  /* We only store one timeslice's worth of connection : */
  /* Even though we set the first substep's conn  in step_ch(), we set it here so that it 
     is dumped in the first gdump file. */
  ncurr = n_beg;
  set_general_conn(t_beg);

#else   /*  IT IS NOT DYNAMIC_SPACETIME */
  ncurr = n_beg = n_mid = n_end = 0;
  t_beg = t_mid = t_end = startx[0];
  set_general_geometry(n_beg,t_beg);
  set_general_conn(t_beg);

#endif

  TRACE_END;
  return;
}


/********************************************************************************/
/********************************************************************************
 get_bbh_traj_data():
 -------------
  -- acts as a translator for the various definitions of "nz_params";
     
********************************************************************************/
void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) 
{

#if( METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_35PN_NZ  ||\
      METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ   )
  bbh_traj->tt      = bh1_traj[ncurr][TT] ;
  bbh_traj->xi1x    = bh1_traj[ncurr][XX] ;
  bbh_traj->xi1y    = bh1_traj[ncurr][YY] ;
  bbh_traj->xi1z    = bh1_traj[ncurr][ZZ] ;
  bbh_traj->xi2x    = bh2_traj[ncurr][XX] ;
  bbh_traj->xi2y    = bh2_traj[ncurr][YY] ;
  bbh_traj->xi2z    = bh2_traj[ncurr][ZZ] ;

#elif( METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_FAST||\
       METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ     )
  bbh_traj->xi1x    = nz_params[ 0] ;
  bbh_traj->xi1y    = nz_params[ 1] ;
  bbh_traj->xi2x    = nz_params[ 2] ;
  bbh_traj->xi2y    = nz_params[ 3] ;
  bbh_traj->v1x     = nz_params[ 4] ;
  bbh_traj->v1y     = nz_params[ 5] ;
  bbh_traj->v2x     = nz_params[ 6] ;
  bbh_traj->v2y     = nz_params[ 7] ;
  bbh_traj->v1      = nz_params[ 8] ;
  bbh_traj->v2      = nz_params[ 9] ;
  bbh_traj->v12x    = nz_params[10] ;
  bbh_traj->v12y    = nz_params[11] ;
  bbh_traj->v21x    = nz_params[12] ;
  bbh_traj->v21y    = nz_params[13] ;
  bbh_traj->v12     = nz_params[14] ;
  bbh_traj->v21     = nz_params[15] ;
  bbh_traj->v1v2    = nz_params[16] ;
  bbh_traj->v2v1    = nz_params[17] ;
  bbh_traj->t_c     = nz_params[18] ;
  bbh_traj->phi     = nz_params[19] ;
  bbh_traj->omega   = nz_params[20] ;
  bbh_traj->r12     = nz_params[21] ;
  bbh_traj->r21     = nz_params[22] ;
  bbh_traj->xi1z    = nz_params[23] ;
  bbh_traj->xi2z    = nz_params[24] ;
  bbh_traj->v1z     = nz_params[25] ;
  bbh_traj->v2z     = nz_params[26] ;
  bbh_traj->v12z    = nz_params[27] ;
  bbh_traj->v21z    = nz_params[28] ;
  bbh_traj->tt      = nz_params[36] ;

#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ ||\
       METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NR_FZ      )
  bbh_traj->xi1x    =   nz_params[ 0];
  bbh_traj->xi1y    =   nz_params[ 1];
  bbh_traj->xi2x    =   nz_params[ 2];
  bbh_traj->xi2y    =   nz_params[ 3];
  bbh_traj->v1x	    =   nz_params[ 4];
  bbh_traj->v1y	    =   nz_params[ 5];
  bbh_traj->v2x	    =   nz_params[ 6];
  bbh_traj->v2y	    =   nz_params[ 7];
  bbh_traj->v1	    =   nz_params[ 8];
  bbh_traj->v2	    =   nz_params[ 9];
  bbh_traj->v12x    =   nz_params[10];
  bbh_traj->v12y    =   nz_params[11];
  bbh_traj->v12	    =   nz_params[12];
  bbh_traj->v1v2    =   nz_params[13];
  bbh_traj->n12x    =   nz_params[14];
  bbh_traj->n12y    =   nz_params[15];
  bbh_traj->n12v12  =   nz_params[16]; 
  bbh_traj->n12v1   =   nz_params[17];
  bbh_traj->n12v2   =   nz_params[18];
  bbh_traj->t_c     =   nz_params[19];
  bbh_traj->phi     =   nz_params[20];
  bbh_traj->omega   =   nz_params[21];
  bbh_traj->r12     =   nz_params[22];
  bbh_traj->xi1z    =   nz_params[23];
  bbh_traj->xi2z    =   nz_params[24];
  bbh_traj->v1z	    =   nz_params[25];
  bbh_traj->v2z	    =   nz_params[26];
  bbh_traj->v12z    =   nz_params[27];
  bbh_traj->n12z    =   nz_params[28];
  bbh_traj->r1T0    =   nz_params[29];
  bbh_traj->w1T0    =   nz_params[30];
  bbh_traj->r2T0    =   nz_params[31];
  bbh_traj->w2T0    =   nz_params[32];
  bbh_traj->xNT0    =   nz_params[33];
  bbh_traj->wNT0    =   nz_params[34];
  bbh_traj->lambda  =   nz_params[35];
  bbh_traj->tt      =   nz_params[36];
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW ||\
       METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_DROP      ||\
       METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND    ||\
       METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND ||\
       METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_SIMPLE_NEWTONIAN )
  bbh_traj->xi1x    =   nz_params[ 0];
  bbh_traj->xi1y    =   nz_params[ 1];
  bbh_traj->xi2x    =   nz_params[ 2];
  bbh_traj->xi2y    =   nz_params[ 3];
  bbh_traj->v1x	    =   nz_params[ 4];
  bbh_traj->v1y	    =   nz_params[ 5];
  bbh_traj->v2x	    =   nz_params[ 6];
  bbh_traj->v2y	    =   nz_params[ 7];
  bbh_traj->v1	    =   nz_params[ 8];
  bbh_traj->v2	    =   nz_params[ 9];
  bbh_traj->v12x    =   nz_params[10];
  bbh_traj->v12y    =   nz_params[11];
  bbh_traj->v12	    =   nz_params[12];
  bbh_traj->v1v2    =   nz_params[13];
  bbh_traj->n12x    =   nz_params[14];
  bbh_traj->n12y    =   nz_params[15];
  bbh_traj->n12v12  =   nz_params[16]; 
  bbh_traj->n12v1   =   nz_params[17];
  bbh_traj->n12v2   =   nz_params[18];
  bbh_traj->t_c     =   nz_params[19];
  bbh_traj->phi     =   nz_params[20];
  bbh_traj->omega   =   nz_params[21];
  bbh_traj->r12     =   nz_params[22];
  bbh_traj->r12dot  =   nz_params[23];
  bbh_traj->xi1z    =   nz_params[24];
  bbh_traj->xi2z    =   nz_params[25];
  bbh_traj->v1z	    =   nz_params[26];
  bbh_traj->v2z	    =   nz_params[27];
  bbh_traj->v12z    =   nz_params[28];
  bbh_traj->n12z    =   nz_params[29];
  bbh_traj->r1T0    =   nz_params[30];
  bbh_traj->w1T0    =   nz_params[31];
  bbh_traj->r2T0    =   nz_params[32];
  bbh_traj->w2T0    =   nz_params[33];
  bbh_traj->xNT0    =   nz_params[34];
  bbh_traj->wNT0    =   nz_params[35];
  bbh_traj->lambda  =   nz_params[36];
  bbh_traj->tt      =   nz_params[37];

#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )
  bbh_traj->xi1x    =   bbh_params->xi1x;
  bbh_traj->xi1y    =   bbh_params->xi1y;
  bbh_traj->xi2x    =   bbh_params->xi2x;
  bbh_traj->xi2y    =   bbh_params->xi2y;
  bbh_traj->v1x     =   bbh_params->v1x;
  bbh_traj->v1y     =   bbh_params->v1y;
  bbh_traj->v2x     =   bbh_params->v2x;
  bbh_traj->v2y     =   bbh_params->v2y;
  bbh_traj->v1	    =   bbh_params->v1;
  bbh_traj->v2	    =   bbh_params->v2;
  bbh_traj->v12x    =   bbh_params->v12x;
  bbh_traj->v12y    =   bbh_params->v12y;
  bbh_traj->v12     =   bbh_params->v12;
  bbh_traj->v1v2    =   bbh_params->v1v2;
  bbh_traj->n12x    =   bbh_params->n12x;
  bbh_traj->n12y    =   bbh_params->n12y;
  bbh_traj->n12v12  =   bbh_params->n12v12; 
  bbh_traj->n12v1   =   bbh_params->n12v1;
  bbh_traj->n12v2   =   bbh_params->n12v2;
  bbh_traj->t_c     =   bbh_params->t_c;
  bbh_traj->phi     =   bbh_params->phi;
  bbh_traj->omega   =   bbh_params->omega;
  bbh_traj->r12     =   bbh_params->r12;
  bbh_traj->r12dot  =   bbh_params->r12dot;
  bbh_traj->xi1z    =   bbh_params->xi1z;
  bbh_traj->xi2z    =   bbh_params->xi2z;
  bbh_traj->v1z     =   bbh_params->v1z;
  bbh_traj->v2z     =   bbh_params->v2z;
  bbh_traj->v12z    =   bbh_params->v12z;
  bbh_traj->n12z    =   bbh_params->n12z;
  bbh_traj->r1T0    =   bbh_params->r1T0;
  bbh_traj->w1T0    =   bbh_params->w1T0;
  bbh_traj->r2T0    =   bbh_params->r2T0;
  bbh_traj->w2T0    =   bbh_params->w2T0;
  bbh_traj->xNT0    =   bbh_params->xNT0;
  bbh_traj->wNT0    =   bbh_params->wNT0;
  bbh_traj->lambda  =   bbh_params->lambda;
  bbh_traj->tt      =   bbh_params->t;
#else
  bbh_traj->xi1x   = 0. ;
  bbh_traj->xi1y   = 0. ;
  bbh_traj->xi2x   = 0. ;
  bbh_traj->xi2y   = 0. ;
  bbh_traj->v1x    = 0. ;
  bbh_traj->v1y    = 0. ;
  bbh_traj->v2x    = 0. ;
  bbh_traj->v2y    = 0. ;
  bbh_traj->v1     = 0. ;
  bbh_traj->v2     = 0. ;
  bbh_traj->v12x   = 0. ;
  bbh_traj->v12y   = 0. ;
  bbh_traj->v21x   = 0. ;
  bbh_traj->v21y   = 0. ;
  bbh_traj->v12    = 0. ;
  bbh_traj->v21    = 0. ;
  bbh_traj->v1v2   = 0. ; // dot product of v1 and v2
  bbh_traj->v2v1   = 0. ; // dot product of v2 and v1 (different from v1v2 due to PN approximation)
  bbh_traj->t_c    = 0. ; // time to merger from t=0
  bbh_traj->phi    = 0. ; // orbital phase change from t=0
  bbh_traj->omega  = 0. ; // current orbital phase rate of change
  bbh_traj->r12    = 0. ; // current separation
  bbh_traj->r21    = 0. ; // current separation (from BH2's perspective)
  bbh_traj->xi1z   = 0. ; //  
  bbh_traj->xi2z   = 0. ; //  
  bbh_traj->v1z    = 0. ; //  
  bbh_traj->v2z    = 0. ; //  
  bbh_traj->v12z   = 0. ; //  
  bbh_traj->v21z   = 0. ; //  
  bbh_traj->r1T0   = 0. ; // Inner1-Near radius trans. func. parameter 
  bbh_traj->w1T0   = 0. ; // Inner1-Near width trans. func. parameter 
  bbh_traj->r2T0   = 0. ; // Inner2-Near radius trans. func. parameter 
  bbh_traj->w2T0   = 0. ; // Inner2-Near width trans. func. parameter 
  bbh_traj->xNT0   = 0. ; // Near-Inner radius trans. func. parameter 
  bbh_traj->wNT0   = 0. ; // Near-Inner width trans. func. parameter 
  bbh_traj->lambda = 0. ; // Near-Far trans. func. parameter 
  bbh_traj->tt     = 0. ;
  bbh_traj->n12x   = 0. ; 
  bbh_traj->n12y   = 0. ; 
  bbh_traj->n12v12 = 0. ; 
  bbh_traj->n12v1  = 0. ; 
  bbh_traj->n12v2  = 0. ; 
#endif   /*  METRIC_DYNAMIC_TYPE_CHOICE  == ...  */


 return;
}


/********************************************************************************/
/********************************************************************************
 gcov_func():
 -------------
  -- Driver routine for calculating the covariant form of the metric  
     given numerical coordinates and outputted w.r.t. numerical coordinates; 

  -- this function is a top-level function so it is supposed to act with the 
     top-level numerical coordinates specifically and is different in this 
     respect to the other "gcov" routines in this file that use "normal" coordinates;
     
********************************************************************************/
void gcov_func(struct of_coord *coords, double gcov[][NDIM]) 
{
  int i,j;
  double *xg;
  double weight=0.;
  double denom=0.;

  /* Get "normal" coordinates */

  xg = coords->x;

#if(   METRIC_TYPE_CHOICE == METRIC_MINK_CARTESIAN  )
  DLOOP2 gcov[i][j] = mink[i][j] ;

#elif( METRIC_TYPE_CHOICE == METRIC_MINK_SPHERICAL  )
  mink_spherical_gcov_func( xg, gcov );

#elif( METRIC_TYPE_CHOICE == METRIC_KS_SPHERICAL    )
  ks_gcov_func( xg, gcov );

#elif( METRIC_TYPE_CHOICE == METRIC_KS_CARTESIAN    )
  ks_cart_gcov_func( xg, gcov );

#elif( METRIC_TYPE_CHOICE == METRIC_BL_SPHERICAL    )
  bl_gcov_func( xg, gcov );


#elif( METRIC_TYPE_CHOICE == METRIC_GENERAL_STATIC || METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC  )

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
    xg = coords->xcart ; 
# endif

  general_gcov_func(xg,gcov); 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
  transform_rank2cov2(coords->dxc_dxp, gcov);
# endif    


#elif( METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG )
# if( TOP_TYPE_CHOICE == TOP_CARTESIAN )   
    need-to-fix-for-cartesian-case-x-assumed-to-be-spherical-coordinates
# endif

    int ip;
    double gcov_avg[NDIM][NDIM];
    struct of_coord coords_loc;

    copy_coord(coords,&coords_loc);

#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov_avg[i][j] = 0.; } 
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] = 0.;
#endif

    for(ip=0; ip < n_phi_avg; ip++ ) { 
      coords_loc.x[PH] = ip * d_phi_avg; 
      coords_loc.xp[PH] = coords_loc.x[PH];

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      dx_dxp_calc( coords_loc.x, coords_loc.xp,  coords_loc.dx_dxp );
      xcart_of_xspher_special( &coords_loc );
      xg = coords_loc.xcart ; 
#else 
      xg = coords_loc.x;
# endif

      general_gcov_func(xg,gcov); 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      transform_rank2cov2(coords_loc.dxc_dxp, gcov);
# endif

      weight = sqrt(gcov[PH][PH]);
      denom += weight; 

#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov_avg[i][j] += gcov[i][j]*weight ; } 
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] += gcov[0][i]*weight ;
#endif

    }
    denom = 1./denom;
#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov[i][j] = gcov_avg[i][j] * denom ; } 
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov[0][i] = gcov_avg[0][i] * denom;
#endif

#elif( METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG2 )
    /* Here, calculate the metric in r,theta,phi coordinates and not in xp or numerical coordinates */

# if( TOP_TYPE_CHOICE == TOP_CARTESIAN )   
    need-to-fix-for-cartesian-case-x-assumed-to-be-spherical-coordinates-2
# endif

    int ip;
    double gcov_avg[NDIM][NDIM];
    struct of_coord coords_loc;

    copy_coord(coords,&coords_loc);

#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov_avg[i][j] 0.; }
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] = 0.;
#endif

    for(ip=0; ip < n_phi_avg; ip++ ) { 
      coords_loc.x[PH] = ip * d_phi_avg; 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      xcart_of_xspher_special2( &coords_loc );
      xg = coords_loc.xcart ; 
#else 
      xg = coords_loc.x;
# endif

      general_gcov_func(xg,gcov); 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      /* note that here, "dxc_dxp" is really  "dxc_dx" : */ 
      transform_rank2cov2(coords_loc.dxc_dxp, gcov);
# endif

      weight = sqrt(gcov[PH][PH]);
      denom += weight; 

#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov_avg[i][j] += gcov[i][j]*weight ; } 
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] += gcov[0][i]*weight ;
#endif

    }
    denom = 1./denom;
#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov[i][j] = gcov_avg[i][j] * denom ; } 
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov[0][i] = gcov_avg[0][i] * denom;
#endif

#else
  fprintf(stderr,"gcov_func():  Invalid value for METRIC_TYPE_CHOICE : %d \n", 
	  METRIC_TYPE_CHOICE);
  fflush(stderr);
  fail( FAIL_METRIC,0 );
#endif 

  /* Transform to numerical coordinates : */
#if( METRIC_TOP_TYPE == TOP_TYPE_CHOICE )   
  transform_rank2cov2(coords->dx_dxp, gcov);
#endif
  
  return;
}

/********************************************************************************/
/********************************************************************************
 calc_all_geom
 -------------
  -- Driver routine for calculating the metric, inverse metric and connection 
     coefficients 

********************************************************************************/
void calc_all_geom(void) 
{
  int ii,jj, kk, i, j, k, l,pos ;
  struct of_geom *geom;
  struct of_coord *coords;

  /* Set certian coordinate constants : */ 
  ncurr = n_beg = n_mid = n_end = 0;
  asq = a*a;

  DLOOP1 { invdx[i] = 1./dx[i] ; } 

  check_coords(); 

#if(!( DYNAMIC_COORDINATES ) )
  calc_all_coord(ncurr,t);   
#endif

#if(METRIC_TYPE_CHOICE == METRIC_MINK_CARTESIAN)

  GEOM_LOOP {
    get_geometry(i,j,k,pos,ncurr,geom);
    get_coord(   i,j,k,pos,ncurr,coords);
    for(ii=0; ii<NDIM; ii++)  for(jj=0; jj<NDIM; jj++) { 
      geom->gcov[ii][jj] = geom->gcon[ii][jj] = mink[ii][jj] ;
    }
    geom->g = 1.;
    transform_all(coords, geom);
  }
  CONN_LOOP  for(ii=0; ii<NDIM; ii++) for(jj=0; jj<NDIM; jj++) for(kk=0; kk<NDIM; kk++) {  conn[l][ii][jj][kk] = 0.; }


#elif( METRIC_TYPE_CHOICE == METRIC_MINK_SPHERICAL )

  GEOM_LOOP {
    get_geometry(i,j,k,pos,ncurr,geom);
    get_coord(   i,j,k,pos,ncurr,coords);
    mink_spherical_gcov_func( coords->x, geom->gcov );
    mink_spherical_gcon_func( coords->x, geom->gcon );
    mink_spherical_gdet_func( coords->x, &(geom->g));
    transform_all(coords, geom);
  }
  l = 0;
  CONN2_LOOP { 
    get_coord(   i,j,k,CENT,ncurr,coords);
    mink_spherical_conn_func( coords, conn[l] );
    l++;
  }


#elif( METRIC_TYPE_CHOICE == METRIC_KS_SPHERICAL   )
  GEOM_LOOP {
    get_geometry(i,j,k,pos,ncurr,geom);
    get_coord(   i,j,k,pos,ncurr,coords);
    ks_gcov_func( coords->x, geom->gcov );
    ks_gcon_func( coords->x, geom->gcon );
    ks_gdet_func( coords->x, &(geom->g));
    transform_all(coords, geom);
  }
  l = 0;
  CONN2_LOOP { 
    get_coord(   i,j,k,CENT,ncurr,coords);
    ks_conn_func( coords, conn[l] );
    l++;
  }


#elif( METRIC_TYPE_CHOICE == METRIC_KS_CARTESIAN   )
  GEOM_LOOP {
    get_geometry(i,j,k,pos,ncurr,geom);
    get_coord(   i,j,k,pos,ncurr,coords);
    ks_cart_gcov_func( coords->x, geom->gcov );
    ks_cart_gcon_func( coords->x, geom->gcon );
    ks_cart_gdet_func( coords->x, &(geom->g));
    transform_all(coords, geom);
  }
  l = 0;
  CONN2_LOOP { 
    get_coord(   i,j,k,CENT,ncurr,coords);
    ks_cart_conn_func( coords, conn[l] );
    l++;
  }


#elif( METRIC_TYPE_CHOICE == METRIC_BL_SPHERICAL   )
  GEOM_LOOP {
    get_geometry(i,j,k,pos,ncurr,geom);
    get_coord(   i,j,k,pos,ncurr,coords);
    bl_gcov_func( coords->x, geom->gcov );
    bl_gcon_func( coords->x, geom->gcon );
    bl_gdet_func( coords->x, &(geom->g) );
    transform_all(coords, geom);
  }
  l = 0;
  CONN2_LOOP { 
    get_coord(   i,j,k,CENT,ncurr,coords);
    bl_conn_func( coords,conn[l]);
    l++;
  }

#elif( METRIC_TYPE_CHOICE == METRIC_GENERAL_STATIC || METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC || METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG || METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG2 )
  init_general_metric();
#else
  fprintf(stderr,"calc_all_geom():  Invalid value for METRIC_TYPE_CHOICE : %d \n",
	  METRIC_TYPE_CHOICE);
  fflush(stderr);
  fail( FAIL_METRIC,0 );
#endif 

  /* Calculate the lapse and shift from the just calculated metric :  */
#if( METRIC_TYPE_CHOICE != METRIC_GENERAL_STATIC && METRIC_TYPE_CHOICE != METRIC_GENERAL_DYNAMIC && METRIC_TYPE_CHOICE != METRIC_GENERAL_PHI_AVG && METRIC_TYPE_CHOICE != METRIC_GENERAL_PHI_AVG2 )
  GEOM_LOOP {
    get_geometry(i,j,k,pos,ncurr,geom);
    geom->g_inv  = 1./geom->g;
    geom->ncon[0] = sqrt(-geom->gcon[0][0]) ; // 1/alpha
    geom->alpha   =  1./geom->ncon[0];
    geom->ncon[1] = -geom->alpha * geom->gcon[0][1];
    geom->ncon[2] = -geom->alpha * geom->gcon[0][2];
    geom->ncon[3] = -geom->alpha * geom->gcon[0][3];
    geom->beta[0] =  -geom->gcon[0][1]/geom->gcon[0][0];
    geom->beta[1] =  -geom->gcon[0][2]/geom->gcon[0][0];
    geom->beta[2] =  -geom->gcon[0][3]/geom->gcon[0][0];
  }
#endif


  /**************************************************************************************************
    Set global arrays used for 
  **************************************************************************************************/

  /* Set arrays used for transforming to and from x and xp, and regularizing vectors : */
#if( RESCALE_REGULARIZE )
  unsigned int ndim_m1 = (NDIM-1);
  unsigned int npos_m1 = (NPOS-1);

  ALL_LOOP for(pos=0; pos<npos_m1; pos++) for(l=0; l<ndim_m1; l++) { regularize_gf[i][j][k][pos][l] = 1.; }

  ALL_LOOP { 
    get_coord(i,j,k,CENT,ncurr,coords);
    regularize_rank1con( coords->x, regularize_gf[i][j][k][CENT] ); 

    for(pos=FACE1;pos<=FACE3;pos++) { 
      get_coord(i,j,k,pos,ncurr,coords);
      regularize_rank1con( coords->x, regularize_gf[i][j][k][pos] ); 
    }
  }

  ALL_LOOP for(pos=0; pos<npos_m1; pos++) for(l=0; l<ndim_m1; l++) { unregularize_gf[i][j][k][pos][l] = 1./regularize_gf[i][j][k][pos][l]; }
#endif


  dx_global_min = find_min_dx();

  return;
}


/********************************************************************************/
/********************************************************************************
 get_special_geometry():
 ---------------------
  -- Load a different spacetime for temporary use: only sets gcov and gcon for now 
  !! Note that this routine does not set all the elements of the of_geom structure!
       -- I think so far this routine is used just for the metric and its inverse...

********************************************************************************/
void get_special_geometry( struct of_coord *coords, struct of_geom *geom, int geom_type)
{
  int i,j;
  int k,l,ip;
  double *xg;
  double gcov_avg[NDIM][NDIM], tmpout;
  double weight=0.;
  double denom=0.;
  struct of_coord coords_loc;

  xg = coords->x; 

  switch( geom_type ) { 

  case METRIC_MINK_CARTESIAN :
    DLOOP2 geom->gcon[i][j] = geom->gcov[i][j] = mink[i][j] ;  
    geom->g = 1.;
    break;

  case METRIC_MINK_SPHERICAL :
    mink_spherical_gcov_func( xg, geom->gcov ); 
    mink_spherical_gcon_func( xg, geom->gcon ); 
    mink_spherical_gdet_func( xg, &(geom->g) );
    break;

  case METRIC_KS_SPHERICAL   : 
    ks_gcov_func( xg, geom->gcov );  
    ks_gcon_func( xg, geom->gcon );  
    ks_gdet_func( xg, &(geom->g));
    break;

  case METRIC_KS_CARTESIAN   :
    ks_cart_gcov_func( xg, geom->gcov ); 
    ks_cart_gcon_func( xg, geom->gcon ); 
    ks_cart_gdet_func( xg, &(geom->g));
    break;

  case METRIC_BL_SPHERICAL   :
    bl_gcov_func( xg, geom->gcov ); 
    bl_gcon_func( xg, geom->gcon ); 
    bl_gdet_func( xg, &(geom->g) );
    break;

  case METRIC_GENERAL_STATIC        :
  case METRIC_GENERAL_DYNAMIC       :

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
    xg = coords->xcart ; 
# endif

    general_gcov_func(xg,geom->gcov); 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
    transform_rank2cov2(coords->dxc_dxp,geom->gcov);
# endif    

    gcon_func(geom->gcov,geom->gcon,&(geom->g));
    break;


  case METRIC_GENERAL_PHI_AVG :
# if( TOP_TYPE_CHOICE == TOP_CARTESIAN )   
    fprintf(stdout,"get_special_geometry():  Need to fix phi-avg for cartesian coordinates, on line %d  of %s  \n",__LINE__,__FILE__); fflush(stdout); fail(FAIL_BASIC,0);
# endif

    copy_coord(coords,&coords_loc);

#if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 {  gcov_avg[i][j] = 0.; }
#else
    for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] = 0.;
#endif

    for(ip=0; ip < n_phi_avg; ip++ ) { 
      coords_loc.x[PH] = ip * d_phi_avg; 
      coords_loc.xp[PH] = coords_loc.x[PH];

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      dx_dxp_calc( coords_loc.x, coords_loc.xp,  coords_loc.dx_dxp );
      xcart_of_xspher_special( &coords_loc );
      xg = coords_loc.xcart ; 
# endif

      general_gcov_func(xg,geom->gcov); 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      transform_rank2cov2(coords_loc.dxc_dxp, geom->gcov);
# endif 

      weight = sqrt(geom->gcov[PH][PH]);
      denom += weight; 

# if( USE_STRICT_ARRAY_BOUNDS ) 
      DLOOP2 { gcov_avg[i][j] += geom->gcov[i][j]*weight ; }
# else 
      for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] += geom->gcov[0][i]*weight ;
# endif

    }
    
    denom = 1./denom;

# if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 { geom->gcov[i][j] = gcov_avg[i][j] * denom; }
# else
    for(i = 0; i < NDIM*NDIM; i++) geom->gcov[0][i] = gcov_avg[0][i] * denom;
# endif

    gcon_func(geom->gcov,geom->gcon,&(geom->g));

    break;


  case METRIC_GENERAL_PHI_AVG2 :
# if( TOP_TYPE_CHOICE == TOP_CARTESIAN )
    fprintf(stdout,"get_special_geometry():  Need to fix phi-avg2 for cartesian coordinates, on line %d  of %s  \n",__LINE__,__FILE__); fflush(stdout); fail(FAIL_BASIC,0);
# endif

    copy_coord(coords,&coords_loc);

# if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 { gcov_avg[i][j] = 0.; }
# else 
    for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] = 0.;
# endif

    for(ip=0; ip < n_phi_avg; ip++ ) { 
      coords_loc.x[PH] = ip * d_phi_avg; 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      xcart_of_xspher_special2( &coords_loc );
      xg = coords_loc.xcart ; 
# endif

      general_gcov_func(xg,geom->gcov); 

# if( METRIC_TOP_TYPE != TOP_TYPE_CHOICE )   
      transform_rank2cov2(coords_loc.dxc_dxp, geom->gcov);
# endif 

      weight = sqrt(geom->gcov[PH][PH]);
      denom += weight; 

# if( USE_STRICT_ARRAY_BOUNDS ) 
      DLOOP2 { gcov_avg[i][j] += geom->gcov[i][j]*weight ; } 
# else
      for(i = 0; i < NDIM*NDIM; i++) gcov_avg[0][i] += geom->gcov[0][i]*weight ;
# endif

    }
    
    denom = 1./denom;

# if( USE_STRICT_ARRAY_BOUNDS ) 
    DLOOP2 { geom->gcov[i][j] = gcov_avg[i][j] * denom; }
# else
    for(i = 0; i < NDIM*NDIM; i++) geom->gcov[0][i] = gcov_avg[0][i] * denom;
# endif

    gcon_func(geom->gcov,geom->gcon,&(geom->g));

    break;


  default                    :
    fprintf(stderr,
	    "get_special_geometry():  Invalid value for METRIC_TYPE_CHOICE : %d \n", 
	    geom_type); 
    fflush(stderr);
    fail( FAIL_METRIC ,0);
    break;

  }


  /* Set auxiliary, derived metric functions : */ 
  geom->g_inv  = 1./geom->g;
  geom->ncon[0] = sqrt(-geom->gcon[0][0]) ; // 1/alpha
  geom->alpha   =  1./geom->ncon[0];
  geom->ncon[1] = -geom->alpha * geom->gcon[0][1];
  geom->ncon[2] = -geom->alpha * geom->gcon[0][2];
  geom->ncon[3] = -geom->alpha * geom->gcon[0][3];
  geom->beta[0] =  -geom->gcon[0][1]/geom->gcon[0][0];
  geom->beta[1] =  -geom->gcon[0][2]/geom->gcon[0][0];
  geom->beta[2] =  -geom->gcon[0][3]/geom->gcon[0][0];


  return;
}



/********************************************************************************/
/********************************************************************************
 test_geom
 -------------
  -- Calculates various identites of the metric and connection to check 
     for errors; 

********************************************************************************/
void test_geom(void) 
{
  int i,j,k,l, ii, jj, kk, ll;
  double f1, f2, favg, reldiff;
  double sin_th, cos_th;
  double *x, *xp, xspher[NDIM];
  double dx_dxp_dxp[NDIM][NDIM][NDIM]; 
  double dxc_dxs[NDIM][NDIM],dxs_dxc[NDIM][NDIM];
  double dxdxpdxp[NDIM][NDIM][NDIM];
  double ***conn, ***conn2;
  double mk_gcov[NDIM][NDIM], mk_gcon[NDIM][NDIM];
  double ks_gcov[NDIM][NDIM], ks_gcon[NDIM][NDIM];
  double bl_gcov[NDIM][NDIM], bl_gcon[NDIM][NDIM];
  double ks_cart_gcov[NDIM][NDIM], ks_cart_gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  double vcon3[NDIM], vcon2[NDIM], vcon[NDIM];
  double vcov3[NDIM], vcov2[NDIM], vcov[NDIM];
  double del_xp, xp_h[NDIM], xp_l[NDIM], x_h[NDIM], x_l[NDIM];
  double xp_hh[NDIM], xp_lh[NDIM], xp_hl[NDIM], xp_ll[NDIM];
  double x_hh[NDIM], x_lh[NDIM], x_hl[NDIM], x_ll[NDIM];
  struct of_geom *geom; 
  struct of_coord *coords; 

  double  ranc(int iseed); 

  ALLOC_3D_ARRAY(conn ,NDIM,NDIM,NDIM);
  ALLOC_3D_ARRAY(conn2,NDIM,NDIM,NDIM);

  LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom);
    get_coord(   i,j,k,CENT,ncurr,coords);
    dx_dxp_dxp_calc( coords->x, coords->xp, dx_dxp_dxp );
    x = coords->x;
    xp = coords->xp;
    ks_conn_func( coords, conn );
    //    ks_cart_conn_func( x, xp, conn );
    ks_gcov_func( x, ks_gcov ); 
    ks_gcon_func( x, ks_gcon ); 
    ks_cart_gcov_func( x, ks_cart_gcov ); 
    ks_cart_gcon_func( x, ks_cart_gcon ); 
    bl_gcov_func( x, bl_gcov ); 
    bl_gcon_func( x, bl_gcon ); 
    mink_spherical_gcov_func( x, mk_gcov);
    mink_spherical_gcon_func( x, mk_gcon);
    trigvalues( x[TH], &cos_th, &sin_th );


    /* dx_dxp() and dxp_dx() (back and forth) */ 
    vcon[TT] = 1.231;
    vcon[RR] = 1.34242;
    vcon[TH] = 0.1231 ; 
    vcon[PH] = 34.213;
    for( ii = 0 ; ii < NDIM ; ii++ ) {  vcon2[ii] = vcon[ii] ; } 
    transform_rank1con( x, xp, vcon2 );   // transform to xp 
    for( ii = 0 ; ii < NDIM ; ii++ ) { 
      vcon3[ii] = 0.;
      for( jj = 0 ; jj < NDIM ; jj++ ) {  
	vcon3[ii] += vcon2[jj] * coords->dx_dxp[ii][jj] ; 
      }
    }
    for( ii = 0 ; ii < NDIM ; ii++ ) { vcon2[ii] = REL_DIFF_FUNC(vcon[ii],vcon3[ii]); }
    for( ii = 0 ; ii < NDIM ; ii++ ) {
      fprintf(stdout,"test-geom:dx_dxp: (%d,%d,%d) %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,vcon[ii],vcon3[ii],vcon2[ii]); 
    }
    fflush(stdout);
    

    /* Connection symmetry:  */
    for( ii = 0 ; ii < NDIM ; ii++ ) { 
      for( jj = 0 ; jj < NDIM ; jj++ ) {  
	for( kk = jj+1 ; kk < NDIM ; kk++ ) {  
	  f1 = conn[ii][jj][kk]; f2 = conn[ii][kk][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
	  fprintf(stdout,"test-geom:consym: (%d,%d,%d) %d %d %d : %28.18e %28.18e %28.18e\n", 
		  i,j,k,ii,jj,kk,f1,f2,reldiff);
	}
      }
    }
    fflush(stdout);

    /* Connection analytic vs. numerical :  */
    conn_func( coords, geom, conn2 ); 
    for(ii=0; ii<NDIM; ii++) for(jj=0; jj<NDIM; jj++) for(kk=0; kk<NDIM; kk++) {  
      f1 = conn[ii][jj][kk] ;  f2 = conn2[ii][jj][kk] ; reldiff = REL_DIFF_FUNC(f1,f2);
	  fprintf(stdout,"test-geom:concmp: (%d,%d,%d) %d %d %d : %28.18e %28.18e %28.18e\n", 
		  i,j,k,ii,jj,kk,f1,f2,reldiff);
    }
    fflush(stdout);

//    /* KS Covariant Metric  Symmetry */
//    for( jj = 0 ; jj < NDIM ; jj++ ) {  
//      for( kk = jj+1 ; kk < NDIM ; kk++ ) {  
//	f1 = ks_gcov[jj][kk]; f2 = ks_gcov[kk][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//	fprintf(stdout,"test-geom:ksvsym: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//		i,j,k,jj,kk,f1,f2,reldiff);
//      }
//    }
//    fflush(stdout);
//
//    /* KS Contravariant Metric  Symmetry */
//    for( jj = 0 ; jj < NDIM ; jj++ ) {  
//      for( kk = jj+1 ; kk < NDIM ; kk++ ) {  
//	f1 = ks_gcon[jj][kk]; f2 = ks_gcon[kk][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//	fprintf(stdout,"test-geom:ksnsym: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//		i,j,k,jj,kk,f1,f2,reldiff);
//      }
//    }
//    fflush(stdout);
//
//    /* KS Covariant Metric  Symmetry */
//    for( jj = 0 ; jj < NDIM ; jj++ ) {  
//      for( kk = jj+1 ; kk < NDIM ; kk++ ) {  
//	f1 = bl_gcov[jj][kk]; f2 = bl_gcov[kk][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//	fprintf(stdout,"test-geom:blvsym: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//		i,j,k,jj,kk,f1,f2,reldiff);
//      }
//    }
//    fflush(stdout);
//
//    /* BL Contravariant Metric  Symmetry */
//    for( jj = 0 ; jj < NDIM ; jj++ ) {  
//      for( kk = jj+1 ; kk < NDIM ; kk++ ) {  
//	f1 = bl_gcon[jj][kk]; f2 = bl_gcon[kk][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//	fprintf(stdout,"test-geom:blnsym: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//		i,j,k,jj,kk,f1,f2,reldiff);
//      }
//    }
//    fflush(stdout);
//
//    /* Compare conn_func()  to spherical MINK  connection */
//    mink_spherical_conn_func( x, conn); 
//    conn_func( xp, mk_gcon, conn2 ); 
//    for( ii = 0 ; ii < NDIM ; ii++ ) { 
//      for( jj = 0 ; jj < NDIM ; jj++ ) {  
//	for( kk = 0 ; kk < NDIM ; kk++ ) {  
//	  f1 = conn[ii][jj][kk]; f2 = conn2[ii][jj][kk]; reldiff = REL_DIFF_FUNC(f1,f2);
//	  fprintf(stdout,"test-geom:mkconn: (%d,%d,%d) %d %d %d : %28.18e %28.18e %28.18e\n", 
//		  i,j,k,ii,jj,kk,f1,f2,reldiff);
//	}
//      }
//    }
//    fflush(stdout);
//
//    /* Transformation from BL to KS (covariant metric) */
//    //    for( ii = 0 ; ii < NDIM*NDIM; ii++ ) { gcov[0][ii] = bl_gcov[0][ii] ; } 
//    for( ii = 0 ; ii < NDIM; ii++ ) {  bl_to_ks_cov( x, bl_gcov[ii], gcov[ii] ); }
//    for( jj = 0 ; jj < NDIM; jj++ ) {
//      for( ii = 0 ; ii < NDIM; ii++ )  { vcov2[ii]    =  gcov[ii][jj]; }
//      bl_to_ks_cov( x, vcov2, vcov );
//      for( ii = 0 ; ii < NDIM; ii++ )  { gcov[ii][jj] = vcov[ii]; } 
//    }
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = ks_gcov[ii][jj]; f2 = gcov[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:b2kgcv: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//	
//
//    /* Transformation from BL to KS (contravariant metric) */
//    for( ii = 0 ; ii < NDIM; ii++ ) {  bl_to_ks_con( x, bl_gcon[ii], gcon[ii] ); }
//    for( jj = 0 ; jj < NDIM; jj++ ) {
//      for( ii = 0 ; ii < NDIM; ii++ )  { vcon2[ii]    =  gcon[ii][jj]; }
//      bl_to_ks_con( x, vcon2, vcon );
//      for( ii = 0 ; ii < NDIM; ii++ )  { gcon[ii][jj] = vcon[ii]; } 
//    }
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = ks_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:b2kgcn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//
//    /* Transformation from KS to BL (covariant metric) */
//    for( ii = 0 ; ii < NDIM; ii++ ) {  ks_to_bl_cov( x, ks_gcov[ii], gcov[ii] ); }
//    for( jj = 0 ; jj < NDIM; jj++ ) {
//      for( ii = 0 ; ii < NDIM; ii++ )  { vcov2[ii]    =  gcov[ii][jj]; }
//      ks_to_bl_cov( x, vcov2, vcov );
//      for( ii = 0 ; ii < NDIM; ii++ )  { gcov[ii][jj] = vcov[ii]; } 
//    }
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = bl_gcov[ii][jj]; f2 = gcov[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:k2bgcv: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//
//    /* Transformation from KS to BL (contravariant metric) */
//    for( ii = 0 ; ii < NDIM; ii++ ) {  ks_to_bl_con( x, ks_gcon[ii], gcon[ii] ); }
//    for( jj = 0 ; jj < NDIM; jj++ ) {
//      for( ii = 0 ; ii < NDIM; ii++ )  { vcon2[ii]    =  gcon[ii][jj]; }
//      ks_to_bl_con( x, vcon2, vcon );
//      for( ii = 0 ; ii < NDIM; ii++ )  { gcon[ii][jj] = vcon[ii]; } 
//    }
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = bl_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:k2bgcn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//
//    /* Transformation from BL to KS to BL (covariant vector) */
//    vcov[TT] = 4.4897;  vcov[RR] = 0.156;    vcov[TH] = -0.747831 ;  vcov[PH] = -10.285;
//    bl_to_ks_cov( x, vcov, vcov2 ); 
//    ks_to_bl_cov( x, vcov2, vcov3); 
//    for( ii = 0 ; ii < NDIM; ii++ ) { 
//      f1 = vcov[ii]; f2 = vcov3[ii]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:b2k2bv: (%d,%d,%d) %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//
//    /* Transformation from BL to KS to BL (contravariant vector) */
//    vcon[TT] = 3.4897;  vcon[RR] = -2.156;    vcon[TH] = 3.747831 ;  vcon[PH] = -1.285;
//    bl_to_ks_con( x, vcon, vcon2 ); 
//    ks_to_bl_con( x, vcon2, vcon3); 
//    for( ii = 0 ; ii < NDIM; ii++ ) { 
//      f1 = vcon[ii]; f2 = vcon3[ii]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:b2k2bn: (%d,%d,%d) %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//
//    /* for KS : Set gcov, invert, then compare to gcov */
//    gcon_func( ks_gcov, gcon) ; 
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = ks_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:kinvgn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//    /* for BL : Set gcov, invert, then compare to gcov */
//    gcon_func( bl_gcov, gcon) ; 
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = bl_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:binvgn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//    /* for MK-SPH : Set gcov, invert, then compare to gcov */
//    gcon_func( mk_gcov, gcon) ; 
//    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
//      f1 = mk_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:minvgn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);
//
//    /* bl_gdet test (compare to numerical determinant)*/
//    bl_gdet_func( x, &f1 ); f2 = gdet_func( bl_gcov ); reldiff = REL_DIFF_FUNC(f1,f2);
//    fprintf(stdout,"test-geom:blgdet: (%d,%d,%d) : %28.18e %28.18e %28.18e\n", 
//	    i,j,k,f1,f2,reldiff);
//
//    /* ks_gdet test (compare to numerical determinant)*/
//    ks_gdet_func( x, &f1 ); f2 = gdet_func( ks_gcov ); reldiff = REL_DIFF_FUNC(f1,f2);
//    fprintf(stdout,"test-geom:ksgdet: (%d,%d,%d) : %28.18e %28.18e %28.18e\n", 
//	    i,j,k,f1,f2,reldiff);
//
//    /* mink_gdet test (tests the determinant routine*/
//    f1 = x[RR]*x[RR]*fabs(sin_th); f2 = gdet_func( mk_gcov ); reldiff = REL_DIFF_FUNC(f1,f2);
//    fprintf(stdout,"test-geom:mkgdet: (%d,%d,%d) : %28.18e %28.18e %28.18e\n", 
//	    i,j,k,f1,f2,reldiff);


//    /* gdet tranformation test (using ks) */
//    for( ii=0; ii<NDIM*NDIM; ii++ ) { gcov[0][ii] = ks_gcov[0][ii] ; } 
//    transform_rank2cov( x, xp, gcov );
//    ks_gdet_func( x, &f1 );
//    f1 *= fabs(det_dx_dxp_calc(x,xp)); f2 = gdet_func(gcov); reldiff = REL_DIFF_FUNC(f1,f2);
//    fprintf(stdout,"test-geom:trgdet: (%d,%d,%d) : %28.18e %28.18e %28.18e\n", 
//	    i,j,k,f1,f2,reldiff);


    del_xp = 1.e-4;
    /* x_of_xp test ?? */

    /* dx_dxp with numerical value */
    for( jj=0; jj<NDIM; jj++ ) {
      for( ii=0; ii<NDIM; ii++ ) { xp_h[ii] = xp_l[ii] = xp[ii]; } 
      xp_h[jj] += del_xp;     x_of_xp( x_h, xp_h );
      xp_l[jj] -= del_xp;     x_of_xp( x_l, xp_l );
      for( ii=0; ii<NDIM; ii++ ) { gcov[ii][jj] = 0.5 * ( x_h[ii] - x_l[ii] ) / del_xp; }
    }
    for( ii=0; ii<NDIM; ii++ ) for( jj=0; jj<NDIM; jj++ )  { 
      f1 = coords->dx_dxp[ii][jj]; f2 = gcov[ii][jj];  reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:dxdxpn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);

    /* dx_dxp_dxp with numerical value */
    for( jj=0; jj<NDIM; jj++ ) {
      for( kk=jj+1; kk<NDIM; kk++ ) {
	for( ii=0; ii<NDIM; ii++ ) { xp_hh[ii] = xp[ii]; } 
	for( ii=0; ii<NDIM; ii++ ) { xp_hl[ii] = xp[ii]; } 
	for( ii=0; ii<NDIM; ii++ ) { xp_lh[ii] = xp[ii]; } 
	for( ii=0; ii<NDIM; ii++ ) { xp_ll[ii] = xp[ii]; } 
	xp_hh[jj] += del_xp;   xp_hh[kk] += del_xp;  x_of_xp( x_hh, xp_hh );
	xp_hl[jj] += del_xp;   xp_hl[kk] -= del_xp;  x_of_xp( x_hl, xp_hl );
	xp_lh[jj] -= del_xp;   xp_lh[kk] += del_xp;  x_of_xp( x_lh, xp_lh );
	xp_ll[jj] -= del_xp;   xp_ll[kk] -= del_xp;  x_of_xp( x_ll, xp_ll );
	for( ii=0; ii<NDIM; ii++ ) { 
	  dxdxpdxp[ii][jj][kk] = 0.25*(x_hh[ii]-x_hl[ii]-x_lh[ii]+x_ll[ii])/(del_xp*del_xp); 
	}
      }
      kk = jj;
      for( ii=0; ii<NDIM; ii++ ) { xp_h[ii] = xp[ii]; } 
      for( ii=0; ii<NDIM; ii++ ) { xp_l[ii] = xp[ii]; } 
      xp_h[jj] += del_xp;  x_of_xp( x_h, xp_h );
      xp_l[jj] -= del_xp;  x_of_xp( x_l, xp_l );
      for( ii=0; ii<NDIM; ii++ ) { 
	dxdxpdxp[ii][jj][kk] = ( x_h[ii] - 2*x[ii] + x_l[ii] ) / (del_xp*del_xp); 
      }
    }
    for( ii=0; ii<NDIM; ii++ ) for( jj=0; jj<NDIM; jj++ ) for( kk=jj; kk<NDIM; kk++ ) { 
      f1 = dx_dxp_dxp[ii][jj][kk]; f2 = dxdxpdxp[ii][jj][kk];  reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:dxdxp2: (%d,%d,%d) %d %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,kk,f1,f2,reldiff);
    }
    fflush(stdout);
    

    /* dxp_dx test using inversion identity of dx_dxp */
    for( ii=0; ii<NDIM; ii++ ) {
      for( jj=0; jj<NDIM; jj++ ) {
	gcov[ii][jj] = 0.; 
	for( kk=0; kk<NDIM; kk++ ) {
	  gcov[ii][jj] += coords->dx_dxp[ii][kk] * coords->dxp_dx[kk][jj];
	}
      }
    }
    for( ii=0; ii<NDIM; ii++ ) for( jj=0; jj<NDIM; jj++ )  { 
      f1 = DELTA(ii,jj); f2 = gcov[ii][jj];  reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:idxdxp: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);
    
//    /* transform_rank2cov/con  */
//    for( ii=0; ii<NDIM*NDIM; ii++ ) { gcov[0][ii] = ks_gcov[0][ii] ; } 
//    for( ii=0; ii<NDIM*NDIM; ii++ ) { gcon[0][ii] = ks_gcon[0][ii] ; } 
//    transform_rank2cov( x, xp, gcov );
//    transform_rank2con( x, xp, gcon );
//    for( ii=0; ii<NDIM; ii++ )  for( jj=0; jj<NDIM; jj++ ) {
//      f2 = 0.; 
//      for( kk=0; kk<NDIM; kk++ ) {
//	f2 += gcov[ii][kk] * gcon[kk][jj] ; 
//      }
//      f1 = DELTA(ii,jj);  reldiff = REL_DIFF_FUNC(f1,f2);
//      fprintf(stdout,"test-geom:invtrn: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
//	      i,j,k,ii,jj,f1,f2,reldiff);
//    }
//    fflush(stdout);


    /* KS pos  Cartesian -> Spherical -> Cartesian test : */ 
    ks_cart_to_ks_spher_pos(x, xspher);
    ks_spher_to_ks_cart_pos(xspher, xp);
    for( ii = 0 ; ii < NDIM; ii++ ) { 
      f1 = x[ii]; f2 = xp[ii]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:ks-cart-spher-pos: (%d,%d,%d) %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,f1,f2,reldiff);
      fprintf(stdout,"test-geom:ks-cart-spher-pos2: (%d,%d,%d) %d : %28.18e \n", 
	      i,j,k,ii,xspher[ii]);
    }
    fflush(stdout);
    

    /* dxs_dxc test using inversion identity of dx_dxp */
    ks_dxs_dxc_calc(xspher, x, dxs_dxc); 
    ks_dxc_dxs_calc(x, xspher, dxc_dxs);
    for( ii=0; ii<NDIM; ii++ ) {
      for( jj=0; jj<NDIM; jj++ ) {
	gcov[ii][jj] = 0.; 
	for( kk=0; kk<NDIM; kk++ ) {
	  gcov[ii][jj] += dxs_dxc[ii][kk] * dxc_dxs[kk][jj];
	}
      }
    }
    for( ii=0; ii<NDIM; ii++ ) for( jj=0; jj<NDIM; jj++ )  { 
      f1 = DELTA(ii,jj); f2 = gcov[ii][jj];  reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:idxcdxs: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);
    
    
    /* KS vcon  Cartesian -> Spherical -> Cartesian test : */ 
    vcon[RR] = ranc(0); 
    vcon[TH] = ranc(0); 
    vcon[PH] = ranc(0); 
    setutcon(vcon,geom->gcov); 
    ks_cart_to_ks_spher_con(x, xspher, vcon, vcon2);
    ks_spher_to_ks_cart_con(xspher, x, vcon2,vcon3);
    for( ii = 0 ; ii < NDIM; ii++ ) { 
      f1 = vcon[ii]; f2 = vcon3[ii]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:ks-cart-spher-vcon: (%d,%d,%d) %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,f1,f2,reldiff);
    }
    fflush(stdout);
    
    
    /* KS vcov  Cartesian -> Spherical -> Cartesian test : */ 
    lower(vcon,geom,vcov); 
    ks_cart_to_ks_spher_cov(x, xspher, vcov, vcov2);
    ks_spher_to_ks_cart_cov(xspher, x, vcov2,vcov3);
    for( ii = 0 ; ii < NDIM; ii++ ) { 
      f1 = vcov[ii]; f2 = vcov3[ii]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:ks-cart-spher-vcov: (%d,%d,%d) %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,f1,f2,reldiff);
    }
    fflush(stdout);

    /* KS Cartesian gcov/gcon inverse test : */ 
    for( ii=0; ii<NDIM; ii++ ) {
      for( jj=0; jj<NDIM; jj++ ) {
	gcov[ii][jj] = 0.; 
	for( kk=0; kk<NDIM; kk++ ) {
	  gcov[ii][jj] += ks_cart_gcon[ii][kk] * ks_cart_gcov[kk][jj];
	}
      }
    }
    for( ii=0; ii<NDIM; ii++ ) for( jj=0; jj<NDIM; jj++ )  { 
      f1 = DELTA(ii,jj); f2 = gcov[ii][jj];  reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:ks-cart-gcov-gcon-id: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);
    
    /* Transformation from KS-cart to KS-spher (covariant metric) */
    for( ii = 0 ; ii < NDIM; ii++ ) for( jj = 0 ; jj < NDIM; jj++ ) { 
      gcov[ii][jj] = 0.;
      for( kk = 0 ; kk < NDIM; kk++ ) for( ll = 0 ; ll < NDIM; ll++ ) {  
	gcov[ii][jj] += dxc_dxs[kk][ii] * dxc_dxs[ll][jj] * ks_cart_gcov[kk][ll]; 
      }
    }
    ks_gcov_func( xspher, ks_gcov ); 
    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
      f1 = ks_gcov[ii][jj]; f2 = gcov[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:KS-cart-spher-gcov: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);


    /* Transformation from KS-cart to KS-spher (contravariant metric) */
    for( ii = 0 ; ii < NDIM; ii++ ) for( jj = 0 ; jj < NDIM; jj++ ) { 
      gcon[ii][jj] = 0.;
      for( kk = 0 ; kk < NDIM; kk++ ) for( ll = 0 ; ll < NDIM; ll++ ) {  
	gcon[ii][jj] += dxs_dxc[ii][kk] * dxs_dxc[jj][ll] * ks_cart_gcon[kk][ll]; 
      }
    }
    ks_gcon_func( xspher, ks_gcon ); 
    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
      f1 = ks_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:KS-cart-spher-gcon: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);

    /* Transformation from KS-spher to KS-cart (covariant metric) */
    for( ii = 0 ; ii < NDIM; ii++ ) for( jj = 0 ; jj < NDIM; jj++ ) { 
      gcov[ii][jj] = 0.;
      for( kk = 0 ; kk < NDIM; kk++ ) for( ll = 0 ; ll < NDIM; ll++ ) {  
	gcov[ii][jj] += dxs_dxc[kk][ii] * dxs_dxc[ll][jj] * ks_gcov[kk][ll]; 
      }
    }
    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
      f1 = ks_cart_gcov[ii][jj]; f2 = gcov[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:KS-spher-cart-gcov: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);


    /* Transformation from KS-spher to KS-cart (contravariant metric) */
    for( ii = 0 ; ii < NDIM; ii++ ) for( jj = 0 ; jj < NDIM; jj++ ) { 
      gcon[ii][jj] = 0.;
      for( kk = 0 ; kk < NDIM; kk++ ) for( ll = 0 ; ll < NDIM; ll++ ) {  
	gcon[ii][jj] += dxc_dxs[ii][kk] * dxc_dxs[jj][ll] * ks_gcon[kk][ll]; 
      }
    }
    for( ii = 0 ; ii < NDIM; ii++ )  for( jj = 0 ; jj < NDIM; jj++ ) { 
      f1 = ks_cart_gcon[ii][jj]; f2 = gcon[ii][jj]; reldiff = REL_DIFF_FUNC(f1,f2);
      fprintf(stdout,"test-geom:KS-spher-cart-gcov: (%d,%d,%d) %d %d : %28.18e %28.18e %28.18e\n", 
	      i,j,k,ii,jj,f1,f2,reldiff);
    }
    fflush(stdout);

  }

  DEALLOC_3D_ARRAY(conn ,NDIM,NDIM,NDIM);
  DEALLOC_3D_ARRAY(conn2,NDIM,NDIM,NDIM);

  return;
}  

/********************************************************************************/
/********************************************************************************
 check_proximity_to_pt_from_cell_pos()
 ---------------------
  -- given the position and cell index, returns the distance from position to a 
     given point; 
  -- handles various exceptions to prevent out of bounds issues;
********************************************************************************/
int check_proximity_to_pt_from_cell_pos(int i, int j, int k, int pos, double *xpt, double rsq )
{

  double *xc;
  double delx[NDIM-1];
  struct of_coord *coords;

  if( (i >= N1TOT) || (j >= N2TOT) || (k >= N3TOT) ) { 
    double xp2[NDIM], x2[NDIM];
    coord(i,j,k,pos,xp2);    

#if( TOP_TYPE_CHOICE == TOP_CARTESIAN )  
    x_of_xp( x2, xp2 ) ; 
#elif( TOP_TYPE_CHOICE == TOP_SPHERICAL )
    double xspher2[NDIM];
    x_of_xp( xspher2, xp2 ) ; 
    xcart_of_xspher_only( x2 , xspher2);
#else
    fprintf(stdout,"check_proximity_to_pt_from_cell_pos() TOP_CYLINDRICAL choice not implemented yet!\n"); fflush(stdout);
    fail(FAIL_BASIC,0);
#endif
    xc = x2; 
  }
  else { 
    get_coord(i,j,k,pos,ncurr,coords);
    xc = coords->xcart; 
  }

  delx[0] = xpt[XX] - xc[XX];
  delx[1] = xpt[YY] - xc[YY];
  delx[2] = xpt[ZZ] - xc[ZZ];

  double dist_sq = delx[0]*delx[0] + delx[1]*delx[1] + delx[2]*delx[2] ; 

  return( (dist_sq <= rsq) );

}

/********************************************************************************/
/********************************************************************************
 set_excision_mask():
 ---------------------
  -- sets evol_mask[] (our mask array) to indicate the excision volume;
  -- currently only built for general metrics (read binary black hole spacetimes);
  -- currently setup to excise a set of nuclei or spherical-like regions based
      on horizon radius of each set;
  -- to excise the whole black-hole horizon, set EXCISION_FUZZ to 2.0 below,
     otherwise keep its default 1.1; and EXCISION_BIGFUZZ to some number slightly
     larger than it;
  -- Excision should be done within the horizon in order to ensure that nothing 
     errant propagates out from that inner surface;  
********************************************************************************/
#define EXCISION_FUZZ    (1.1)
#define EXCISION_BIGFUZZ (2.)
//#define EXCISION_FUZZ    (2.0)
//#define EXCISION_BIGFUZZ (2.2)
#define MAX_N_NUCLEI     (2)

void set_excision_mask( short unsigned int ***evol_mask_loc )
{
  int i,j,k,n, ii,jj,kk ;
  int n_nuclei;
  double rsq_nuclei[MAX_N_NUCLEI];
  double r_nuclei[MAX_N_NUCLEI];
  double rsq_excised[MAX_N_NUCLEI];
  double rsq_close[MAX_N_NUCLEI];
  double loc_nuclei[MAX_N_NUCLEI][NDIM];
  double bbox_nuclei[MAX_N_NUCLEI][2][NDIM];

  unsigned short int val, overlap;
  unsigned short int n_overlap_dims[MAX_N_NUCLEI];

  double *x, x2[NDIM], xp2[NDIM], delx[NDIM],xspher2[NDIM];
  double *xg;
  double dist_sq;
  struct of_coord *coords;
  struct of_bbh_traj bbh_traj;

  static unsigned short int local_first_time = 1;

  TRACE_BEG;

  /* Setup the parameters for the max routine specific to our problem */
  n_nuclei = 2;
  if( n_nuclei > MAX_N_NUCLEI ) { 
    fprintf(stderr,"set_excision_mask():  n_nuclei > MAX_N_NUCLEI  :  %d  %d \n", n_nuclei, MAX_N_NUCLEI); 
    fflush(stderr); 
    fail(FAIL_BASIC,0); 
  }

  i = 0;
#if( EXCISION_TYPE_CHOICE == EXCISE_HORIZON )
  r_nuclei[i++] = 0.3*r_horizon1; 
  r_nuclei[i++] = 0.3*r_horizon2; 
#else
  r_nuclei[i++] = 1.5; 
  r_nuclei[i++] = 1.5; 
#endif  
  i = 0;
  rsq_nuclei[i] = r_nuclei[i]*r_nuclei[i] ; i++;
  rsq_nuclei[i] = r_nuclei[i]*r_nuclei[i] ; 

  get_bbh_traj_data(&bbh_traj) ;
  i = j = 0;
  loc_nuclei[i][j++] = bbh_traj.tt  ;
  loc_nuclei[i][j++] = bbh_traj.xi1x;
  loc_nuclei[i][j++] = bbh_traj.xi1y;
  loc_nuclei[i][j++] = bbh_traj.xi1z;

  i++;  j=0; 
  loc_nuclei[i][j++] = bbh_traj.tt  ;
  loc_nuclei[i][j++] = bbh_traj.xi2x;
  loc_nuclei[i][j++] = bbh_traj.xi2y;
  loc_nuclei[i][j++] = bbh_traj.xi2z;

  for(i=0;i<n_nuclei;i++) { 
    rsq_excised[i] = EXCISION_FUZZ    * EXCISION_FUZZ    * rsq_nuclei[i]; 
    rsq_close[i]   = EXCISION_BIGFUZZ * EXCISION_BIGFUZZ * rsq_nuclei[i];
  }

#if( TOP_TYPE_CHOICE == TOP_CARTESIAN )
  for(i=0;i<n_nuclei;i++) { 
    bbox_nuclei[i][0][0] = bbox_nuclei[i][1][0] = 0.;
    for(j=1;j<NDIM;j++) { 
      bbox_nuclei[i][0][j] = loc_nuclei[i][j] - EXCISION_FUZZ * r_nuclei[i];
    }
    for(j=1;j<NDIM;j++) { 
      bbox_nuclei[i][1][j] = loc_nuclei[i][j] + EXCISION_FUZZ * r_nuclei[i];
    }
  }
#elif( TOP_TYPE_CHOICE == TOP_SPHERICAL )
  extern void xspher_of_xcart_only(double *xspher, double *xcart);
  double drtmp;
  for(i=0;i<n_nuclei;i++) { 
    bbox_nuclei[i][0][0] = bbox_nuclei[i][1][0] = 0.;
    xspher_of_xcart_only( x2, loc_nuclei[i] ); 
    drtmp = EXCISION_FUZZ * r_nuclei[i];
    if( x2[RR] > SMALL ) { 
      bbox_nuclei[i][0][RR] = x2[RR] - drtmp;
      bbox_nuclei[i][1][RR] = x2[RR] + drtmp;
      bbox_nuclei[i][0][TH] = x2[TH] - drtmp/x2[RR];
      bbox_nuclei[i][1][TH] = x2[TH] + drtmp/x2[RR];
      bbox_nuclei[i][0][PH] = x2[PH] - drtmp/x2[RR];
      bbox_nuclei[i][1][PH] = x2[PH] + drtmp/x2[RR];
      /* forcing phi range to be between [0,2pi[ */
      while( bbox_nuclei[i][0][PH] >= 2.*M_PI ) bbox_nuclei[i][0][PH] -= 2.*M_PI ;
      while( bbox_nuclei[i][0][PH] < 0.       ) bbox_nuclei[i][0][PH] += 2.*M_PI ;
      while( bbox_nuclei[i][1][PH] >= 2.*M_PI ) bbox_nuclei[i][1][PH] -= 2.*M_PI ;
      while( bbox_nuclei[i][1][PH] < 0.       ) bbox_nuclei[i][1][PH] += 2.*M_PI ;
    }
    else { 
      bbox_nuclei[i][0][RR] = -SMALL;
      bbox_nuclei[i][1][RR] = drtmp;
      bbox_nuclei[i][0][TH] = -SMALL;
      bbox_nuclei[i][1][TH] = M_PI + SMALL;
      bbox_nuclei[i][0][PH] = -SMALL;
      bbox_nuclei[i][1][PH] = 2.*M_PI + SMALL;
    }
  }
#else
  fprintf(stdout,"set_excision_mask() TOP_CYLINDRICAL choice not implemented yet!\n"); fflush(stdout);
  fail(FAIL_BASIC,0);
#endif
  

  /* Output the local static parameters at the first call : */
  if( (myid == printer_pid) && local_first_time ) { 
    fprintf(stdout,"set_excision_mask(): n_nuclei             = %d \n",n_nuclei);
    for(i=0;i<n_nuclei;i++) { 
      fprintf(stdout,"set_excision_mask(): rsq_nuclei[%03d]   = %26.16e \n",i,rsq_nuclei[i]);
    }
    for(i=0;i<n_nuclei;i++) { 
      fprintf(stdout,"set_excision_mask(): rsq_excised[%03d]  = %26.16e \n",i,rsq_excised[i]);
    }
    for(i=0;i<n_nuclei;i++) { 
      fprintf(stdout,"set_excision_mask(): rsq_close[%03d]    = %26.16e \n",i,rsq_close[i]);
    }
    for(i=0;i<n_nuclei;i++) { 
      for(j=0;j<NDIM;j++) { 
	fprintf(stdout,"set_excision_mask(): loc_nuclei[%03d][%1d] = %26.16e \n",i,j,loc_nuclei[i][j]);
      }
    }
    for(i=0;i<n_nuclei;i++) { 
      for(j=0;j<2;j++) { 
	for(k=0;k<NDIM;k++) { 
	  fprintf(stdout,"set_excision_mask(): bbox_nuclei[%03d][%1d][%1d] = %26.16e \n",i,j,k,bbox_nuclei[i][j][k]);
	}
      }
    }
    local_first_time = 0; 
  }

  /**************************************************************
    First see if the horizons overlap with this node's grid: 
     -- note that the boundary box comparison *should* be done in xp[1-3] coordinates, but 
        sometimes  the x[1-3] ->  xp[1-3]  transformation is difficult; 
        -- knowing this and recognizing that x[1-3] is often suitable if we set the boundary box to the 
             extremal values measured on the domain when the coordinates are calculated in calc_all_coord(); 
            -- thus, we need to transform the horizon's bbox's to the appropriate coordinate system 
                  depending on topology choice:
  *****************************************************************/
  for(i=0;i<n_nuclei;i++) { 
    n_overlap_dims[i] = 0; 
    /* fewest comparisons to verify overlap I believe: */
    j=1; 
    if( (x1_min <= bbox_nuclei[i][1][j]) && (bbox_nuclei[i][0][j] <= x1_max) ) {  n_overlap_dims[i]++;  }
    j++;
    if( (x2_min <= bbox_nuclei[i][1][j]) && (bbox_nuclei[i][0][j] <= x2_max) ) {  n_overlap_dims[i]++;  }
    j++;
    /* due to the discontinuity of the phi coordinate we need to take care if we're in the 1st or 4th quadrant */
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
    if( x3_min < 0.5*M_PI ) { bbox_nuclei[i][0][j] -= 2*M_PI ; }
    if( x3_max > 1.5*M_PI ) { bbox_nuclei[i][1][j] += 2*M_PI ; }
#endif
    if( (x3_min <= bbox_nuclei[i][1][j]) && (bbox_nuclei[i][0][j] <= x3_max) ) {  n_overlap_dims[i]++;  }
  }

  overlap = 0; 
  for(i=0;i<n_nuclei;i++) { 
    if( n_overlap_dims[i] == (NDIM-1) ) { overlap = 1; }
  }

#if( 0 )
  fprintf(stdout,"mask: overlap = %d \n", overlap); 
  fprintf(stdout,"mask: bbox lower = %26.16e %26.16e %26.16e \n", x1_min,x2_min,x3_min);
  fprintf(stdout,"mask: bbox upper = %26.16e %26.16e %26.16e \n", x1_max,x2_max,x3_max);
  fprintf(stdout,"mask: n_overlap_dims = %d %d \n", n_overlap_dims[0],n_overlap_dims[1]);
  fflush(stdout);
#endif

  /* Return now if no intersection with local domain: */
  if( !overlap ) { 
    ALL_LOOP {evol_mask_loc[i][j][k] = MASK_NORMAL;} 
    TRACE_END;
    return;
  }


  /************************************************************************************************
   All cell loop to set excised region: 
      -- need to go into ghosts to avoid the flux calculation on the excised surface
  *************************************************************************************************/
  overlap = 0;

  ALL_LOOP { 
    val = MASK_NORMAL;  /* assume we are evolving */

    get_coord(i,j,k,CENT,ncurr,coords);
    xg = coords->xcart; 

    for(n=0;n<n_nuclei;n++) { 
      delx[XX] = xg[XX] - loc_nuclei[n][XX];
      delx[YY] = xg[YY] - loc_nuclei[n][YY];
      delx[ZZ] = xg[ZZ] - loc_nuclei[n][ZZ];
    
      /* Only check the nearest corner if we are close to the excision region: */
      dist_sq = delx[XX]*delx[XX] + delx[YY]*delx[YY] + delx[ZZ]*delx[ZZ] ; 

      if( dist_sq <= rsq_excised[n] ) { 
	val = MASK_EXCISED; 
	overlap = 1; 	
	break;
      }

      /* Continue with stricter more precise test if we are at least close to the excision volume: 
	 -- we do not have to do better than test each position in the cell as no other point in the cell is used by reconstruction, etc. */
      if( dist_sq <= rsq_close[n] ) { 
	/* 6 faces and 8 corners: */
	if( check_proximity_to_pt_from_cell_pos(i  ,j  ,k  ,FACE1,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j  ,k  ,FACE2,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j  ,k  ,FACE3,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j  ,k  ,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i+1,j  ,k  ,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i+1,j  ,k  ,FACE1,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j+1,k  ,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j+1,k  ,FACE2,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j  ,k+1,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j  ,k+1,FACE3,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i+1,j+1,k  ,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i  ,j+1,k+1,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i+1,j  ,k+1,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
	if( check_proximity_to_pt_from_cell_pos(i+1,j+1,k+1,CORN ,loc_nuclei[n],rsq_excised[n]) ) {  val = MASK_EXCISED;  overlap = 1; break; }
      }
    }  /* n_nuclei  loop */

    evol_mask_loc[i][j][k] = val;

  }

  if( !overlap ) { 
    TRACE_END;
    return; 
  }
   
  /************************************************************************************************
   All cell loop to set buffer region: 
      -- defined as those cells sharing faces to excised cells;
      -- we can safely neglect setting the outermost cells on the domain to be buffer cells;
!!!-- can we optimize this further???
  *************************************************************************************************/
  /* Rules for excision: 
     1) Excised  points lie somewhere within the event horizon in order to mask out singular, divergent, or 
     "hard to invert" points; 
     
     2) Buffer zones surround the Excised zones so that Evolved (aka "Normal") zones do not use any fluxes from 
        Excised zones;

     3) Buffer zones should also buffer enough to prevent an Evolved cell from using the cell-centered values for its 
        reconstruction process;
     
     Implications:
     -- Rule 1 implies that we excise within some radius that is a fraction of the event horizon's radius;
     -- Rule 2 implies that we have to have enough buffer zones so that the FluxCT stencil of a Evolved zone does
        not intersect that of a Excised zone; 
	The FluxCT stencils (for linear and parabolic CT schemes) are given in comments in flux_ct.c; 
     -- Rule 3 implies that NG zones per dimenssion on each side of a Excised cell should be buffered;
	 
  */
  
  int itmp,jtmp,ktmp;
  short unsigned int *evtmp;  

    ALL_LOOP { 
      if( evol_mask_loc[i][j][k] == MASK_EXCISED ) { 
	for(ii=-BUFFER_WIDTH;ii<=BUFFER_WIDTH;ii++) {
	  itmp = i+ii;
	  if( itmp <  0     ) { continue; }
	  if( itmp >= N1TOT ) { continue; }

	  for(jj=-BUFFER_WIDTH;jj<=BUFFER_WIDTH;jj++) { 
	    jtmp = j+jj;
	    if( jtmp <  0     ) { continue; }
	    if( jtmp >= N2TOT ) { continue; }

	    for(kk=-BUFFER_WIDTH;kk<=BUFFER_WIDTH;kk++) { 
	      ktmp = k+kk;
	      if( ktmp <  0     ) { continue; }
	      if( ktmp >= N3TOT ) { continue; }
	      
	      evtmp = &(evol_mask_loc[itmp][jtmp][ktmp]);
	      if( *evtmp == MASK_NORMAL ) { 
		*evtmp = MASK_BUFFER; 
	      }
	    }
	  }
	}
      }
    }
	
  TRACE_END;

  return;
}

#undef EXCISION_FUZZ 
#undef EXCISION_BIGFUZZ 
#undef MAX_N_NUCLEI  


#undef M

#undef DEL_X_CONN
