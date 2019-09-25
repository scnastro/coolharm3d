
#include "decs.h"
#include "metric.h"

#define VERBOSE  (0) 

static int levi_civita[NDIM][NDIM][NDIM][NDIM];

/* Global predicate arrays used to do reconstruction and find fluxes : */ 
static const int idir[] = {1,0,0};
static const int jdir[] = {0,1,0};
static const int kdir[] = {0,0,1};


static void misc_source(double *ph, struct of_state *q, struct of_geom *geom,
			int ii, int jj, int kk, double *dU, double Dt);
static void source_cooling_func(double *ph,struct of_state *q,int ii,int jj,int kk,double *dU, int bound );

/** more Physics **/
/***********************************************************************************************/
/***********************************************************************************************
  primtoflux():
  ---------
   --  calculate fluxes in direction dir, 
        
***********************************************************************************************/

void primtoflux(double *pr, struct of_state *q, int dir, struct of_geom *geom, double *flux) 
{
	int j,k,l ;
	register double eta_ucon,bcon_dir;
	double mhd[NDIM] ;

	//	mhd_calc(pr, dir, q, mhd) ;

	eta_ucon = q->eta * q->ucon[dir];
	bcon_dir =  q->bcon[dir];
	mhd[0] = eta_ucon * q->ucov[0] - bcon_dir * q->bcov[0] ;
	mhd[1] = eta_ucon * q->ucov[1] - bcon_dir * q->bcov[1] ;
	mhd[2] = eta_ucon * q->ucov[2] - bcon_dir * q->bcov[2] ;
	mhd[3] = eta_ucon * q->ucov[3] - bcon_dir * q->bcov[3] ;

	mhd[dir] += q->ptot;

	/* particle number flux */
	flux[RHO] = pr[RHO]*q->ucon[dir] ;

	/* MHD stress-energy tensor w/ first index up, 
	 * second index down. */
	flux[UU] = mhd[0] + flux[RHO] ;
	flux[U1] = mhd[1] ;
	flux[U2] = mhd[2] ;
	flux[U3] = mhd[3] ;

	/* dual of Maxwell tensor */
	flux[B1]  = q->bcon[1]*q->ucon[dir] - q->bcon[dir]*q->ucon[1] ;
	flux[B2]  = q->bcon[2]*q->ucon[dir] - q->bcon[dir]*q->ucon[2] ;
	flux[B3]  = q->bcon[3]*q->ucon[dir] - q->bcon[dir]*q->ucon[3] ;

	PLOOP flux[l] *= geom->g;

	return; 
}

/* calculate "conserved" quantities; provided strictly for
 * historical reasons */
void primtoU(double *pr, struct of_state *q, struct of_geom *geom, double *U)
{

	primtoflux(pr,q,0,geom, U) ;
	return ;
}

/* calculate magnetic field four-vector */
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) 
{
  register double inv_ut;
  //  register unsigned int j;

  inv_ut = 1./ucon[TT];
  bcon[TT] = pr[B1]*ucov[1] + pr[B2]*ucov[2] + pr[B3]*ucov[3] ;
  bcon[RR] = (pr[B1] + bcon[TT]*ucon[RR])*inv_ut ;
  bcon[TH] = (pr[B2] + bcon[TT]*ucon[TH])*inv_ut ;
  bcon[PH] = (pr[B3] + bcon[TT]*ucon[PH])*inv_ut ;
	
//	for(j=1;j<4;j++)
//		bcon[j] = (pr[B1-1+j] + bcon[TT]*ucon[j])/ucon[TT] ;

	return ;
}

/* MHD stress tensor, with first index up, second index down */
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd) 
{
  register double eta_ucon,bcon_dir;
  //  int i ;

//  P = (gam - 1.)*pr[UU] ;
//  eta = P + pr[RHO] + pr[UU] + q->bsq ;
//  ptot = P + 0.5*q->bsq;

  /* single row of mhd stress tensor, 
   * first index up, second index down */
//  DLOOP1 mhd[i] = eta * q->ucon[dir] * q->ucov[i] 
//                      - q->bcon[dir] * q->bcov[i] ;
  
  eta_ucon = q->eta * q->ucon[dir];
  bcon_dir =  q->bcon[dir];
  mhd[0] = eta_ucon * q->ucov[0] - bcon_dir * q->bcov[0] ;
  mhd[1] = eta_ucon * q->ucov[1] - bcon_dir * q->bcov[1] ;
  mhd[2] = eta_ucon * q->ucov[2] - bcon_dir * q->bcov[2] ;
  mhd[3] = eta_ucon * q->ucov[3] - bcon_dir * q->bcov[3] ;
  
  mhd[dir] += q->ptot;

}

/* add in source terms to equations of motion */
/* Note that for stationary phi-indep. metrics and since we are using the densitized 
   stress tensor, the T^t_t and T^t_\phi sources are zero.  */
void source(double *ph, struct of_geom *geom, int ii, int jj, int kk, double *dU, double Dt)
{
  int i,j,k,l,bound ;
	double mhd[NDIM][NDIM];
	struct of_state q ;

	//-fast 	PLOOP dU[l] = 0. ;
	dU[0] = 0. ;	dU[1] = 0. ;	dU[2] = 0. ;	dU[3] = 0. ;
	dU[4] = 0. ;	dU[5] = 0. ;	dU[6] = 0. ;	dU[7] = 0. ;

#if( METRIC_DIM > 0 )
	l = CONN_ID(ii,jj,kk);
	get_state(ph,geom,&q) ;

	mhd_calc(ph, 0, &q, mhd[0]) ;	mhd_calc(ph, 1, &q, mhd[1]) ;
	mhd_calc(ph, 2, &q, mhd[2]) ;	mhd_calc(ph, 3, &q, mhd[3]) ;

#if( USE_PRESSURE_FLUX_FIX )
	SDLOOP1 {  mhd[i][i] -= q.ptot; }
#endif

	/* contract mhd stress tensor with connection */
	DLOOP2 {

		dU[U1] += mhd[i][j]*conn[l][j][1][i] ;
		dU[U2] += mhd[i][j]*conn[l][j][2][i] ;

#if( METRIC_DIM > 2 )
		dU[UU] += mhd[i][j]*conn[l][j][0][i] ;
		dU[U3] += mhd[i][j]*conn[l][j][3][i] ;
#endif

	}

#endif 

	//misc_source(ph, ii, jj, kk, geom, &q, dU, Dt) ;



#if( USE_COOLING_FUNCTION ) 
#if( METRIC_DIM == 0 )
	get_state(ph,geom,&q) ;  /* calculate state since we haven't done so already */
#endif

	/* only cool for bounded matter:  */
	bound = (  (q.ucov[0]*( q.p + ph[UU] + ph[RHO] + q.bsq ))  >  (-ph[RHO]) ) ?  1 : 0 ; 
	  
	source_cooling_func(ph, &q,ii,jj,kk,dU,bound);
#endif 

#if( USE_KINEMATIC_VISCOSITY )
#if( METRIC_DIM == 0 )
	get_state(ph,geom,&q) ;  /* calculate state since we haven't done so already */
#endif
	set_visc_stress_functions(n_substep, ii, jj, kk, ph, &q);
#endif 

	
	dU[UU] *= geom->g ;
	dU[U1] *= geom->g ;
	dU[U2] *= geom->g ;
	dU[U3] *= geom->g ;

	/* done! */
	return;
}

/* returns b^2 (i.e., twice magnetic pressure) */
double bsq_calc(double *pr, struct of_geom *geom)
{
	struct of_state q ;

	get_state(pr,geom,&q) ;
	return( q.bsq ) ;
}

/* find ucon, ucov, bcon, bcov from primitive variables */
void get_state(double *pr, struct of_geom *geom, struct of_state *q)
{

	/* get ucon */
	ucon_calc(pr, geom, q->ucon) ;
	lower(q->ucon, geom, q->ucov) ;
	bcon_calc(pr, q->ucon, q->ucov, q->bcon) ;
	lower(q->bcon, geom, q->bcov) ;
	q->bsq = DOT(q->bcon,q->bcov);
	q->p   = (gam - 1.)*pr[UU];
	q->eta = q->p + pr[RHO] + pr[UU] + q->bsq ;
	q->ptot = q->p + 0.5*q->bsq;

	return ;
}

/* find contravariant four-velocity */
void ucon_calc(double *pr, struct of_geom *geom, double *ucon)
{
	double gamma ;
	//	double alpha,beta[NDIM] ;
	//	int i ;

	//	alpha = 1./sqrt(-geom->gcon[TT][TT]) ;
	//	SDLOOP1 beta[i] = geom->gcon[TT][i]*alpha*alpha ;

	if( gamma_calc(pr,geom,&gamma) ) { 
	  fflush(stderr);
	  fprintf(stderr,"\nucon_calc(): gamma failure \n");
	  fflush(stderr);
	  fail(FAIL_GAMMA_CALC,0);
	}

	//	ucon[TT] = gamma/alpha ;
	//	SDLOOP1 ucon[i] = pr[U1+i-1] - ucon[TT]*beta[i];

	ucon[TT] = gamma*geom->ncon[0] ;
	ucon[1] = pr[U1] - ucon[TT]*geom->beta[0];
	ucon[2] = pr[U2] - ucon[TT]*geom->beta[1];
	ucon[3] = pr[U3] - ucon[TT]*geom->beta[2];

	return ;
}

/* find gamma-factor wrt normal observer */
int gamma_calc(double *pr, struct of_geom *geom, double *gamma)
{
  double qsq ;

  qsq =     geom->gcov[1][1]*pr[U1]*pr[U1]
    + geom->gcov[2][2]*pr[U2]*pr[U2]
    + geom->gcov[3][3]*pr[U3]*pr[U3]
    + 2*(geom->gcov[1][2]*pr[U1]*pr[U2]
	 + geom->gcov[1][3]*pr[U1]*pr[U3]
	 + geom->gcov[2][3]*pr[U2]*pr[U3]) ;

  if( qsq < 0. ){
    if( fabs(qsq) > 1.E-10 ){ // then assume not just machine precision
      fprintf(stderr,"gamma_calc():  failed: i,j,qsq = %d %d %28.18e \n", icurr,jcurr,qsq);
      fprintf(stderr,"v[1-3] = %28.18e %28.18e %28.18e  \n",pr[U1],pr[U2],pr[U3]);
      *gamma = 1.;
      return (1);
    }
    else qsq=1.E-10; // set floor
  }

  *gamma = sqrt(1. + qsq);

  return(0) ;
}

/********************************************************************
 gamma_calc2(): 
 --------------
   -- Alternate way of calculating gamma, using Taylor expansion of 
        sqrt(1+x) instead, for large x and for small x ;
*********************************************************************/
int gamma_calc2(double *pr, struct of_geom *geom, double *gamma)
{
  double qsq ;
  double *gcovtmp;

  qsq =  
    pr[U1] * ( geom->gcov[1][1]*pr[U1] + geom->gcov[1][2]*pr[U2] + geom->gcov[1][3]*pr[U3] ) + 
    pr[U2] * ( geom->gcov[1][2]*pr[U1] + geom->gcov[2][2]*pr[U2] + geom->gcov[2][3]*pr[U3] ) + 
    pr[U3] * ( geom->gcov[1][3]*pr[U1] + geom->gcov[2][3]*pr[U2] + geom->gcov[3][3]*pr[U3] ) ;

  if( qsq < 0. ){
    if( fabs(qsq) > 1.E-10 ){ // then assume not just machine precision
      fprintf(stderr,"gamma_calc():  failed: i,j,qsq = %d %d %28.18e \n", icurr,jcurr,qsq);
      fprintf(stderr,"v[1-3] = %28.18e %28.18e %28.18e  \n",pr[U1],pr[U2],pr[U3]);
      *gamma = 1.;
      return (1);
    }
    else qsq=1.E-10; // set floor
  }

  *gamma = MY_SQRT_ONE_PLUS_EPS(qsq);

  return(0) ;
}

/* finds the primitive velocity components from the 4-velocity */
void ucon2pr( double *pr, double ucon[NDIM], double gcon[][NDIM] )
{
  double f_tmp ;

  f_tmp = ucon[TT] / gcon[TT][TT] ; 
  pr[ U1] = ucon[RR] - gcon[TT][RR] * f_tmp ;
  pr[ U2] = ucon[TH] - gcon[TT][TH] * f_tmp ;
  pr[ U3] = ucon[PH] - gcon[TT][PH] * f_tmp ;

  return ;
}


/*  
 * VCHAR():
 * 
 * calculate components of magnetosonic velocity 
 * corresponding to primitive variables p 
 *
 * cfg 7-10-01
 * 
 */

void vchar(double *pr, struct of_state *q, struct of_geom *geom, int js, 
		double *vmax, double *vmin)
{
	double discr,vp,vm,bsq,EE,EF,va2,cs2,cms2,rho,u ;
	double Acov[NDIM],Bcov[NDIM],Acon[NDIM],Bcon[NDIM] ;
	double Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,C ;
	int i,j ;


	DLOOP1 Acov[i] = 0. ;
	Acov[js] = 1. ;
	raise(Acov,geom,Acon) ;

	DLOOP1 Bcov[i] = 0. ;
	Bcov[TT] = 1. ;
	raise(Bcov,geom,Bcon) ;

	/* find fast magnetosonic speed */
	bsq = q->bsq;
	rho = pr[RHO] ;
	u = pr[UU] ;
	EF = rho + gam*u ;
	EE = bsq + EF ;
	va2 = bsq/EE ;
	cs2 = gam*(gam - 1.)*u/EF ;

//	if(cs2 < 0.) cs2 = SMALL ;
//	if(cs2 > 1.) cs2 = 1. ;
//	if(va2 < 0.) va2 = SMALL ;
//	if(va2 > 1.) va2 = 1. ;

	cms2 = cs2 + va2 - cs2*va2 ;	/* and there it is... */

	//cms2 *= 1.1 ;

	/* check on it! */
	if(cms2 < 0.) {
	  fprintf(stdout,"vcharneg: %g %g %g %g %g %g %g \n", bsq, rho,u,EF,EE,va2,cs2);
	  fflush(stdout);
	  fail(FAIL_VCHAR_NEG,0) ;
	  cms2 = SMALL ;
	}
	if(cms2 > 1.) {
	  fprintf(stdout,"vcharsup: %g %g %g %g %g %g %g \n", bsq, rho,u,EF,EE,va2,cs2);
	  fflush(stdout);
	  fail(FAIL_VCHAR_SUPER,0) ; 
	  cms2 = 1. ;
	}

	/* now require that speed of wave measured by observer 
	   q->ucon is cms2 */
	Asq = DOT(Acon,Acov) ;
	Bsq = DOT(Bcon,Bcov) ;
	Au =  DOT(Acov,q->ucon) ;
	Bu =  DOT(Bcov,q->ucon) ;
	AB =  DOT(Acon,Bcov) ;
	Au2 = Au*Au ;
	Bu2 = Bu*Bu ;
	AuBu = Au*Bu ;

	A =      Bu2  - (Bsq + Bu2)*cms2 ;
	B = 2.*( AuBu - (AB + AuBu)*cms2 ) ;
	C =      Au2  - (Asq + Au2)*cms2 ;

	discr = B*B - 4.*A*C ;
	if((discr<0.0)&&(discr>-1.e-10)) discr=0.0;
	else if(discr < -1.e-10) {
		fprintf(stderr,"\n\t %g %g %g %g %g\n",A,B,C,discr,cms2) ;
		fprintf(stderr,"\n\t q->ucon: %g %g %g %g\n",q->ucon[0],q->ucon[1],
				q->ucon[2],q->ucon[3]) ;
		fprintf(stderr,"\n\t q->bcon: %g %g %g %g\n",q->bcon[0],q->bcon[1],
				q->bcon[2],q->bcon[3]) ;
		fprintf(stderr,"\n\t Acon: %g %g %g %g\n",Acon[0],Acon[1],
				Acon[2],Acon[3]) ;
		fprintf(stderr,"\n\t Bcon: %g %g %g %g\n",Bcon[0],Bcon[1],
				Bcon[2],Bcon[3]) ;
		fail(FAIL_VCHAR_DISCR,0) ;
		discr = 0. ;
	}

	discr = sqrt(discr) ;
	vp = -(-B + discr)/(2.*A) ;
	vm = -(-B - discr)/(2.*A) ;

	if(vp > vm) {
		*vmax = vp ;
		*vmin = vm ;
	}
	else {
		*vmax = vm ;
		*vmin = vp ;
	}

	return ;
}


/*  
 * VCHAR_fast():
 *  -- faster version of vchar();
 * calculate components of magnetosonic velocity 
 * corresponding to primitive variables p 
 *
 * SCN 2009
 * 
 */

void vchar_fast(struct of_state *q, struct of_geom *geom, int js, 
		double *vmax, double *vmin)
{
  register double discr, cms2;
  register double Au2,Bu2,AuBu,A,B,C;
  register double vp,vm, denom ;

	cms2 = (gam * q->p  + q->bsq ) / q->eta ;

	//cms2 *= 1.1 ;

	/* check on it! */
	if(cms2 < 0.) {
	  fprintf(stdout,"vcharneg: %g %g %g \n", q->p, q->bsq, q->eta);
	  fflush(stdout);
	  fail(FAIL_VCHAR_NEG,0) ;
	  cms2 = SMALL ;
	}
	else if(cms2 > 1.) {
	  fprintf(stdout,"vcharsup: %g %g %g \n", q->p, q->bsq, q->eta);
	  fflush(stdout);
	  fail(FAIL_VCHAR_SUPER,0) ; 
	  cms2 = 1. ;
	}

	/* now require that speed of wave measured by observer 
	   q->ucon is cms2 */
	//	Asq = Acon[js];  // = gcon[js][js]
	//	Bsq = Bcon[TT];  // = gcon[TT][TT]
	//	Au =  q->ucon[js];
	//	Bu =  q->ucon[TT];
	//	AB =  Acon[TT];  // = gcon[TT][js] 
	Au2 = q->ucon[js]*q->ucon[js];
	Bu2 = q->ucon[TT]*q->ucon[TT] ;
	AuBu = q->ucon[js]*q->ucon[TT] ;

	A =      Bu2  - (geom->gcon[TT][TT] + Bu2 )*cms2 ;
	B =      AuBu - (geom->gcon[TT][js] + AuBu)*cms2 ;
	C =      Au2  - (geom->gcon[js][js] + Au2 )*cms2 ;

	discr = B*B - A*C ;
	if( discr < 0. ) { 
	  if( discr > -1.e-10 ) {  discr=0.0; }
	  else { 
	    fprintf(stderr,"\n\t %g %g %g %g %g\n",A,B,C,discr,cms2) ;
	    fprintf(stderr,"\n\t q->ucon: %g %g %g %g\n",q->ucon[0],q->ucon[1],
		    q->ucon[2],q->ucon[3]) ;
	    fprintf(stderr,"\n\t q->bcon: %g %g %g %g\n",q->bcon[0],q->bcon[1],
		    q->bcon[2],q->bcon[3]) ;
	    fail(FAIL_VCHAR_DISCR,0) ;
	    discr = 0. ;
	  }
	}

	discr = sqrt(discr) ;
	denom = 1. / A; 
	vp = -(-B + discr)*denom ;
	vm = -(-B - discr)*denom ;

	if(vp > vm) {
		*vmax = vp ;
		*vmin = vm ;
	}
	else {
		*vmax = vm ;
		*vmin = vp ;
	}

	return ;
}
/********************************************************************************/
/********************************************************************************
 vchar_fast_luminal():
 --------------------
********************************************************************************/
//void vchar_fast_luminal(struct of_state *q, struct of_geom *geom, int js, 
//			double *vmax, double *vmin)
//
//lmin1 = (gcon01 + sqrt(gcon01*gcon01-gcon00*gcon11))/gcon00
//lmax1 = (gcon01 - sqrt(gcon01*gcon01-gcon00*gcon11))/gcon00
//
//lmin2 = (gcon02 + sqrt(gcon02*gcon02-gcon00*gcon22))/gcon00
//lmax2 = (gcon02 - sqrt(gcon02*gcon02-gcon00*gcon22))/gcon00
//
//lmin3 = (gcon03 + sqrt(gcon03*gcon03-gcon00*gcon33))/gcon00
//lmax3 = (gcon03 - sqrt(gcon03*gcon03-gcon00*gcon33))/gcon00
//
//lall1 = [abs(lmin1),abs(lmax1)]
//lall2 = [abs(lmin2),abs(lmax2)]
//lall3 = [abs(lmin3),abs(lmax3)]
//
//lgt1 = abs(lmin1) > abs(lmax1)
//     lgt2 = abs(lmin2) > abs(lmax2)
//     lgt3 = abs(lmin3) > abs(lmax3)
//

/******************************************************************************/
/******************************************************************************
 source_cooling_func():
 ---------------------------
   -- general set of source terms used when employing a cooling function that is 
      isotropic in the fluid frame; 

   --   T^\mu_{\nu ; \mu} = -f(\rho,u) u_\nu

        where u is the internal internal energy density, \rho is the rest-mass 
        density and f() is the cooling rate in the fluid rest frame. 

******************************************************************************/
static void source_cooling_func(double *ph, struct of_state *q,
				int ii, int jj, int kk, double *dU, int bound )
{
  unsigned int ind;
  double f_cool=0., dU_cool[NDIM];

//   /* Since we assume axisymmetric near-Keplerian flow, check for correct coordinates */
// #if((METRIC_TYPE_CHOICE!=METRIC_KS_SPHERICAL) && (METRIC_TYPE_CHOICE!=METRIC_BL_SPHERICAL))
//   fprintf(stderr,
// 	  "cooling_func_hr_disk(): Coordinate specific cooling function bad METRIC_TYPE_CHOICE");
//   fflush(stderr);
//   fail(FAIL_BASIC,0); 
// #endif 

  /* Set scalar (local fluid frame) cooling function : */
#if( USE_COOLING_FUNCTION == 1 ) 
  if( bound ) { f_cool = cooling_func_hr_disk( ii, jj, kk, ph ); } 
#elif( USE_COOLING_FUNCTION == 2 || USE_COOLING_FUNCTION == 3 ) 
  if( bound ) { f_cool = cooling_func_isentropic_disk( ii, jj, kk, ph, q ); } 
#elif( USE_COOLING_FUNCTION >= 4 )
  ind = kk-N3S + N3*((jj-N2S) + N2*(ii-N1S));
  f_cool = coolflux[0][ind] ;
#endif


  /*  dU_a  = -f u_a  */ 
   dU[UU ]  -= f_cool * q->ucov[0] ;
   dU[U1 ]  -= f_cool * q->ucov[1] ;
   dU[U2 ]  -= f_cool * q->ucov[2] ;
   dU[U3 ]  -= f_cool * q->ucov[3] ;

#if( MAKE_RADFLUX ) 
# if( USE_COOLING_FUNCTION < 4 )
  if( n_substep == (N0-1) ) {
    ind = kk-N3S + N3*((jj-N2S) + N2*(ii-N1S));
    coolflux[0][ind] = f_cool;
  }
# endif
#endif

  return;
}


/******************************************************************************/
/******************************************************************************
 cooling_func_hr_disk():
 ---------------------------
   -- cooling function to equilibrate a thin disk so that it tends to a 
       near steady constant H/R; 
   -- cooling function was derived by Julian Krolik to turn off if the disk 
       gets too thin;
   -- after performing a run at a=0.9, we found that the disk became
      unexpectedly thinner at small radius.  This was because we were 
      neglecting a relativistic factor in the temperature equation. 
      See Section 7.3.2 "Relativistic Effects" of Krolik's AGN book. 
   -- Also, we feel that the density scale height, 
        \int d\theta d\phi  \sqrt{-g}  \sqrt{g_{\theta\theta}} |\theta - \pi/2| \rho 
        ----------------------------------------------------------------------------
        \int d\theta d\phi  \sqrt{-g}  \rho 

      is the correct measure for the height of the disk.  The canonical height
      (i.e. the one used in cooling_func_hr_disk_old()) is really the standard
      deviation of the resulting Gaussian profile.  If you integrate a Gaussian
      in the scale height equation you will get that 
   
          H_{scaleheight} =  \sqrt{\pi/2} \sigma 
  
       where \sigma is standard deviation from the Gaussian profile:  
                   exp[ - ( z/(2 \sigma) )^2 ]

       This factor is included in our equation here. 
 
    -- Another correction from cooling_func_hr_disk_old() is the temperature
       profile we were using was off by a factor of (gam-1).  This factor is included 
       now. 

   -- The temperature profile, T(r) = P(r,z)/\rho(r,z), is derived from the 
        pressure balance equation: 

           dP/dz =  - \rho z  G_z(r)

       After assuming that P(r,z) = T(r) \rho(r,z), one can integrate this 
       equation to get  \rho(r,z) --- a Gaussian --- with 
          \sigma = \sqrt{ T / G_z }
       This density profile yields a scale height: 
            H = \sqrt{ \pi T / ( 2 G_z ) } 
       This then gives us the temperature profile: 

           Pi G_z H^2    Pi G_z r^2  
       T = ----------- = ----------  h_o_r^2 
             2              2        

   --  The switch 
       E = u(g-1)/(rho T)

      f(rho,u) = s u Omega [ (E - 1) + |E - 1| ]^q

      where s \lesssim 1 and q \gtrsim 1.    
      The expression (E-1) + |E-1| > 0 if E > 1, but = 0 if E < 1.   
      That way, there is cooling only when the temperature exceeds 
      the target temperature, and the cooling rate goes continuously 
      to zero as E approaches unity from above.   The cooling rate is 
      of order \Omega, but slower if s < 1.    If we want very rapid 
      equilibration, we can choose q > 1.


******************************************************************************/
double cooling_func_hr_disk( int i, int jj, int kk, double *ph)
{
  unsigned int j,k;
  double E_ratio, T_r, dE, Omega_r;
  double ucov[NDIM],ucon[NDIM],ucon_isco[NDIM],ucov_isco[NDIM];

  const double h_o_r  = 5.e-2;  /* Value of A=H/R at which to aim */ 
  const double q_exp  = 0.5;     /* Exponent to control cooling rate */
  const double s_norm = 1.;     /* Overall O(1) normalization factor */ 
  
  struct of_coord coords_loc;
  struct of_coord *coords;
  static int n1tot_old = 0;
  static double *Omega_gf=NULL; 
  static double   *T_r_gf=NULL;
  static unsigned int local_first_time = 1; 

#if( (USE_COOLING_FUNCTION != 1)  && \
     (USE_COOLING_FUNCTION != 4)  && \
     (USE_COOLING_FUNCTION != 5) )
    fprintf(stderr,"cooling_func_hr_disk(): Bad value USE_COOLING_FUNCTION  should not be herer  %d   !! \n",USE_COOLING_FUNCTION); 
    fflush(stderr);    fail(FAIL_BASIC,0);
#endif

  if( local_first_time ) { 
    if( myid == printer_pid ) { 
      fprintf(stdout,"cooling_func_hr_disk(): \n##################################################\n");
      fprintf(stdout,"cooling_func_hr_disk(): PARAMETERS \n------------------------------------\n");
      fprintf(stdout,"cooling_func_hr_disk(): \t  s_norm         =  %28.18e \n",s_norm     ); 
      fprintf(stdout,"cooling_func_hr_disk(): \t  q_exp          =  %28.18e \n",q_exp     ); 
      fprintf(stdout,"cooling_func_hr_disk(): \t  h_o_r          =  %28.18e \n",h_o_r ); 
      fprintf(stdout,"cooling_func_hr_disk(): \t  rate method    =  %d      \n",USE_COOLING_FUNCTION   ); 
      fprintf(stdout,"cooling_func_hr_disk(): \n------------------------------------\n");
      fflush(stdout);
    }
  }

#if((COORD_TYPE_CHOICE==COORD_IDENTITY)||(COORD_TYPE_CHOICE==COORD_DIAGONAL)||(COORD_TYPE_CHOICE==COORD_MIXED)||(COORD_TYPE_CHOICE==COORD_DIAGONAL2)||(COORD_TYPE_CHOICE==COORD_DIAGONAL3)) 
  if( local_first_time || (n1tot_old != N1TOT) ) { 

    if( Omega_gf != NULL ) { FREE(Omega_gf);  }
    if(   T_r_gf != NULL ) { FREE(  T_r_gf);  }
    
    n1tot_old = N1TOT;
    ALLOC_ARRAY(Omega_gf,N1TOT);
    ALLOC_ARRAY(  T_r_gf,N1TOT);
    
     /* radius */
    for(j=0;j<N1TOT;j++) { 
      get_coord(j,N2S,N3S,CENT,ncurr,coords);      

      /* Orbital frequency */
      Omega_gf[j] = omega_circular_equatorial(coords->x[RR]); 

      /* The T_r_gf[] here is really T/(gam-1) */
      T_r_gf[j] = target_temperature(coords->x[RR],h_o_r); 

      /* Dump the cooling temperature profile : */
#if( VERBOSE )
      fprintf(stdout,"cooling omega temperature (%4d,%5d,%28.18e) =  %28.18e %28.18e \n", 
	      myid,j,coords->x[RR],Omega_gf[j],T_r_gf[j]);     fflush(stdout);
#endif
    }
    local_first_time = 0;
  }
  
  Omega_r = Omega_gf[i];
  T_r     =   T_r_gf[i];

#elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
  static int nlast = -10;

  if( (nlast != ncurr) || (n1tot_old != N1TOT) || local_first_time ) {

    if( n1tot_old != N1TOT ) { 
      if( Omega_gf != NULL ) { FREE(Omega_gf);  }
      if(   T_r_gf != NULL ) { FREE(  T_r_gf);  }
    
      n1tot_old = N1TOT;
      ALLOC_ARRAY(Omega_gf,N1TOT);
      ALLOC_ARRAY(  T_r_gf,N1TOT);
    }

    /* radius */
    for(j=0;j<N1TOT;j++) { 
      get_coord(j,N2S,N3S,CENT,ncurr,coords);      

      /* Orbital frequency */
      Omega_gf[j] = omega_circular_equatorial(coords->x[RR]); 

      /* The T_r_gf[] here is really T/(gam-1) */
      T_r_gf[j] = target_temperature(coords->x[RR],h_o_r); 

      /* Dump the cooling temperature profile : */
    }
    nlast = ncurr; 
  }
  if( local_first_time ) { 
    /* Dump the cooling temperature profile : */
#if( VERBOSE )
    for(j=0;j<N1TOT;j++) { 
      fprintf(stdout,"cooling omega temperature (%4d,%5d,%28.18e) =  %28.18e %28.18e \n", 
	      myid,j,coords->x[RR],Omega_gf[j],T_r_gf[j]); 
    }
    fflush(stdout);
#endif 
    local_first_time = 0;
  }

  Omega_r = Omega_gf[i];
  T_r     =   T_r_gf[i];


# elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )

  get_coord(i,jj,kk,CENT,ncurr,coords);      
  Omega_r = omega_circular_equatorial(coords->x[RR]); 
  T_r     = target_temperature(coords->x[RR],h_o_r); 

#else 
  
  fprintf(stderr,"cooling_func_hr_disk():  need to implement this for COORD_TYPE_CHOICE = %d \n", 
	  COORD_TYPE_CHOICE); 
  fflush(stderr);  fail(FAIL_BASIC,0);
#endif


  E_ratio = (gam-1.)*ph[UU] / (ph[RHO] * T_r) ; 
  dE = E_ratio - 1.;

  return( 
	 s_norm * ph[UU] * Omega_r
	 * pow( (dE + fabs(dE)), q_exp ) 
	 );

}

/******************************************************************************/
/******************************************************************************
 target_temperature():
 ---------------------------
   -- Returns with the value of the target temperature at a specific radial 
      index; 

           Pi G_z H^2    Pi G_z r^2  
       T = ----------- = ----------  h_o_r^2 
             2              2        
 
     where G_z = R_z/r^3 ,  R_z is given in Krolik's book except there's 
     a typographical error carried over from Abramowicz, Lanza, Percival (1997), 
     whose equation for \mathcal{L}_\star  should be (immediately after eq. (27)): 
           \mathcal{L}_\star = \mathcal{L}^2 - a^2(E^2 - 1)

   -- Use ISCO values of u_t and u_phi  for  r < r_isco 

******************************************************************************/
double target_temperature( double r_loc, double h_o_r )
{
  double r, x, rm1, xm3, R_z, term1;
  double C_f, F_f, G_f;  /* NT functions */
  double ut_sq, uphi_sq;
  static double ut_sq_isco, uphi_sq_isco;
  static int local_first_time = 1;

  /* Set ISCO values for future use if we need them */
  if( local_first_time ) { 
    r = r_isco;
    rm1 = 1./r;     x = sqrt(r);    xm3 = 1./(x*x*x);    term1 = a*a*rm1;

    C_f = 1. - 3*rm1 + 2*a*xm3;
    F_f = 1. - 2*a*xm3 + term1*rm1;
    G_f = 1. - 2*rm1 + a*xm3 ;

    ut_sq_isco   = G_f * G_f     / C_f ;
    uphi_sq_isco = F_f * F_f * r / C_f;
    local_first_time = 0;
  }

  r = r_loc;

#if( USE_ISOTROPIC_COORDS )
  r += M;  /* Transform to BL coordinates */
#endif

  rm1 = 1./r;     x = sqrt(r);    xm3 = 1./(x*x*x);    term1 = a*a*rm1;

  if( r > r_isco ) { 
    C_f = 1. - 3*rm1 + 2*a*xm3;
    F_f = 1. - 2*a*xm3 + term1*rm1;
    G_f = 1. - 2*rm1 + a*xm3 ;
    ut_sq   = G_f * G_f     / C_f ;
    uphi_sq = F_f * F_f * r / C_f;
    R_z = uphi_sq*rm1 - term1*(ut_sq - 1.) ;
    term1 = 0.5*M_PI*R_z*h_o_r*h_o_r*rm1;
  }
  else { 
    R_z = uphi_sq_isco*rm1 - term1*(ut_sq_isco - 1.) ;
    term1 = 0.5*M_PI*R_z*h_o_r*h_o_r*rm1;
  }
  
  return( term1 ) ;
}

/******************************************************************************/
/******************************************************************************
 omega_circular_equatorial():
 ---------------------------
  -- orbital frequency of a circular orbit in the equator of the black hole 
     given the radial index of the location;
  -- assumes the orbit's energy and angular momentum are the ISCO values 
     for radii within the ISCO;
******************************************************************************/
double omega_circular_equatorial(double r_loc)
{
  double omega,vcon[NDIM],r;
  struct of_geom geom;
  struct of_coord coords;

  static double vcov_isco[NDIM];
  static int local_first_time=1;

  /* First find the ISCO values : */
  if( local_first_time ) { 
    if( REL_DIFF_FUNC(r_isco,0.) > SMALL ) { 
#if( USE_ISOTROPIC_COORDS )
      fprintf(stdout,"omega_circular_equatorial(): should not be here1 \n"); fflush(stdout);
      fail( FAIL_BASIC, 0 ) ; 
#endif    
      coord_of_r(r_isco, &coords);
      //      fprintf(stdout,"omega_circular_equatorial(): th diff =  %26.16e   %26.16e \n",coords.x[2],0.5*M_PI);  fflush(stdout);
      get_special_geometry( &coords, &geom, METRIC_TYPE_CHOICE ); 
      vel_circular_equatorial(r_isco, vcon);
      lower(vcon,&geom,vcov_isco);
    }
      local_first_time = 0;
  }
   
  r = r_loc;

  if( r >= r_isco ) { 
#if( USE_ISOTROPIC_COORDS )
    r += M;  /* Transform to BL coordinates */
#endif
    omega = sqrt(M)/(pow(r,1.5) + a); 
  }
  else { 
#if( USE_ISOTROPIC_COORDS )
    fprintf(stdout,"omega_circular_equatorial(): should not be here2 \n"); fflush(stdout);
    fail( FAIL_BASIC, 0 ) ; 
#endif    

    coord_of_r(r, &coords);
    get_special_geometry( &coords, &geom, METRIC_TYPE_CHOICE ); 
    raise(vcov_isco,&geom,vcon);
    omega = vcon[3]/vcon[0];
  }

  return(omega);
}

//Dennis new function here
/******************************************************************************/
/******************************************************************************
 omega_circular_equatorial_general():
 ---------------------------
  -- orbital frequency of a circular orbit in the equator of the black hole 
     given the radial index of the location;
  -- assumes the orbit's energy and angular momentum are the ISCO values 
     for radii within the ISCO;
******************************************************************************/
double omega_circular_equatorial_general(double r_loc, double spin, double mass)
{
  double omega,vcon[NDIM],r, isco;
  struct of_geom geom;
  struct of_coord coords;

  static double vcov_isco[NDIM];
  static int local_first_time=1;

  /* First find the ISCO values : */
  isco = risco_calc_general(1, spin, mass);
  r = r_loc;
/* #if( USE_ISENTROPIC_COORDS )  */
/*   r += mass; */
/* #endif */
  if( r <= isco ) { r = isco; }
  
  omega = sqrt(mass)/(pow(r,1.5) + spin*sqrt(mass));

  //  fprintf(stdout,"omega-circ(r=%g,a=%g) = %28.18e \n", r,a,omega); fflush(stdout);
  
  return(omega);
}
//Dennis end new function

/******************************************************************************/
/******************************************************************************
 vel_circular_equatorial():
 ---------------------------
  -- returns with the 4-velocity of a circular orbit in the equator of the black hole 
     given the radial index of the location;
  -- valid for MKS, KS and BL coordinates; 
  -- verifies that the 4-velocity is finite (checks for coordinate singularities)
  -- assumes  central gravitating mass M=1. 
******************************************************************************/
void vel_circular_equatorial(double r, double vcon[NDIM])
{
  double uphi,r3halves; 

  r3halves = pow(r,1.5);
  uphi = 2*a*r3halves + r*r*(r-3.) ;
  if( uphi <= 0. ) {
    fprintf(stdout,"vel_circular_equatorial(): Invalid radius , r = %28.18e \n",r);
    fprintf(stdout,"vel_circular_equatorial(): Returning velocity at the ISCO ... \n");
    fflush(stdout);
    r = r_isco; 
    r3halves = pow(r,1.5);
    uphi = 2*a*r3halves + r*r*(r-3.) ;
  }
  
  uphi = 1. / sqrt( uphi ); 
  vcon[1] = vcon[2] = 0.;
  vcon[3] = uphi ; 
  vcon[0] = uphi * ( r3halves + a );

  return;
}


/******************************************************************************/
/******************************************************************************
 cooling_func_isentropic_disk():
 ---------------------------
   -- originally written for GRHydro 

   -- Calculates the local, fluid-frame cooling function using local quantities. 
 
   -- this is a spherically symmetric version assuming spherical coordinates ;
       -- also assumes that the "PH" or 3rd coordinate is not different from \phi, 
            i.e. it is uniform and without a scaling factor; 

   -- This method implicitly assumes isotropy (independence of emission angle) 
      as it is a scalar quantity.   Multiply it by the 4-velocity of the fluid 
      and you get a flux. 

   -- Unlike the cooling function that has been used in HARM3d in the past, this one 
      cools toward a constant value of entropy, thereby eliminating the gauge 
      dependency of the constant H/R profile.  The target entropy should be 
      that used in the initial disk in order to track dissipation. 

   -- cooling function was derived by Julian Krolik to turn off if the disk 
       gets too thin;

   --  The switch 
       E =  S/S_t   (S is local entropy and S_t is target entropy)

      f(rho,u) = s u f [ (E - 1) + |E - 1| ]^q

      where s \lesssim 1 and q \gtrsim 1    

              and 

             f = Omega / 2 pi   (Omega = "orbital  velocity" ) 

      The expression (E-1) + |E-1| > 0 if E > 1, but = 0 if E < 1.   
      That way, there is cooling only when the entropy  exceeds 
      the target entropy, and the cooling rate goes continuously 
      to zero as E approaches unity from above.   The cooling rate is 
      of order \Omega, but slower if s < 1.    If we want very rapid 
      equilibration, we can choose q > 1.

******************************************************************************/
double cooling_func_isentropic_disk( int i,  int jj, int kk, double *ph, struct of_state *q)
{
  int j;
  double E_ratio, dE, entropy_loc, r, cooling_rate=0.;
  struct of_coord *coords;

  const double target_entropy  = 1.e-2;  /* Entropy to cool towards */
  const double q_exp  = 0.5;             /* Exponent to control cooling rate */
  const double s_norm = 1.;              /* Overall O(1) normalization factor */ 

  static int n1tot_old = 0;
  static double *freq_orbit_gf=NULL;
  static double inv_2pi = (0.5/M_PI) ; 
  static unsigned int local_first_time = 1; 

#if( BBH_SPACETIME )
  struct of_bbh_traj bbh_traj;
  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  get_bbh_traj_data(&bbh_traj);
  double xx, xx1, xx2, yy, yy1, yy2, zz, zz1, zz2;
  double r1, r2, omega;
  double r_min_orbit =  1.5 * initial_bbh_separation ;
#endif 

#if((USE_COOLING_FUNCTION != 2) && \
    (USE_COOLING_FUNCTION != 3) )
  fprintf(stderr,"cooling_func_isentropic_disk(): Bad value USE_COOLING_FUNCTION  should not be herer  %d   !! \n",USE_COOLING_FUNCTION); 
  fflush(stderr);    fail(FAIL_BASIC,0);
#endif

  if( local_first_time ) { 
    if( myid == printer_pid ) { 
      fprintf(stdout,"cooling_func_isentropic_disk(): \n##################################################\n");
      fprintf(stdout,"cooling_func_isentropic_disk(): PARAMETERS \n------------------------------------\n");
      fprintf(stdout,"cooling_func_isentropic_disk(): \t  s_norm         =  %28.18e \n",s_norm     ); 
      fprintf(stdout,"cooling_func_isentropic_disk(): \t  q_exp          =  %28.18e \n",q_exp     ); 
      fprintf(stdout,"cooling_func_isentropic_disk(): \t  target_entropy =  %28.18e \n",target_entropy ); 
      fprintf(stdout,"cooling_func_isentropic_disk(): \t  rate method    =  %d      \n",USE_COOLING_FUNCTION   ); 
      fprintf(stdout,"cooling_func_isentropic_disk(): \n------------------------------------\n");
      fflush(stdout);
    }
  }

#if( USE_COOLING_FUNCTION == 2 ) 

# if((COORD_TYPE_CHOICE==COORD_IDENTITY)||(COORD_TYPE_CHOICE==COORD_DIAGONAL)||(COORD_TYPE_CHOICE==COORD_MIXED)||(COORD_TYPE_CHOICE==COORD_DIAGONAL2)||(COORD_TYPE_CHOICE==COORD_DIAGONAL3)) 
  if( local_first_time || (n1tot_old != N1TOT) ) { 
    if( freq_orbit_gf != NULL ) { FREE(freq_orbit_gf);  }
    n1tot_old = N1TOT;
    ALLOC_ARRAY(freq_orbit_gf,N1TOT);

    /* Orbital frequency */
    for(j=0;j<N1TOT;j++) { 
      get_coord(j,N2S,N3S,CENT,ncurr,coords);      
      freq_orbit_gf[j] = inv_2pi * omega_circular_equatorial(coords->x[RR]); 
#if( VERBOSE )
      fprintf(stdout,"cooling rate (%4d,%5d,%28.18e) =  %28.18e \n", myid,j,coords->x[RR],freq_orbit_gf[j]);    fflush(stdout);
#endif
    }
    local_first_time = 0;
  }

  cooling_rate = freq_orbit_gf[i] ;

# elif( COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD )
  static int nlast = -10;

  if( (nlast != ncurr) || (n1tot_old != N1TOT) || local_first_time ) {
    if( n1tot_old != N1TOT ) { 
      if( freq_orbit_gf != NULL ) { FREE(freq_orbit_gf);  }
      n1tot_old = N1TOT;
      ALLOC_ARRAY(freq_orbit_gf,N1TOT);
    }

    /* radius */
    for(j=0;j<N1TOT;j++) { 
      get_coord(j,N2S,N3S,CENT,ncurr,coords);      
      freq_orbit_gf[j] = inv_2pi * omega_circular_equatorial(coords->x[RR]); 
    }
    nlast = ncurr; 
  }
  if( local_first_time ) { 
    /* Dump the cooling temperature profile : */
#if( VERBOSE )
    for(j=0;j<N1TOT;j++) { 
      fprintf(stdout,"cooling rate (%4d,%5d,%28.18e) =  %28.18e \n", myid,j,coords->x[RR],freq_orbit_gf[j]); 
    }
    fflush(stdout);
#endif
    local_first_time = 0;
  }

  cooling_rate = freq_orbit_gf[i] ;

# elif( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )

  get_coord(i,jj,kk,CENT,ncurr,coords);      
  r = coords->r ; 

#  if( BBH_SPACETIME )
  /*Find distance to BHs and see if part of mini-disk*/
  xx  = coords->xcart[XX]; yy  = coords->xcart[YY]; zz  = coords->xcart[ZZ];
  xx1 = bbh_traj.xi1x;     yy1 = bbh_traj.xi1y;     zz1 = bbh_traj.xi1z;
  xx2 = bbh_traj.xi2x;     yy2 = bbh_traj.xi2y;     zz2 = bbh_traj.xi2z;

  r1 = sqrt( (xx-xx1)*(xx-xx1) + (yy-yy1)*(yy-yy1) + (zz-zz1)*(zz-zz1) );
  r2 = sqrt( (xx-xx2)*(xx-xx2) + (yy-yy2)*(yy-yy2) + (zz-zz2)*(zz-zz2) );
  /* Calculate rdisk properly from coord transforms */
  double xNZ[NDIM], xIZ[NDIM], xBL[NDIM];
  xNZ[TT] = t ; xNZ[XX] = xx;
  xNZ[YY] = yy; xNZ[ZZ] = zz;
  
  if( r1 <= 0.45 * bbh_traj.r12 ) { 
    coordsIZ_of_coordsNZ( MBH1, xNZ, xIZ );
    coordsBL_of_coordsIZ( m_bh1, 0., xIZ, xBL);
    omega = omega_circular_equatorial_general(xBL[RR],0,m_bh1); 
  }
  else if( r2 <= 0.45 * bbh_traj.r12 ) {
    coordsIZ_of_coordsNZ( MBH2, xNZ, xIZ );
    coordsBL_of_coordsIZ( m_bh2, 0., xIZ, xBL );
    omega = omega_circular_equatorial_general(xBL[RR],0,m_bh2);
  }
  else if( r < r_min_orbit ) { omega = omega_circular_equatorial(r_min_orbit); }
  else { omega = omega_circular_equatorial(r); }
  //if( r < r_min_orbit ) { r = r_min_orbit; }
  cooling_rate = inv_2pi * omega;
#  else
  cooling_rate = inv_2pi * omega_circular_equatorial(r); 
#  endif

# else
    fprintf(stderr,"cooling_func_isentropic_disk(): We assume that r[i] does not change with j and t !! \n"); 
    fflush(stderr);    fail(FAIL_BASIC,0);
# endif

#endif

    local_first_time = 0;

#if( USE_COOLING_FUNCTION == 3 ) 
      cooling_rate = inv_2pi * q->ucon[PH] / q->ucon[TT] ;    /*  u^\phi / u^t  / 2pi  */ 
#endif

    entropy_loc = (gam-1.)*ph[UU] / pow(ph[RHO],gam);
    E_ratio = entropy_loc / target_entropy;
    dE = E_ratio - 1.;

  return( 
	 s_norm * ph[UU] * cooling_rate
	 * pow( (dE + fabs(dE)), q_exp ) 
	 );

}


/******************************************************************************/
/******************************************************************************
 misc_source():
 ---------------------------
  -- example source function;
******************************************************************************/
static void misc_source(double *ph, struct of_state *q, struct of_geom *geom,
			int ii, int jj, int kk, double *dU, double Dt)
{
  
  /* This is merely an example and does not represent any physical source term that */
  /* I can think of.  Place your calculation for the extra source terms here.  */
  dU[RHO]  += ph[RHO] ;
  dU[UU ]  += ph[UU ] ;
  dU[U1 ]  += ph[U1 ] ;
  dU[U2 ]  += ph[U2 ] ;
  dU[U3 ]  += ph[U3 ] ;

}

/**************************************************************************************/
/***************************************************************************************
  set_levi_civita():
  --------------
     -- calculates the 4d anti-symmetric array once so that we don't have to 
        do hundreds of boolean operations per cell per time step to calculate 
        the Faraday tensor;
***************************************************************************************/
void set_levi_civita(void)
{
  int i,j,k,l,n,index[NDIM], do_sort, n_perm,val,n_swap;

  for(i=0;i<NDIM;i++)  for(j=0;j<NDIM;j++)  for(k=0;k<NDIM;k++)  for(l=0;l<NDIM;l++) { 
    /* The easy ones: */
    if( i==j || i==k || i==l || j==k || j==l || k==l ) { 
      levi_civita[i][j][k][l] = 0;
    }
    else{
      index[0] = i; index[1] = j; index[2] = k; index[3] = l; 
      do_sort = 1; 
      n_perm = 0; 
      while( do_sort ) { 
	n_swap = 0;
	for(n=0;n<NDIM-1;n++) { 
	  if( index[n] > index[n+1] ) { 
	    n_perm++;   
	    n_swap++;
	    val        = index[n  ];	    
	    index[n  ] = index[n+1];	    
	    index[n+1] = val;
	  }
	}
	do_sort = n_swap; 
      }
      levi_civita[i][j][k][l] = ( n_perm % 2 )  ? -1 : 1 ;
    }
  }

  /* Test levi-civita : */
//  for(i=0;i<NDIM;i++)  for(j=0;j<NDIM;j++)  for(k=0;k<NDIM;k++)  for(l=0;l<NDIM;l++) { 
//    fprintf(stdout,"levi-civita[%d%d%d%d] = %d \n", i,j,k,l,levi_civita[i][j][k][l] ); 
//    fflush(stdout);
//  }

}

/****************************************************************************************/
/*****************************************************************************************
  faraday_calc():
  --------------
     -- calculates the contravariant components of the Faraday tensor times sqrt{-g} ;
         or rather just: 

      fcon[j] =  gdet * Faraday[i][j] =  [ijkl] b_k u_l 
                
                where  i=dim here 
                  and [ijkl] is the anti-symmetric tensor;
                    
*****************************************************************************************/
void faraday_calc(int dim, struct of_state *q, double fcon[])
{
  int i,j,k,l ;

  for(i=0;i<NDIM;i++) fcon[i] = 0.;

  /* Calculate the tensor, but skip most of the terms that are components that are zero to minimize loop
     iterations */
  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++) { 
    if( (j != k) && (j != dim) && (dim != k) ) {
      for(l=0;l<NDIM;l++) { 
	fcon[j] += levi_civita[dim][j][k][l] * q->bcov[k] * q->ucov[l];
      }
    }
  }

  return; 
}

/****************************************************************************************/
/*****************************************************************************************
  current_calc():
  --------------
     -- uses faraday_calc() to calculate faraday tensor at the necessary locations;
     -- in order to reduce memory usage, we eliminate faraday calculations that are 
        never used; in order to do this, certain components of faraday tensor 
        exist only at certain locations;  it turns out that, 
        faraday[i][j][k][l][d] is located at the l-face;  
        the derivative in the l-direction should use the component that exists at the 
        l-face so that we need to use the following equation for the current:
 
       jcon[i] =  (1/gdet) partial_j ( gdet * F[i][j] ) 
                 -(1/gdet) partial_j ( gdet * F[j][i] )
                 -(1/gdet) partial_j ( faraday[j][i] )
*****************************************************************************************/
void current_calc( void ) 
{
  int i,j,k,l,d,dim ;
  int i_e,j_e,k_e,nsave;
  
  double pface[NP], inv_gdet_ij, inv_dt_old;
  struct of_state qface;
  struct of_geom   *geom;

  /* Initializing nt_curr[] this way ensures that the first call just sets faraday[] */
  static int nt_curr[2] = { -100, -100 };

  static int n          = 0; 
  static int np1        = 0; 
  static int first_time = 1;

  //  void current_calc2( void ) ;

  if( myid == printer_pid ) { fprintf(stdout,"Calculating current at n,t = %d %28.18e \n", nstep, t);  fflush(stdout); }

#if( CALC_CURRENT ) 

  /****************************************************************************************
     Set faraday[] using t=t^n time step's reconstructed values:
  ****************************************************************************************/
  /* Spatial faces : */
  FACE_LOOP {
    dim = d + 1;  // Position in a 4-vector to which this direction corresponds
    i_e = N1E + idir[d]; j_e = N2E + jdir[d]; k_e = N3E + kdir[d]; 

    for(i=N1S;i<=i_e;i++) for(j=N2S;j<=j_e;j++) for(k=N3S;k<=k_e;k++) {
      get_geometry(i,j,k,d,n_mid,geom);
      PLOOP  pface[l] = 0.5*(p_L[i][j][k][d][l] + p_R[i][j][k][d][l]);
      get_state(pface, geom, &qface);
      faraday_calc( dim, &qface, faraday[np1][i][j][k][dim] );
    }
  }
  
  /* Temporal faces: t=t^n face (cell centered) */
  dim = 0 ;
  LOOP { 
    get_geometry(i,j,k,CENT,n_beg,geom);
    get_state(p[i][j][k], geom, &qface);
    faraday_calc( dim, &qface, faraday[np1][i][j][k][dim] );
  }

  nt_curr[np1] = nstep; 

  /****************************************************************************************
    Make sure that the two Faraday tensors are at adjacent time steps, else we need to 
    recalculate the old one : 
      -- also take into account of roundoff error;
  ****************************************************************************************/
  if( (nstep - nt_curr[n])  > 1 ) {
    n   = (n  +1)%2;
    np1 = (np1+1)%2;
    return;  /* postpone current calc until we have adjacent-t faraday[]'s */
  }

  
  /****************************************************************************************
    If we have calculated both times, then we can calculate the current;
  ****************************************************************************************/
  inv_dt_old = 1./dt_old;
  LOOP {
    get_geometry(i,j,k,CENT,n_mid,geom);
    inv_gdet_ij = geom->g_inv;

    for(l=0; l<NDIM; l++) { 
      jcon[i][j][k][l] = -inv_gdet_ij * ( 
					 /* time part : */
					 (faraday[np1][i][j][k][0][l] - faraday[n][i][j][k][0][l])*inv_dt_old +

					 /* spatial part: (need to time average here) */
					 0.5 * ((faraday[n][i+1][j  ][k  ][1][l] + faraday[np1][i+1][j  ][k  ][1][l] -
						 faraday[n][i  ][j  ][k  ][1][l] - faraday[np1][i  ][j  ][k  ][1][l] ) * invdx[1] + 
						(faraday[n][i  ][j+1][k  ][2][l] + faraday[np1][i  ][j+1][k  ][2][l] -
						 faraday[n][i  ][j  ][k  ][2][l] - faraday[np1][i  ][j  ][k  ][2][l] ) * invdx[2] + 
						(faraday[n][i  ][j  ][k+1][3][l] + faraday[np1][i  ][j  ][k+1][3][l] -
						 faraday[n][i  ][j  ][k  ][3][l] - faraday[np1][i  ][j  ][k  ][3][l] ) * invdx[3] )
					  );
    }
  }

  n   = (n  +1)%2;
  np1 = (np1+1)%2;

#endif
  //  current_calc2(); 

  return; 
}

/****************************************************************************************/
/*****************************************************************************************
  calc_visc_stress_source():
  --------------
     -- calculates the source terms of the stress-energy equations that arise when there is a 
        dynamic stress tensor.  

     -- the total stress-energy tensor becomes 

        T_{a b} =  T_{a b}[MHD]  + T_{a b}[Visc]  

                where 
 
               T_{a b}[Visc] =  - 2 eta sigma_{a b}

               eta  =  rho nu 

               nu   = kinematic viscosity parameter
       
               sigma_{a b} =  shear tensor
               sigma_{a b} =  (1/2)*( u_{a ;b} + u_{b ;a} + u_a u^c u_{b ;c}+ u_b u^c u_{a ;c} )
                            - (1/3)*( u^c_{;c} * ( g_{a b} + u_a u_b ) )

     -- the source term we want to calculate is : 

         V_a  =  - sqrt(-g) \nabla_c T^c_a[Visc] =  - sqrt(-g) g^{c d} \nabla_c T_{d a}[Visc] 
              = sqrt(-g) g^{c d} \nabla_c ( 2 eta sigma_{d a} )
              = sqrt(-g) g^{c d} 2 [ sigma_{d a} \partial_c ( eta ) + eta \nabla_c ( sigma_{d a} ) ]

     --  we use 2nd-order centered differencing, so we only need one neighbor on each side; 

*****************************************************************************************/
void  calc_visc_stress_source( void )
{
#if( USE_KINEMATIC_VISCOSITY )
#define NGV      (2*NDC)
#define NV1_LOOP for(i=N1S-NGV;i<=N1E+NGV;i++) 
#define NV2_LOOP for(j=N2S-NGV;j<=N2E+NGV;j++) 
#define NV3_LOOP for(k=N3S-NGV;k<=N3E+NGV;k++) 

  int i, j, k, l;
  int n, nm1, np1;
  int nu,mu,kp,lm;

  static const unsigned short int shift[NDIM][NDIM] = {{ 1,0,0,0},
						       { 0,1,0,0},
						       { 0,0,1,0},
						       { 0,0,0,1}};

  double nu_loc, factor;
  struct of_viscosity *fmm, *f00, *fpp, *fmp, *fpm;
  struct of_geom    *geom;

  double stress[NDIM][NDIM]; 
  double ducov[NDIM][NDIM]; 
  double *ucov, *ucon;
  double *ucov1, *ucon1;
  double *ucov2, *ucon2;
  double acov[NDIM];   /* covariant 4-acceleration */
  double divu, eta ;
  const double two_thirds = (2./3.);

  double dfactors[NDIM];

  /* Mnemonics for time steps:  */
  nm1 = 0;
  n   = nm1+1;
  np1 = n+1;
  invdx[TT] = 1./dx[TT];

  for( mu = 0; mu < NDIM; mu++ ) {   dfactors[mu] = 0.5*invdx[mu];   }

  /*********************************************************************
   Loop over all physical cells and beyond to calculate viscous stress
  *********************************************************************/
  for(i=N1S-1;i<=N1E+1;i++) {
    for(j=N2S-1;j<=N2E+1;j++) { 
      for(k=N3S-1;k<=N3E+1;k++)  {

	f00 = &(visc_funcs[n][i][j][k]);
	ucov = f00->ucov;
	ucon = f00->ucon;
	eta  = f00->eta_visc;


	/* Calculate the 1st-order partial deriviatives  */
	for( mu = 0; mu < NDIM; mu++ ) { 
	  ucov2 = visc_funcs[n+shift[0][mu]][i+shift[1][mu]][j+shift[2][mu]][k+shift[3][mu]].ucov;
	  ucov1 = visc_funcs[n-shift[0][mu]][i-shift[1][mu]][j-shift[2][mu]][k-shift[3][mu]].ucov;
	  factor =  dfactors[mu];
	  for( nu = 0; nu < NDIM; nu++ ) { 
	    ducov[mu][nu] = factor * (ucov2[nu] - ucov1[nu]);
	  }
	}
	
	l = CONN_ID(i,j,k);
	get_geometry(i,j,k,CENT,n_mid,geom);

	/* Calculate the 1st-order covariant deriviatives  */
	for( mu = 0; mu < NDIM; mu++ ) { 
	  for( nu = 0; nu < NDIM; nu++ ) { 
	    divu = 0. ; 
	    for( lm = 0; lm < NDIM; lm++ ) { 
	      divu += conn[l][lm][mu][nu] * f00->ucov[lm];
	    }
	    ducov[mu][nu] -= divu;
	  }
	}

	for( mu = 0; mu < NDIM; mu++ ) { acov[mu] = 0.; } 
	for( mu = 0; mu < NDIM; mu++ ) for( nu = 0; nu < NDIM; nu++ ) { acov[mu] += ucon[nu] * ducov[nu][mu] ; }

	divu = 0.;
	for( mu = 0; mu < NDIM; mu++ ) { 
	  for( nu = 0; nu < NDIM; nu++ ) { 
	    divu += geom->gcon[mu][nu] * ducov[mu][nu];
	  }
	}

	divu *= two_thirds;

	for( mu = 0; mu < NDIM; mu++ ) { 
	  for( nu = 0; nu < NDIM; nu++ ) { 
	    stress[mu][nu] = eta * (
				    ducov[mu][nu] 
				    + ducov[nu][mu] 
				    + ucov[mu] * acov[nu] 
				    + ucov[nu] * acov[mu] 
				    - divu * ( geom->gcov[mu][nu] + ucov[mu]*ucov[nu] )
				    );
	  }
	}

	/* Raise stress and multiply by gdet: */
	eta = geom->g;
	for( mu = 0; mu < NDIM; mu++ ) { 
	  for( nu = 0; nu < NDIM; nu++ ) { 
	    divu = 0.;
	    for( lm = 0; lm < NDIM; lm++ ) { 
	      divu += geom->gcon[mu][lm] * stress[lm][nu] ; 
	    }
	    f00->stress[mu][nu] = divu * eta;
	  }
	}

      }
    }
  }


  /*********************************************************************
   Loop over all physical cells to calculate viscous SOURCE
        -- stress[][] is now   =   2 \sqrt{-g} eta sigma^\mu_\nu  
  *********************************************************************/
  LOOP { 

    f00 = &(visc_funcs[n][i][j][k]);

    /* acov_j =  \partial_i ( stress^i_j )   here :  */
    for( mu = 0; mu < NDIM; mu++ ) { acov[mu] = 0.; }

    for( mu = 0; mu < NDIM; mu++ ) { 
      fpp = &(visc_funcs[n+shift[0][mu]][i+shift[1][mu]][j+shift[2][mu]][k+shift[3][mu]]);
      fmm = &(visc_funcs[n-shift[0][mu]][i-shift[1][mu]][j-shift[2][mu]][k-shift[3][mu]]);
      factor =  dfactors[mu];
      for( nu = 0; nu < NDIM; nu++ ) { 
	acov[nu] += factor * ( fpp->stress[mu][nu] - fmm->stress[mu][nu] );
      }
    }

    l = CONN_ID(i,j,k);

    for( nu = 0; nu < NDIM; nu++ ) { 
      divu = 0.;
      for( mu = 0; mu < NDIM; mu++ ) { 
	for( lm = 0; lm < NDIM; lm++ ) { 
	  divu += conn[l][lm][mu][nu] * f00->stress[mu][lm] ; 
	}
      }
      visc_source[i][j][k][nu] = acov[nu] - divu;
    }

  }


#endif
  return;
}

/****************************************************************************************/
/*****************************************************************************************
  set_visc_stress_functions():
  --------------
     -- saves the 4-velocity and kinematic viscosity at the current time level to be used in the 
        future for finite differencing; 
*****************************************************************************************/
void  set_visc_stress_functions(int nn, int ii, int jj, int kk, double *prim, struct of_state *q )
{
#if( USE_KINEMATIC_VISCOSITY )
  struct of_viscosity *f00;
  
  f00 = &(visc_funcs[nn][ii][jj][kk]); 
  f00->ucon[0] = q->ucon[0];
  f00->ucon[1] = q->ucon[1];
  f00->ucon[2] = q->ucon[2];
  f00->ucon[3] = q->ucon[3];
  f00->ucov[0] = q->ucov[0];
  f00->ucov[1] = q->ucov[1];
  f00->ucov[2] = q->ucov[2];
  f00->ucov[3] = q->ucov[3];
  
  f00->eta_visc = prim[RHO] * nu_visc[ii]; 

#endif
  return;
}

/****************************************************************************************/
/*****************************************************************************************
  set_all_visc_stress_functions():
  --------------
     -- loops over all necessary cells and calculates the functions needed for 
        future for finite differencing; 
*****************************************************************************************/
void  set_all_visc_stress_functions(int nn, double ****prim)
{
#if( USE_KINEMATIC_VISCOSITY )

  int i, j, k;
  double nu_loc;
  struct of_viscosity *f00;
  struct of_geom   *geom;

#define BOUND_BODY(ii,jj,kk)    {			 \
    f00 = &(visc_funcs[nn][(ii)][(jj)][(kk)]);		 \
    get_geometry((ii),(jj),(kk),CENT,ncurr,geom);	 \
    ucon_calc(prim[(ii)][(jj)][(kk)], geom, f00->ucon) ; \
    lower(f00->ucon, geom, f00->ucov) ;                  \
    f00->eta_visc = prim[(ii)][(jj)][(kk)][RHO] * nu_loc; }
      

  /* nn is for  visc_funcs   and   ncurr is for geometry, they are not necessary equal */ 

  NV1_LOOP { 
    nu_loc = nu_visc[i]; 
    NV2_LOOP  NV3_LOOP {  BOUND_BODY(i,j,k);   }  
  }

#endif
  return;
}

/****************************************************************************************/
/*****************************************************************************************
  bound_visc_stress_functions():
  --------------
     -- sets the boundary values for those functions used to calculate the viscous source terms;
     -- since we use 2nd-order centered finite differencing, we only need one neighbor on each side, 
         or only need to set one ghost cell's worth of data per surface of the simulation volume;
*****************************************************************************************/
void  bound_visc_stress_functions(int nn, double ****prim)
{
#if( USE_KINEMATIC_VISCOSITY )
#define NGV_LOOP  for(l=1;l<=NGV;l++) 

  int i, j, k,l;
  double nu_loc;
  struct of_viscosity *f00;
  struct of_geom   *geom;

  /* nn is for  visc_funcs   and   ncurr is for geometry, they are not necessary equal */ 

  NGV_LOOP { 
    /* Lower X1-faces  */
    i = N1S-l; 
    nu_loc = nu_visc[i]; 
    NV2_LOOP NV3_LOOP {  BOUND_BODY(i,j,k);   }  

    /* Upper X1-faces  */
    i = N1E+l; 
    nu_loc = nu_visc[i]; 
    NV2_LOOP NV3_LOOP {  BOUND_BODY(i,j,k);   }  

    /* Lower X2-face  */
    j = N2S-l; 
    NV1_LOOP { 
      nu_loc = nu_visc[i]; 
      NV3_LOOP  {  BOUND_BODY(i,j,k);   }  
    }
  
    /* Upper X2-face  */
    j = N2E+l; 
    NV1_LOOP { 
      nu_loc = nu_visc[i]; 
      NV3_LOOP  {  BOUND_BODY(i,j,k);   }  
    }
  
    /* Lower X3-face  */
    k = N3S-l; 
    NV1_LOOP {
      nu_loc = nu_visc[i]; 
      NV2_LOOP {  BOUND_BODY(i,j,k);   }  
    }
  
    /* Upper X3-face  */
    k = N3E+l; 
    NV1_LOOP {
      nu_loc = nu_visc[i]; 
      NV2_LOOP {  BOUND_BODY(i,j,k);   }  
    }

  }  

#undef BOUND_BODY 
#undef NV1_LOOP 
#undef NV2_LOOP 
#undef NV3_LOOP 
#undef NGV
#endif
  return;
}

#undef VERBOSE
