/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1dvsq2fix1.c: 
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS;
  

    Uses the 1D^*_{v^2} method: 
       -- solves for one independent variable (v^2, or vsq) via a 1D 
          Newton-Raphson method 
       -- like the 1D_{v^2} method, except can be used (in principle)
          with a general equation of state. The main difference is how
          it calculates the intermediary value of W: by performing a 
          nested set of Newton-Rapshon iterations to get the best W value
          for the current v^2 value.  

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#include "u2p_util.h"

#undef NDIM
#undef RHO
#undef UU 

#include "decs.h"

#define NEWT_DIM (1)

#define LTRACE (0)

/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D,alpha,K_atm, inv_gm1, gfactor;


// Declarations: 
static int Utoprim_new_body(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);

static void func_1d_gnr(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df);

static int general_newton_raphson( FTYPE x[], int n, 
				   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						  FTYPE [][NEWT_DIM], FTYPE *, 
						  FTYPE *) );

static FTYPE dWdvsq_calc(FTYPE vsq, FTYPE rho, FTYPE p);
static FTYPE W_of_vsq(FTYPE vsq, FTYPE *p, FTYPE *rho, FTYPE *u, FTYPE *sqrt_gtmp);

/**********************************************************************/
/******************************************************************

  Utoprim_1dvsq2fix1():
  
  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  

     It assumes that on input/output:

              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))  
              |  T^t_\mu           |  
              \   B^i              /


             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

     ala HARM. 

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

   -- different from Utoprim_1dvsq2() in that it also assumes an EOS of 
            P = K rho^\Gamma  or   P = K rho    
        depending on the value of G_ATM flag in u2p_defs.h


******************************************************************/
int Utoprim_1dvsq2fix1(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq, FTYPE K )
{

  FTYPE U_tmp[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE geomfactor,inv_gdet;


  if( U[0] <= 0. ) { 
    return(-100);
  }

  /* Set the geometry variables: */
  inv_gdet = geom->g_inv; 
  alpha = geom->alpha;
  geomfactor = alpha * inv_gdet;

  K_atm = K ; 
  inv_gm1 = 1./(GAMMA-1.);
  gfactor = GAMMA*(2.-G_ATM)*inv_gm1;

  /* First update the primitive B-fields */
  prim[BCON1] = U[BCON1] * inv_gdet ;
  prim[BCON2] = U[BCON2] * inv_gdet ;
  prim[BCON3] = U[BCON3] * inv_gdet ;

  
  /* Transform the CONSERVED variables into the new system */
  U_tmp[RHO]    = geomfactor * U[RHO] ;
  U_tmp[UU]     = geomfactor * (U[UU] - U[RHO]) ;
  U_tmp[UTCON1] = geomfactor * U[UTCON1] ;
  U_tmp[UTCON2] = geomfactor * U[UTCON2] ;
  U_tmp[UTCON3] = geomfactor * U[UTCON3] ;
  U_tmp[BCON1 ] = geomfactor * U[BCON1] ;
  U_tmp[BCON2 ] = geomfactor * U[BCON2] ;
  U_tmp[BCON3 ] = geomfactor * U[BCON3] ;

  
  /* Transform the PRIMITIVE variables into the new system */

  prim_tmp[RHO   ] = prim[RHO   ];
  prim_tmp[UU    ] = prim[UU    ];
  prim_tmp[UTCON1] = prim[UTCON1];
  prim_tmp[UTCON2] = prim[UTCON2];
  prim_tmp[UTCON3] = prim[UTCON3];

  prim_tmp[BCON1] =  U_tmp[BCON1 ] ;
  prim_tmp[BCON2] =  U_tmp[BCON2 ] ;
  prim_tmp[BCON3] =  U_tmp[BCON3 ] ;

  ret = Utoprim_new_body(U_tmp, geom, prim_tmp, gamma, bsq);

  /* Use new primitive variables if there was no problem : */ 
  if( ret == 0 ) {
    prim[RHO   ] = prim_tmp[RHO   ];
    prim[UU    ] = prim_tmp[UU    ];
    prim[UTCON1] = prim_tmp[UTCON1];
    prim[UTCON2] = prim_tmp[UTCON2];
    prim[UTCON3] = prim_tmp[UTCON3];
  }

  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu    |
             \  alpha B^i        /



               /    rho        \
	prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/
static int Utoprim_new_body(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma_out, FTYPE *bsq)
{
  FTYPE x_1d[1];
  FTYPE QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],Qsq,Qtcon[NDIM];
  FTYPE rho0,u,p,w,gammasq,gamma,sqrt_gtmp,W_last,W,utsq,vsq;
  FTYPE g_o_WBsq, QdB_o_W;
  int i,j, n, retval, i_increase ;



  // Assume ok initially:
  retval = 0;
  
  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  Bcon[1] = U[BCON1];
  Bcon[2] = U[BCON2];
  Bcon[3] = U[BCON3];
  //-faster   lower(Bcon,geom,Bcov);
  //-faster   Bsq = Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3] ;
  Bsq =   geom->gcov[1][1]*Bcon[1]*Bcon[1]
        + geom->gcov[2][2]*Bcon[2]*Bcon[2]
        + geom->gcov[3][3]*Bcon[3]*Bcon[3]
     + 2*(geom->gcov[1][2]*Bcon[1]*Bcon[2]
	+ geom->gcov[1][3]*Bcon[1]*Bcon[3]
	+ geom->gcov[2][3]*Bcon[2]*Bcon[3]) ;


  Qcov[0] = U[QCOV0] ;  
  Qcov[1] = U[QCOV1] ;  
  Qcov[2] = U[QCOV2] ;  
  Qcov[3] = U[QCOV3] ;  
  raise(Qcov,geom,Qcon) ;

  QdotB = Qcov[1]*Bcon[1] + Qcov[2]*Bcon[2] + Qcov[3]*Bcon[3] ;
  QdotBsq = QdotB*QdotB ;

  //-fast  ncov_calc(gcon,ncov) ;
//-faster  ncov[0] = -geom->alpha;
//-faster  ncov[1] = ncov[2] = ncov[3] = 0.;

  //-fast   raise_g(ncov,gcon,ncon);
//-faster   ncon[0] = -alpha * geom->gcon[0][0]; 
//-faster   ncon[1] = -alpha * geom->gcon[0][1]; 
//-faster   ncon[2] = -alpha * geom->gcon[0][2]; 
//-faster   ncon[3] = -alpha * geom->gcon[0][3]; 
  
  //-faster  Qdotn = Qcon[0]*ncov[0] ;
  Qdotn = -alpha*Qcon[0];

  //-fast Qsq = 0. ;
  //-fast  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;
  Qsq = DOT(Qcov,Qcon);

  Qtsq = Qsq + Qdotn*Qdotn ;

  D = U[RHO] ;

  /* calculate W from last timestep and use for guess */
  //-fast   utsq = 0. ;
  //-fast   for(i=1;i<4;i++)
  //-fast     for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;

  utsq =  geom->gcov[1][1]*prim[U1]*prim[U1]
        + geom->gcov[2][2]*prim[U2]*prim[U2]
        + geom->gcov[3][3]*prim[U3]*prim[U3]
     + 2*(geom->gcov[1][2]*prim[U1]*prim[U2]
	+ geom->gcov[1][3]*prim[U1]*prim[U3]
	+ geom->gcov[2][3]*prim[U2]*prim[U3]) ;

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  W_last = w*gammasq ;


  // Initialize independent variables for Newton-Raphson:
  x_1d[0] = utsq / gammasq;


  // Find vsq via Newton-Raphson:
  retval = general_newton_raphson( x_1d, 1, func_1d_gnr) ; 

  /* Problem with solver, so return denoting error before doing anything further */
  if( retval != 0 ) { 
    retval = retval*100+1;
    return(retval);
  }

  // Calculate v^2 :
  vsq = x_1d[0];
  if( (vsq >= 1.) || (vsq < 0.) ) {
    retval = 4;
    return(retval) ;
  }

  // Find W from this vsq:
  W = W_of_vsq(vsq, &p, &rho0, &u, &sqrt_gtmp);


  // Recover the primitive variables from the scalars and conserved variables:
  gamma = 1./sqrt_gtmp;
  *gamma_out = gamma;

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( treat_floor_as_failure && ((rho0 <= 0.) || (u <= 0.)) ) { 
    retval = 5;
    return(retval) ;
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;


  //-fast  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;

  //-fast  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;
  g_o_WBsq = gamma/(W+Bsq);
  QdB_o_W  = QdotB / W; 
  *bsq = Bsq * (1.-vsq) + QdB_o_W*QdB_o_W;
  prim[UTCON1] = g_o_WBsq * ( Qcon[1] + geom->ncon[1] * Qdotn + QdB_o_W*Bcon[1] ) ;
  prim[UTCON2] = g_o_WBsq * ( Qcon[2] + geom->ncon[2] * Qdotn + QdB_o_W*Bcon[2] ) ;
  prim[UTCON3] = g_o_WBsq * ( Qcon[3] + geom->ncon[3] * Qdotn + QdB_o_W*Bcon[3] ) ;
	
  /* set field components */
  //-fast  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* done! */
  return(retval) ;
    
}
  
/**********************************************************************/
/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0] is physical, based upon its definition;
    -- "corrects" it if it is not physical;

*********************************************************************/

static void validate_x(FTYPE x[1], FTYPE x0[1] ) 
{
  
  FTYPE small = 1.e-10;

  x[0] = (x[0] >= 1.0)    ?  ( 0.5*(x0[0] + 1.) )    : x[0];
  x[0] = (x[0] <  -small) ?  ( 0.5*x0[0] )           : x[0];
  x[0] = fabs(x[0]);

  return;

}


/**********************************************************************/
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, 
				   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						  FTYPE [][NEWT_DIM], FTYPE *, 
						  FTYPE *) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old, rho,p,u;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df =  f = 1.;
  i_extra = doing_extra = 0;

  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  vsq_old = vsq = W = W_old = 0.;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    // METHOD specific
    validate_x( x, x_old );


    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific
    W_old = W;
    W = W_of_vsq( x[0], &p, &rho, &u, &errx);
    errx  = (W==0.) ?  fabs(W-W_old) : fabs((W-W_old)/W);
    errx += (x[0]==0.) ?  fabs(x[0]-x_old[0]) : fabs((x[0]-x_old[0])/x[0]);


//    fprintf(stderr,"NR: %26.20e  %26.20e  %26.20e  %26.20e  %26.20e  \n", x_old[0], x[0], resid[0], jac[0][0], errx);
//    fflush(stderr);

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*   before stopping                                                         */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    // See if we've done the extra iterations, or have done too many iterations:
    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
	|| (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) ) {
    return(2);
  }

  // Return in different ways depending on whether a solution was found:
  if( fabs(errx) > MIN_NEWT_TOL){
#if(LTRACE)
    fprintf(stderr," totalcount = %d   0   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr);                      
#endif
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    //fprintf(stderr," totalcount = %d   1   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr);
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    //fprintf(stderr," totalcount = %d   2   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr); 
    return(0);
  }

  return(0);

}



/********************************************************************************/
/********************************************************************** 
   func_1d_gnr(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/

static void func_1d_gnr(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df)
{
  FTYPE vsq,W,Wsq,W3,dWdvsq , fact_tmp, rho, p, u  ;
  int retval, iters; 


  vsq = x[0];

  // Calculate best value for W given current guess for vsq: 
  W = W_of_vsq(vsq, &p, &rho, &u, &Wsq);  /* Wsq is used here as a temp var */
  Wsq = W*W;
  W3 = W*Wsq;

  // Doing this assuming  P = (G-1) u :

  dWdvsq = dWdvsq_calc(vsq, rho, p);

  fact_tmp = (Bsq + W) ;

  resid[0] = Qtsq  -  vsq * fact_tmp * fact_tmp  +  QdotBsq * ( Bsq + 2.*W ) / Wsq ; 
  jac[0][0] =  -fact_tmp * ( fact_tmp +   2. * dWdvsq * ( vsq + QdotBsq/W3 ) ) ; 

  dx[0] = -resid[0]/jac[0][0];

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

}
/********************************************************************** 
 ********************************************************************** 
   
 The following routines specify the equation of state.  All routines 
  above here should be indpendent of EOS.  If the user wishes 
  to use another equation of state, the below functions must be replaced 
  by equivalent routines based upon the new EOS. 

 **********************************************************************
**********************************************************************/

/* 
W as a function of v^2
*/
static FTYPE W_of_vsq(FTYPE vsq, FTYPE *p, FTYPE *rho, FTYPE *u, FTYPE *sqrt_gtmp)
{
  FTYPE gtmp;
  gtmp = (1. - vsq);
  *sqrt_gtmp = sqrt(gtmp);
  *rho = D * (*sqrt_gtmp);

#if( USE_ISENTROPIC ) 
  *p = K_atm * pow( *rho, G_ATM );
#else
  *p = K_atm * (*rho);
#endif

  *u = *p * inv_gm1; 
  
  return( (*rho + *u + *p ) / (gtmp)  );

}

/* 
dW/dvsq as a function of v^2, rho, p
*/
static FTYPE dWdvsq_calc(FTYPE vsq, FTYPE rho, FTYPE p)
{

  return(  ( gfactor*p   + rho ) / ( 2.*(1.-vsq)*(1.-vsq) )   ) ;
  
}


/****************************************************************************** 
             END   OF   UTOPRIM_1DVSQ2FIX.C
 ******************************************************************************/
