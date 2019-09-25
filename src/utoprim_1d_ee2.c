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

utoprim_1d_ee2.c: 
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS, which 
        is  

             P = Sc rho^(GAMMA-1) / gamma
  
    Uses a method similiar to  1D_W method: 
       -- solves for one independent variable (rho) via a 1D
          Newton-Raphson method 
       -- by substituting 
          W = Dc ( Dc + GAMMA Sc rho^(GAMMA-1) / (GAMMA-1) ) / rho 
           into Qtsq equation, one can get one equation for 
           one unknown (rho)

       -- can be used (in principle) with a general equation of state. 

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/


#include "u2p_util.h"

#define NEWT_DIM 1

#define LTRACE 0

#undef NDIM
#undef RHO
#undef UU 

#include "decs.h"


/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D, Sc, rho_gm1;
FTYPE W_for_gnr2, rho_for_gnr2, W_for_gnr2_old, rho_for_gnr2_old, drho_dW;
FTYPE g_o_gm1;


// Declarations: 
static int Utoprim_new_body(FTYPE U[], struct of_geom *geom, FTYPE prim[], FTYPE *gamma);

static int general_newton_raphson( FTYPE x[], int n, 
				   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						  FTYPE [][NEWT_DIM], FTYPE *, 
						  FTYPE *, int) );

static void func_rho(FTYPE x[], FTYPE dx[], FTYPE resid[], 
		     FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);

/**********************************************************************/
/******************************************************************

  Utoprim_1d_ee2():

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


         S = sqrt(-det(g_{a b})) u^t P / rho^(GAMMA-1)

     ala HARM. 

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
         S       = entropy density  = sqrt(-det(g_{a b})) u^t P / rho^(GAMMA-1)
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

int Utoprim_1d_ee2(FTYPE U[NPR], FTYPE S, struct of_geom *geom, FTYPE prim[NPR], FTYPE *gamma )
{
 
  FTYPE U_tmp[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE alpha,geomfactor,inv_gdet;


  if( U[0] <= 0. ) { 
    return(-100);
  }

  g_o_gm1 = GAMMA/(GAMMA-1.);
  inv_gdet = 1./geom->g; 

  /* First update the primitive B-fields */
  prim[BCON1] = U[BCON1] * inv_gdet ;
  prim[BCON2] = U[BCON2] * inv_gdet ;
  prim[BCON3] = U[BCON3] * inv_gdet ;

  /* Set the geometry variables: */
  alpha = geom->alpha;
  geomfactor = alpha*inv_gdet;
  
  /* Transform the CONSERVED variables into the new system */
  Sc = geomfactor * S; 
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

  ret = Utoprim_new_body(U_tmp, geom, prim_tmp, gamma);

  /* Transform new primitive variables back if there was no problem : */ 
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

static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *geom, FTYPE prim[NPR], FTYPE *gamma_out)
{

  FTYPE x_1d[1];
  FTYPE QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qsq,Qtcon[NDIM];
  FTYPE rho0,u,p,w,gammasq,gamma,gamma_sq,W_last,W,utsq,tmpdiff ;
  FTYPE g_o_WBsq, QdB_o_W;
  int i,j, retval, i_increase,ntries ;

  // Assume ok initially:
  retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  Bcon[1] = U[BCON1];
  Bcon[2] = U[BCON2];
  Bcon[3] = U[BCON3];
  lower(Bcon,geom,Bcov);

  Qcov[0] = U[QCOV0] ;  
  Qcov[1] = U[QCOV1] ;  
  Qcov[2] = U[QCOV2] ;  
  Qcov[3] = U[QCOV3] ;  
  raise(Qcov,geom,Qcon) ;

  Bsq = Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3] ;

  QdotB = Qcov[1]*Bcon[1] + Qcov[2]*Bcon[2] + Qcov[3]*Bcon[3] ;
  QdotBsq = QdotB*QdotB ;

  ncov[0] = -geom->alpha;
  ncov[1] = ncov[2] = ncov[3] = 0.;
  
  ncon[0] = -geom->alpha * geom->gcon[0][0]; 
  ncon[1] = -geom->alpha * geom->gcon[0][1]; 
  ncon[2] = -geom->alpha * geom->gcon[0][2]; 
  ncon[3] = -geom->alpha * geom->gcon[0][3]; 

  Qdotn = Qcon[0]*ncov[0] ;

  Qsq = DOT(Qcov,Qcon);

  Qtsq = Qsq + Qdotn*Qdotn ;

  D = U[RHO] ;

  /* calculate W from last timestep and use for guess */
//  utsq = 0. ;
//  for(i=1;i<4;i++)
//    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;

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
	
  rho0 = D / gamma ;
  rho_gm1 = pow(rho0,(GAMMA-1.));
  p = Sc * rho_gm1 / gamma; 
  u = p / (GAMMA-1.);
  w = rho0 + u + p ;


#if(LTRACE)
  if( rho0 <= 0. ) { 
    fprintf(stderr,"beg neg rho fix1, rho,D,gamma = %26.20e %26.20e %26.20e \n", rho0, D, gamma); fflush(stderr);
    rho0 = fabs(rho0);
  }
#endif 
    
  x_1d[0] = rho0;
  retval = general_newton_raphson( x_1d, 1, func_rho);
  rho0 = x_1d[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (rho0 == FAIL_VAL) ) {
    retval = retval*100+1;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, rho = %d %26.20e \n", retval,rho0);fflush(stderr);
#endif
    return(retval);
  }
  else{
    if(rho0 > W_TOO_BIG) {
      retval = 3;
#if(LTRACE)
      fprintf(stderr,"fix1: retval, rho = %d %26.20e \n", retval,rho0);fflush(stderr);
#endif
      return(retval) ;
    }
  }

  // Calculate v^2 : 
  utsq = (D-rho0)*(D+rho0)/(rho0*rho0);
  gamma_sq = 1.+utsq;
  gamma = sqrt(gamma_sq);
  *gamma_out = gamma;
  
  if( utsq < 0. ) {
    retval = 4;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, rho, utsq = %d %26.20e %26.20e \n", retval,rho0,utsq);fflush(stderr);
#endif
    return(retval) ;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  rho_gm1 = pow(rho0,(GAMMA-1.));
  p = Sc * rho_gm1 / gamma; 
  u = p / (GAMMA-1.);
  w = rho0 + u + p ;
  W = w * gamma_sq;


  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( treat_floor_as_failure && ((rho0 <= 0.) || (u <= 0.)) ) { 
    retval = 5;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho,utsq,u = %d %26.20e %26.20e %26.20e %26.20e \n", retval,W, rho0,utsq,u);fflush(stderr);
#endif
    return(retval) ;
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;

  g_o_WBsq = gamma/(W+Bsq);
  QdB_o_W  = QdotB / W; 
  prim[UTCON1] = g_o_WBsq * ( Qcon[1] + ncon[1] * Qdotn + QdB_o_W*Bcon[1] ) ;
  prim[UTCON2] = g_o_WBsq * ( Qcon[2] + ncon[2] * Qdotn + QdB_o_W*Bcon[2] ) ;
  prim[UTCON3] = g_o_WBsq * ( Qcon[3] + ncon[3] * Qdotn + QdB_o_W*Bcon[3] ) ;


  /* done! */
  return(retval) ;

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
						  FTYPE *, int) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,W,W_old;

  int   keep_iterating, i_increase;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  W = W_old = 0.;

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

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = fabs(x[0]);


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) || (finite(x[0])==0)  ) {
#if(LTRACE)
    fprintf(stderr,"\ngnr not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
	    f,df,x[0],x_old[0],W_for_gnr2_old,W_for_gnr2,rho_for_gnr2_old,rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }


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

/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho Qtsq equation with 
            the definition of W 
        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho
              substituted in. 

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void func_rho(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			 FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{

  FTYPE t1,t2;
  FTYPE t12;
  FTYPE t17;
  FTYPE t3, t100,rhosq,t200,rho,W,dWdrho,dvsqdrho,vsq;

  rho = x[0]; 
  rhosq = rho*rho;
  t200 = 1./(D*D);
  t100 = g_o_gm1*Sc*pow(rho,(GAMMA-1.));
  W = D*( D + t100 ) / rho ; 
  dWdrho = D * ( -D + t100*(GAMMA-2.) ) / rhosq; 
  t1 = W*W;
  t2 = Bsq+W;
  //    t3 = pow(Bsq+W,2.0);
  t3 = t2*t2;
  vsq = (D-rho)*(D+rho)*t200;
  dvsqdrho = -2*rho*t200;
  resid[0] = t1*(Qtsq-vsq*t3)+QdotBsq*(t2+W);
  t12 = Bsq*Bsq;
  t17 = dWdrho*vsq;
  jac[0][0] = 2*QdotBsq*dWdrho
    +((Qtsq-vsq*t12)*2*dWdrho+(-6*t17*Bsq-dvsqdrho*t12
			       +(-2*dvsqdrho*Bsq-4*t17-dvsqdrho*W)*W)*W)*W;

  dx[0] = -resid[0]/jac[0][0];
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

  return;

}

/****************************************************************************** 
             END   OF   UTOPRIM_1D_EE2.C
 ******************************************************************************/

