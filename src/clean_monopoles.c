
#include <time.h>

#include <limits.h>

#include <hdf5.h>

#include "decs.h"
#include "metric.h"

double div_calc_CORN( double ****func, int i, int j, int k ) ;
double div_calc_CENT( double ****func, int i, int j, int k ) ;
void grad_calc_CORN ( double *gradf, double ***func, int i, int j, int k ) ;
void grad_calc_CENT ( double *gradf, double ***func, int i, int j, int k ) ;

void get_dpsi(double *ftmp, int i, int j, int k) ;
void get_func(double *ftmp, int i, int j, int k) ;

double get_fact_CORN( int i, int j, int k) ;
double get_fact_CENT( int i, int j, int k) ;

void SOR(const double omegaSOR, const double tol, const int itmax,
         double *divb,
         void (*grad_calc)( double *, double ***, int, int, int),
         double (*div_calc)( double ****, int, int, int ), 
         double (*get_fact)( int, int, int) ) ;

void bounds_func(double ****f) ;

double ****func;
double *divb_CORN;
double *divb_CENT;
double *res;
double ***psi;
double ****dpsi;


double norm2(double *v, int N)
{
  int i;
  double result[N], v2[N] ;
  double norm2 = 0;

  for(i=0; i<N; i++) {
    v2[i] = v[i]*v[i] ;
  }

  mpi_sum_to_head(v2, result, N) ;

  for(i=0; i<N; i++) {
    norm2 += result[i] ;
  }

  return sqrt(norm2) ;
}


/* for testing purposes only */
#if(0)
void get_fake_B(void)
{
  int i,j,k ;

  /* initialize random seed: */
  //  srand ( time(NULL) );
  srand(0);          // let us have it deterministic for now...
  double amp   = 1.e-6;
  //  double amp   = 0.;
  double pert ;

  ALL_LOOP {
    /* random number between 0 and amp */
    pert  = ((double) rand()) / RAND_MAX * amp ;

    p[i][j][k][B1] *= 1. + pert ;
    p[i][j][k][B2] *= 1. + pert ;
    p[i][j][k][B3] *= 1. + pert ;
  }
  //  fixup(p) ;
  bounds(p,0) ;

  return;
}
#endif


void get_initial_values(void)
{
  int i,j,k, ind;
  double ****btmp, det_g ;

  ALLOC_4D_ARRAY(btmp,N1TOT,N2TOT,N3TOT,3);

  struct of_geom  *geom;

  ALL_LOOP {
    get_geometry(i,j,k,CENT,ncurr,geom);  det_g = geom->g;

    /* initial guess for psi */
    psi[i][j][k]     = 0. ;
    //    psi[i][j][k]     =  ((double) rand()) / RAND_MAX * amp ; 

    btmp[i][j][k][0] = p[i][j][k][B1] * det_g ;
    btmp[i][j][k][1] = p[i][j][k][B2] * det_g ;
    btmp[i][j][k][2] = p[i][j][k][B3] * det_g ;
  }

  /* we store the RHS of the equation (divb) */
  ind = 0;
  LOOP {
    //    divb_CORN[ind]      = divb_calc( i, j, k ) ; // do not use this 
    divb_CORN[ind]      = div_calc_CORN( btmp, i, j, k ) ;
    divb_CENT[ind]      = div_calc_CENT( btmp, i, j, k ) ;
    ind++;
  }

  /* use the primitive variables array p to store (locally) psi and dpsi/dxi in
  the place for rho and Bi */
  LOOP {
    (grad_calc_CORN)( &dpsi[i][j][k][0], psi, i, j, k ) ;
    p[i][j][k][RHO] = psi[i][j][k] ;
    p[i][j][k][B1]  = dpsi[i][j][k][0] ;
    p[i][j][k][B2]  = dpsi[i][j][k][1] ;
    p[i][j][k][B3]  = dpsi[i][j][k][2] ;
  }
  bounds(p,0) ;

  /* filling in the temporary function, to take the divergence afterwards */
  ALL_LOOP {
    get_func(&func[i][j][k][0], i,j,k) ; 
  }

  DEALLOC_4D_ARRAY(btmp,N1TOT,N2TOT,N3TOT,3);

  return;
}



double clean_monopoles(double tol_CORN, int itmax_CORN)
{

  TRACE_BEG ;

  int i,j,k, ind;

  struct of_geom  *geom;
  double ***rho0;
  double ****B0;

  /* 
     our goal is to solve the equation laplacian psi = div B for psi. in this
     process, we need to impose physical boundary conditions in both psi and its
     first spatial derivatives. to help the process, we will, at least for now,
     use the primitive variables array p to store (locally) psi and dpsi/dxi in
     the place for rho and Bi, respectively. this way, when calling bounds, the
     boundary conditions for rho are imposed on psi and those for Bi are imposed
     on its derivatives. in the end, we copy back the original rho and Bi to the
     p array. note that we should not use the routine divb_calc (in diag.c)
     directly, since this returns the absolute value of the divergence.
   */

  ALLOC_4D_ARRAY(func,N1TOT,N2TOT,N3TOT,3);
  ALLOC_ARRAY(divb_CORN,NCELLS);
  ALLOC_ARRAY(divb_CENT,NCELLS);
  ALLOC_ARRAY(res,NCELLS);
  ALLOC_3D_ARRAY(psi,N1TOT,N2TOT,N3TOT);
  ALLOC_4D_ARRAY(dpsi,N1TOT,N2TOT,N3TOT,3);
  ALLOC_3D_ARRAY(rho0,N1TOT,N2TOT,N3TOT)   ;
  ALLOC_4D_ARRAY(B0  ,N1TOT,N2TOT,N3TOT,3) ;

  /* create some fictitious, divb != 0, B field. for testing purposes only.

  static unsigned int local_first_time = 1;
  if (local_first_time ) {  
    get_fake_B()         ;
    local_first_time = 0 ;
  }
  */

  ALL_LOOP {
    /* we copy the original rho and Bi to be able to restore later */
    rho0[i][j][k]    = p[i][j][k][RHO] ;
    B0 [i][j][k][0]  = p[i][j][k][B1] ;
    B0 [i][j][k][1]  = p[i][j][k][B2] ;
    B0 [i][j][k][2]  = p[i][j][k][B3] ;
  }

  /* initialize variables */
  get_initial_values() ;


  // it seems that using the cell-centred stencil does not improve things at all.
  /* 
  const double omegaSOR_CENT = 1.0 - 0.5 ;
  const double tol_CENT      = 1.e-9 ;
  const int itmax_CENT       = 10000 ;

  SOR(omegaSOR_CENT, tol_CENT, itmax_CENT,
      divb_CENT,
      &grad_calc_CENT, &div_calc_CENT, &get_fact_CENT ) ;
  */



  //  const double omegaSOR_CORN = 1./3. ; // stable?
  const double omegaSOR_CORN = 1./3 - 0.1 ; // stable



  SOR(omegaSOR_CORN, tol_CORN, itmax_CORN,
      divb_CORN,
      &grad_calc_CORN, &div_calc_CORN, &get_fact_CORN ) ;
  // it's not necessary to explicitly pass the address of the function, but anyway...



  /* we restore back the original primitive variables, adjust B to its new
     value, and move merrily along */
  double hdpsi[3] ;

  ind = 0;
  ALL_LOOP {
    get_geometry(i,j,k,CENT,ncurr,geom);

    p[i][j][k][RHO] = rho0[i][j][k]    ;

    get_dpsi(hdpsi, i,j,k) ;
    p[i][j][k][B1]  =   B0[i][j][k][0] - hdpsi[0] / geom->alpha;
    p[i][j][k][B2]  =   B0[i][j][k][1] - hdpsi[1] / geom->alpha;
    p[i][j][k][B3]  =   B0[i][j][k][2] - hdpsi[2] / geom->alpha;
    ind++ ;
  }
  //  bounds(p,0) ;


  ind = 0;
  LOOP {
    res[ind]      = divb_calc( i, j, k ) ; // here it's fine since we'll be squaring the value
    ind++ ;
  }

  double restot ;
  restot = norm2(res,NCELLS);
    //    sync_val(&restot) ;
  if( myid == printer_pid ) {
    fprintf(stdout, "final residual    = %28.18e \n", restot); fflush(stdout) ;
  }

  DEALLOC_4D_ARRAY(func,N1TOT,N2TOT,N3TOT,3);
  DEALLOC_ARRAY(divb_CORN,NCELLS);
  DEALLOC_ARRAY(divb_CENT,NCELLS);
  DEALLOC_ARRAY(res,NCELLS);
  DEALLOC_3D_ARRAY(psi,N1TOT,N2TOT,N3TOT);
  DEALLOC_4D_ARRAY(dpsi,N1TOT,N2TOT,N3TOT,3);
  DEALLOC_3D_ARRAY(rho0,N1TOT,N2TOT,N3TOT)   ;
  DEALLOC_4D_ARRAY(B0  ,N1TOT,N2TOT,N3TOT,3) ;

  TRACE_END ;

  return restot;
}


void SOR(const double omegaSOR, const double tol, const int itmax,
         double *divb,
         void (*grad_calc)( double *, double ***, int, int, int),
         double (*div_calc)( double ****, int, int, int ), 
         double (*get_fact)( int, int, int) ) 
{
  int i,j,k, ind, it;
  double fact, psi0, restot;



  ind = 0;
  LOOP {
    /* current residual */
    res[ind]      = (*div_calc)( func, i, j, k ) - divb[ind] ;
    ind++ ;
  }
  restot = norm2(res,NCELLS) ;
  sync_val(&restot) ;

  if( myid == printer_pid ) {
    fprintf(stdout, "initial residual  = %28.18e \n",  restot); fflush(stdout) ;
  }
  it = 0 ;
  while( restot > tol && it < itmax ) {

    restot = 0; 
    ind    = 0;
    LOOP {
      /* current residual */
      res[ind]        = (*div_calc)( func, i, j, k ) - divb[ind] ;

      //      fprintf(stdout, "res[%d] =  %28.18e  \n",ind, res[ind]  ); fflush(stdout) ;

      /* we get the correcting factor */
      fact            = (*get_fact)(i,j,k) ;

      psi0            = psi[i][j][k] ;
      /* we update the function psi to a more accurate value */
      psi[i][j][k]   += fact * res[ind] ;                                  // Gauss-Seidel relaxation
      psi[i][j][k]    = omegaSOR * psi[i][j][k] + (1. - omegaSOR) * psi0 ; // SOR relaxation

      /* we fill in the derivatives */
      (*grad_calc)( &dpsi[i][j][k][0], psi, i, j, k ) ;

      p[i][j][k][RHO] = psi[i][j][k] ;
      p[i][j][k][B1]  = dpsi[i][j][k][0] ;
      p[i][j][k][B2]  = dpsi[i][j][k][1] ;
      p[i][j][k][B3]  = dpsi[i][j][k][2] ;

      /* filling in the temporary function, to take the divergence afterwards */
      get_func(&func[i][j][k][0], i,j,k) ;

      ind++ ;
    }
    /* we apply the boundary conditions */
    bounds(p,0) ;
    /* need to fill ghostzones of the func array */
#if(0)
    ALL_LOOP {
      get_func(&func[i][j][k][0], i,j,k) ; 
    }
#else
    bounds_func(func) ;
#endif

    restot = norm2(res,NCELLS) ;
    sync_val(&restot) ;

    if( (myid == printer_pid ) && (it % 100 == 0) ) {
      fprintf(stdout, "it = %7d, res = %28.18e \n", it, restot); fflush(stdout) ;
    }

    it++ ;
  }


  if( myid == printer_pid ) {
    if( it >= itmax ) {
      fprintf(stdout, "reached maximum number of iterations, giving up.\nit = %7d, res = %28.18e \n", 
              it, restot); fflush(stdout) ;
    }
    else {
      fprintf(stdout, "it = %7d, res = %28.18e \n", it, restot); fflush(stdout) ;
    }
  }

  ALL_LOOP {
    psi[i][j][k]     = p[i][j][k][RHO] ;
    dpsi[i][j][k][0] = p[i][j][k][B1] ;
    dpsi[i][j][k][1] = p[i][j][k][B2] ;
    dpsi[i][j][k][2] = p[i][j][k][B3] ;
  }

  return ;
}


void bounds_func(double ****f)
{
  int i,j,k,g ;

  /* 1-lower face */
  GLOOP     N2_LOOP   N3_LOOP {
    get_func(&f[N1S-1-g][j][k][0], N1S-1-g,j,k) ; 
  }

  /* 1-upper face */
  GLOOP     N2_LOOP   N3_LOOP {
    get_func(&f[N1E+1+g][j][k][0], N1E+1+g,j,k) ; 
  }

  /* 2-lower face */
  N1ALL_LOOP      GLOOP    N3_LOOP  {
    get_func(&f[i][N2S-1-g][k][0], i,N2S-1-g,k) ; 
  }

  /* 2-upper face */
  N1ALL_LOOP      GLOOP    N3_LOOP  {
    get_func(&f[i][N2E+1+g][k][0], i,N2E+1+g,k) ; 
  }

  /* 3-lower face  */
  N1ALL_LOOP      N2ALL_LOOP     GLOOP    {
    get_func(&f[i][j][N3S-1-g][0], i,j,N3S-1-g) ; 
  }

  /* 3-upper face */
  N1ALL_LOOP      N2ALL_LOOP     GLOOP    {
    get_func(&f[i][j][N3E+1+g][0], i,j,N3E+1+g) ; 
  }

  return ;
}



/*
  the following function takes a CELL-centred 3-d vector and returns its
  CORNER-centred (at the ijk cell) divergence.
 */
double div_calc_CORN( double ****f, int i, int j, int k ) 
{
  unsigned int id;
  double div, ftmp1[8], ftmp2[8], ftmp3[8], ftmp;

  id = 0;
  ftmp1[id++] = -f[i-1][j-1][k-1][0];
  ftmp1[id++] = -f[i-1][j-1][k  ][0];
  ftmp1[id++] = -f[i-1][j  ][k-1][0];
  ftmp1[id++] = -f[i-1][j  ][k  ][0];
  ftmp1[id++] =  f[i  ][j-1][k-1][0];
  ftmp1[id++] =  f[i  ][j-1][k  ][0];
  ftmp1[id++] =  f[i  ][j  ][k-1][0];
  ftmp1[id  ] =  f[i  ][j  ][k  ][0];

  id = 0;
  ftmp2[id++] = -f[i-1][j-1][k-1][1];
  ftmp2[id++] = -f[i-1][j-1][k  ][1];
  ftmp2[id++] =  f[i-1][j  ][k-1][1];
  ftmp2[id++] =  f[i-1][j  ][k  ][1];
  ftmp2[id++] = -f[i  ][j-1][k-1][1];
  ftmp2[id++] = -f[i  ][j-1][k  ][1];
  ftmp2[id++] =  f[i  ][j  ][k-1][1];
  ftmp2[id  ] =  f[i  ][j  ][k  ][1];

  id = 0;
  ftmp3[id++] = -f[i-1][j-1][k-1][2];
  ftmp3[id++] =  f[i-1][j-1][k  ][2];
  ftmp3[id++] = -f[i-1][j  ][k-1][2];
  ftmp3[id++] =  f[i-1][j  ][k  ][2];
  ftmp3[id++] = -f[i  ][j-1][k-1][2];
  ftmp3[id++] =  f[i  ][j-1][k  ][2];
  ftmp3[id++] = -f[i  ][j  ][k-1][2];
  ftmp3[id  ] =  f[i  ][j  ][k  ][2];

  div = 0.;

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += ftmp1[id]; } 
  div += ftmp*invdx[1];

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += ftmp2[id]; } 
  div += ftmp*invdx[2];

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += ftmp3[id]; } 
  div += ftmp*invdx[3];

  div *= 0.25;

  /*
  if( i==45 && j==64 && k==3) {
    for(id=0;id<8;id++) {
      fprintf(stdout, "ftmp1[%d] = %28.18e \n",id,ftmp1[id]); fflush(stdout) ;
      fprintf(stdout, "ftmp2[%d] = %28.18e \n",id,ftmp2[id]); fflush(stdout) ;
      fprintf(stdout, "ftmp3[%d] = %28.18e \n",id,ftmp3[id]); fflush(stdout) ;
    }
  }
  */

  return div ;
}



/*
  the following function takes a CORNER-centred (scalar) function and returns
  its CELL-centred (at the ijk cell) first derivatives in all directions.
 */
void grad_calc_CORN( double *gradf, double ***f, int i, int j, int k ) 
{
  unsigned int id;
  double  ftmp1[8], ftmp2[8], ftmp3[8], ftmp;

  id = 0;
  ftmp1[id++] = -f[i  ][j  ][k  ];
  ftmp1[id++] = -f[i  ][j  ][k+1];
  ftmp1[id++] = -f[i  ][j+1][k  ];
  ftmp1[id++] = -f[i  ][j+1][k+1];
  ftmp1[id++] =  f[i+1][j  ][k  ];
  ftmp1[id++] =  f[i+1][j  ][k+1];
  ftmp1[id++] =  f[i+1][j+1][k  ];
  ftmp1[id  ] =  f[i+1][j+1][k+1];

  id = 0;
  ftmp2[id++] = -f[i  ][j  ][k  ];
  ftmp2[id++] = -f[i  ][j  ][k+1];
  ftmp2[id++] =  f[i  ][j+1][k  ];
  ftmp2[id++] =  f[i  ][j+1][k+1];
  ftmp2[id++] = -f[i+1][j  ][k  ];
  ftmp2[id++] = -f[i+1][j  ][k+1];
  ftmp2[id++] =  f[i+1][j+1][k  ];
  ftmp2[id  ] =  f[i+1][j+1][k+1];

  id = 0;
  ftmp3[id++] = -f[i  ][j  ][k  ];
  ftmp3[id++] =  f[i  ][j  ][k+1];
  ftmp3[id++] = -f[i  ][j+1][k  ];
  ftmp3[id++] =  f[i  ][j+1][k+1];
  ftmp3[id++] = -f[i+1][j  ][k  ];
  ftmp3[id++] =  f[i+1][j  ][k+1];
  ftmp3[id++] = -f[i+1][j+1][k  ];
  ftmp3[id  ] =  f[i+1][j+1][k+1];


  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += ftmp1[id]; } 
  gradf[0] = ftmp*invdx[1]*0.25;

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += ftmp2[id]; } 
  gradf[1] = ftmp*invdx[2]*0.25;

  ftmp = 0.;
  for(id=0;id<8;id++) { ftmp += ftmp3[id]; } 
  gradf[2] = ftmp*invdx[3]*0.25;


  return ;
}



/*
  the following ftion takes a CELL-centred (scalar) ftion and returns
  its CELL-centred (at the ijk cell) first derivatives in all directions.
 */
void grad_calc_CENT( double *gradf, double ***f, int i, int j, int k ) 
{

  gradf[0] = invdx[1]*0.5 * ( f[i+1][j  ][k  ] - f[i-1][j  ][k  ] ) ;
  gradf[1] = invdx[2]*0.5 * ( f[i  ][j+1][k  ] - f[i  ][j-1][k  ] ) ;
  gradf[2] = invdx[3]*0.5 * ( f[i  ][j  ][k+1] - f[i  ][j  ][k-1] ) ;

  return ;
}


/*
  the following function takes a CELL-centred 3-d vector and returns its
  CELL-centred (at the ijk cell) divergence.
 */
double div_calc_CENT( double ****f, int i, int j, int k ) 
{
  double div ;

  div  =  invdx[1]*0.5 * ( f[i+1][j  ][k  ][0] - f[i-1][j  ][k  ][0] )
        + invdx[2]*0.5 * ( f[i  ][j+1][k  ][1] - f[i  ][j-1][k  ][1] ) 
        + invdx[3]*0.5 * ( f[i  ][j  ][k+1][2] - f[i  ][j  ][k-1][2] ) ;

  return div;
}


void get_dpsi(double *ftmp, int i, int j, int k) 
{
  double psi_x, psi_y, psi_z ;
  double huxx, huxy, huxz, huyy, huyz, huzz; 

  struct of_geom  *geom;

  get_geometry(i,j,k,CENT,ncurr,geom);

  /* gamma^ij (contravariant 3-metric) */
  huxx = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;

  /* psi spacial derivatives. note that this assumes that we're storing psi's
     derivatives in the p[i][j][k][Bi] array. if this changes, then we must
     change here accordingly. */ 
  psi_x = p[i][j][k][B1] ;
  psi_y = p[i][j][k][B2] ;
  psi_z = p[i][j][k][B3] ;

  ftmp[0] =  huxx * psi_x + huxy * psi_y + huxz * psi_z  ;
  ftmp[1] =  huxy * psi_x + huyy * psi_y + huyz * psi_z  ;
  ftmp[2] =  huxz * psi_x + huyz * psi_y + huzz * psi_z  ;

  return ;
}


void get_func(double *ftmp, int i, int j, int k) 
{
  double hdet ;
  double hdpsi[3] ;

  struct of_geom  *geom;
  get_geometry(i,j,k,CENT,ncurr,geom);
  /* sqrt(gamma) */
  hdet = geom->g / geom->alpha ;

  get_dpsi(hdpsi, i,j,k) ; 

  ftmp[0] = hdet * hdpsi[0] ;
  ftmp[1] = hdet * hdpsi[1] ;
  ftmp[2] = hdet * hdpsi[2] ;

  return ;
}


double get_fact_CORN( int i, int j, int k) 
{
  unsigned int id;
  double huxx[8], huyy[8], huzz[8]; 
  double huxy[8], huxz[8], huyz[8]; 
  double hdet[8], ftmp[8], fact ;

  struct of_geom  *geom;

  get_geometry(i-1,j-1,k-1,CENT,ncurr,geom);   id = 0 ;
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i-1,j-1,k  ,CENT,ncurr,geom);   id++ ; // id = 1
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i-1,j  ,k-1,CENT,ncurr,geom);   id++ ; // id = 2
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i-1,j  ,k  ,CENT,ncurr,geom);   id++ ; // id = 3
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i  ,j-1,k-1,CENT,ncurr,geom);   id++ ; // id = 4
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i  ,j-1,k  ,CENT,ncurr,geom);   id++ ; // id = 5
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i  ,j  ,k-1,CENT,ncurr,geom);   id++ ; // id = 6
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i  ,j  ,k  ,CENT,ncurr,geom);   id++ ; // id = 7
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  huxy[id] = geom->gcon[XX][YY] + geom->beta[0] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huxz[id] = geom->gcon[XX][ZZ] + geom->beta[0] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  huyz[id] = geom->gcon[YY][ZZ] + geom->beta[1] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;


  for(id=0;id<8;id++) { ftmp[id] = huxx[id]*hdet[id]; } 

  fact  = 0.0625      * invdx[1]*invdx[1] * ( - ftmp[0] - ftmp[1]
                                              - ftmp[2] - ftmp[3]
                                              - ftmp[4] - ftmp[5]
                                              - ftmp[6] - ftmp[7] ) ;


  for(id=0;id<8;id++) { ftmp[id] = huxy[id]*hdet[id]; } 

  fact += 2. * 0.0625 * invdx[1]*invdx[2] * ( - ftmp[0] - ftmp[1]
                                              + ftmp[2] + ftmp[3]
                                              + ftmp[4] + ftmp[5]
                                              - ftmp[6] - ftmp[7] ) ;


  for(id=0;id<8;id++) { ftmp[id] = huxz[id]*hdet[id]; } 

  fact += 2. * 0.0625 * invdx[1]*invdx[3] * ( - ftmp[0] + ftmp[1]
                                              - ftmp[2] + ftmp[3]
                                              + ftmp[4] - ftmp[5]
                                              + ftmp[6] - ftmp[7] ) ;


  for(id=0;id<8;id++) { ftmp[id] = huyy[id]*hdet[id]; } 

  fact += 0.0625       * invdx[2]*invdx[2] * ( - ftmp[0] - ftmp[1]
                                               - ftmp[2] - ftmp[3]
                                               - ftmp[4] - ftmp[5]
                                               - ftmp[6] - ftmp[7] ) ;


  for(id=0;id<8;id++) { ftmp[id] = huyz[id]*hdet[id]; } 

  fact += 2. * 0.0625 * invdx[2]*invdx[3] * ( - ftmp[0] + ftmp[1]
                                              + ftmp[2] - ftmp[3]
                                              - ftmp[4] + ftmp[5]
                                              + ftmp[6] - ftmp[7] ) ;


  for(id=0;id<8;id++) { ftmp[id] = huzz[id]*hdet[id]; } 

  fact += 0.0625       * invdx[3]*invdx[3] * ( - ftmp[0] - ftmp[1]
                                               - ftmp[2] - ftmp[3]
                                               - ftmp[4] - ftmp[5]
                                               - ftmp[6] - ftmp[7] ) ;

  return -1./fact ;
}


double get_fact_CENT( int i, int j, int k) 
{
  unsigned int id;
  double huxx[2], huyy[2], huzz[2]; 
  double hdet[2], fact ;

  struct of_geom  *geom;

  get_geometry(i-1,j  ,k  ,CENT,ncurr,geom);   id = 0 ;
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i+1,j  ,k  ,CENT,ncurr,geom);   id++ ; // id = 1
  /* gamma^ij (contravariant 3-metric) */
  huxx[id] = geom->gcon[XX][XX] + geom->beta[0] * geom->beta[0]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  fact  = 0.25 * invdx[1]*invdx[1] * ( hdet[0] * huxx[0] + hdet[1] * huxx[1] ) ;


  get_geometry(i  ,j-1,k  ,CENT,ncurr,geom);   id = 0 ;
  /* gamma^ij (contravariant 3-metric) */
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i  ,j+1,k  ,CENT,ncurr,geom);   id++ ;
  /* gamma^ij (contravariant 3-metric) */
  huyy[id] = geom->gcon[YY][YY] + geom->beta[1] * geom->beta[1]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  fact  += 0.25 * invdx[2]*invdx[2] * ( hdet[0] * huyy[0] + hdet[1] * huyy[1] ) ;


  get_geometry(i  ,j  ,k-1,CENT,ncurr,geom);   id = 0 ;
  /* gamma^ij (contravariant 3-metric) */
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  get_geometry(i  ,j  ,k+1,CENT,ncurr,geom);   id++ ;
  /* gamma^ij (contravariant 3-metric) */
  huzz[id] = geom->gcon[ZZ][ZZ] + geom->beta[2] * geom->beta[2]/(geom->alpha * geom->alpha) ;
  /* sqrt(gamma) */
  hdet[id] = geom->g / geom->alpha ;

  fact  += 0.25 * invdx[3]*invdx[3] * ( hdet[0] * huzz[0] + hdet[1] * huzz[1] ) ;


  return 1./fact ;
}
