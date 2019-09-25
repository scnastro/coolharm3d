/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************

  This file contains the functions define the functional relationship 
   between the "x" coordinates (e.g. r,th,phi)  to the "xp" or "x^\prime"  
   coordinates (e.g. x1,x2,x3).  One needs to specify : 
      x(xp) :   x_to_(, x'(xp), x''(

****************************************************************************************
****************************************************************************************/

#include "decs.h"
#include "metric.h"

static double norm_diag2,beg_diag2,end_diag2,len_diag2;
static double *ai_diag2,*bi_diag2;

static double x2_of_xp2_diag2(double xp2);
static double dx2_dxp2_diag2(double xp2);
static double dx2_dxp2_dxp2_diag2(double xp2);
static void calc_ai_bi_di_diag2(void);
static void calc_ai_bi_diag2(void);
static double theoretical_norm_diag2(void);
static void int_step(double x,double x0,double s,double *intstep, double *intstepx);
static double seesaw(double x, double x0, double xn);

static double stepfunc(double x, double x1, double hh);
static double square(double x, double x1, double h, double delta);
static double square_per(double x, double x1, double h, double delta);
static double stepfunc_int(double x, double x1, double h);
static double square_int(double x, double x1, double h, double delta);
static double square_int_per(double x, double x1, double h, double delta);
static double dsquare(double x, double x1, double h, double delta) ;
static double dsquare_per(double x, double x1, double h, double delta) ;
static double d2square(double x, double x1, double h, double delta) ;
static double d2square_per( double x, double x1, double h, double delta) ;

static void rfuncs(double *rf, double *drdy, double *d2rdy2, double *ar_i,
                   double y, double ybh_i, double Rout_i, double Rin_i, double br_i, double s_i) ;
static double sinh_terms(double s, double y, double y0 );

static void dxp_dx_calc_default( double *x, double *xp, double dxp_dx[][NDIM]);
static void dxp_dx_calc2_default( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM]);
static void dxp_dx_calc_dyn_rad( double *x, double *xp, double dxp_dx[][NDIM]);
static void dxp_dx_calc2_dyn_rad( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM]);
static void dxp_dx_calc_general( double *x, double *xp, double dxp_dx[][NDIM]);
static void dxp_dx_calc2_general( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM]);


#define NEWT_DIM 1
 /* your choice of floating-point data type */
#define FTYPE double 

static void newt_raphs_func(FTYPE x[], FTYPE dx[], FTYPE resid[], 
                            FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE df[], int n, 
                            double Rout_i, double Rin_i, double br_i, double s_i, double rbh_i) ;



#define MAX_NEWT_ITER 100     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-15    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER 2       /* Keep this at ZERO in order to prevent the solver from climbing back out of the steep residual hill, which often happens for this equation when trying extra iterations */

// from utoprim_1d.c, with adaptations
static int newt_raphs_warped_spherical( FTYPE x[], FTYPE dfdx[], int n, 
                                        double Rout_i, double Rin_i, double br_i, double s_i, double rbh_i,
                                        void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
                                                       FTYPE [][NEWT_DIM], FTYPE *, 
                                                       FTYPE [], int, 
                                                       double Rout_i, double Rin_i, double br_i, double s_i, double rbh_i)
                                        ) ;


#if(COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD)
# define dxp_dx_calc            dxp_dx_calc_dyn_rad
# define dxp_dx_calc2          dxp_dx_calc2_dyn_rad
# define det_dx_dxp_calc    det_dx_dxp_calc_default
# define det_dx_dxp_calc2  det_dx_dxp_calc2_default
#elif(COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL||COORD_TYPE_CHOICE==COORD_WARPED_CARTESIAN)
# define dxp_dx_calc            dxp_dx_calc_general
# define dxp_dx_calc2          dxp_dx_calc2_general
# define det_dx_dxp_calc    det_dx_dxp_calc_general
# define det_dx_dxp_calc2  det_dx_dxp_calc2_general
#else 
# define dxp_dx_calc            dxp_dx_calc_default
# define dxp_dx_calc2          dxp_dx_calc2_default
# define det_dx_dxp_calc    det_dx_dxp_calc_default
# define det_dx_dxp_calc2  det_dx_dxp_calc2_default
#endif

#define sech2(x)  (1./( (cosh((x))) * (cosh((x))) ))

// we define shorthand notation for some parameters in the
// of_coord_warped_spherical structure
#if( COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL )
# define ar_ coord_params->ar_
# define rbh_ coord_params->rbh_
# define delta_z1 coord_params->delta_z1 
# define delta_z2 coord_params->delta_z2 
# define a_z10    coord_params->a_z10 
# define h_z1     coord_params->h_z1 
# define h_z2     coord_params->h_z2 
# define s_      coord_params->s_
# define br_     coord_params->br_
# define Rin_    coord_params->Rin_
# define Rout_   coord_params->Rout_
#endif

#if( COORD_TYPE_CHOICE==COORD_WARPED_CARTESIAN ) 
#define xmax coord_params->xmax
#define xmin coord_params->xmin
#define ymax coord_params->ymax
#define ymin coord_params->ymin
#endif

#if( (COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL) || (COORD_TYPE_CHOICE==COORD_WARPED_CARTESIAN) )
# define ybh_ coord_params->ybh_
# define xbh_ coord_params->xbh_
# define zbh_ coord_params->zbh_
# define dxbh1_dt coord_params->dxbh1_dt 
# define dxbh2_dt coord_params->dxbh2_dt
# define dybh1_dt coord_params->dybh1_dt
# define dybh2_dt coord_params->dybh2_dt
# define dybh3_dt coord_params->dybh3_dt
# define d2xbh1_dt2 coord_params->d2xbh1_dt2
# define d2xbh2_dt2 coord_params->d2xbh2_dt2 
# define d2ybh1_dt2 coord_params->d2ybh1_dt2 
# define d2ybh2_dt2 coord_params->d2ybh2_dt2 
# define d2ybh3_dt2 coord_params->d2ybh3_dt2 
# define delta_x1 coord_params->delta_x1 
# define delta_x2 coord_params->delta_x2 
# define delta_x3 coord_params->delta_x3 
# define delta_x4 coord_params->delta_x4 
# define delta_y1 coord_params->delta_y1 
# define delta_y2 coord_params->delta_y2 
# define delta_y3 coord_params->delta_y3 
# define delta_y4 coord_params->delta_y4 
# define a_x10    coord_params->a_x10 
# define a_x20    coord_params->a_x20 
# define a_y10    coord_params->a_y10 
# define a_y20    coord_params->a_y20 
# define h_x1     coord_params->h_x1 
# define h_x2     coord_params->h_x2  
# define h_x3     coord_params->h_x3 
# define h_x4     coord_params->h_x4 
# define h_y1     coord_params->h_y1   
# define h_y2     coord_params->h_y2 
# define h_y3     coord_params->h_y3 
# define h_y4     coord_params->h_y4 
#endif

/***************************************************************************/
/***************************************************************************
    coord():
    -------
       -- given the indices i,j,k and location in the cell, return with 
          the values of X1,X2,X3 there;  
       -- the locations are defined by : 

            -----------------------       
            |                     |   X2 ^    
            |                     |      |       
            |FACE1   CENT         |      |    
            |                     |      |    
            |CORN    FACE2        |      ---> 
            ----------------------          X1
   

            -----------------------   
            |                     |   X2 ^    
            |                     |      |    
            |FACE3   CENT         |      |    
            |                     |      |    
            |CORN    FACE2        |      ---> 
            ----------------------          X3
        

***************************************************************************/
void coord(int i, int j, int k, int loc, double *xp)
{
  int ig,jg,kg;

  xp[0] = t; 

  get_global_ijk(i,j,k,&ig,&jg,&kg);

  switch( loc ) { 
    
  case FACE1 : 
    xp[1] = startx[1] +  ig       *dx[1] ;
    xp[2] = startx[2] + (jg + 0.5)*dx[2] ;
    xp[3] = startx[3] + (kg + 0.5)*dx[3] ;
    return ;

  case FACE2 : 
    xp[1] = startx[1] + (ig + 0.5)*dx[1] ;
    xp[2] = startx[2] +  jg       *dx[2] ;
    xp[3] = startx[3] + (kg + 0.5)*dx[3] ;
    return ;
    
  case FACE3 : 
    xp[1] = startx[1] + (ig + 0.5)*dx[1] ;
    xp[2] = startx[2] + (jg + 0.5)*dx[2] ;
    xp[3] = startx[3] +  kg       *dx[3] ;
    return ;

  case CENT :
    xp[1] = startx[1] + (ig + 0.5)*dx[1] ;
    xp[2] = startx[2] + (jg + 0.5)*dx[2] ;
    xp[3] = startx[3] + (kg + 0.5)*dx[3] ;
    return ;

  case CORN :
    xp[1] = startx[1] + ig*dx[1] ;
    xp[2] = startx[2] + jg*dx[2] ;
    xp[3] = startx[3] + kg*dx[3] ;
    return ;

  default : 
    fprintf(stderr,"coord(): Bad location value :  %d \n", loc); 
    fflush(stderr);
    fail(FAIL_BASIC,0);
    return ;
  }
}

/****************************************************************************************

 x_of_xp():
 ----------
   -- finds "normal" coordinates from numerical coordinates;
   -- x  = e.g.  "t,r,th,phi", "t,x,y,z"
   -- xp = X[0-3]
   -- this is a little counter intuitive since we use X[] for xp all the tme;

****************************************************************************************/
void x_of_xp( double *x, double *xp ) 
{
  double ftmp;
  double xp2;
 
#if(   COORD_TYPE_CHOICE  == COORD_DIAGONAL  ) 
  xp2 = seesaw(xp[2],0.,1.);
  x[TT]  = xp[0];   
  x[RR]  = R0 + exp(xp[1]) ;
  //-newcoords  x[TH]  = M_PI * xp[2]  +  h_slope * sin(2 * M_PI * xp[2]); 
  //  x[TH]  = th_length * xp[2] + th_beg +  h_slope * sin(2 * M_PI * xp[2]); 
  x[TH]  = th_length * xp2 + th_beg +  h_slope * sin(2 * M_PI * xp2); 
  x[PH]  = xp[3];   


#elif( COORD_TYPE_CHOICE  == COORD_MIXED     ) 
  x[TT]  = xp[0];   
  x[RR]  = R0 + exp(xp[1]) ;
//-newcoords  x[TH]  = M_PI * xp[2] + h_slope * sin(2 * M_PI * xp[2]) 
//-newcoords                          * 4 * atan(X1_slope*(X1_0 - xp[1])) / M_PI  ;
  x[TH]  = th_length * xp[2] + th_beg + h_slope * sin(2 * M_PI * xp[2]) 
                                       * 4 * atan(X1_slope*(X1_0 - xp[1])) / M_PI  ;
  x[PH]  = xp[3];   

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL2 ) 
  x[TT]  = xp[0];   
  x[RR]  = R0 + exp(xp[1]) ;
  x[TH]  = (x2_of_xp2_diag2(xp[2]) - beg_diag2) * norm_diag2;
  x[PH]  = xp[3];   

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3 ) 
  //  ftmp = 2 * ( PERIODIC( (xp[2]) , (0.), (1.) ) )  - 1.;
  ftmp = 2 * seesaw(xp[2],0.,1.) - 1.;
  x[TT]  = xp[0];   
  x[RR]  = R0 + exp(xp[1]);
  x[TH]  = 0.5*M_PI*( 1. + h_slope*ftmp  + diag3_factor*pow(ftmp,diag3_exponent) );
  x[PH]  = xp[3];   

#elif( COORD_TYPE_CHOICE  == COORD_IDENTITY  ) 
  x[0] = xp[0] ; 
  x[1] = xp[1] ; 
  x[2] = xp[2] ; 
  x[3] = xp[3] ; 


#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3_DYN_RAD )
  ftmp = 2 * seesaw(xp[2],0.,1.) - 1.;
  x[TT]  = xp[0];   
  x[RR]  = coord_params->Rin_of_t * ( coord_params->f0_r  + (1. - coord_params->f0_r) * exp(xp[1]) ); 
  x[TH]  = 0.5*M_PI*( 1. + h_slope*ftmp  + diag3_factor*pow(ftmp,diag3_exponent) );
  x[PH]  = xp[3];   


#elif( COORD_TYPE_CHOICE  == COORD_WARPED_SPHERICAL )

  int n;
  double r_of_yt_[3] ;

  for(n=0; n<3; n++) {
//--orig    r_of_yt_[n] =  Rin_[n] + ( br_[n] - ar_[n]*s_[n] ) * xp[1]
//--orig                 + ar_[n] * (   sinh(s_[n] * (xp[1] - ybh_[n]))
//--orig                              - sinh(-s_[n] * ybh_[n]) ) ;

    r_of_yt_[n] =  Rin_[n] + ( br_[n] ) * xp[1]
      + ar_[n] * (  sinh_terms( s_[n], xp[1], ybh_[n] ) );

  }

  x[TT] = xp[0];   

  x[RR] = r_of_yt_[0] + square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                (r_of_yt_[1] - r_of_yt_[0]) * square_per(xp[3],xbh_[1],h_x3,delta_x3)
                                                              + (r_of_yt_[2] - r_of_yt_[0]) * square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                              ) ;

  x[TH] = M_PI * ( xp[2] - a_z10 * (   square_int_per(xp[2],zbh_[1],h_z1,delta_z1) 
                                     - square_int_per(zbh_[1],zbh_[1],h_z1,delta_z1) 
                                     - (xp[2]-zbh_[1]) * (  square_int_per(1,zbh_[1],h_z1,delta_z1) 
                                                          - square_int_per(0,zbh_[1],h_z1,delta_z1)) 
                                    )
                   );

  x[PH] = 2.*M_PI * (   xp[3] 
                      - square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                    a_x10 * square(xp[1],ybh_[1],h_y3,delta_y3) * (  square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                                                                                                   - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                                                                                                   - (xp[3]-xbh_[1])*(  square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                                                                                                      - square_int_per(0,xbh_[1],h_x1,delta_x1)) )
                                                                  + a_x20 * square(xp[1],ybh_[2],h_y4,delta_y4) * (  square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                                                                                                   - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                                                                                                   - (xp[3]-xbh_[2])*(  square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                                                                                                      - square_int_per(0,xbh_[2],h_x2,delta_x2)) )  
                                                                   )
                        ) ;


#elif( COORD_TYPE_CHOICE  == COORD_WARPED_CARTESIAN )

  x[TT] = xp[0] ; 

  x[XX] = xp[1] 
                 -  a_x10 * square_per(xp[2],ybh_[1],h_y3,delta_y3) * (
                                                                     square_int_per(xp[1],xbh_[1],h_x1,delta_x1)
                                                                   - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1)
                                                                   - (xp[1]-xbh_[1])*(  square_int_per(1,xbh_[1],h_x1,delta_x1)
                                                                                      - square_int_per(0,xbh_[1],h_x1,delta_x1))
                                                                        )
                 -  a_x20 * square_per(xp[2],ybh_[2],h_y4,delta_y4) * (
                                                                  square_int_per(xp[1],xbh_[2],h_x2,delta_x2)
                                                                - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2)
                                                                - (xp[1]-xbh_[2])*(  square_int_per(1,xbh_[2],h_x2,delta_x2)
                                                                                   - square_int_per(0,xbh_[2],h_x2,delta_x2))
                                                                        ) ;
  x[XX] = xmin + (xmax-xmin)*x[XX] ;

  x[YY] = xp[2] 
                 -  a_y10 * square_per(xp[1],xbh_[1],h_x3,delta_x3) * ( 
                                                                   square_int_per(xp[2],ybh_[1],h_y1,delta_y1)
                                                                 - square_int_per(ybh_[1],ybh_[1],h_y1,delta_y1)
                                                                 - (xp[2]-ybh_[1])*(  square_int_per(1,ybh_[1],h_y1,delta_y1)
                                                                                    - square_int_per(0,ybh_[1],h_y1,delta_y1)) 
                                                                        )
                 -  a_y20 * square_per(xp[1],xbh_[2],h_x4,delta_x4) * ( 
                                                                  square_int_per(xp[2],ybh_[2],h_y2,delta_y2)
                                                                - square_int_per(ybh_[2],ybh_[2],h_y2,delta_y2)
                                                                - (xp[2]-ybh_[2])*(  square_int_per(1,ybh_[2],h_y2,delta_y2)
                                                                                   - square_int_per(0,ybh_[2],h_y2,delta_y2)) 
                                                                        )  ;
  x[YY] = ymin + (ymax-ymin)*x[YY] ;

  x[ZZ] = xp[3] ;





#else 
  fprintf(stderr,"x_of_xp(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif

  return;
}

/****************************************************************************************

 xp_at_eq():
 ----------
   -- returns xp for a point at radius r at the equator
****************************************************************************************/
void xp_at_eq( double r, double *xp ) 
{
  double x[NDIM];

  xp[TT] = xp[PH] = 0.;
  xp[TH] = 0.5;

#if(COORD_TYPE_CHOICE==COORD_DIAGONAL || COORD_TYPE_CHOICE==COORD_MIXED || COORD_TYPE_CHOICE==COORD_DIAGONAL2 || COORD_TYPE_CHOICE==COORD_DIAGONAL3 )

  xp[RR] = log(r-R0);

#elif(COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD)

  xp[RR]  = log( ( r/coord_params->Rin_of_t - coord_params->f0_r ) / ( 1. - coord_params->f0_r ) ); 


#elif( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL ) 
  xp[RR] = log(r-R0);

#else 
  fprintf(stdout,"xp_at_eq(): should not be here \n"); fflush(stdout);
  fail(FAIL_BASIC,0);
#endif

  x_of_xp(x,xp);
  if( REL_DIFF_FUNC( x[TH] , (0.5*M_PI) ) > 10.*SMALL ) { 
    fprintf(stdout,"xp_at_eq(): th equator mismatch :  %26.16e  %26.16e \n", x[TH], (0.5*M_PI)); fflush(stdout);
  }


  return;
}


/***************************************************************

   dx_dxp_calc():
   ------
          -- finds the transformation array dx^\mu / dx^{\mu'} 
             that is used to transform x' contravariant 
             vectors to x contravariant vectors;

          -- Here, dx_dxp[i][j] =  dx^i / dx^{j'} 
 
***************************************************************/
void dx_dxp_calc( double *x, double *xp, double dx_dxp[][NDIM] )
{

  int i, j;
  double htmp1,htmp2;
  double cth, sth;

  // Start from delta function:

  DLOOP2   dx_dxp[i][j] = DELTA(i,j);


  /* Now set those elements that are different from identity: */

#if(   COORD_TYPE_CHOICE  == COORD_DIAGONAL ) 
    dx_dxp[RR][1] = x[RR] - R0;
    //-newcoords    dx_dxp[TH][2] = M_PI * ( 1.  +  2 * h_slope * cos(2 * M_PI * xp[2]) );
    dx_dxp[TH][2] = th_length  +  2 * M_PI * h_slope * cos(2 * M_PI * xp[2]) ;

#elif( COORD_TYPE_CHOICE  == COORD_MIXED    ) 
  htmp2 = X1_slope*(X1_0 - xp[1]);
  htmp1 = 2 * M_PI * xp[2]; 
  sincos(htmp1, &sth, &cth); 

  dx_dxp[RR][1] = x[RR] - R0;
  dx_dxp[TH][1] =  -4 * X1_slope * h_slope * sth / (M_PI * (1. + htmp2*htmp2));
  //-newcoords  dx_dxp[TH][2] =  M_PI  +  8 * h_slope * cth * atan(htmp2)  ; 
  dx_dxp[TH][2] =  th_length  +  8 * h_slope * cth * atan(htmp2)  ; 

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL2 ) 
    dx_dxp[RR][1] = x[RR] - R0;
    dx_dxp[TH][2] = dx2_dxp2_diag2(xp[2]) * norm_diag2;

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3 ) 
    //    htmp1 = 2 * ( PERIODIC( (xp[2]) , (0.), (1.) ) )  - 1.;
    htmp1 = 2 * seesaw(xp[2],0.,1.) - 1.;

    dx_dxp[RR][1] = x[RR] - R0;
    dx_dxp[TH][2] = M_PI*( h_slope + diag3_exponent*diag3_factor*pow(htmp1,(diag3_exponent-1)) );

#elif( COORD_TYPE_CHOICE  == COORD_IDENTITY  ) 
  /* Leave unchanged */


#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3_DYN_RAD )
  htmp1 = 2 * seesaw(xp[2],0.,1.) - 1.;

  dx_dxp[RR][0] = coord_params->dln_Rin_of_t_dt * x[RR] ;
  dx_dxp[RR][1] = x[RR] -  coord_params->Rin_of_t * coord_params->f0_r ;
  dx_dxp[TH][2] = M_PI*( h_slope + diag3_exponent*diag3_factor*pow(htmp1,(diag3_exponent-1)) );

#elif( COORD_TYPE_CHOICE  == COORD_WARPED_SPHERICAL )


  double dr_dxbh1, dr_dxbh2, dr_dybh1, dr_dybh2, dr_dybh3;
  double dphi_dxbh1, dphi_dxbh2, dphi_dybh1, dphi_dybh2, dphi_dybh3;

  int n;

  double r_of_yt_[3], drdy_of_yt_[3] ; 

  for(n=0; n<3; n++) {
//--orig     r_of_yt_[n] =  Rin_[n] + ( br_[n] - ar_[n]*s_[n] ) * xp[1]
//--orig               + ar_[n] * (   sinh(s_[n] * (xp[1] - ybh_[n]))
//--orig                        - sinh(-s_[n] * ybh_[n]) ) ;

    r_of_yt_[n] =  Rin_[n] + ( br_[n] ) * xp[1]
      + ar_[n] * (  sinh_terms( s_[n], xp[1], ybh_[n] ) );

    drdy_of_yt_[n] =  br_[n] 
        + s_[n] * ar_[n] * ( cosh( s_[n]*(xp[1] - ybh_[n]) ) - 1. ) ;
  }


  double dr1_dybh1, dr2_dybh2, dr3_dybh3 ;

  dr1_dybh1 =  ar_[1]*ar_[1]/(Rin_[1] - Rout_[1] + br_[1]) * s_[1] * (
                                                                          sinh(s_[1]*xp[1]) 
                                                                        + sinh(s_[1]*(1. - xp[1])) 
                                                                        - sinh(s_[1])
                                                                        + s_[1] * (
                                                                                   - cosh( s_[1]*(xp[1] - ybh_[1]) )  
                                                                                   + (1. - xp[1]) * cosh(s_[1]*ybh_[1]) 
                                                                                   + xp[1]*cosh( s_[1]*(1. - ybh_[1]) )
                                                                                   )
                                                                      ) ;

  dr2_dybh2 =  ar_[2]*ar_[2]/(Rin_[2] - Rout_[2] + br_[2]) * s_[2] * (
                                                                          sinh(s_[2]*xp[1]) 
                                                                        + sinh(s_[2]*(1. - xp[1])) 
                                                                        - sinh(s_[2])
                                                                        + s_[2] * (
                                                                                   - cosh( s_[2]*(xp[1] - ybh_[2]) ) 
                                                                                   + (1. - xp[1]) * cosh(s_[2]*ybh_[2]) 
                                                                                   + xp[1]*cosh( s_[2]*( 1. - ybh_[2]) ) 
                                                                                   )
                                                                      ) ;

  dr3_dybh3 =  ar_[0]*ar_[0]/(Rin_[0] - Rout_[0] + br_[0]) * s_[0] * (
                                                                          sinh(s_[0]*xp[1]) 
                                                                        + sinh(s_[0]*(1. - xp[1])) 
                                                                        - sinh(s_[0])
                                                                        + s_[0] * (
                                                                                   - cosh( s_[0]*(xp[1] - ybh_[0]) ) 
                                                                                   + (1. - xp[1]) * cosh(s_[0]*ybh_[0]) 
                                                                                   + xp[1]*cosh( s_[0]*( 1. - ybh_[0]) ) 
                                                                                   )
                                                                      ) ;


  dr_dxbh1 = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) *
             (r_of_yt_[1] - r_of_yt_[0]);

  dr_dxbh2 = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) *
             (r_of_yt_[2] - r_of_yt_[0]);

  dr_dybh1 = dr1_dybh1 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  dr_dybh2 = dr2_dybh2 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  dr_dybh3 = dr3_dybh3 * ( 1. - square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                              square_per(xp[3],xbh_[1],h_x3,delta_x3) 
                                                                            + square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                           )
                           ) ;



  dphi_dxbh1 = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                 * ( - square_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                     + square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[3] - xbh_[1]) * ( square_per(1,xbh_[1],h_x1,delta_x1) 
                                                           - square_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  dphi_dxbh2 = - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                 * ( - square_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                     + square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                     - square_int_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[3] - xbh_[2]) * ( square_per(1,xbh_[2],h_x2,delta_x2) 
                                                           - square_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  dphi_dybh1 = 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                 * (  square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                    - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                    - (xp[3] - xbh_[1]) * ( square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                          - square_int_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  dphi_dybh2 = 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                 * (  square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                    - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                    - (xp[3] - xbh_[2]) * ( square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                          - square_int_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  dx_dxp[RR][0] = dxbh1_dt * dr_dxbh1 + dxbh2_dt * dr_dxbh2 + dybh1_dt * dr_dybh1 + dybh2_dt * dr_dybh2 + dybh3_dt * dr_dybh3 ;

  dx_dxp[RR][2] = dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                (r_of_yt_[1] - r_of_yt_[0]) * square_per(xp[3],xbh_[1],h_x3,delta_x3)
                                                              + (r_of_yt_[2] - r_of_yt_[0]) * square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                              ) ;

  dx_dxp[RR][3] = - dr_dxbh1 - dr_dxbh2 ;

  dx_dxp[RR][1] = drdy_of_yt_[0] + square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                (drdy_of_yt_[1] - drdy_of_yt_[0]) * square_per(xp[3],xbh_[1],h_x3,delta_x3)
                                                              + (drdy_of_yt_[2] - drdy_of_yt_[0]) * square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                              ) ;


  dx_dxp[PH][3] = 2.*M_PI - 2.*M_PI * a_x10 *  square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                            * (square_per(xp[3],xbh_[1],h_x1,delta_x1) - square_int_per(1,xbh_[1],h_x1,delta_x1) + square_int_per(0,xbh_[1],h_x1,delta_x1))
                          - 2.*M_PI * a_x20 *  square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                            * (square_per(xp[3],xbh_[2],h_x2,delta_x2) - square_int_per(1,xbh_[2],h_x2,delta_x2) + square_int_per(0,xbh_[2],h_x2,delta_x2)) ;

  dx_dxp[PH][1] = - 2.*M_PI * a_x10 *   square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                    * ( 
                                          square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                        - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1)
                                        - (xp[3] - xbh_[1]) * (   square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                                - square_int_per(0,xbh_[1],h_x1,delta_x1) ) 
                                        )
                  - 2.*M_PI * a_x20 *   square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                    * ( 
                                         square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                       - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2)
                                       - (xp[3] - xbh_[2]) * (   square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                               - square_int_per(0,xbh_[2],h_x2,delta_x2) ) 
                                        ) ;

  dx_dxp[PH][2] = - 2.*M_PI * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                          a_x10 * square(xp[1],ybh_[1],h_y3,delta_y3) * (   square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                                                                                                          - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                                                                                                          - (xp[3]-xbh_[1])*(  square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                                                                                                             - square_int_per(0,xbh_[1],h_x1,delta_x1)) )
                                                                        + a_x20 * square(xp[1],ybh_[2],h_y4,delta_y4) * (   square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                                                                                                          - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                                                                                                          - (xp[3]-xbh_[2])*(   square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                                                                                                              - square_int_per(0,xbh_[2],h_x2,delta_x2)) )  
                                                                          ) ;

  dx_dxp[PH][0] =  dxbh1_dt * dphi_dxbh1 + dxbh2_dt * dphi_dxbh2 + dybh1_dt * dphi_dybh1 + dybh2_dt * dphi_dybh2 ;


  dx_dxp[TH][2] = M_PI - M_PI * a_z10 * (   square_per(xp[2],zbh_[1],h_z1,delta_z1) 
                                          - square_int_per(1,zbh_[1],h_z1,delta_z1)
                                          + square_int_per(0,zbh_[1],h_z1,delta_z1) ) ;

  dx_dxp[TH][0] = 0. ;
  dx_dxp[TH][1] = 0. ;
  dx_dxp[TH][3] = 0. ;


#elif( COORD_TYPE_CHOICE  == COORD_WARPED_CARTESIAN )

  // dX/dx
  dx_dxp[XX][1] = 1.
                     -  a_x10 * square_per(xp[2],ybh_[1],h_y3,delta_y3) * ( 
                                                                            square_per(xp[1],xbh_[1],h_x1,delta_x1)
                                                                           - (  square_int_per(1,xbh_[1],h_x1,delta_x1)
                                                                              - square_int_per(0,xbh_[1],h_x1,delta_x1)) 
                                                                            )
                     -  a_x20 * square_per(xp[2],ybh_[2],h_y4,delta_y4) * ( 
                                                                           square_per(xp[1],xbh_[2],h_x2,delta_x2)
                                                                          - (  square_int_per(1,xbh_[2],h_x2,delta_x2)
                                                                             - square_int_per(0,xbh_[2],h_x2,delta_x2)) 
                                                                            ) ;  
  dx_dxp[XX][1] = (xmax-xmin)*dx_dxp[XX][1] ;

  // dX/dy
  dx_dxp[XX][2] = -  a_x10 * dsquare_per(xp[2],ybh_[1],h_y3,delta_y3) * ( 
                                                                   square_int_per(xp[1],xbh_[1],h_x1,delta_x1)
                                                                 - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1)
                                                                 - (xp[1]-xbh_[1])*(  square_int_per(1,xbh_[1],h_x1,delta_x1)
                                                                                    - square_int_per(0,xbh_[1],h_x1,delta_x1)) 
                                                                        )
                  -  a_x20 * dsquare_per(xp[2],ybh_[2],h_y4,delta_y4) * ( 
                                                                  square_int_per(xp[1],xbh_[2],h_x2,delta_x2)
                                                                - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2)
                                                                - (xp[1]-xbh_[2])*(  square_int_per(1,xbh_[2],h_x2,delta_x2)
                                                                                   - square_int_per(0,xbh_[2],h_x2,delta_x2)) 
                                                                        ) ;  
  dx_dxp[XX][2] = (xmax-xmin)*dx_dxp[XX][2] ;

  // dY/dx
  dx_dxp[YY][1] = -  a_y10 * dsquare_per(xp[1],xbh_[1],h_x3,delta_x3) * ( 
                                                                   square_int_per(xp[2],ybh_[1],h_y1,delta_y1)
                                                                 - square_int_per(ybh_[1],ybh_[1],h_y1,delta_y1)
                                                                 - (xp[2]-ybh_[1])*(  square_int_per(1,ybh_[1],h_y1,delta_y1)
                                                                                    - square_int_per(0,ybh_[1],h_y1,delta_y1)) 
                                                                        )
                  -  a_y20 * dsquare_per(xp[1],xbh_[2],h_x4,delta_x4) * ( 
                                                                  square_int_per(xp[2],ybh_[2],h_y2,delta_y2)
                                                                - square_int_per(ybh_[2],ybh_[2],h_y2,delta_y2)
                                                                - (xp[2]-ybh_[2])*(  square_int_per(1,ybh_[2],h_y2,delta_y2)
                                                                                   - square_int_per(0,ybh_[2],h_y2,delta_y2)) 
                                                                        )  ;
  dx_dxp[YY][1] = (ymax-ymin)*dx_dxp[YY][1] ;

  // dY/dy
  dx_dxp[YY][2] = 1.
                     -  a_y10 * square_per(xp[1],xbh_[1],h_x3,delta_x3) * ( 
                                                                   square_per(xp[2],ybh_[1],h_y1,delta_y1)
                                                                 - (  square_int_per(1,ybh_[1],h_y1,delta_y1)
                                                                    - square_int_per(0,ybh_[1],h_y1,delta_y1)) 
                                                                            )
                     -  a_y20 * square_per(xp[1],xbh_[2],h_x4,delta_x4) * ( 
                                                                  square_per(xp[2],ybh_[2],h_y2,delta_y2)
                                                                - (  square_int_per(1,ybh_[2],h_y2,delta_y2)
                                                                   - square_int_per(0,ybh_[2],h_y2,delta_y2)) 
                                                                            )  ;
  dx_dxp[YY][2] = (ymax-ymin)*dx_dxp[YY][2] ;


  double dx_dxbh1, dx_dxbh2, dx_dybh1, dx_dybh2; 
  double dy_dxbh1, dy_dxbh2, dy_dybh1, dy_dybh2; 
  
  dx_dxbh1 = - a_x10 * square_per(xp[2],ybh_[1],h_y3,delta_y3)
                                 * ( - square_per(xp[1],xbh_[1],h_x1,delta_x1) 
                                     + square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[1] - xbh_[1]) * ( square_per(1,xbh_[1],h_x1,delta_x1) 
                                                           - square_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  dx_dxbh2 = - a_x20 * square_per(xp[2],ybh_[2],h_y4,delta_y4)
                                 * ( - square_per(xp[1],xbh_[2],h_x2,delta_x2) 
                                     + square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                     - square_int_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[1] - xbh_[2]) * ( square_per(1,xbh_[2],h_x2,delta_x2) 
                                                           - square_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  dx_dybh1 =  a_x10 * dsquare_per(xp[2],ybh_[1],h_y3,delta_y3) * ( 
                                                                   square_int_per(xp[1],xbh_[1],h_x1,delta_x1)
                                                                 - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1)
                                                                 - (xp[1]-xbh_[1])*(  square_int_per(1,xbh_[1],h_x1,delta_x1)
                                                                                    - square_int_per(0,xbh_[1],h_x1,delta_x1)) 
                                                                                     ) ;

  dx_dybh2 =  a_x20 * dsquare_per(xp[2],ybh_[2],h_y4,delta_y4) * ( 
                                                                  square_int_per(xp[1],xbh_[2],h_x2,delta_x2)
                                                                - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2)
                                                                - (xp[1]-xbh_[2])*(  square_int_per(1,xbh_[2],h_x2,delta_x2)
                                                                                   - square_int_per(0,xbh_[2],h_x2,delta_x2)) 
                                                                        ) ;  


  dy_dybh1 = - a_y10 * square_per(xp[1],xbh_[1],h_x3,delta_x3)
                                 * ( - square_per(xp[2],ybh_[1],h_y1,delta_y1) 
                                     + square_int_per(1,ybh_[1],h_y1,delta_y1) 
                                     - square_int_per(0,ybh_[1],h_y1,delta_y1)
                                     + (xp[2] - ybh_[1]) * ( square_per(1,ybh_[1],h_y1,delta_y1) 
                                                           - square_per(0,ybh_[1],h_y1,delta_y1) )
                                     ) ;

  dy_dybh2 = - a_y20 * square_per(xp[1],xbh_[2],h_x4,delta_x4)
                                 * ( - square_per(xp[2],ybh_[2],h_y2,delta_y2) 
                                     + square_int_per(1,ybh_[2],h_y2,delta_y2) 
                                     - square_int_per(0,ybh_[2],h_y2,delta_y2)
                                     + (xp[2] - ybh_[2]) * ( square_per(1,ybh_[2],h_y2,delta_y2) 
                                                           - square_per(0,ybh_[2],h_y2,delta_y2) )
                                     ) ;
  

  dy_dxbh1 = a_y10 * dsquare_per(xp[1],xbh_[1],h_x3,delta_x3) * ( 
                                                                   square_int_per(xp[2],ybh_[1],h_y1,delta_y1)
                                                                 - square_int_per(ybh_[1],ybh_[1],h_y1,delta_y1)
                                                                 - (xp[2]-ybh_[1])*(  square_int_per(1,ybh_[1],h_y1,delta_y1)
                                                                                    - square_int_per(0,ybh_[1],h_y1,delta_y1)) 
                                                                  ) ;

  dy_dxbh2 = a_y20 * dsquare_per(xp[1],xbh_[2],h_x4,delta_x4) * ( 
                                                                  square_int_per(xp[2],ybh_[2],h_y2,delta_y2)
                                                                - square_int_per(ybh_[2],ybh_[2],h_y2,delta_y2)
                                                                - (xp[2]-ybh_[2])*(  square_int_per(1,ybh_[2],h_y2,delta_y2)
                                                                                   - square_int_per(0,ybh_[2],h_y2,delta_y2)) 
                                                                  )  ;

  // dX/dt
  dx_dxp[XX][0] = dxbh1_dt * dx_dxbh1 + dxbh2_dt * dx_dxbh2 + dybh1_dt * dx_dybh1 + dybh2_dt * dx_dybh2 ;
  dx_dxp[XX][0] = (xmax-xmin)*dx_dxp[XX][0] ;

  // dY/dt
  dx_dxp[YY][0] = dxbh1_dt * dy_dxbh1 + dxbh2_dt * dy_dxbh2 + dybh1_dt * dy_dybh1 + dybh2_dt * dy_dybh2 ;

  dx_dxp[YY][0] = (ymax-ymin)*dx_dxp[YY][0] ;



#else 
  fprintf(stderr,"dx_dxp_calc(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif 

  return; 

}
  

/***************************************************************

   det_dx_dxp_calc_default():
   ------
          -- finds the determinant of the 
             transformation matrix dx^\mu / dx^{\mu'},  
             which is used to transform x' contravariant 
             vectors to x contravariant vectors;
          -- the square of the determinant is useful for 
             transforming the determinant of the metric in x^i
             coordinates;
          -- not very efficient since it recalculates dx_dxp() 
              that you probably have already calculated...
               oh well, it's neater this way and we will only 
               do this at t=0 anyway;
 
          -- assumes the only non-diagonal terms are 
                dxdxp[1][2], dxdxp[2][1]

***************************************************************/
double det_dx_dxp_calc_default( double *x, double *xp )
{
  
  double b[NDIM][NDIM];

  dx_dxp_calc( x, xp,  b ); 

  return( b[0][0]*b[3][3]*( b[1][1]*b[2][2] - b[1][2]*b[2][1] ) );

}
  

/***************************************************************

   det_dx_dxp_calc2_default():
   ------
          -- like det_dx_dxp_calc() but requires dx_dxp[]

***************************************************************/
double det_dx_dxp_calc2_default( double dx_dxp[][NDIM] )
{
  
  return( dx_dxp[0][0]*dx_dxp[3][3]*( dx_dxp[1][1]*dx_dxp[2][2] - dx_dxp[1][2]*dx_dxp[2][1] ) );

}
  
/***************************************************************

   det_dx_dxp_calc_general():
   ------
          -- finds the determinant of the 
             transformation matrix dx^\mu / dx^{\mu'},  
             which is used to transform x' contravariant 
             vectors to x contravariant vectors;
          -- the square of the determinant is useful for 
             transforming the determinant of the metric in x^i
             coordinates;
          -- not very efficient since it recalculates dx_dxp() 
              that you probably have already calculated...
               oh well, it's neater this way and we will only 
               do this at t=0 anyway;
 
          -- assumes the only non-diagonal terms are 
                dxdxp[1][2], dxdxp[2][1]

***************************************************************/
double det_dx_dxp_calc_general( double *x, double *xp )
{
  
  fprintf(stdout,"det_dx_dxp_calc_general(): not implemented yet!! \n"); fflush(stdout); 
  fail( FAIL_BASIC,0 );

  double b[NDIM][NDIM];

  dx_dxp_calc( x, xp,  b ); 

  return( b[0][0]*b[3][3]*( b[1][1]*b[2][2] - b[1][2]*b[2][1] ) );

}
  

/***************************************************************

   det_dx_dxp_calc2_general():
   ------
          -- like det_dx_dxp_calc() but requires dx_dxp[]

***************************************************************/
double det_dx_dxp_calc2_general( double dx_dxp[][NDIM] )
{

  fprintf(stdout,"det_dx_dxp_calc2_general(): not implemented yet!! \n"); fflush(stdout); 
  fail( FAIL_BASIC,0 );

  return( dx_dxp[0][0]*dx_dxp[3][3]*( dx_dxp[1][1]*dx_dxp[2][2] - dx_dxp[1][2]*dx_dxp[2][1] ) );

}
  

/***************************************************************

   dx_dxp_dxp_calc():
   ------
          -- finds the quantitiesy d^2 x^\mu / dx^{\mu'} dx^{\nu'} 
             that is used in finding the connection coefficients;

          -- Here, dx_dxp_dxp[i][j][k] =  d^2 x^i / dx^{j'} dx^{k'} 
 
***************************************************************/
void dx_dxp_dxp_calc( double *x, double *xp, double dx_dxp_dxp[][NDIM][NDIM] )
{

  int i, j, k;
  double htmp1,htmp2,htmp3,hprime;
  double sth, cth;

  // Start from nothing:
#if( USE_STRICT_ARRAY_BOUNDS )
  DLOOP3  { dx_dxp_dxp[i][j][k] = 0.; }
#else
  for( i=0 ; i < NDIM*NDIM*NDIM; i++ ) { dx_dxp_dxp[0][0][i] = 0.; }
#endif


#if(   COORD_TYPE_CHOICE  == COORD_DIAGONAL ) 
  dx_dxp_dxp[RR][1][1] = x[RR] - R0;
  dx_dxp_dxp[TH][2][2] = -4 * M_PI * M_PI * h_slope * sin( 2 * M_PI * xp[2] ); 


#elif( COORD_TYPE_CHOICE  == COORD_MIXED    ) 
  htmp2 = X1_slope*(X1_0 - xp[1]);
  htmp1 = 2 * M_PI * xp[2]; 
  hprime =  X1_slope  / (1. + htmp2*htmp2);
  sincos(htmp1, &sth, &cth); 

  dx_dxp_dxp[RR][1][1] = x[RR] - R0;
  dx_dxp_dxp[TH][1][1] =  -8 * h_slope * sth * htmp2 * hprime * hprime / M_PI ;
  dx_dxp_dxp[TH][1][2] =  dx_dxp_dxp[TH][2][1]  =  -8 * h_slope * cth * hprime ;
  dx_dxp_dxp[TH][2][2] =  -16 * M_PI * h_slope * sth * atan(htmp2)  ; 

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL2 ) 
  dx_dxp_dxp[RR][1][1] = x[RR] - R0;
  dx_dxp_dxp[TH][2][2] = dx2_dxp2_dxp2_diag2(xp[2]) * norm_diag2;

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3 ) 
  //  htmp1 = 2 * ( PERIODIC( (xp[2]) , (0.), (1.) ) )  - 1.;
  htmp1 = 2 * seesaw(xp[2],0.,1.) - 1.;

  dx_dxp_dxp[RR][1][1] = x[RR] - R0;
  dx_dxp_dxp[TH][2][2] = 2*M_PI*diag3_exponent*(diag3_exponent-1)*diag3_factor*pow(htmp1,(diag3_exponent-2));

#elif( COORD_TYPE_CHOICE  == COORD_IDENTITY  ) 
  /* Leave unchanged */

#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3_DYN_RAD )
  fprintf(stderr,"dx_dxp_dxp_calc(): Not implemented yet for  COORD_TYPE_CHOICE =  COORD_DIAGONAL3_DYN_RAD = : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);


#elif( COORD_TYPE_CHOICE  == COORD_WARPED_SPHERICAL )

  fprintf(stderr,"\nWARNING! \nAs it stands, the call to this function (dx_dxp_dxp_calc)\nwith COORD_WARPED_SPHERICAL only works when both r12_dot_dot and omega_dot are zero. \nSee TODO under update_coord_time_funcs.\n\n"); 
  fflush(stderr);
  int n ;
  double d2ri_dybhi2_[3] ;

  double dr_dxbh1, dr_dxbh2, dr_dybh1, dr_dybh2, dr_dybh3;
  double dphi_dxbh1, dphi_dxbh2, dphi_dybh1, dphi_dybh2, dphi_dybh3;

  double r_of_yt_[3], drdy_of_yt_[3], d2rdy2_of_yt_[3] ; 

  for(n=0; n<3; n++) {
//--orig     r_of_yt_[n] =  Rin_[n] + ( br_[n] - ar_[n]*s_[n] ) * xp[1]
//--orig               + ar_[n] * (   sinh(s_[n] * (xp[1] - ybh_[n]))
//--orig                        - sinh(-s_[n] * ybh_[n]) ) ;

    r_of_yt_[n] =  Rin_[n] + ( br_[n] ) * xp[1]
      + ar_[n] * (  sinh_terms( s_[n], xp[1], ybh_[n] ) );

    drdy_of_yt_[n] =  br_[n] 
        + s_[n] * ar_[n] * ( cosh( s_[n]*(xp[1] - ybh_[n]) ) - 1. ) ;

    d2rdy2_of_yt_[n] = s_[n]*s_[n] * ar_[n] * sinh( s_[n]*(xp[1] - ybh_[n]) ) ;  
  }

  double dr1_dybh1, dr2_dybh2, dr3_dybh3 ;

  dr1_dybh1 =  ar_[1]*ar_[1]/(Rin_[1] - Rout_[1] + br_[1]) * s_[1] * (
                                                                          sinh(s_[1]*xp[1]) 
                                                                        + sinh(s_[1]*(1. - xp[1])) 
                                                                        - sinh(s_[1])
                                                                        + s_[1] * (
                                                                                   - cosh( s_[1]*(xp[1] - ybh_[1]) )  
                                                                                   + (1. - xp[1]) * cosh(s_[1]*ybh_[1]) 
                                                                                   + xp[1]*cosh( s_[1]*(1. - ybh_[1]) )
                                                                                   )
                                                                      ) ;

  dr2_dybh2 =  ar_[2]*ar_[2]/(Rin_[2] - Rout_[2] + br_[2]) * s_[2] * (
                                                                          sinh(s_[2]*xp[1]) 
                                                                        + sinh(s_[2]*(1. - xp[1])) 
                                                                        - sinh(s_[2])
                                                                        + s_[2] * (
                                                                                   - cosh( s_[2]*(xp[1] - ybh_[2]) ) 
                                                                                   + (1. - xp[1]) * cosh(s_[2]*ybh_[2]) 
                                                                                   + xp[1]*cosh( s_[2]*( 1. - ybh_[2]) ) 
                                                                                   )
                                                                      ) ;

  dr3_dybh3 =  ar_[0]*ar_[0]/(Rin_[0] - Rout_[0] + br_[0]) * s_[0] * (
                                                                          sinh(s_[0]*xp[1]) 
                                                                        + sinh(s_[0]*(1. - xp[1])) 
                                                                        - sinh(s_[0])
                                                                        + s_[0] * (
                                                                                   - cosh( s_[0]*(xp[1] - ybh_[0]) ) 
                                                                                   + (1. - xp[1]) * cosh(s_[0]*ybh_[0]) 
                                                                                   + xp[1]*cosh( s_[0]*( 1. - ybh_[0]) ) 
                                                                                   )
                                                                      ) ;



  dr_dxbh1 = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) *
             (r_of_yt_[1] - r_of_yt_[0]);

  dr_dxbh2 = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) *
             (r_of_yt_[2] - r_of_yt_[0]);

  dr_dybh1 = dr1_dybh1 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  dr_dybh2 = dr2_dybh2 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  dr_dybh3 = dr3_dybh3 * ( 1. - square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                              square_per(xp[3],xbh_[1],h_x3,delta_x3) 
                                                                            + square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                           )
                           ) ;



  dphi_dxbh1 = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                 * ( - square_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                     + square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[3] - xbh_[1]) * ( square_per(1,xbh_[1],h_x1,delta_x1) 
                                                           - square_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  dphi_dxbh2 = - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                 * ( - square_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                     + square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                     - square_int_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[3] - xbh_[2]) * ( square_per(1,xbh_[2],h_x2,delta_x2) 
                                                           - square_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  dphi_dybh1 = 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                 * (  square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                    - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                    - (xp[3] - xbh_[1]) * ( square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                          - square_int_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  dphi_dybh2 = 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                 * (  square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                    - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                    - (xp[3] - xbh_[2]) * ( square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                          - square_int_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  for(n=0; n<3; n++) {

    d2ri_dybhi2_[n] =   3*s_[n] - 6*s_[n]*xp[1] 
                      + 4*sinh(s_[n]*(xp[1] - ybh_[n])) 
                      - 2*sinh(s_[n]*(1 + xp[1] - ybh_[n])) 
                      + 2*sinh(s_[n]*ybh_[n]) 
                      + 2*sinh(s_[n]*(1 - xp[1] + ybh_[n]))
                      + 2*sinh(s_[n]*(-2 + xp[1] + ybh_[n]))
                      - 4*sinh(s_[n]*(-1 + xp[1] + ybh_[n]))
                      + 2*sinh(s_[n]*(xp[1] + ybh_[n]))
                      - 2*sinh(s_[n] - s_[n]*ybh_[n])
                      - s_[n]*(  
                                 (3 - 6*xp[1])*cosh(s_[n])
                                + 3*cosh(s_[n]*xp[1])
                                - 3*cosh(s_[n] - s_[n]*xp[1])
                                + cosh(s_[n]*(xp[1] - 2*ybh_[n]))
                                - cosh(s_[n]*(1 + xp[1] - 2*ybh_[n]))
                                - cosh(2*s_[n]*ybh_[n]) 
                                + cosh(s_[n] - 2*s_[n]*ybh_[n])
                                - 2*s_[n]*(  sinh(s_[n]*(xp[1] - ybh_[n]))
                                           + sinh(s_[n]*ybh_[n]) )
                                + xp[1]*(   
                                              cosh(2*s_[n]*(-1 + ybh_[n]))
                                            + cosh(2*s_[n]*ybh_[n]) 
                                            - 2*cosh(s_[n] - 2*s_[n]*ybh_[n])
                                            + 2*s_[n]*(   sinh(s_[n]*ybh_[n])
                                                        + sinh(s_[n] - s_[n]*ybh_[n])) 
                                            )
                                 ) 
                        + 2*sinh(2*s_[n] - s_[n]*ybh_[n])
                        - 2*sinh(s_[n] + s_[n]*ybh_[n]) ;

      d2ri_dybhi2_[n] =   d2ri_dybhi2_[n] * 0.5*s_[n]*s_[n] * ar_[n]*ar_[n]*ar_[n]  
                        / ( (Rin_[n] - Rout_[n] + br_[n])*(Rin_[n] - Rout_[n] + br_[n]) ) ;

  }

  double d2r_dxbh1_dxbh1, d2r_dxbh2_dxbh2, d2r_dybh1_dybh1, d2r_dybh2_dybh2, d2r_dybh3_dybh3 ;

  d2r_dxbh1_dxbh1 = square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square_per(xp[3],xbh_[1],h_x3,delta_x3) *
             (r_of_yt_[1] - r_of_yt_[0]);

  d2r_dxbh2_dxbh2 = square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square_per(xp[3],xbh_[2],h_x4,delta_x4) *
             (r_of_yt_[2] - r_of_yt_[0]);


  d2r_dybh1_dybh1 = d2ri_dybhi2_[1] * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  d2r_dybh2_dybh2 = d2ri_dybhi2_[2] * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  d2r_dybh3_dybh3 = d2ri_dybhi2_[0] * ( 1. - square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                                          square_per(xp[3],xbh_[1],h_x3,delta_x3) 
                                                                                        + square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                                        )
                                        ) ;

  double d2r_dxbh1_dybh1, d2r_dxbh2_dybh2, d2r_dxbh1_dybh3, d2r_dxbh2_dybh3 ;

  d2r_dxbh1_dybh1 = - dr1_dybh1 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  d2r_dxbh2_dybh2 = - dr2_dybh2 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  d2r_dxbh1_dybh3 = dr3_dybh3 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  d2r_dxbh2_dybh3 = dr3_dybh3 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) ;


  // d^2r/dt^2:
  dx_dxp_dxp[RR][0][0] =  d2xbh1_dt2 * dr_dxbh1 + d2xbh2_dt2 * dr_dxbh2
                        + d2ybh1_dt2 * dr_dybh1 + d2ybh2_dt2 * dr_dybh2 + d2ybh3_dt2 * dr_dybh3 
                        + d2r_dxbh1_dxbh1 * dxbh1_dt * dxbh1_dt 
                        + d2r_dxbh2_dxbh2 * dxbh2_dt * dxbh2_dt
                        + d2r_dybh1_dybh1 * dybh1_dt * dybh1_dt
                        + d2r_dybh2_dybh2 * dybh2_dt * dybh2_dt
                        + d2r_dybh3_dybh3 * dybh3_dt * dybh3_dt
                        + 2 * d2r_dxbh1_dybh1 * dxbh1_dt * dybh1_dt
                        + 2 * d2r_dxbh2_dybh2 * dxbh2_dt * dybh2_dt
                        + 2 * d2r_dxbh1_dybh3 * dxbh1_dt * dybh3_dt
                        + 2 * d2r_dxbh2_dybh3 * dxbh2_dt * dybh3_dt ;


  // d^2r/dtdx:
  double d2r_dxbh1_dx, d2r_dxbh2_dx ;
  double d2r_dybh1_dx, d2r_dybh2_dx, d2r_dybh3_dx;

  d2r_dxbh1_dx = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square_per(xp[3],xbh_[1],h_x3,delta_x3) *
                   (r_of_yt_[1] - r_of_yt_[0]);

  d2r_dxbh2_dx = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square_per(xp[3],xbh_[2],h_x4,delta_x4) *
                   (r_of_yt_[2] - r_of_yt_[0]);

  d2r_dybh1_dx = dr1_dybh1 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  d2r_dybh2_dx = dr2_dybh2 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  d2r_dybh3_dx = - dr3_dybh3 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                            dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) 
                                                                          + dsquare_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                          ) ;


  dx_dxp_dxp[RR][0][3] =   d2r_dxbh1_dx * dxbh1_dt + d2r_dxbh2_dx * dxbh2_dt 
                         + d2r_dybh1_dx * dybh1_dt + d2r_dybh2_dx * dybh2_dt + d2r_dybh3_dx * dybh3_dt ;



  // d^2r/dtdy:
  double d2r1_dybh1_dy, d2r2_dybh2_dy, d2r3_dybh3_dy ;

  d2r1_dybh1_dy = - ar_[1]*ar_[1]/(Rin_[1] - Rout_[1] + br_[1]) * s_[1] * s_[1]
                       * ( 
                          - cosh(s_[1]*xp[1]) + cosh(s_[1] - s_[1]*xp[1])
                          + cosh(s_[1]*ybh_[1]) - cosh(s_[1] - s_[1]*ybh_[1]) 
                          + s_[1]*sinh(s_[1]*(xp[1] - ybh_[1]))
                           ) ;

  d2r2_dybh2_dy = - ar_[2]*ar_[2]/(Rin_[2] - Rout_[2] + br_[2]) * s_[2] * s_[2]
                       * ( 
                          - cosh(s_[2]*xp[1]) + cosh(s_[2] - s_[2]*xp[1])
                          + cosh(s_[2]*ybh_[2]) - cosh(s_[2] - s_[2]*ybh_[2]) 
                          + s_[2]*sinh(s_[2]*(xp[1] - ybh_[2]))
                           ) ;

  d2r3_dybh3_dy = - ar_[0]*ar_[0]/(Rin_[0] - Rout_[0] + br_[0]) * s_[0] * s_[0]
                       * ( 
                          - cosh(s_[0]*xp[1]) + cosh(s_[0] - s_[0]*xp[1])
                          + cosh(s_[0]*ybh_[0]) - cosh(s_[0] - s_[0]*ybh_[0]) 
                          + s_[0]*sinh(s_[0]*(xp[1] - ybh_[0]))
                           ) ;

  double d2r_dxbh1_dy, d2r_dxbh2_dy ;
  double d2r_dbyh1_dy, d2r_dybh2_dy, d2r_dybh3_dy ;

  d2r_dxbh1_dy = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) *
                   (drdy_of_yt_[1] - drdy_of_yt_[0]);

  d2r_dxbh2_dy = - square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) *
                   (drdy_of_yt_[2] - drdy_of_yt_[0]);

  d2r_dbyh1_dy = d2r1_dybh1_dy * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  d2r_dybh2_dy = d2r2_dybh2_dy * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  d2r_dybh3_dy = d2r3_dybh3_dy * ( 1. - square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                                        square_per(xp[3],xbh_[1],h_x3,delta_x3) 
                                                                                      + square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                                   )
                                   ) ;


  dx_dxp_dxp[RR][0][1] =   d2r_dxbh1_dy * dxbh1_dt + d2r_dxbh2_dy * dxbh2_dt 
                         + d2r_dbyh1_dy * dybh1_dt + d2r_dybh2_dy * dybh2_dt + d2r_dybh3_dy * dybh3_dt ;



  // d^2r/dtdz:
  double d2r_dxbh1_dz, d2r_dxbh2_dz ;
  double d2r_dybh1_dz, d2r_dybh2_dz, d2r_dybh3_dz;

  d2r_dxbh1_dz = - dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) *
                   (r_of_yt_[1] - r_of_yt_[0]);

  d2r_dxbh2_dz = - dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) *
                   (r_of_yt_[2] - r_of_yt_[0]);

  d2r_dybh1_dz = dr1_dybh1 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[1],h_x3,delta_x3) ;

  d2r_dybh2_dz = dr2_dybh2 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * square_per(xp[3],xbh_[2],h_x4,delta_x4) ;

  d2r_dybh3_dz = - dr3_dybh3 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                            square_per(xp[3],xbh_[1],h_x3,delta_x3) 
                                                                          + square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                          ) ;

  dx_dxp_dxp[RR][0][2] =   d2r_dxbh1_dz * dxbh1_dt + d2r_dxbh2_dz * dxbh2_dt 
                         + d2r_dybh1_dz * dybh1_dt + d2r_dybh2_dz * dybh2_dt + d2r_dybh3_dz * dybh3_dt ;


  // d2r/dxdt
  dx_dxp_dxp[RR][3][0] = dx_dxp_dxp[RR][0][3] ;

  // d2r/dxdx
  dx_dxp_dxp[RR][3][3] =   square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square_per(xp[3],xbh_[1],h_x3,delta_x3) *
                                                                     (r_of_yt_[1] - r_of_yt_[0])
                         + square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square_per(xp[3],xbh_[2],h_x4,delta_x4) *
                                                                     (r_of_yt_[2] - r_of_yt_[0]);

  // d2r/dxdy
  dx_dxp_dxp[RR][3][1] =   square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) *
                                                                     (drdy_of_yt_[1] - drdy_of_yt_[0])
                         + square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) *
                                                                     (drdy_of_yt_[2] - drdy_of_yt_[0]);

  // d2r/dxdz
  dx_dxp_dxp[RR][3][2] =   dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) *
                                                                      (r_of_yt_[1] - r_of_yt_[0])
                         + dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) *
                                                                      (r_of_yt_[2] - r_of_yt_[0]);


  // d2r/dydt
  dx_dxp_dxp[RR][1][0] = dx_dxp_dxp[RR][0][1] ;

  // d2r/dydx
  dx_dxp_dxp[RR][1][3] = dx_dxp_dxp[RR][3][1] ;

  // d2r/dydy
  dx_dxp_dxp[RR][1][1] = d2rdy2_of_yt_[0] + square_per(xp[2],zbh_[2],h_z2,delta_z2) * ( 
                                                                                         (d2rdy2_of_yt_[1] - d2rdy2_of_yt_[0]) * square_per(xp[3],xbh_[1],h_x3,delta_x3)
                                                                                       + (d2rdy2_of_yt_[2] - d2rdy2_of_yt_[0]) * square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                                        ) ;

  // d2r/dydz
  dx_dxp_dxp[RR][1][2] =  dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                        (drdy_of_yt_[1] - drdy_of_yt_[0]) * square_per(xp[3],xbh_[1],h_x3,delta_x3)
                                                                      + (drdy_of_yt_[2] - drdy_of_yt_[0]) * square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                      ) ;


  // d2r/dzdt
  dx_dxp_dxp[RR][2][0] = dx_dxp_dxp[RR][0][2] ;

  // d2r/dzdx
  dx_dxp_dxp[RR][2][3] = dx_dxp_dxp[RR][3][2] ;

  // d2r/dzdy
  dx_dxp_dxp[RR][2][1] = dx_dxp_dxp[RR][1][2] ;

  // d2r/dzdz
  dx_dxp_dxp[RR][2][2] = d2square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                        (r_of_yt_[1] - r_of_yt_[0]) * square_per(xp[3],xbh_[1],h_x3,delta_x3)
                                                                      + (r_of_yt_[2] - r_of_yt_[0]) * square_per(xp[3],xbh_[2],h_x4,delta_x4)
                                                                      ) ;



  // d^2phi/dt^2:
  double d2phi_dxbh1_dxbh1, d2phi_dxbh2_dxbh2, d2phi_dxbh1_dybh1 ;
  double d2phi_dxbh2_dybh2, d2phi_dybh1_dybh1, d2phi_dybh2_dybh2 ;

  d2phi_dxbh1_dxbh1 = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                * (    dsquare_per(xp[3],xbh_[1],h_x1,delta_x1)
                                     - square_per(1,xbh_[1],h_x1,delta_x1)
                                     + square_per(0,xbh_[1],h_x1,delta_x1)
                                     - square_per(1,xbh_[1],h_x1,delta_x1) 
                                     + square_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[3] - xbh_[1]) * ( - dsquare_per(1,xbh_[1],h_x1,delta_x1)
                                                             + dsquare_per(0,xbh_[1],h_x1,delta_x1) )
                                    ) ;

  d2phi_dxbh2_dxbh2 = - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                * (    dsquare_per(xp[3],xbh_[2],h_x2,delta_x2)
                                     - square_per(1,xbh_[2],h_x2,delta_x2)
                                     + square_per(0,xbh_[2],h_x2,delta_x2)
                                     - square_per(1,xbh_[2],h_x2,delta_x2)
                                     + square_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[3] - xbh_[2]) * ( - dsquare_per(1,xbh_[2],h_x2,delta_x2)
                                                             + dsquare_per(0,xbh_[2],h_x2,delta_x2) )
                                    ) ;

  d2phi_dxbh1_dybh1 = 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                              * (    - square_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                     + square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[3] - xbh_[1]) * ( square_per(1,xbh_[1],h_x1,delta_x1) 
                                                           - square_per(0,xbh_[1],h_x1,delta_x1) )
                                  ) ;
  
  d2phi_dxbh2_dybh2 = 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                              * (    - square_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                     + square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                     - square_int_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[3] - xbh_[2]) * ( square_per(1,xbh_[2],h_x2,delta_x2) 
                                                           - square_per(0,xbh_[2],h_x2,delta_x2) )
                                  ) ;

  d2phi_dybh1_dybh1 = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square(xp[1],ybh_[1],h_y3,delta_y3)
                                * (  square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                   - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                   - (xp[3] - xbh_[1]) * ( square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                         - square_int_per(0,xbh_[1],h_x1,delta_x1) )
                                   ) ;

  d2phi_dybh2_dybh2 = - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square(xp[1],ybh_[2],h_y4,delta_y4)
                                * (   square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                    - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                    - (xp[3] - xbh_[2]) * ( square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                          - square_int_per(0,xbh_[2],h_x2,delta_x2) )
                                   ) ;


  dx_dxp_dxp[PH][0][0] =  d2xbh1_dt2 * dphi_dxbh1 + d2xbh2_dt2 * dphi_dxbh2 
                        + d2ybh1_dt2 * dphi_dybh1 + d2ybh2_dt2 * dphi_dybh2
                        + d2phi_dxbh1_dxbh1 * dxbh1_dt * dxbh1_dt 
                        + d2phi_dxbh2_dxbh2 * dxbh2_dt * dxbh2_dt
                        + 2 * d2phi_dxbh1_dybh1 * dxbh1_dt * dybh1_dt
                        + 2 * d2phi_dxbh2_dybh2 * dxbh2_dt * dybh2_dt
                        + d2phi_dybh1_dybh1 * dybh1_dt * dybh1_dt
                        + d2phi_dybh2_dybh2 * dybh2_dt * dybh2_dt ;


  // d^2phi/dtdx:
  double d2phi_dxbh1_dx, d2phi_dxbh2_dx, d2phi_dybh1_dx, d2phi_dybh2_dx ;

  d2phi_dxbh1_dx = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                 * ( - dsquare_per(xp[3],xbh_[1],h_x1,delta_x1)
                                     + square_per(1,xbh_[1],h_x1,delta_x1)
                                     - square_per(0,xbh_[1],h_x1,delta_x1)
                                     ) ;

  d2phi_dxbh2_dx = - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                 * ( - dsquare_per(xp[3],xbh_[2],h_x2,delta_x2)
                                     + square_per(1,xbh_[2],h_x2,delta_x2)
                                     - square_per(0,xbh_[2],h_x2,delta_x2)
                                     ) ;

  d2phi_dybh1_dx =  2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                 * (  square_per(xp[3],xbh_[1],h_x1,delta_x1)
                                    - square_int_per(1,xbh_[1],h_x1,delta_x1)
                                    + square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     ) ;

  d2phi_dybh2_dx = 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                 * (  square_per(xp[3],xbh_[2],h_x2,delta_x2)
                                    - square_int_per(1,xbh_[2],h_x2,delta_x2)
                                    + square_int_per(0,xbh_[2],h_x2,delta_x2)
                                    ) ;

  dx_dxp_dxp[PH][0][3] =   d2phi_dxbh1_dx * dxbh1_dt + d2phi_dxbh2_dx * dxbh2_dt 
                         + d2phi_dybh1_dx * dybh1_dt + d2phi_dybh2_dx * dybh2_dt ;


  // d^2phi/dtdy:
  double d2phi_dxbh1_dy, d2phi_dxbh2_dy, d2phi_dybh1_dy, d2phi_dybh2_dy ;

  d2phi_dxbh1_dy = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                 * ( - square_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                     + square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[3] - xbh_[1]) * ( square_per(1,xbh_[1],h_x1,delta_x1) 
                                                           - square_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  d2phi_dxbh2_dy = - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                 * ( - square_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                     + square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                     - square_int_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[3] - xbh_[2]) * ( square_per(1,xbh_[2],h_x2,delta_x2) 
                                                           - square_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  d2phi_dybh1_dy = 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square(xp[1],ybh_[1],h_y3,delta_y3)
                                 * (  square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                    - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                    - (xp[3] - xbh_[1]) * ( square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                          - square_int_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  d2phi_dybh2_dy = 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square(xp[1],ybh_[2],h_y4,delta_y4)
                                 * (  square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                    - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                    - (xp[3] - xbh_[2]) * ( square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                          - square_int_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  dx_dxp_dxp[PH][0][1] =   d2phi_dxbh1_dy * dxbh1_dt + d2phi_dxbh2_dy * dxbh2_dt 
                         + d2phi_dybh1_dy * dybh1_dt + d2phi_dybh2_dy * dybh2_dt ;




  // d^2phi/dtdz:
  double d2phi_dxbh1_dz, d2phi_dxbh2_dz, d2phi_dybh1_dz, d2phi_dybh2_dz ;

  d2phi_dxbh1_dz = - 2.*M_PI * a_x10 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                 * ( - square_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                     + square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)
                                     + (xp[3] - xbh_[1]) * ( square_per(1,xbh_[1],h_x1,delta_x1) 
                                                           - square_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  d2phi_dxbh2_dz = - 2.*M_PI * a_x20 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                 * ( - square_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                     + square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                     - square_int_per(0,xbh_[2],h_x2,delta_x2)
                                     + (xp[3] - xbh_[2]) * ( square_per(1,xbh_[2],h_x2,delta_x2) 
                                                           - square_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  d2phi_dybh1_dz = 2.*M_PI * a_x10 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                 * (  square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                    - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                    - (xp[3] - xbh_[1]) * ( square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                          - square_int_per(0,xbh_[1],h_x1,delta_x1) )
                                     ) ;

  d2phi_dybh2_dz = 2.*M_PI * a_x20 * dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                 * (  square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                    - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                    - (xp[3] - xbh_[2]) * ( square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                          - square_int_per(0,xbh_[2],h_x2,delta_x2) )
                                     ) ;

  dx_dxp_dxp[PH][0][2] =  d2phi_dxbh1_dz * dxbh1_dt + d2phi_dxbh2_dz * dxbh2_dt 
                        + d2phi_dybh1_dz * dybh1_dt + d2phi_dybh2_dz * dybh2_dt ;
  


  // d2phi/dxdt            
  dx_dxp_dxp[PH][3][0] = dx_dxp_dxp[PH][0][3] ;
                         
  // d2phi/dxdx            
  dx_dxp_dxp[PH][3][3] = - 2.*M_PI * a_x10 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                           * dsquare_per(xp[3],xbh_[1],h_x1,delta_x1)
                         - 2.*M_PI * a_x20 * square_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                           * dsquare_per(xp[3],xbh_[2],h_x2,delta_x2) ;

  // d2phi/dxdy            
  dx_dxp_dxp[PH][3][1] =  - 2.*M_PI * a_x10 *  square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                            * (square_per(xp[3],xbh_[1],h_x1,delta_x1) - square_int_per(1,xbh_[1],h_x1,delta_x1) + square_int_per(0,xbh_[1],h_x1,delta_x1))
                          - 2.*M_PI * a_x20 *  square_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                            * (square_per(xp[3],xbh_[2],h_x2,delta_x2) - square_int_per(1,xbh_[2],h_x2,delta_x2) + square_int_per(0,xbh_[2],h_x2,delta_x2)) ;
                         
  // d2phi/dxdz            
  dx_dxp_dxp[PH][3][2] =  - 2.*M_PI * a_x10 *  dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[1],h_y3,delta_y3)
                                            * (square_per(xp[3],xbh_[1],h_x1,delta_x1) - square_int_per(1,xbh_[1],h_x1,delta_x1) + square_int_per(0,xbh_[1],h_x1,delta_x1))
                          - 2.*M_PI * a_x20 *  dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * square(xp[1],ybh_[2],h_y4,delta_y4)
                                            * (square_per(xp[3],xbh_[2],h_x2,delta_x2) - square_int_per(1,xbh_[2],h_x2,delta_x2) + square_int_per(0,xbh_[2],h_x2,delta_x2)) ;

  // d2phi/dydt            
  dx_dxp_dxp[PH][1][0] = dx_dxp_dxp[PH][0][1] ;
                         
  // d2phi/dydx            
  dx_dxp_dxp[PH][1][3] = dx_dxp_dxp[PH][3][1] ;
                         
  // d2phi/dydy            
  dx_dxp_dxp[PH][1][1] =  - 2.*M_PI * a_x10 *   square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square(xp[1],ybh_[1],h_y3,delta_y3)
                                    * ( 
                                          square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                        - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1)
                                        - (xp[3] - xbh_[1]) * (   square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                                - square_int_per(0,xbh_[1],h_x1,delta_x1) ) 
                                        )
                          - 2.*M_PI * a_x20 *   square_per(xp[2],zbh_[2],h_z2,delta_z2) * d2square(xp[1],ybh_[2],h_y4,delta_y4)
                                    * ( 
                                         square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                       - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2)
                                       - (xp[3] - xbh_[2]) * (   square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                               - square_int_per(0,xbh_[2],h_x2,delta_x2) ) 
                                        ) ;

  // d2phi/dydz
  dx_dxp_dxp[PH][1][2] = - 2.*M_PI * a_x10 *  dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[1],h_y3,delta_y3)
                                    * ( 
                                          square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                        - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1)
                                        - (xp[3] - xbh_[1]) * (   square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                                - square_int_per(0,xbh_[1],h_x1,delta_x1) ) 
                                        )
                          - 2.*M_PI * a_x20 *   dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) * dsquare(xp[1],ybh_[2],h_y4,delta_y4)
                                    * ( 
                                         square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                       - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2)
                                       - (xp[3] - xbh_[2]) * (   square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                               - square_int_per(0,xbh_[2],h_x2,delta_x2) ) 
                                        ) ;

  // d2phi/dzdt            
  dx_dxp_dxp[PH][2][0] = dx_dxp_dxp[PH][0][2] ;
                         
  // d2phi/dzdx            
  dx_dxp_dxp[PH][2][3] = dx_dxp_dxp[PH][3][2] ;
                         
  // d2phi/dzdy            
  dx_dxp_dxp[PH][2][1] = dx_dxp_dxp[PH][1][2] ;
                        
  // d2phi/dzdz
  dx_dxp_dxp[PH][2][2] = - 2.*M_PI * d2square_per(xp[2],zbh_[2],h_z2,delta_z2) * (
                                                                                  a_x10 * square(xp[1],ybh_[1],h_y3,delta_y3) * (   square_int_per(xp[3],xbh_[1],h_x1,delta_x1) 
                                                                                                                                  - square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) 
                                                                                                                                  - (xp[3]-xbh_[1])*(  square_int_per(1,xbh_[1],h_x1,delta_x1) 
                                                                                                                                                     - square_int_per(0,xbh_[1],h_x1,delta_x1)) )
                                                                                + a_x20 * square(xp[1],ybh_[2],h_y4,delta_y4) * (   square_int_per(xp[3],xbh_[2],h_x2,delta_x2) 
                                                                                                                                  - square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) 
                                                                                                                                  - (xp[3]-xbh_[2])*(   square_int_per(1,xbh_[2],h_x2,delta_x2) 
                                                                                                                                                      - square_int_per(0,xbh_[2],h_x2,delta_x2)) )  
                                                                                  ) ;

  // d2theta/dzdz (only non-null component)
  dx_dxp_dxp[TH][2][2] = - M_PI * a_z10 * dsquare_per(xp[2],zbh_[1],h_z1,delta_z1) ;


#else 
  fprintf(stderr,"dx_dxp_dxp_calc(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif

  return;

}
  
/***************************************************************

   dxp_dx_calc_default():
   -------------------
    -- finds the transformation array dx^{\mu'} / dx^\mu
       that is used to transform x contravariant 
       vectors to x' contravariant vectors;

    -- Here, dxp_dx[i][j] =  dx^{i'} / dx^{j} 

    -- pretty much calculate the elements from dx^{i} / dx^{j'} 
       as given by dx_dxp_calc(), using the following 
       identity:

 (dx^{i'} / dx^{k}) (dx^{k} / dx^{j'})  = DELTA(i',j')

    -- currently implemented to just take the inverse 
       of the inner 2x2 matrix ([1][RR],[1][TH],[2][RR],[2][TH])
       and leave the rest as the identity transformation; 
       SO IT ASSUMES THAT ONLY X1,X2 ARE COUPLED !!
 
***************************************************************/
void dxp_dx_calc_default( double *x, double *xp, double dxp_dx[][NDIM])
{

  int i, j;
  double inv_det, a11, a12, a21, a22;

  /* Start from inverse matrix:  */
  dx_dxp_calc( x, xp, dxp_dx );
  
  a11 = dxp_dx[RR][1] ;
  a12 = dxp_dx[RR][2] ; 
  a21 = dxp_dx[TH][1] ;
  a22 = dxp_dx[TH][2] ; 

  inv_det = 1. / ( a11 * a22 - a12 * a21 ) ; 

  
  /* Start from the identity matrix: */
  DLOOP2   dxp_dx[i][j] = DELTA(i,j);

  dxp_dx[1][RR] =  a22 * inv_det; 
  dxp_dx[1][TH] = -a12 * inv_det; 
  dxp_dx[2][RR] = -a21 * inv_det; 
  dxp_dx[2][TH] =  a11 * inv_det; 

  return;

}


/***************************************************************

   dxp_dx_calc2_default():
   --------------
    -- like dxp_dx_calc_default() but uses dx_dxp() as an argument.
 
***************************************************************/
void dxp_dx_calc2_default( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM])
{

  int i, j;
  double inv_det, a11, a12, a21, a22;

  /* Start from inverse matrix:  */
  a11 = dx_dxp[RR][1] ;
  a12 = dx_dxp[RR][2] ; 
  a21 = dx_dxp[TH][1] ;
  a22 = dx_dxp[TH][2] ; 

  inv_det = 1. / ( a11 * a22 - a12 * a21 ) ; 
  
  /* Start from the identity matrix: */
  DLOOP2   dxp_dx[i][j] = DELTA(i,j);

  dxp_dx[1][RR] =  a22 * inv_det; 
  dxp_dx[1][TH] = -a12 * inv_det; 
  dxp_dx[2][RR] = -a21 * inv_det; 
  dxp_dx[2][TH] =  a11 * inv_det; 

  return;

}


/***************************************************************

   dxp_dx_calc_dyn_rad():
   --------------
    -- finds the transformation array dx^{\mu'} / dx^\mu
       that is used to transform x contravariant 
       vectors to x' contravariant vectors;

    -- Here, dxp_dx[i][j] =  dx^{i'} / dx^{j} 

    -- pretty much calculate the elements from dx^{i} / dx^{j'} 
       as given by dx_dxp_calc(), using the following 
       identity:

 (dx^{i'} / dx^{k}) (dx^{k} / dx^{j'})  = DELTA(i',j')

    -- assumes that only coupling between  r,th and x1,x2  and r=r(t)
 
***************************************************************/
void dxp_dx_calc_dyn_rad( double *x, double *xp, double dxp_dx[][NDIM])
{

  int i, j;
  double a10, a11, a22;

  /* Start from inverse matrix:  */
  dx_dxp_calc( x, xp, dxp_dx );
  
  a10 = dxp_dx[RR][0] ;
  a11 = dxp_dx[RR][1] ;
  a22 = dxp_dx[TH][2] ; 

  /* Start from the identity matrix: */
  DLOOP2   dxp_dx[i][j] = DELTA(i,j);

 /* note we assume a00 = 1 */
  dxp_dx[1][TT] = -a10 / a11;
  dxp_dx[1][RR] =  1./a11;
  dxp_dx[2][TH] =  1./a22;

  return;

}

/***************************************************************

   dxp_dx_calc2_dyn_rad():
   --------------
    -- like dxp_dx_calc_dyn_rad() but uses dx_dxp() as an argument.
 
***************************************************************/
void dxp_dx_calc2_dyn_rad( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM])
{

  int i, j;
  double a10, a11, a22;

  a10 = dx_dxp[RR][0] ;
  a11 = dx_dxp[RR][1] ;
  a22 = dx_dxp[TH][2] ; 

  /* Start from the identity matrix: */
  DLOOP2   dxp_dx[i][j] = DELTA(i,j);

 /* note we assume a00 = 1 */
  dxp_dx[1][TT] = -a10 / a11;
  dxp_dx[1][RR] =  1./a11;
  dxp_dx[2][TH] =  1./a22;

  return;
}


/***************************************************************

   dxp_dx_calc_general():
   --------------
    -- finds the transformation array dx^{\mu'} / dx^\mu
       that is used to transform x contravariant 
       vectors to x' contravariant vectors;

    -- Here, dxp_dx[i][j] =  dx^{i'} / dx^{j} 

    -- pretty much calculate the elements from dx^{i} / dx^{j'} 
       as given by dx_dxp_calc(), using the following 
       identity:

 (dx^{i'} / dx^{k}) (dx^{k} / dx^{j'})  = DELTA(i',j')

    -- assumes nothing about the dx_dxp[][]
 
***************************************************************/
void dxp_dx_calc_general( double *x, double *xp, double dxp_dx[][NDIM])
{

  int i, j;
  double det_tmp, dx_dxp[NDIM][NDIM];

  dx_dxp_calc( x, xp, dx_dxp );
  
  if( invert_matrix2( dx_dxp, dxp_dx, &det_tmp )  ) {
    fprintf(stdout,"dxp_dx_calc_general(): singular coordinate transformation detected!!  \n"); fflush(stdout); 
    fail(FAIL_BASIC,0); 
  }
  
  return;

}

/***************************************************************

   dxp_dx_calc2_general():
   --------------
    -- like dxp_dx_calc_general() but uses dx_dxp() as an argument.
 
***************************************************************/
void dxp_dx_calc2_general( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM])
{

  int i, j;
  double det_tmp;

  if( invert_matrix2( dx_dxp, dxp_dx, &det_tmp )  ) {
    fprintf(stdout,"dxp_dx_calc2_general(): singular coordinate transformation detected!!  \n"); fflush(stdout); 
    fail(FAIL_BASIC,0); 
  }
  
  return;

}

/***************************************************************

   dxp_dx_calc3_general():
   --------------
    -- like dxp_dx_calc_general() but uses dx_dxp() as an argument.
 
***************************************************************/
void dxp_dx_calc3_general( double *x, double *xp, double dx_dxp[][NDIM], double dxp_dx[][NDIM], double *det)
{

  int i, j;
  double det_tmp;

  if( invert_matrix2( dx_dxp, dxp_dx, &det_tmp )  ) {
    fprintf(stdout,"dxp_dx_calc2_general(): singular coordinate transformation detected!!  \n"); fflush(stdout); 
    fail(FAIL_BASIC,0); 
  }
  *det = det_tmp*det_tmp;
  
  return;

}

/***************************************************************************/
/***************************************************************************
   transform_rank2cov():
  ---------------
   -- transforms any rank-2 covariant (indices down) tensor from standard 
      coordinates  (x) to numerical coordinates (xp); 
   -- see inverse_transform_rank2cov() to  transform any rank-2 covariant tensor from xp to x; 
***************************************************************************/
void transform_rank2cov(double *x, double *xp,  double vcov[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, dx_dxp[NDIM][NDIM], gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
  dx_dxp_calc( x, xp, dx_dxp );

#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcov[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcov[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dx_dxp[i][ip] * dx_dxp[j][jp];
    }
    vcov[ip][jp] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank2cov2():
  ---------------
   -- like transform_rank2cov()  but takes dx_dxp[] as an argument;
   -- transforms any rank-2 covariant (indices down) tensor from standard 
      coordinates  (x) to numerical coordinates (xp); 
   -- see inverse_transform_rank2cov() to  transform any rank-2 covariant tensor from xp to x; 
***************************************************************************/
void transform_rank2cov2( double dx_dxp[][NDIM], double vcov[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcov[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcov[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dx_dxp[i][ip] * dx_dxp[j][jp];
    }
    vcov[ip][jp] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank2cov():
  ---------------
   -- transforms any rank-2 covariant (indices down) tensor from 
      numerical coordinates (xp) to standard  coordinates  (x)
***************************************************************************/
void inverse_transform_rank2cov(double *x, double *xp,  double vcov[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, dxp_dx[NDIM][NDIM], gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
  dxp_dx_calc( x, xp, dxp_dx );

#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcov[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcov[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dxp_dx[i][ip] * dxp_dx[j][jp];
    }
    vcov[ip][jp] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank2cov2():
  ---------------
   -- like inverse_transform_rank2cov()  but uses existing dxp_dx[]
***************************************************************************/
void inverse_transform_rank2cov2(double dxp_dx[][NDIM], double vcov[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcov[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcov[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dxp_dx[i][ip] * dxp_dx[j][jp];
    }
    vcov[ip][jp] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank2con():
  ---------------
   -- transforms any rank-2 contravariant (indices up) tensor from standard 
      coordinates (x) to numerical coordinates (xp); 
   -- see inverse_transform_rank2con() to  transform any rank-2 contravariant tensor from xp to x; 
***************************************************************************/
void transform_rank2con( double *x, double *xp, double vcon[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, dxp_dx[NDIM][NDIM], gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
  dxp_dx_calc( x, xp, dxp_dx );

#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcon[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcon[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dxp_dx[ip][i] * dxp_dx[jp][j];
    }
    vcon[ip][jp] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank2con2():
  ---------------
   -- like transform_rank2con()  but takes dxp_dx[] as an argument;
***************************************************************************/
void transform_rank2con2( double dxp_dx[][NDIM], double vcon[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcon[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcon[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dxp_dx[ip][i] * dxp_dx[jp][j];
    }
    vcon[ip][jp] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank2con():
  ---------------
   -- transforms any rank-2 contravariant (indices up) tensor from 
      numerical coordinates (xp) to standard coordinates (x) ;
***************************************************************************/
void inverse_transform_rank2con( double *x, double *xp, double vcon[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, dx_dxp[NDIM][NDIM], gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
  dx_dxp_calc( x, xp, dx_dxp );

#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcon[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcon[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dx_dxp[ip][i] * dx_dxp[jp][j];
    }
    vcon[ip][jp] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank2con2():
  ---------------
   -- like    inverse_transform_rank2con() but uses existing dx_dxp[]
***************************************************************************/
void inverse_transform_rank2con2( double dx_dxp[][NDIM], double vcon[][NDIM])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM][NDIM];
  
  /* Transform into numerical coordinates : */ 
#if(USE_STRICT_ARRAY_BOUNDS)
  DLOOP2  gtmp[i][j] = vcon[i][j] ;
#else
  for(i = 0; i < NDIM*NDIM; i++) gtmp[0][i] = vcon[0][i] ;
#endif

  for( ip = 0 ; ip < NDIM; ip++ ) for( jp = 0 ; jp < NDIM; jp++ ) { 
    tmpout = 0.;
    DLOOP2 {
      tmpout +=  gtmp[i][j] * dx_dxp[ip][i] * dx_dxp[jp][j];
    }
    vcon[ip][jp] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank1cov():
  ---------------
   -- transforms any rank-1 covariant (indices down) tensor from standard 
      coordinates (x) to numerical coordinates (xp); 
   -- see inverse_transform_rank1cov to transform any rank-1 covariant tensor from xp to x; 
***************************************************************************/
void transform_rank1cov(double *x, double *xp,  double vcov[])
{
  int i,j,ip,jp;
  double tmpout, dx_dxp[NDIM][NDIM], gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  dx_dxp_calc( x, xp, dx_dxp );

  DLOOP1 gtmp[i] = vcov[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dx_dxp[i][ip];
    }
    vcov[ip] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank1cov2():
  ---------------
   -- like    transform_rank1cov() but uses existing dx_dxp[]
***************************************************************************/
void transform_rank1cov2(double dx_dxp[][NDIM], double vcov[])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  DLOOP1 gtmp[i] = vcov[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dx_dxp[i][ip];
    }
    vcov[ip] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank1cov():
  ---------------
   -- transforms any rank-1 covariant (indices down) tensor from 
      numerical coordinates (xp) to standard coordinates (x);
***************************************************************************/
void inverse_transform_rank1cov(double *x, double *xp,  double vcov[])
{
  int i,j,ip,jp;
  double tmpout, dxp_dx[NDIM][NDIM], gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  dxp_dx_calc( x, xp, dxp_dx );

  DLOOP1 gtmp[i] = vcov[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dxp_dx[i][ip];
    }
    vcov[ip] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank1cov2():
  ---------------
   -- like    inverse_transform_rank1cov() but uses existing dxp_dx[]
***************************************************************************/
void inverse_transform_rank1cov2(double dxp_dx[][NDIM], double vcov[])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  DLOOP1 gtmp[i] = vcov[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dxp_dx[i][ip];
    }
    vcov[ip] = tmpout;
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank1con():
  ---------------
   -- transforms any rank-1 contravariant (indices up) tensor from standard 
      coordinates (x) to numerical coordinates (xp); 
   -- see inverse_transform_rank1con to transform any rank-1 contravariant tensor from xp to x; 
***************************************************************************/
void transform_rank1con( double *x, double *xp, double vcon[])
{
  int i,j,ip,jp;
  double tmpout, dxp_dx[NDIM][NDIM], gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  dxp_dx_calc( x, xp, dxp_dx );

  DLOOP1 gtmp[i] = vcon[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dxp_dx[ip][i] ;
    }
    vcon[ip] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_rank1con2():
  ---------------
   -- like   transform_rank1con() but uses dxp_dx
***************************************************************************/
void transform_rank1con2( double dxp_dx[][NDIM], double vcon[])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  DLOOP1 gtmp[i] = vcon[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dxp_dx[ip][i] ;
    }
    vcon[ip] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank1con():
  ---------------
   -- transforms any rank-1 contravariant (indices up) tensor from 
      numerical coordinates (xp) to standard coordinates (x) ; 
***************************************************************************/
void inverse_transform_rank1con( double *x, double *xp, double vcon[])
{
  int i,j,ip,jp;
  double tmpout, dx_dxp[NDIM][NDIM], gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  dx_dxp_calc( x, xp, dx_dxp );

  DLOOP1 gtmp[i] = vcon[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dx_dxp[ip][i] ;
    }
    vcon[ip] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   inverse_transform_rank1con2():
  ---------------
   -- like    inverse_transform_rank1con() but uses existing dx_dxp[]
***************************************************************************/
void inverse_transform_rank1con2(double dx_dxp[NDIM][NDIM],  double vcon[])
{
  int i,j,ip,jp;
  double tmpout, gtmp[NDIM];
  
  /* Transform into numerical coordinates : */ 
  DLOOP1 gtmp[i] = vcon[i] ;

  for( ip = 0 ; ip < NDIM; ip++ ) {
    tmpout = 0.;
    DLOOP1 {
      tmpout +=  gtmp[i] * dx_dxp[ip][i] ;
    }
    vcon[ip] = tmpout; 
  }

  return;
}

/***************************************************************************/
/***************************************************************************
   transform_connection():
  ---------------
   -- transforms a affine connection (aka  "\Gamma^i_j_k") from non-primed 
      coordinates (x)  to primed/numerical coordinates (xp).
   -- nb.: the connection does not transform like a rank-3 tensor; 
***************************************************************************/
void transform_connection( double *x, double *xp, double conn[][NDIM][NDIM], double ***connp)
{
  int i,j,k, ii,jj,kk;
  double dx_dxp[NDIM][NDIM],dxp_dx[NDIM][NDIM],dx_dxp_dxp[NDIM][NDIM][NDIM];

  dx_dxp_calc( x, xp, dx_dxp); 
  dxp_dx_calc2( x, xp, dx_dxp, dxp_dx); 
  dx_dxp_dxp_calc( x, xp, dx_dxp_dxp) ; 

  DLOOP3 { connp[i][j][k] = 0. ;  } 

  /* This is the extra part that makes it not your standard rank-3 tensor transform : */
  DLOOP3  for(ii=0;ii<NDIM;ii++) { 
    connp[i][j][k] += dxp_dx[i][ii] * dx_dxp_dxp[ii][j][k] ;
  }

  /* This is the equivalent to a rank-3 tensor transformation: */
  DLOOP3  for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
    connp[i][j][k] += dxp_dx[i][ii] * dx_dxp[jj][j] * dx_dxp[kk][k] * conn[ii][jj][kk] ; 
  }
  
  return;
}

/***************************************************************************/
/***************************************************************************
   transform_connection2():
  ---------------
   -- like transform_connection() but uses of_coord structure
   -- nb.: the connection does not transform like a rank-3 tensor; 
***************************************************************************/
void transform_connection2( struct of_coord *coords, double conn[][NDIM][NDIM], double ***connp)
{
  int i,j,k, ii,jj,kk;
  double dx_dxp_dxp[NDIM][NDIM][NDIM];

  dx_dxp_dxp_calc( coords->x, coords->xp, dx_dxp_dxp) ; 

  DLOOP3 { connp[i][j][k] = 0. ;  } 

  /* This is the extra part that makes it not your standard rank-3 tensor transform : */
  DLOOP3  for(ii=0;ii<NDIM;ii++) { 
    connp[i][j][k] += coords->dxp_dx[i][ii] * dx_dxp_dxp[ii][j][k] ;
  }

  /* This is the equivalent to a rank-3 tensor transformation: */
  DLOOP3  for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) for(kk=0;kk<NDIM;kk++) { 
    connp[i][j][k] += coords->dxp_dx[i][ii] * coords->dx_dxp[jj][j] * coords->dx_dxp[kk][k] * conn[ii][jj][kk] ; 
  }
  
  return;
}

/***************************************************************************/
/***************************************************************************
   transform_all():
  ---------------
   -- uses of_coord and of_geom structures to transform all metric functions to xp coordinates;
***************************************************************************/
void transform_all(struct of_coord *coords, struct of_geom *geom)
{

  /* gcov: */
  transform_rank2cov2( coords->dx_dxp, geom->gcov);

  /* gcon: */
  transform_rank2con2( coords->dxp_dx, geom->gcon);

  /* gdet */
  geom->g *= coords->det_dx_dxp;

  return;
}

/***************************************************************************/
/***************************************************************************
   regularize_rank1con:
  ---------------
   -- regularizes a contravariant vector that has already been 
      transformed to "x" coordinates ; 
   -- regularization means normalization of the vector, so that it is now 
      a linear combination of the cartesian components, e.g. factors of 
      r and sin/cos(th);
   -- assumes that only the spatial components need to be regularized
***************************************************************************/
void regularize_rank1con( double *x, double vcon[NDIM-1])
{

#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
  double sth; 

  sth = fabs(sin(x[TH]));
  vcon[TH-1] *= x[RR]/sth; 
  vcon[PH-1] *= x[RR]*sth;

#else
  fprintf(stderr,
	  "regularize_rank1con():  Need to write option for this value of TOP_TYPE_CHOICE = %d \n",
	  TOP_TYPE_CHOICE); 
  fflush(stderr); fail(FAIL_BASIC,0);

#endif

  return;
}

/***************************************************************************/
/***************************************************************************
   unregularize_rank1con:
  ---------------
   -- performs inverse of regularization a contravariant vector that 
      has already been transformed to "x" coordinates ; 
***************************************************************************/
void unregularize_rank1con( double *x, double vcon[NDIM-1])
{
  double inv[NDIM-1]={1.,1.,1.};
  
  regularize_rank1con( x, inv);
  vcon[0] /= inv[0];  vcon[1] /= inv[1];  vcon[2] /= inv[2];

  return;
}


/**************************************************************************************/
/**************************************************************************************
  regularize_prim();
  --------------
   -- for a single set of primitive variables, 
      first transforms vectors to "x" coordinates then 
      scales the vector components of the primitive variables 
      to eliminate any geometric factors;

**************************************************************************************/
void regularize_prim( double prim[NP], double reg[NDIM-1], double dxdxp[NDIM][NDIM] )
{
  register double rf1;

#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  /* Just regularize */
  prim[U1] *= reg[0];   prim[U2] *= reg[1];   prim[U3] *= reg[2]; 
  prim[B1] *= reg[0];   prim[B2] *= reg[1];   prim[B3] *= reg[2]; 

#elif( COORD_TYPE_CHOICE == COORD_DIAGONAL || COORD_TYPE_CHOICE == COORD_DIAGONAL2  || COORD_TYPE_CHOICE == COORD_DIAGONAL3 || COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD)
  /* Transform to x coordinates and then regularize */
  /* Only need to take into account diagonal factors in coordinate transformation: */
  /* We assume also that  regularize_gf[0] = dx_dxp[3][3] = 1. :  */
  rf1 = reg[1] * dxdxp[TH][2];
  prim[U1] *= dxdxp[RR][1];       prim[U2] *= rf1;       prim[U3] *= reg[2]; 
  prim[B1] *= dxdxp[RR][1];       prim[B2] *= rf1;       prim[B3] *= reg[2]; 

#elif( COORD_TYPE_CHOICE == COORD_MIXED )
  /* Transform to x coordinates and then regularize */
  /* Need to consider off-diagonal elements between RR/1 and TH/2 */
  /* We assume also that  regularize_gf[0] = dx_dxp[3][3] = 1. :  */
  prim[U1]  =         dxdxp[RR][1]*prim[U1]  +  dxdxp[RR][2]*prim[U2] ;
  prim[U2]  = reg[1]*(dxdxp[TH][1]*prim[U1]  +  dxdxp[TH][2]*prim[U2]); 
  prim[U3] *= reg[2]; 
  prim[B1]  =         dxdxp[RR][1]*prim[B1]  +  dxdxp[RR][2]*prim[B2] ; 
  prim[B2]  = reg[1]*(dxdxp[TH][1]*prim[B1]  +  dxdxp[TH][2]*prim[B2]); 
  prim[B3] *= reg[2]; 

#else 
  fprintf(stderr,"regularize_prim(): COORD_TYPE_CHOICE = %d  not supported! \n",COORD_TYPE_CHOICE);
  fflush(stderr);  fail(FAIL_BASIC,0); 
#endif


  return;
}


/**************************************************************************************/
/**************************************************************************************
  unregularize_prim();
  --------------
   -- preforms inverse operation of regularize_prim();
**************************************************************************************/
void unregularize_prim( double prim[NP], double unreg[NDIM-1], double dxpdx[NDIM][NDIM] )
{
  register double rf1;


#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  prim[U1] *= unreg[0];   prim[U2] *= unreg[1];   prim[U3] *= unreg[2]; 
  prim[B1] *= unreg[0];   prim[B2] *= unreg[1];   prim[B3] *= unreg[2]; 

#elif( COORD_TYPE_CHOICE == COORD_DIAGONAL || COORD_TYPE_CHOICE == COORD_DIAGONAL2  || COORD_TYPE_CHOICE == COORD_DIAGONAL3 || COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD)
  rf1 = unreg[1] * dxpdx[2][TH];
  prim[U1] *= dxpdx[1][RR];       prim[U2] *= rf1;       prim[U3] *= unreg[2]; 
  prim[B1] *= dxpdx[1][RR];       prim[B2] *= rf1;       prim[B3] *= unreg[2]; 

#elif( COORD_TYPE_CHOICE == COORD_MIXED )
  prim[U1]  =           dxpdx[1][RR]*prim[U1]  +  dxpdx[1][TH]*prim[U2] ;
  prim[U2]  = unreg[1]*(dxpdx[2][RR]*prim[U1]  +  dxpdx[2][TH]*prim[U2]); 
  prim[U3] *= unreg[2]; 
  prim[B1]  =           dxpdx[1][RR]*prim[B1]  +  dxpdx[1][TH]*prim[B2] ; 
  prim[B2]  = unreg[1]*(dxpdx[2][RR]*prim[B1]  +  dxpdx[2][TH]*prim[B2]); 
  prim[B3] *= unreg[2]; 

#else 
  fprintf(stderr,"unregularize_prim(): COORD_TYPE_CHOICE = %d  not supported! \n",COORD_TYPE_CHOICE);
  fflush(stderr);  fail(FAIL_BASIC,0); 
#endif


  return;
}


/**************************************************************************************/
/**************************************************************************************
  unregularize_prim_LR();
  --------------
   -- like unregularize_prim() but for left and right states of the primitive variables; 
**************************************************************************************/
void unregularize_prim_LR( double prim_L[NP], double prim_R[NP], 
			   double unreg[NDIM-1], double dxpdx[NDIM][NDIM] )
{
  register double rf1;


#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
  prim_L[U1] *= unreg[0];   prim_L[U2] *= unreg[1];   prim_L[U3] *= unreg[2]; 
  prim_L[B1] *= unreg[0];   prim_L[B2] *= unreg[1];   prim_L[B3] *= unreg[2]; 
  prim_R[U1] *= unreg[0];   prim_R[U2] *= unreg[1];   prim_R[U3] *= unreg[2]; 
  prim_R[B1] *= unreg[0];   prim_R[B2] *= unreg[1];   prim_R[B3] *= unreg[2]; 

#elif( COORD_TYPE_CHOICE == COORD_DIAGONAL || COORD_TYPE_CHOICE == COORD_DIAGONAL2  || COORD_TYPE_CHOICE == COORD_DIAGONAL3 || COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD)
  rf1 = unreg[1] * dxpdx[2][TH];
  prim_L[U1] *= dxpdx[1][RR];       prim_L[U2] *= rf1;       prim_L[U3] *= unreg[2]; 
  prim_L[B1] *= dxpdx[1][RR];       prim_L[B2] *= rf1;       prim_L[B3] *= unreg[2]; 
  prim_R[U1] *= dxpdx[1][RR];       prim_R[U2] *= rf1;       prim_R[U3] *= unreg[2]; 
  prim_R[B1] *= dxpdx[1][RR];       prim_R[B2] *= rf1;       prim_R[B3] *= unreg[2]; 

#elif( COORD_TYPE_CHOICE == COORD_MIXED )
  prim_L[U1]  =           dxpdx[1][RR]*prim_L[U1]  +  dxpdx[1][TH]*prim_L[U2] ;
  prim_L[U2]  = unreg[1]*(dxpdx[2][RR]*prim_L[U1]  +  dxpdx[2][TH]*prim_L[U2]); 
  prim_L[U3] *= unreg[2]; 
  prim_L[B1]  =           dxpdx[1][RR]*prim_L[B1]  +  dxpdx[1][TH]*prim_L[B2] ; 
  prim_L[B2]  = unreg[1]*(dxpdx[2][RR]*prim_L[B1]  +  dxpdx[2][TH]*prim_L[B2]); 
  prim_L[B3] *= unreg[2]; 

  prim_R[U1]  =           dxpdx[1][RR]*prim_R[U1]  +  dxpdx[1][TH]*prim_R[U2] ;
  prim_R[U2]  = unreg[1]*(dxpdx[2][RR]*prim_R[U1]  +  dxpdx[2][TH]*prim_R[U2]); 
  prim_R[U3] *= unreg[2]; 
  prim_R[B1]  =           dxpdx[1][RR]*prim_R[B1]  +  dxpdx[1][TH]*prim_R[B2] ; 
  prim_R[B2]  = unreg[1]*(dxpdx[2][RR]*prim_R[B1]  +  dxpdx[2][TH]*prim_R[B2]); 
  prim_R[B3] *= unreg[2]; 

#else 
  fprintf(stderr,"unregularize_prim_LR(): COORD_TYPE_CHOICE = %d  not supported! \n",COORD_TYPE_CHOICE);
  fflush(stderr);  fail(FAIL_BASIC,0); 
#endif


  return;
}



/***********************************************************************/
/***********************************************************************
  Rin_calc(): 
  --------------------
   -- returns the innermost radial coordinate assuming that we want 
       "n_within_horizon" cells within the horizon (excluding ghost cells)
       and assuming that the cell boundary cannot be on the horizon;

   -- (OLD WAY) make first super-horizon cell's boundary be within SMALL of horizon; 
   -- (NEW WAY) make sure horizon falls in between a cell border and its center 
       so that gcon[0][0] is always non-zero ; 
   -- assumes that we are using logarithmically spaced radial coordinates;
   -- using similar HARM2D method;
***********************************************************************/
double Rin_calc(void)
{
  int n ; 
  double x0,fn;
  
//-old-way   n = n_within_horizon;
//-old-way   x0 = ( totalsize[1]*log(r_horizon - SMALL - R0) - n*log(Rout-R0) ) / (totalsize[1]-n) ; 

  fn = ((double) n_within_horizon) - 0.25 ;
  x0 = ( totalsize[1]*log(r_horizon - SMALL - R0) - fn*log(Rout-R0) ) / (totalsize[1]-fn) ; 
  
//  fprintf(stdout,"x0 = %28.18e \n", x0); 
//  fflush(stdout);

  return( R0 + exp(x0) );
}


/***********************************************************************/
/***********************************************************************
  find_min_dx():
  --------------------
   -- returns the smallest real length scale in non-numerical coordinates;
***********************************************************************/
double find_min_dx(void)
{
  int i,j,k,l;
  double *xp_l, *xp_h, *x_l, *x_h;
  double dxcell[NDIM], dxmin;
  struct of_coord *coords_l, *coords_h;


  /* Now determine the smallest value of dx : */ 
  DLOOP1 dxcell[i] = 0.;
  dxmin   = 1.e200;

#if(USE_MASK)
  LOOP if( evol_mask[ncurr][i][j][k] != MASK_EXCISED ) { 
#else 
  LOOP { 
#endif

    get_coord(i  ,j  ,k  ,CORN,ncurr,coords_l);
    get_coord(i+1,j+1,k+1,CORN,ncurr,coords_h);

    xp_l = coords_l->xp;
     x_l = coords_l->x;
    xp_h = coords_h->xp;
     x_h = coords_h->x;

#if( TOP_TYPE_CHOICE == TOP_CARTESIAN )
    dxcell[1] = x_h[1] - x_l[1]; 
    dxcell[2] = x_h[2] - x_l[2]; 
    dxcell[3] = x_h[3] - x_l[3]; 
#elif( TOP_TYPE_CHOICE == TOP_SPHERICAL )
    dxcell[1] = x_h[1] - x_l[1];
    dxcell[2] = 0.5*(x_h[1]+x_l[1]) * (x_h[2] - x_l[2]);
    dxcell[3] = 0.5*( (x_h[1]*sin(x_h[2])) + (x_l[1]*sin(x_l[2])) ) * (x_h[3] - x_l[3]);
#else
    dxcell[1] = x_h[1] - x_l[1];
    dxcell[2] = 0.5*(x_h[1]+x_l[1]) * (x_h[2] - x_l[2]);
    dxcell[3] = x_h[3] - x_l[3];
#endif

    dxmin = MIN( dxmin, MIN( dxcell[1], MIN(dxcell[2], dxcell[3])) ) ; 

  }
  
  /* Get the smallest dx over all domains: */
  mpi_global_min( &dxmin );

  if( myid == printer_pid ) {  fprintf(stdout,"\ndxmin = %28.18e \n", dxmin); fflush(stdout); }

  return( dxmin ) ; 
}

/***********************************************************************/
/***********************************************************************
  find_min_dt():
  --------------------
   -- returns the light crossing time of the smallest cell dimension
      in the grid; 
   -- also checks to see if th_cutout,th_beg,th_end are ok; 
***********************************************************************/
double find_min_dt(void)
{
  int i,j,k,l;
  double dtmin[NDIM], f1, f2, r_hor,inv_fmax,inv_fmin;
  struct of_geom *geom;

  TRACE_BEG;

  /* Now determine the smallest value of dx : */ 
  dtmin[0] = 1.e200; 

#if(USE_MASK)
  LOOP if( evol_mask[ncurr][i][j][k] != MASK_EXCISED ) { 
#else 
  LOOP { 
#endif

    get_geometry(i,j,k,CENT,ncurr,geom);

    for( l = 1; l < NDIM; l++ ) { 
      f2 = sqrt(geom->gcon[0][l]*geom->gcon[0][l] - geom->gcon[0][0]*geom->gcon[l][l]);
      inv_fmin = fabs(geom->gcon[0][0]/(geom->gcon[0][l] + f2)) ;
      inv_fmax = fabs(geom->gcon[0][0]/(geom->gcon[0][l] - f2)) ;
      dtmin[l] = MIN( inv_fmin , inv_fmax ); 
      dtmin[l] *= dx[l];
      if( dtmin[l] < dtmin[0] ) { dtmin[0] = dtmin[l] ; }
    }
  }
  
  /* Get the smallest dx over all domains: */
  dtmin_glob = dtmin[0] ; 
  mpi_global_min( &dtmin_glob );

  if( myid == printer_pid ) { 
    fprintf(stdout,"\ndtmin = %28.18e \n", dtmin_glob); fflush(stdout);
  }

  TRACE_END;

  return( dtmin_glob ) ; 
}

/***********************************************************************/
/***********************************************************************
  dump_min_dt():
  --------------------
   -- wrapper for find_min_dt();
***********************************************************************/
void dump_min_dt(void)
{
  TRACE_BEG;

  find_min_dt();

  TRACE_END;

  return;
}

/***********************************************************************/
/***********************************************************************
  check_cutout():
  --------------------
   -- after all the grid specifications have been set, this routine is called
      to insure that the ghost zones are within regular bounds on theta; 
***********************************************************************/
void check_cutout(void)
{
  int i,j,k,l;
  struct of_coord *coords;

  if( th_beg <= 1.e2*SMALL ) { return; }

  /* lower bounds check: */ 
  j = 0; 
  N1ALL_LOOP N3ALL_LOOP { 
    get_coord(i,j,k,CORN,ncurr,coords);
    if( coords->x[2] < 0. ) {
      fprintf(stderr,"check_cutout(): need to increase th_beg = %g : th = %28.18e \n", 
	      th_beg, coords->x[2] );
      fflush(stderr);
      myexit(102);  /* not using fail() on purpose here */
    }
  }

  /* upper bounds check: */ 
  j = N2TOT-1; 
  N1ALL_LOOP N3ALL_LOOP { 
    get_coord(i,j,k,CORN,ncurr,coords);
    if( coords->x[2] > M_PI ) {
      fprintf(stderr,"check_cutout(): Need to decrease th_end  = %g : th = %28.18e  xp2 = %28.18e    pi = %28.18e , loc = [%6d,%6d,%6d,%6d]\n", 
	      th_end, coords->x[2], coords->xp[2], M_PI, ncurr,i,j,k );
      fflush(stderr);
      myexit(102);    /* not using fail() on purpose here */
    }
  }

  return;
}

/****************************************************************************************

 check_coords():
 ----------
   -- verifies that the parameters used to setup a coordinate system are valid; 

****************************************************************************************/
void check_coords( void )
{

#if(   COORD_TYPE_CHOICE  == COORD_DIAGONAL          ) 
  /* no tests */
#elif( COORD_TYPE_CHOICE  == COORD_MIXED             ) 
  /* no tests */
#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL2         ) 
  /* no tests */
#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3         ) 
  /* no tests */
#elif( COORD_TYPE_CHOICE  == COORD_IDENTITY          ) 
  /* no tests */
#elif( COORD_TYPE_CHOICE  == COORD_DIAGONAL3_DYN_RAD )

  if( (coord_params->fa_r < 0.) || (coord_params->fa_r > 1.) ) { fprintf(stderr,"check_coords(): fa_r must be within [0,1] : %26.16e  \n", coord_params->fa_r);  fflush(stderr); fail(FAIL_BASIC,0); }
  if( (coord_params->f0_r < 0.) || (coord_params->f0_r > 1.) ) { fprintf(stderr,"check_coords(): f0_r must be within [0,1] : %26.16e  \n", coord_params->f0_r);  fflush(stderr); fail(FAIL_BASIC,0); }
  if( coord_params->upsilon_r < 1.      ) { fprintf(stderr,"check_coords(): upsilon_r must be greater than 1 : %26.16e  \n", coord_params->upsilon_r);  fflush(stderr); fail(FAIL_BASIC,0); }

#elif( COORD_TYPE_CHOICE  == COORD_WARPED_SPHERICAL )
  // TODO

#elif( COORD_TYPE_CHOICE  == COORD_WARPED_CARTESIAN )
  // TODO

#else 
  fprintf(stderr,"check_coords(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif


#if(TOP_TYPE_CHOICE==TOP_SPHERICAL)
  check_cutout();
#endif

  return;
}


/***********************************************************************/
/***********************************************************************
  xcart_of_xspher():
  --------------------
   -- calculates Cartesian coordinates from spherical coordinates;
   -- also calculates the covariant transformation matrix;
***********************************************************************/
void xcart_of_xspher(double *xcart, double *xspher, double dxc_dxs[][NDIM] )
{
  int i; 
  double r, th, ph;
  double x,y;
  double sth,cth,sph,cph;
  
  
#if(USE_STRICT_ARRAY_BOUNDS)
  int j; 
  DLOOP2  {  dxc_dxs[i][j] = 0. ; } 
#else
  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  dxc_dxs[0][i] = 0. ; } 
#endif
  
  r  = xspher[RR]; 
  th = xspher[TH];
  ph = xspher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  x = r * sth * cph; 
  y = r * sth * sph; 

  xcart[TT] = xspher[TT];
  xcart[XX] = x;
  xcart[YY] = y;
  xcart[ZZ] = r * cth;

  dxc_dxs[TT][TT] = 1.                    ;  // dt/dt
  dxc_dxs[XX][RR] =      cph * sth        ;  // dx/dr
  dxc_dxs[XX][TH] =  r * cph * cth        ;  // dx/dtheta
  dxc_dxs[XX][PH] = -y                    ;  // dx/dphi
  dxc_dxs[YY][RR] =      sph * sth        ;  // dy/dr
  dxc_dxs[YY][TH] =  r * sph * cth        ;  // dy/dtheta
  dxc_dxs[YY][PH] = x                     ;  // dy/dphi
  dxc_dxs[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dxs[ZZ][TH] = -r*sth                ;  // dz/dtheta
  //  dxc_dxs[ZZ][PH] = 0.                    ;  // dz/dphi
  


  return;
}


/***********************************************************************/
/***********************************************************************
  xcart_of_xspher_only():
  --------------------
   -- calculates Cartesian coordinates from spherical coordinates;
   -- does not return dxc_dxs;
***********************************************************************/
void xcart_of_xspher_only(double *xcart, double *xspher)
{
  double r, th, ph;
  double x,y;
  double sth,cth,sph,cph;
  
  r  = xspher[RR]; 
  th = xspher[TH];
  ph = xspher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  x = r * sth * cph; 
  y = r * sth * sph; 

  xcart[TT] = xspher[TT];
  xcart[XX] = x;
  xcart[YY] = y;
  xcart[ZZ] = r * cth;

  return;
}


/***********************************************************************/
/***********************************************************************
  xspher_of_xcart_only():
  --------------------
   -- calculates spherical coordinates from Cartesian coordinates;
   -- does not return dxs_dxc;
***********************************************************************/
void xspher_of_xcart_only(double *xspher, double *xcart)
{

  double x = xcart[XX];
  double y = xcart[YY];
  double z = xcart[ZZ];
  
  double r  = sqrt(x*x + y*y + z*z);
  double th =  acos( z / r);
  double ph = atan2( y , x );
  if( ph < 0. ) { ph += 2.*M_PI; }

  xspher[TT] = xcart[TT];
  xspher[RR] = r;
  xspher[TH] = th;
  xspher[PH] = ph;

  return;
}


/***********************************************************************/
/***********************************************************************
  xspher_of_xcart():
  --------------------
   -- calculates spherical coordinates from Cartesian coordinates;
   -- returns with dxs_dxc too;
***********************************************************************/
void xspher_of_xcart(double *xspher, double *xcart, double dxs_dxc[][NDIM] )
{
  int i,j; 

  double x = xcart[XX];
  double y = xcart[YY];
  double z = xcart[ZZ];
  
  double rho_sq     = x*x + y*y;
  double rho        = sqrt(rho_sq); 
  double rho_inv    = 1./rho;
  double rho_inv_sq = rho_inv * rho_inv;
  double r_sq       = rho_sq + z*z;
  double r          = sqrt(r_sq);
  double r_inv      = 1./r;
  double r_inv_sq   = r_inv * r_inv;
  
  double th =  acos( z * r_inv );
  double ph = atan2( y , x );
  if( ph < 0. ) { ph += 2.*M_PI; }

  xspher[TT] = xcart[TT];
  xspher[RR] = r;
  xspher[TH] = th;
  xspher[PH] = ph;

  double dth_term = z * r_inv_sq * rho_inv ;

  /* Note time-dependence is implicitly left to the  x(xp) transformation, so all time derivatives here are zero:  */
  dxs_dxc[TT][TT] = 1.               ;  // dt/dt
  dxs_dxc[TT][XX] = 0.               ;  // dt/dx
  dxs_dxc[TT][YY] = 0.               ;  // dt/dy
  dxs_dxc[TT][ZZ] = 0.               ;  // dt/dz
  dxs_dxc[RR][TT] = 0.               ;  // dr/dt
  dxs_dxc[RR][XX] = x * r_inv        ;  // dr/dx
  dxs_dxc[RR][YY] = y * r_inv        ;  // dr/dy
  dxs_dxc[RR][ZZ] = z * r_inv        ;  // dr    /dz
  dxs_dxc[TH][TT] = 0.               ;  // dth/dt
  dxs_dxc[TH][XX] = x * dth_term     ;  // dtheta/dx
  dxs_dxc[TH][YY] = y * dth_term     ;  // dtheta/dy
  dxs_dxc[TH][ZZ] = -rho * r_inv_sq  ;  // dtheta/dz  
  dxs_dxc[PH][TT] = 0.               ;  // dphi/dt
  dxs_dxc[PH][XX] = -y*rho_inv_sq    ;  // dphi/dx
  dxs_dxc[PH][YY] =  x*rho_inv_sq    ;  // dphi/dy
  dxs_dxc[PH][ZZ] = 0.               ;  // dphi/dz

  return;
}


/***********************************************************************/
/***********************************************************************
  xspher_of_xcart_special3():
  --------------------
   -- named "special3" to match similar name of xcart_of_xspher_special3();
   -- calculates spherical coordinates from Cartesian coordinates and the 
        transformation matrix from numerical-Cartesian coordinates to spherical coordinates;
   -- returns with dxs_dxp too;
***********************************************************************/
void xspher_of_xcart_special3(double *xspher, double *xcart, double *xp, double dxs_dxp[][NDIM] )
{
  int i,j; 

  double x = xcart[XX];
  double y = xcart[YY];
  double z = xcart[ZZ];
  
  double rho_sq     = x*x + y*y;
  double rho        = sqrt(rho_sq); 
  double rho_inv    = 1./rho;
  double rho_inv_sq = rho_inv * rho_inv;
  double r_sq       = rho_sq + z*z;
  double r          = sqrt(r_sq);
  double r_inv      = 1./r;
  double r_inv_sq   = r_inv * r_inv;
  
  double th =  acos( z * r_inv );
  double ph = atan2( y , x );
  if( ph < 0. ) { ph += 2.*M_PI; }

  xspher[TT] = xcart[TT];
  xspher[RR] = r;
  xspher[TH] = th;
  xspher[PH] = ph;

  double dth_term = z * r_inv_sq * rho_inv ;

#if( COORD_TYPE_CHOICE==COORD_IDENTITY  )

  /* Note time-dependence is implicitly left to the  x(xp) transformation, so all time derivatives here are zero:  */
  dxs_dxp[TT][TT] = 1.               ;  // dt/dt
  dxs_dxp[TT][XX] = 0.               ;  // dt/dx
  dxs_dxp[TT][YY] = 0.               ;  // dt/dy
  dxs_dxp[TT][ZZ] = 0.               ;  // dt/dz
  dxs_dxp[RR][TT] = 0.               ;  // dr/dt
  dxs_dxp[RR][XX] = x * r_inv        ;  // dr/dx
  dxs_dxp[RR][YY] = y * r_inv        ;  // dr/dy
  dxs_dxp[RR][ZZ] = z * r_inv        ;  // dr    /dz
  dxs_dxp[TH][TT] = 0.               ;  // dth/dt
  dxs_dxp[TH][XX] = x * dth_term     ;  // dtheta/dx
  dxs_dxp[TH][YY] = y * dth_term     ;  // dtheta/dy
  dxs_dxp[TH][ZZ] = -rho * r_inv_sq  ;  // dtheta/dz  
  dxs_dxp[PH][TT] = 0.               ;  // dphi/dt
  dxs_dxp[PH][XX] = -y*rho_inv_sq    ;  // dphi/dx
  dxs_dxp[PH][YY] =  x*rho_inv_sq    ;  // dphi/dy
  dxs_dxp[PH][ZZ] = 0.               ;  // dphi/dz

#else

  fprintf(stderr,"xspher_of_xcart_special3(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
  //  double dx_dxp[NDIM][NDIM];
  //  dx_dxp_calc( xcart, xp,  dx_dxp );

#endif

  return;
}


/***********************************************************************/
/***********************************************************************
  xcart_of_xspher_special():
  --------------------
   -- calculates Cartesian coordinates from spherical coordinates;
   -- also calculates the covariant transformation matrix to xp ;
***********************************************************************/
void xcart_of_xspher_special(struct of_coord *coords )
{
  int i; 
  double r, th, ph;
  double x,y;
  double sth,cth,sph,cph;
  
#if( COORD_TYPE_CHOICE==COORD_DIAGONAL          ||\
     COORD_TYPE_CHOICE==COORD_IDENTITY          ||\
     COORD_TYPE_CHOICE==COORD_DIAGONAL2         ||\
     COORD_TYPE_CHOICE==COORD_DIAGONAL3         ||\
     COORD_TYPE_CHOICE==COORD_MIXED             ||\
     COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD ||\
     COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL   )

  r  = coords->x[RR]; 
  th = coords->x[TH];
  ph = coords->x[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  x = r * sth * cph; 
  y = r * sth * sph; 

  coords->xcart[TT] = coords->x[TT];
  coords->xcart[XX] = x;
  coords->xcart[YY] = y;
  coords->xcart[ZZ] = r * cth;

# if( COORD_TYPE_CHOICE == COORD_MIXED )
  should-test-this-somehow

  coords->dxc_dxp[TT][ 0] = 1.                                                                               ;  // dt/dt
  coords->dxc_dxp[TT][ 1] = 0.                                                                               ;  // dt/dx1
  coords->dxc_dxp[TT][ 2] = 0.                                                                               ;  // dt/dx2
  coords->dxc_dxp[TT][ 3] = 0.                                                                               ;  // dt/dx3
  coords->dxc_dxp[XX][ 0] = 0.                                                                               ;  // dx/dt
  coords->dxc_dxp[XX][ 1] =      cph * sth  * coords->dx_dxp[RR][1] + r * cph * cth  * coords->dx_dxp[TH][1] ;  // dx/dx1
  coords->dxc_dxp[XX][ 2] =  r * cph * cth  * coords->dx_dxp[TH][2] +     cph * sth  * coords->dx_dxp[RR][2] ;  // dx/dx2
  coords->dxc_dxp[XX][ 3] = -y                                                                               ;  // dx/dx3
  coords->dxc_dxp[YY][ 0] = 0.                                                                               ;  // dy/dt
  coords->dxc_dxp[YY][ 1] =      sph * sth  * coords->dx_dxp[RR][1] + r * sph * cth  * coords->dx_dxp[TH][1] ;  // dy/dx1
  coords->dxc_dxp[YY][ 2] =  r * sph * cth  * coords->dx_dxp[TH][2] +     sph * sth  * coords->dx_dxp[RR][2] ;  // dy/dx2
  coords->dxc_dxp[YY][ 3] = x                                                                                ;  // dy/dx3
  coords->dxc_dxp[ZZ][ 0] = 0.                                                                               ;  // dz/dt
  coords->dxc_dxp[ZZ][ 1] = cth             * coords->dx_dxp[RR][1] - r*sth * coords->dx_dxp[TH][1]          ;  // dz/dx1
  coords->dxc_dxp[ZZ][ 2] = -r*sth          * coords->dx_dxp[TH][2] +   cth * coords->dx_dxp[RR][2]          ;  // dz/dx2
  coords->dxc_dxp[ZZ][ 3] = 0.                                                                               ;  // dz/dx3

# elif( COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD )
  coords->dxc_dxp[TT][ 0] = 1.                                        ;  // dt/dt
  coords->dxc_dxp[TT][ 1] = 0.                                        ;  // dt/dx1
  coords->dxc_dxp[TT][ 2] = 0.                                        ;  // dt/dx2
  coords->dxc_dxp[TT][ 3] = 0.                                        ;  // dt/dx3
  coords->dxc_dxp[XX][ 0] =      cph * sth  * coords->dx_dxp[RR][0]   ;  // dx/dt
  coords->dxc_dxp[XX][ 1] =      cph * sth  * coords->dx_dxp[RR][1]   ;  // dx/dx1
  coords->dxc_dxp[XX][ 2] =  r * cph * cth  * coords->dx_dxp[TH][2]   ;  // dx/dx2
  coords->dxc_dxp[XX][ 3] = -y                                        ;  // dx/dx3
  coords->dxc_dxp[YY][ 0] =      sph * sth  * coords->dx_dxp[RR][0]   ;  // dy/dt
  coords->dxc_dxp[YY][ 1] =      sph * sth  * coords->dx_dxp[RR][1]   ;  // dy/dx1
  coords->dxc_dxp[YY][ 2] =  r * sph * cth  * coords->dx_dxp[TH][2]   ;  // dy/dx2
  coords->dxc_dxp[YY][ 3] = x                                         ;  // dy/dx3
  coords->dxc_dxp[ZZ][ 0] = cth             * coords->dx_dxp[RR][0]   ;  // dz/dt
  coords->dxc_dxp[ZZ][ 1] = cth             * coords->dx_dxp[RR][1]   ;  // dz/dx1
  coords->dxc_dxp[ZZ][ 2] = -r*sth          * coords->dx_dxp[TH][2]   ;  // dz/dx2
  coords->dxc_dxp[ZZ][ 3] = 0.                                        ;  // dz/dx3

# elif( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL )

  int j,k;
  double dxc_dx[NDIM][NDIM];
  dxc_dx[TT][TT] = 1.                    ;  // dt/dt
  dxc_dx[TT][RR] = 0.                    ;  // dt/dt
  dxc_dx[TT][TH] = 0.                    ;  // dt/dt
  dxc_dx[TT][PH] = 0.                    ;  // dt/dt
  dxc_dx[XX][TT] = 0.                    ;  // dx/dt
  dxc_dx[XX][RR] =      cph * sth        ;  // dx/dr
  dxc_dx[XX][TH] =  r * cph * cth        ;  // dx/dtheta
  dxc_dx[XX][PH] = -y                    ;  // dx/dphi
  dxc_dx[YY][TT] = 0.                    ;  // dy/dt
  dxc_dx[YY][RR] =      sph * sth        ;  // dy/dr
  dxc_dx[YY][TH] =  r * sph * cth        ;  // dy/dtheta
  dxc_dx[YY][PH] = x                     ;  // dy/dphi
  dxc_dx[ZZ][TT] = 0.                    ;  // dz/dt
  dxc_dx[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dx[ZZ][TH] = -r*sth                ;  // dz/dtheta
  dxc_dx[ZZ][PH] = 0.                    ;  // dz/dphi

  for(i=0; i<NDIM; i++) for(j=0; j<NDIM; j++) { 
      coords->dxc_dxp[i][j] = 0. ;
      for(k=0; k<NDIM; k++) {
	coords->dxc_dxp[i][j] += dxc_dx[i][k] * coords->dx_dxp[k][j] ;
      }
    }

# else 
  coords->dxc_dxp[TT][ 0] = 1.                                        ;  // dt/dt
  coords->dxc_dxp[TT][ 1] = 0.                                        ;  // dt/dx1
  coords->dxc_dxp[TT][ 2] = 0.                                        ;  // dt/dx2
  coords->dxc_dxp[TT][ 3] = 0.                                        ;  // dt/dx3
  coords->dxc_dxp[XX][ 0] = 0.                                        ;  // dx/dt
  coords->dxc_dxp[XX][ 1] =      cph * sth  * coords->dx_dxp[RR][1]   ;  // dx/dx1
  coords->dxc_dxp[XX][ 2] =  r * cph * cth  * coords->dx_dxp[TH][2]   ;  // dx/dx2
  coords->dxc_dxp[XX][ 3] = -y                                        ;  // dx/dx3
  coords->dxc_dxp[YY][ 0] = 0.                                        ;  // dy/dt
  coords->dxc_dxp[YY][ 1] =      sph * sth  * coords->dx_dxp[RR][1]   ;  // dy/dx1
  coords->dxc_dxp[YY][ 2] =  r * sph * cth  * coords->dx_dxp[TH][2]   ;  // dy/dx2
  coords->dxc_dxp[YY][ 3] = x                                         ;  // dy/dx3
  coords->dxc_dxp[ZZ][ 0] = 0.                                        ;  // dz/dt
  coords->dxc_dxp[ZZ][ 1] = cth             * coords->dx_dxp[RR][1]   ;  // dz/dx1
  coords->dxc_dxp[ZZ][ 2] = -r*sth          * coords->dx_dxp[TH][2]   ;  // dz/dx2
  coords->dxc_dxp[ZZ][ 3] = 0.                                        ;  // dz/dx3
# endif

#else
  fprintf(stderr,"xcart_of_xspher_special(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif

  return;
}


/***********************************************************************/
/***********************************************************************
  xcart_of_xspher_special2():
  --------------------
   -- calculates Cartesian coordinates from spherical coordinates;
   -- also calculates the covariant transformation matrix to x, not xp (this is what 
      makes this routine different from xcart_of_xspher_special()  ;
   -- independent of COORD_TYPE_CHOICE; 
***********************************************************************/
void xcart_of_xspher_special2(struct of_coord *coords )
{
  int i; 
  double r, th, ph;
  double x,y;
  double sth,cth,sph,cph;
  
  r  = coords->x[RR]; 
  th = coords->x[TH];
  ph = coords->x[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  x = r * sth * cph; 
  y = r * sth * sph; 

  coords->xcart[TT] = coords->x[TT];
  coords->xcart[XX] = x;
  coords->xcart[YY] = y;
  coords->xcart[ZZ] = r * cth;

  /* we load dxc_dxp[][]  with dxc_dx[][] values instead of allocating more memory:  */
  coords->dxc_dxp[TT][TT] = 1.                    ;  // dt/dt
  coords->dxc_dxp[TT][RR] = 0.                    ;  // dt/dt
  coords->dxc_dxp[TT][TH] = 0.                    ;  // dt/dt
  coords->dxc_dxp[TT][PH] = 0.                    ;  // dt/dt
  coords->dxc_dxp[XX][TT] = 0.                    ;  // dx/dt
  coords->dxc_dxp[XX][RR] =      cph * sth        ;  // dx/dr
  coords->dxc_dxp[XX][TH] =  r * cph * cth        ;  // dx/dtheta
  coords->dxc_dxp[XX][PH] = -y                    ;  // dx/dphi
  coords->dxc_dxp[YY][TT] = 0.                    ;  // dy/dt
  coords->dxc_dxp[YY][RR] =      sph * sth        ;  // dy/dr
  coords->dxc_dxp[YY][TH] =  r * sph * cth        ;  // dy/dtheta
  coords->dxc_dxp[YY][PH] = x                     ;  // dy/dphi
  coords->dxc_dxp[ZZ][TT] = 0.                    ;  // dz/dt
  coords->dxc_dxp[ZZ][RR] = cth                   ;  // dz/dr
  coords->dxc_dxp[ZZ][TH] = -r*sth                ;  // dz/dtheta
  coords->dxc_dxp[ZZ][PH] = 0.                    ;  // dz/dphi

  return;
}

/***********************************************************************/
/***********************************************************************
  xcart_of_xspher_special():
  --------------------
   -- calculates Cartesian coordinates from spherical coordinates;
   -- also calculates the covariant transformation matrix to xp ;
***********************************************************************/
void xcart_of_xspher_special3(double *xcart, double *xspher, double *xp, double dxc_dxp[][NDIM])
{
  int i; 
  double r, th, ph;
  double x,y;
  double sth,cth,sph,cph;
  double dx_dxp[NDIM][NDIM];

  dx_dxp_calc( xspher, xp,  dx_dxp );
  
#if( COORD_TYPE_CHOICE==COORD_DIAGONAL          ||\
     COORD_TYPE_CHOICE==COORD_IDENTITY          ||\
     COORD_TYPE_CHOICE==COORD_DIAGONAL2         ||\
     COORD_TYPE_CHOICE==COORD_DIAGONAL3         ||\
     COORD_TYPE_CHOICE==COORD_MIXED             ||\
     COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD ||\
     COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL   )

  r  = xspher[RR]; 
  th = xspher[TH];
  ph = xspher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  x = r * sth * cph; 
  y = r * sth * sph; 

  xcart[TT] = xspher[TT];
  xcart[XX] = x;
  xcart[YY] = y;
  xcart[ZZ] = r * cth;

# if( COORD_TYPE_CHOICE == COORD_MIXED )
  should-test-this-somehow

  dxc_dxp[TT][ 0] = 1.                                                                               ;  // dt/dt
  dxc_dxp[TT][ 1] = 0.                                                                               ;  // dt/dx1
  dxc_dxp[TT][ 2] = 0.                                                                               ;  // dt/dx2
  dxc_dxp[TT][ 3] = 0.                                                                               ;  // dt/dx3
  dxc_dxp[XX][ 0] = 0.                                                                               ;  // dx/dt
  dxc_dxp[XX][ 1] =      cph * sth  * dx_dxp[RR][1] + r * cph * cth  * dx_dxp[TH][1] ;  // dx/dx1
  dxc_dxp[XX][ 2] =  r * cph * cth  * dx_dxp[TH][2] +     cph * sth  * dx_dxp[RR][2] ;  // dx/dx2
  dxc_dxp[XX][ 3] = -y                                                                               ;  // dx/dx3
  dxc_dxp[YY][ 0] = 0.                                                                               ;  // dy/dt
  dxc_dxp[YY][ 1] =      sph * sth  * dx_dxp[RR][1] + r * sph * cth  * dx_dxp[TH][1] ;  // dy/dx1
  dxc_dxp[YY][ 2] =  r * sph * cth  * dx_dxp[TH][2] +     sph * sth  * dx_dxp[RR][2] ;  // dy/dx2
  dxc_dxp[YY][ 3] = x                                                                                ;  // dy/dx3
  dxc_dxp[ZZ][ 0] = 0.                                                                               ;  // dz/dt
  dxc_dxp[ZZ][ 1] = cth             * dx_dxp[RR][1] - r*sth * dx_dxp[TH][1]          ;  // dz/dx1
  dxc_dxp[ZZ][ 2] = -r*sth          * dx_dxp[TH][2] +   cth * dx_dxp[RR][2]          ;  // dz/dx2
  dxc_dxp[ZZ][ 3] = 0.                                                                               ;  // dz/dx3

# elif( COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD )
  dxc_dxp[TT][ 0] = 1.                                        ;  // dt/dt
  dxc_dxp[TT][ 1] = 0.                                        ;  // dt/dx1
  dxc_dxp[TT][ 2] = 0.                                        ;  // dt/dx2
  dxc_dxp[TT][ 3] = 0.                                        ;  // dt/dx3
  dxc_dxp[XX][ 0] =      cph * sth  * dx_dxp[RR][0]   ;  // dx/dt
  dxc_dxp[XX][ 1] =      cph * sth  * dx_dxp[RR][1]   ;  // dx/dx1
  dxc_dxp[XX][ 2] =  r * cph * cth  * dx_dxp[TH][2]   ;  // dx/dx2
  dxc_dxp[XX][ 3] = -y                                        ;  // dx/dx3
  dxc_dxp[YY][ 0] =      sph * sth  * dx_dxp[RR][0]   ;  // dy/dt
  dxc_dxp[YY][ 1] =      sph * sth  * dx_dxp[RR][1]   ;  // dy/dx1
  dxc_dxp[YY][ 2] =  r * sph * cth  * dx_dxp[TH][2]   ;  // dy/dx2
  dxc_dxp[YY][ 3] = x                                         ;  // dy/dx3
  dxc_dxp[ZZ][ 0] = cth             * dx_dxp[RR][0]   ;  // dz/dt
  dxc_dxp[ZZ][ 1] = cth             * dx_dxp[RR][1]   ;  // dz/dx1
  dxc_dxp[ZZ][ 2] = -r*sth          * dx_dxp[TH][2]   ;  // dz/dx2
  dxc_dxp[ZZ][ 3] = 0.                                        ;  // dz/dx3

# elif( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL )

  int j,k;
  double dxc_dx[NDIM][NDIM];
  dxc_dx[TT][TT] = 1.                    ;  // dt/dt
  dxc_dx[TT][RR] = 0.                    ;  // dt/dt
  dxc_dx[TT][TH] = 0.                    ;  // dt/dt
  dxc_dx[TT][PH] = 0.                    ;  // dt/dt
  dxc_dx[XX][TT] = 0.                    ;  // dx/dt
  dxc_dx[XX][RR] =      cph * sth        ;  // dx/dr
  dxc_dx[XX][TH] =  r * cph * cth        ;  // dx/dtheta
  dxc_dx[XX][PH] = -y                    ;  // dx/dphi
  dxc_dx[YY][TT] = 0.                    ;  // dy/dt
  dxc_dx[YY][RR] =      sph * sth        ;  // dy/dr
  dxc_dx[YY][TH] =  r * sph * cth        ;  // dy/dtheta
  dxc_dx[YY][PH] = x                     ;  // dy/dphi
  dxc_dx[ZZ][TT] = 0.                    ;  // dz/dt
  dxc_dx[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dx[ZZ][TH] = -r*sth                ;  // dz/dtheta
  dxc_dx[ZZ][PH] = 0.                    ;  // dz/dphi

  for(i=0; i<NDIM; i++) for(j=0; j<NDIM; j++) { 
      dxc_dxp[i][j] = 0. ;
      for(k=0; k<NDIM; k++) {
	dxc_dxp[i][j] += dxc_dx[i][k] * dx_dxp[k][j] ;
      }
    }

# else 
  dxc_dxp[TT][ 0] = 1.                                        ;  // dt/dt
  dxc_dxp[TT][ 1] = 0.                                        ;  // dt/dx1
  dxc_dxp[TT][ 2] = 0.                                        ;  // dt/dx2
  dxc_dxp[TT][ 3] = 0.                                        ;  // dt/dx3
  dxc_dxp[XX][ 0] = 0.                                        ;  // dx/dt
  dxc_dxp[XX][ 1] =      cph * sth  * dx_dxp[RR][1]   ;  // dx/dx1
  dxc_dxp[XX][ 2] =  r * cph * cth  * dx_dxp[TH][2]   ;  // dx/dx2
  dxc_dxp[XX][ 3] = -y                                        ;  // dx/dx3
  dxc_dxp[YY][ 0] = 0.                                        ;  // dy/dt
  dxc_dxp[YY][ 1] =      sph * sth  * dx_dxp[RR][1]   ;  // dy/dx1
  dxc_dxp[YY][ 2] =  r * sph * cth  * dx_dxp[TH][2]   ;  // dy/dx2
  dxc_dxp[YY][ 3] = x                                         ;  // dy/dx3
  dxc_dxp[ZZ][ 0] = 0.                                        ;  // dz/dt
  dxc_dxp[ZZ][ 1] = cth             * dx_dxp[RR][1]   ;  // dz/dx1
  dxc_dxp[ZZ][ 2] = -r*sth          * dx_dxp[TH][2]   ;  // dz/dx2
  dxc_dxp[ZZ][ 3] = 0.                                        ;  // dz/dx3
# endif

#else
  fprintf(stderr,"xcart_of_xspher_special3(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif

  return;
}




/****************************************************************************************

 all_x_of_xp(): 
 -------------
 --  assumes update_coord_time_funcs has been called

****************************************************************************************/

void  all_x_of_xp( struct of_coord *coords ) {

  int i,j,k,n,m,pos;
  double xp[NDIM];
  for(n = 0; n < NDIM; ++n)  xp[n] = coords->xp[n] ;
  i   = coords->i ;
  j   = coords->j ;
  k   = coords->k ;
  pos = coords->pos ;

#if( COORD_TYPE_CHOICE  == COORD_WARPED_SPHERICAL )

  coords->x[TT] = xp[0];   

  coords->x[RR] = coord_params->r_of_yt_[i][pos][0] + coord_params->square_per_xp2_2[j][pos] * (
                                                                                   (coord_params->r_of_yt_[i][pos][1] - coord_params->r_of_yt_[i][pos][0]) * coord_params->square_per_xp3_3[k][pos]
                                                                                 + (coord_params->r_of_yt_[i][pos][2] - coord_params->r_of_yt_[i][pos][0]) * coord_params->square_per_xp3_4[k][pos]
                                                                                                ) ;

  coords->x[TH] = M_PI * ( xp[2] - a_z10 * (   coord_params->square_int_per_xp2_1[j][pos] 
                                             - coord_params->square_int_per_zbh1 
                                             - (xp[2]-zbh_[1]) * ( 
                                                                    coord_params->square_int_per_z1 
                                                                  - coord_params->square_int_per_z0
                                                                   ) 
                                               )
                           );

  coords->x[PH] = 2.*M_PI * (   xp[3] 
                              - coord_params->square_per_xp2_2[j][pos] * (
                                                                          a_x10 * coord_params->square_xp1_3[i][pos] * (  coord_params->square_int_per_xp3_1[k][pos] 
                                                                                                                        - coord_params->square_int_per_xbh1 
                                                                                                                        - (xp[3]-xbh_[1])*(  coord_params->square_int_per_x11 
                                                                                                                                           - coord_params->square_int_per_x01) )
                                                                        + a_x20 * coord_params->square_xp1_4[i][pos] * (  coord_params->square_int_per_xp3_2[k][pos]  
                                                                                                                        - coord_params->square_int_per_xbh2 
                                                                                                                        - (xp[3]-xbh_[2])*(  coord_params->square_int_per_x12 
                                                                                                                                           - coord_params->square_int_per_x02) )  
                                                                          )
                                ) ;


  // derivatives

  // Start from delta function:
  for(n=0  ;n<NDIM ;n++) for(m=0  ;m<NDIM ;m++)  coords->dx_dxp[n][m] = DELTA(n,m) ;



  double dr1_dybh1, dr2_dybh2, dr3_dybh3 ;

  dr1_dybh1 =  ar_[1]*ar_[1]/(Rin_[1] - Rout_[1] + br_[1]) * s_[1] * (
                                                                          sinh(s_[1]*xp[1]) 
                                                                        + sinh(s_[1]*(1. - xp[1])) 
                                                                        - sinh(s_[1])
                                                                        + s_[1] * (
                                                                                   - cosh( s_[1]*(xp[1] - ybh_[1]) )  
                                                                                   + (1. - xp[1]) * cosh(s_[1]*ybh_[1]) 
                                                                                   + xp[1]*cosh( s_[1]*(1. - ybh_[1]) )
                                                                                   )
                                                                      ) ;

  dr2_dybh2 =  ar_[2]*ar_[2]/(Rin_[2] - Rout_[2] + br_[2]) * s_[2] * (
                                                                          sinh(s_[2]*xp[1]) 
                                                                        + sinh(s_[2]*(1. - xp[1])) 
                                                                        - sinh(s_[2])
                                                                        + s_[2] * (
                                                                                   - cosh( s_[2]*(xp[1] - ybh_[2]) ) 
                                                                                   + (1. - xp[1]) * cosh(s_[2]*ybh_[2]) 
                                                                                   + xp[1]*cosh( s_[2]*( 1. - ybh_[2]) ) 
                                                                                   )
                                                                      ) ;

  dr3_dybh3 =  ar_[0]*ar_[0]/(Rin_[0] - Rout_[0] + br_[0]) * s_[0] * (
                                                                          sinh(s_[0]*xp[1]) 
                                                                        + sinh(s_[0]*(1. - xp[1])) 
                                                                        - sinh(s_[0])
                                                                        + s_[0] * (
                                                                                   - cosh( s_[0]*(xp[1] - ybh_[0]) ) 
                                                                                   + (1. - xp[1]) * cosh(s_[0]*ybh_[0]) 
                                                                                   + xp[1]*cosh( s_[0]*( 1. - ybh_[0]) ) 
                                                                                   )
                                                                      ) ;


  double dr_dxbh1, dr_dxbh2, dr_dybh1, dr_dybh2, dr_dybh3;
  double dphi_dxbh1, dphi_dxbh2, dphi_dybh1, dphi_dybh2, dphi_dybh3;

  dr_dxbh1 = - coord_params->square_per_xp2_2[j][pos] * coord_params->dsquare_per_xp3_3[k][pos] *
               (coord_params->r_of_yt_[i][pos][1] - coord_params->r_of_yt_[i][pos][0]);

  dr_dxbh2 = - coord_params->square_per_xp2_2[j][pos] * coord_params->dsquare_per_xp3_4[k][pos] *
               (coord_params->r_of_yt_[i][pos][2] - coord_params->r_of_yt_[i][pos][0]);

  dr_dybh1 = dr1_dybh1 * coord_params->square_per_xp2_2[j][pos] * coord_params->square_per_xp3_3[k][pos] ;

  dr_dybh2 = dr2_dybh2 * coord_params->square_per_xp2_2[j][pos] * coord_params->square_per_xp3_4[k][pos] ;

  dr_dybh3 = dr3_dybh3 * ( 1. - coord_params->square_per_xp2_2[j][pos] * (
                                                                           coord_params->square_per_xp3_3[k][pos]
                                                                         + coord_params->square_per_xp3_4[k][pos]
                                                                          )
                           ) ;


  dphi_dxbh1 = - 2.*M_PI * a_x10 * coord_params->square_per_xp2_2[j][pos] * coord_params->square_xp1_3[i][pos]
                                 * ( - coord_params->square_per_xp3_1[k][pos]
                                     + coord_params->square_int_per_x11 
                                     - coord_params->square_int_per_x01
                                     + (xp[3] - xbh_[1]) * ( coord_params->square_per_x11
                                                           - coord_params->square_per_x01 )
                                     ) ;

  dphi_dxbh2 = - 2.*M_PI * a_x20 * coord_params->square_per_xp2_2[j][pos] *  coord_params->square_xp1_4[i][pos]
                                 * ( - coord_params->square_per_xp3_2[k][pos]
                                     + coord_params->square_int_per_x12
                                     - coord_params->square_int_per_x02
                                     + (xp[3] - xbh_[2]) * ( coord_params->square_per_x12
                                                           - coord_params->square_per_x02 )
                                     ) ;

  dphi_dybh1 = 2.*M_PI * a_x10 * coord_params->square_per_xp2_2[j][pos] * coord_params->dsquare_xp1_3[i][pos]
                                 * (  coord_params->square_int_per_xp3_1[k][pos]
                                    - coord_params->square_int_per_xbh1
                                    - (xp[3] - xbh_[1]) * (  coord_params->square_int_per_x11 
                                                           - coord_params->square_int_per_x01 )
                                      ) ;

  dphi_dybh2 = 2.*M_PI * a_x20 * coord_params->square_per_xp2_2[j][pos] * coord_params->dsquare_xp1_4[i][pos]
                                 * (  coord_params->square_int_per_xp3_2[k][pos]
                                    - coord_params->square_int_per_xbh2
                                    - (xp[3] - xbh_[2]) * (  coord_params->square_int_per_x12 
                                                          -  coord_params->square_int_per_x02  )
                                     ) ;






  coords->dx_dxp[RR][0] = dxbh1_dt * dr_dxbh1 + dxbh2_dt * dr_dxbh2 + dybh1_dt * dr_dybh1 + dybh2_dt * dr_dybh2 + dybh3_dt * dr_dybh3 ;

  coords->dx_dxp[RR][2] = coord_params->dsquare_per_xp2_2[j][pos] * (
                                                                (coord_params->r_of_yt_[i][pos][1] - coord_params->r_of_yt_[i][pos][0]) * coord_params->square_per_xp3_3[k][pos]
                                                              + (coord_params->r_of_yt_[i][pos][2] - coord_params->r_of_yt_[i][pos][0]) * coord_params->square_per_xp3_4[k][pos]
                                                              ) ;

  coords->dx_dxp[RR][3] = - dr_dxbh1 - dr_dxbh2 ;

  coords->dx_dxp[RR][1] =   coord_params->drdy_of_yt_[i][pos][0] 
                  + coord_params->square_per_xp2_2[j][pos] * (
                                                                (coord_params->drdy_of_yt_[i][pos][1] - coord_params->drdy_of_yt_[i][pos][0]) * coord_params->square_per_xp3_3[k][pos]
                                                              + (coord_params->drdy_of_yt_[i][pos][2] - coord_params->drdy_of_yt_[i][pos][0]) * coord_params->square_per_xp3_4[k][pos]
                                                              ) ;

  coords->dx_dxp[PH][3] = 2.*M_PI - 2.*M_PI * a_x10 * coord_params->square_per_xp2_2[j][pos] * coord_params->square_xp1_3[i][pos]
                                            * ( coord_params->square_per_xp3_1[k][pos] - coord_params->square_int_per_x11 + coord_params->square_int_per_x01 )
                          - 2.*M_PI * a_x20 * coord_params->square_per_xp2_2[j][pos] * coord_params->square_xp1_4[i][pos] 
                                            * ( coord_params->square_per_xp3_2[k][pos] - coord_params->square_int_per_x12 + coord_params->square_int_per_x02 ) ;

  coords->dx_dxp[PH][1] = - 2.*M_PI * a_x10 *  coord_params->square_per_xp2_2[j][pos] * coord_params->dsquare_xp1_3[i][pos] 
                                    * ( 
                                          coord_params->square_int_per_xp3_1[k][pos]
                                        - coord_params->square_int_per_xbh1
                                        - (xp[3] - xbh_[1]) * (   coord_params->square_int_per_x11
                                                                - coord_params->square_int_per_x01 )
                                        )
                          - 2.*M_PI * a_x20 *  coord_params->square_per_xp2_2[j][pos] * coord_params->dsquare_xp1_4[i][pos]
                                    * ( 
                                         coord_params->square_int_per_xp3_2[k][pos]
                                       - coord_params->square_int_per_xbh2
                                       - (xp[3] - xbh_[2]) * (   coord_params->square_int_per_x12
                                                               - coord_params->square_int_per_x02 ) 
                                        ) ;

  coords->dx_dxp[PH][2] = - 2.*M_PI * coord_params->dsquare_per_xp2_2[j][pos] * (
                                                                          a_x10 * coord_params->square_xp1_3[i][pos] * (   coord_params->square_int_per_xp3_1[k][pos] 
                                                                                                                         - coord_params->square_int_per_xbh1
                                                                                                                         - (xp[3]-xbh_[1])*(  coord_params->square_int_per_x11
                                                                                                                                            - coord_params->square_int_per_x01) )
                                                                        + a_x20 * coord_params->square_xp1_4[i][pos] * (   coord_params->square_int_per_xp3_2[k][pos] 
                                                                                                                         - coord_params->square_int_per_xbh2
                                                                                                                         - (xp[3]-xbh_[2])*(  coord_params->square_int_per_x12
                                                                                                                                            - coord_params->square_int_per_x02) )  
                                                                          ) ;

  coords->dx_dxp[PH][0] =  dxbh1_dt * dphi_dxbh1 + dxbh2_dt * dphi_dxbh2 + dybh1_dt * dphi_dybh1 + dybh2_dt * dphi_dybh2 ;


  coords->dx_dxp[TH][2] = M_PI - M_PI * a_z10 * ( coord_params->square_per_xp2_1[j][pos]
                                                - coord_params->square_int_per_z1
                                                + coord_params->square_int_per_z0 ) ;
  coords->dx_dxp[TH][0] = 0. ;
  coords->dx_dxp[TH][1] = 0. ;
  coords->dx_dxp[TH][3] = 0. ;


  //  dxp_dx_calc3_general( coords->x, coords->xp,  coords->dx_dxp , coords->dxp_dx, &(coords->det_dx_dxp) );

  double det_tmp;
  if( invert_matrix2( coords->dx_dxp, coords->dxp_dx, &det_tmp )  ) {
    fprintf(stdout,"all_x_of_xp(): singular coordinate transformation detected!!  \n"); fflush(stdout);
    fail(FAIL_BASIC,0);
  }
  coords->det_dx_dxp = det_tmp*det_tmp;



#elif( COORD_TYPE_CHOICE  == COORD_WARPED_CARTESIAN )

  coords->x[TT] = xp[0] ;

  coords->x[XX] = xp[1] 
                 -  a_x10 * coord_params->square_per_xp2_3[j][pos] * ( 
                                                                   coord_params->square_int_per_xp1_1[i][pos]
                                                                 - coord_params->square_int_per_xbh1
                                                                 - (xp[1]-xbh_[1])*( coord_params->square_int_per_x11
                                                                                   - coord_params->square_int_per_x01 ) 
                                                                        )
                 -  a_x20 * coord_params->square_per_xp2_4[j][pos] * ( 
                                                                      coord_params->square_int_per_xp1_2[i][pos]
                                                                    - coord_params->square_int_per_xbh2
                                                                    - (xp[1]-xbh_[2])*( coord_params->square_int_per_x12
                                                                                      - coord_params->square_int_per_x02 )
                                                                       ) ; 
  coords->x[XX] = xmin + (xmax-xmin)*coords->x[XX] ;

  coords->x[YY] = xp[2] 
                 -  a_y10 * coord_params->square_per_xp1_3[i][pos] * ( 
                                                                   coord_params->square_int_per_xp2_1[j][pos]
                                                                 - coord_params->square_int_per_ybh1
                                                                 - (xp[2]-ybh_[1])*( coord_params->square_int_per_y11
                                                                                   - coord_params->square_int_per_y01 ) 
                                                                        )
                 -  a_y20 * coord_params->square_per_xp1_4[i][pos] * ( 
                                                                      coord_params->square_int_per_xp2_2[j][pos]
                                                                    - coord_params->square_int_per_ybh2
                                                                    - (xp[2]-ybh_[2])*( coord_params->square_int_per_y12
                                                                                      - coord_params->square_int_per_y02 )
                                                                       ) ;
  coords->x[YY] = ymin + (ymax-ymin)*coords->x[YY] ;
  
  coords->x[ZZ] = xp[3] ;



  // derivatives

  // Start from delta function:
  for(n=0  ;n<NDIM ;n++) for(m=0  ;m<NDIM ;m++)  coords->dx_dxp[n][m] = DELTA(n,m) ;

  coords->dx_dxp[XX][1] = 1.
                 -  a_x10 * coord_params->square_per_xp2_3[j][pos] * ( 
                                                                   coord_params->square_per_xp1_1[i][pos]
                                                                 - ( coord_params->square_int_per_x11
                                                                   - coord_params->square_int_per_x01 ) 
                                                                        )
                 -  a_x20 * coord_params->square_per_xp2_4[j][pos] * ( 
                                                                      coord_params->square_per_xp1_2[i][pos]
                                                                    - ( coord_params->square_int_per_x12
                                                                      - coord_params->square_int_per_x02 )
                                                                       ) ;  
  coords->dx_dxp[XX][1] = (xmax-xmin)*coords->dx_dxp[XX][1] ;

  coords->dx_dxp[XX][2] = -  a_x10 * coord_params->dsquare_per_xp2_3[j][pos] * ( 
                                                                   coord_params->square_int_per_xp1_1[i][pos]
                                                                 - coord_params->square_int_per_xbh1
                                                                 - (xp[1]-xbh_[1])*( coord_params->square_int_per_x11
                                                                                   - coord_params->square_int_per_x01 ) 
                                                                        )
                          -  a_x20 * coord_params->dsquare_per_xp2_4[j][pos] * ( 
                                                                      coord_params->square_int_per_xp1_2[i][pos]
                                                                    - coord_params->square_int_per_xbh2
                                                                    - (xp[1]-xbh_[2])*( coord_params->square_int_per_x12
                                                                                      - coord_params->square_int_per_x02 )
                                                                       ) ;  
  coords->dx_dxp[XX][2] = (xmax-xmin)*coords->dx_dxp[XX][2] ;

  coords->dx_dxp[YY][1] = -  a_y10 * coord_params->dsquare_per_xp1_3[i][pos] * ( 
                                                                   coord_params->square_int_per_xp2_1[j][pos]
                                                                 - coord_params->square_int_per_ybh1
                                                                 - (xp[2]-ybh_[1])*( coord_params->square_int_per_y11
                                                                                   - coord_params->square_int_per_y01 ) 
                                                                        )
                          -  a_y20 * coord_params->dsquare_per_xp1_4[i][pos] * ( 
                                                                      coord_params->square_int_per_xp2_2[j][pos]
                                                                    - coord_params->square_int_per_ybh2
                                                                    - (xp[2]-ybh_[2])*( coord_params->square_int_per_y12
                                                                                      - coord_params->square_int_per_y02 )
                                                                       ) ;  
  coords->dx_dxp[YY][1] = (ymax-ymin)*coords->dx_dxp[YY][1] ;

  coords->dx_dxp[YY][2] = 1. 
                              -  a_y10 * coord_params->square_per_xp1_3[i][pos] * ( 
                                                                   coord_params->square_per_xp2_1[j][pos]
                                                                 - ( coord_params->square_int_per_y11
                                                                   - coord_params->square_int_per_y01 ) 
                                                                        )
                              -  a_y20 * coord_params->square_per_xp1_4[i][pos] * ( 
                                                                      coord_params->square_per_xp2_2[j][pos]
                                                                    - ( coord_params->square_int_per_y12
                                                                      - coord_params->square_int_per_y02 )
                                                                       ) ;  
  coords->dx_dxp[YY][2] = (ymax-ymin)*coords->dx_dxp[YY][2] ;

  double dx_dxbh1, dx_dxbh2, dx_dybh1, dx_dybh2; 
  double dy_dxbh1, dy_dxbh2, dy_dybh1, dy_dybh2; 

  dx_dxbh1 =   - a_x10 * coord_params->square_per_xp2_3[j][pos] * ( 
                                                                 - coord_params->square_per_xp1_1[i][pos]
                                                                 + coord_params->square_int_per_x11
                                                                 - coord_params->square_int_per_x01
                                                                 + (xp[1]-xbh_[1])*( coord_params->square_per_x11
                                                                                   - coord_params->square_per_x01 ) 
                                                                       ) ;

  dx_dxbh2 =   - a_x20 * coord_params->square_per_xp2_4[j][pos] * ( 
                                                                    - coord_params->square_per_xp1_2[i][pos]
                                                                    + coord_params->square_int_per_x12
                                                                    - coord_params->square_int_per_x02
                                                                    + (xp[1]-xbh_[2])*( coord_params->square_per_x12
                                                                                      - coord_params->square_per_x02 )
                                                                    ) ;

  dx_dybh1 =    a_x10 * coord_params->dsquare_per_xp2_3[j][pos] * ( 
                                                                   coord_params->square_int_per_xp1_1[i][pos]
                                                                 - coord_params->square_int_per_xbh1
                                                                 - (xp[1]-xbh_[1])*( coord_params->square_int_per_x11
                                                                                   - coord_params->square_int_per_x01 ) 
                                                                    ) ;

  dx_dybh2 =    a_x20 * coord_params->dsquare_per_xp2_4[j][pos] * ( 
                                                                      coord_params->square_int_per_xp1_2[i][pos]
                                                                    - coord_params->square_int_per_xbh2
                                                                    - (xp[1]-xbh_[2])*( coord_params->square_int_per_x12
                                                                                      - coord_params->square_int_per_x02 )
                                                                       ) ;  

  dy_dybh1 = - a_y10 * coord_params->square_per_xp1_3[i][pos] * ( 
                                                                 - coord_params->square_per_xp2_1[j][pos]
                                                                 + coord_params->square_int_per_y11
                                                                 - coord_params->square_int_per_y01 
                                                                 + (xp[2]-ybh_[1])*( coord_params->square_per_y11
                                                                                   - coord_params->square_per_y01 ) 
                                                                  ) ;

  dy_dybh2 = - a_y20 * coord_params->square_per_xp1_4[i][pos] * ( 
                                                                    - coord_params->square_per_xp2_2[j][pos]
                                                                    + coord_params->square_int_per_y12
                                                                    - coord_params->square_int_per_y02 
                                                                    + (xp[2]-ybh_[2])*( coord_params->square_per_y12
                                                                                      - coord_params->square_per_y02 )
                                                                  ) ;  

  dy_dxbh1 = a_y10 * coord_params->dsquare_per_xp1_3[i][pos] * ( 
                                                                   coord_params->square_int_per_xp2_1[j][pos]
                                                                 - coord_params->square_int_per_ybh1
                                                                 - (xp[2]-ybh_[1])*( coord_params->square_int_per_y11
                                                                                   - coord_params->square_int_per_y01 ) 
                                                                 ) ;

  dy_dxbh2 = a_y20 * coord_params->dsquare_per_xp1_4[i][pos] * ( 
                                                                      coord_params->square_int_per_xp2_2[j][pos]
                                                                    - coord_params->square_int_per_ybh2
                                                                    - (xp[2]-ybh_[2])*( coord_params->square_int_per_y12
                                                                                      - coord_params->square_int_per_y02 )
                                                                 ) ;  

  // dX/dt
  coords->dx_dxp[XX][0] = dxbh1_dt * dx_dxbh1 + dxbh2_dt * dx_dxbh2 + dybh1_dt * dx_dybh1 + dybh2_dt * dx_dybh2 ;
  coords->dx_dxp[XX][0] = (xmax-xmin)*coords->dx_dxp[XX][0] ;

  // dY/dt
  coords->dx_dxp[YY][0] = dxbh1_dt * dy_dxbh1 + dxbh2_dt * dy_dxbh2 + dybh1_dt * dy_dybh1 + dybh2_dt * dy_dybh2 ;

  coords->dx_dxp[YY][0] = (ymax-ymin)*coords->dx_dxp[YY][0] ;

  double det_tmp;
  if( invert_matrix2( coords->dx_dxp, coords->dxp_dx, &det_tmp )  ) {
    fprintf(stdout,"all_x_of_xp(): singular coordinate transformation detected!!  \n"); fflush(stdout);
    fail(FAIL_BASIC,0);
  }
  coords->det_dx_dxp = det_tmp*det_tmp;



#else 
  fprintf(stdout,"all_x_of_xp(): this routine should be called only for COORD_TYPE_CHOICE = COORD_WARPED_SPHERICAL or COORD_TYPE_CHOICE = COORD_WARPED_CARTESIAN!! : %d \n", 
	  COORD_TYPE_CHOICE); fflush(stdout);
  fail(FAIL_BASIC,0);
#endif
}


/***********************************************************************/
/***********************************************************************
  coord_of_xp():
  --------------------
   -- calculates rest of of_coord structure from xp 
***********************************************************************/
void  coord_of_xp( double *xp, struct of_coord *coords ) 
{
  
  int i,j;

  TRACE_BEG;

  coords->xp[0] = xp[0]; 
  coords->xp[1] = xp[1]; 
  coords->xp[2] = xp[2]; 
  coords->xp[3] = xp[3]; 

#if( (COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL) || (COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN) ) 

  all_x_of_xp( coords );

#else
  x_of_xp( coords->x, coords->xp);
  dx_dxp_calc( coords->x, coords->xp,  coords->dx_dxp );
  dxp_dx_calc2( coords->x, coords->xp,  coords->dx_dxp , coords->dxp_dx );
  coords->det_dx_dxp = det_dx_dxp_calc2( coords->dx_dxp );
#endif

#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
  xcart_of_xspher_special( coords );
  coords->r = coords->x[RR];
#elif( TOP_TYPE_CHOICE == TOP_CARTESIAN )
  DLOOP1 { coords->xcart[i] = coords->x[i]; }
  DLOOP2 { coords->dxc_dxp[i][j] = coords->dx_dxp[i][j]; }
  coords->r = sqrt(coords->x[1]*coords->x[1] + coords->x[2]*coords->x[2] + coords->x[3]*coords->x[3] );
#else
  fprintf(stdout,"coord_of_xp() TOP_CYLINDRICAL choice not implemented yet!\n"); fflush(stdout);
  fail(FAIL_BASIC,0);
#endif


    /* Determine the floor state: */ 
#if( FLOOR_DIM == 0 ) 
  coords->rhoflr = RHOMIN;
  coords->uuflr  = UUMIN;
#else
  coords->rhoflr = RHOMIN*pow(coords->r,RHOPOWER); 
  coords->uuflr  =  UUMIN*pow(coords->r, UUPOWER);
#endif

  TRACE_END;

  return;
}

/***********************************************************************/
/***********************************************************************
  coord_of_xp2():
  --------------------
   -- calculates rest of of_coord structure from xp 
   -- like coord_of_xp() but without "optimizations"
***********************************************************************/
void  coord_of_xp2( double *xp, struct of_coord *coords ) 
{
  
  int i,j;

  coords->xp[0] = xp[0]; 
  coords->xp[1] = xp[1]; 
  coords->xp[2] = xp[2]; 
  coords->xp[3] = xp[3]; 

#if( (COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL) || (COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN) ) 

  x_of_xp( coords->x, coords->xp);
  dx_dxp_calc( coords->x, coords->xp,  coords->dx_dxp );
  dxp_dx_calc3_general( coords->x, coords->xp,  coords->dx_dxp , coords->dxp_dx, &(coords->det_dx_dxp) ); 

#else
  x_of_xp( coords->x, coords->xp);
  dx_dxp_calc( coords->x, coords->xp,  coords->dx_dxp );
  dxp_dx_calc2( coords->x, coords->xp,  coords->dx_dxp , coords->dxp_dx );
  coords->det_dx_dxp = det_dx_dxp_calc2( coords->dx_dxp );
#endif

#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
  xcart_of_xspher_special( coords );
  coords->r = coords->x[RR];
#elif( TOP_TYPE_CHOICE == TOP_CARTESIAN )
  DLOOP1 { coords->xcart[i] = coords->x[i]; }
  DLOOP2 { coords->dxc_dxp[i][j] = coords->dx_dxp[i][j]; }
  coords->r = sqrt(coords->x[1]*coords->x[1] + coords->x[2]*coords->x[2] + coords->x[3]*coords->x[3] );
#else
  fprintf(stdout,"coord_of_xp() TOP_CYLINDRICAL choice not implemented yet!\n"); fflush(stdout);
  fail(FAIL_BASIC,0);
#endif


    /* Determine the floor state: */ 
#if( FLOOR_DIM == 0 ) 
  coords->rhoflr = RHOMIN;
  coords->uuflr  = UUMIN;
#else
  coords->rhoflr = RHOMIN*pow(coords->r,RHOPOWER); 
  coords->uuflr  =  UUMIN*pow(coords->r, UUPOWER);
#endif

  return;
}

/***********************************************************************/
/***********************************************************************
  coord_of_r():
  --------------------
   -- calculates rest of of_coord structure from r, the radius
   -- assume t=phi=0,  th= 0.5*PI
***********************************************************************/
void  coord_of_r( double r, struct of_coord *coords ) 
{
  double xp[NDIM];

  /* Need to set these to valid memory locations in order to work without seg faults for warped coordinates: */
  coords->i   =  N1S;
  coords->j   =  N2S;
  coords->k   =  N3S;
  coords->pos = CENT;

//Dennis
#if( COORD_TYPE_CHOICE == COORD_IDENTITY && TOP_TYPE_CHOICE == TOP_SPHERICAL )
  coords->x[RR] = r;
  coords->x[TH] = 0.5*M_PI;
  coords->x[PH] = 0.;
#else
  xp_at_eq( r, xp );
  coord_of_xp(xp, coords);
#endif

#if( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL ) 
  static unsigned short int local_first_time = 1; 
  if( local_first_time ) { 
    fprintf(stdout,"coord_of_r():  WARNING:  using coord_of_r with  COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL !!!!\n"); 
    fflush(stdout); 
    local_first_time = 0; 
  }
  
  coords->x[RR] = r;
  coords->x[TH] = 0.5*M_PI;
  coords->x[PH] = 0.;
#endif

  if( IS_DIFFERENT_DOUBLE( (r) , (coords->x[RR]) ) ) { 
    fprintf(stdout,"coord_of_r() WARNING:  radii are different: %26.16e %26.16e \n",r,coords->x[RR]); fflush(stdout);
  }


  return;
}

/***********************************************************************/
/***********************************************************************
  calc_coord():
  --------------------
   -- calculates all the parts of a coord structure;
   -- assumes that   update_coord_time_funcs(t_now) has already been called;
   -- "coords" must already be allocated;
***********************************************************************/
void  calc_coord( int i, int j, int k, int pos, struct of_coord *coords)
{
  
  /* Calculate global quantities: */
  coords->i   = i   ;
  coords->j   = j   ;
  coords->k   = k   ;
  coords->pos = pos ;
  coord(i,j,k,pos,coords->xp);
  coords->xp[TT] = coord_params->t;
  coord_of_xp( coords->xp, coords ) ;

  return;
}

/***********************************************************************/
/***********************************************************************
  calc_all_coord():
  --------------------
   -- loads the array of of_coord structures;
***********************************************************************/
void  calc_all_coord(int n, double t_now ) 
{
  
  int i,j,k,pos;
  struct of_coord *coords;

  TRACE_BEG;

#if( USE_MASK) 
    x1_min = x2_min = x3_min =  1e200 ;     
    x1_max = x2_max = x3_max = -1e200 ;   
#endif

  /* Calculate global quantities: */

  ALL_LOOP POSLOOP { 
    get_coord(i,j,k,pos,n,coords);

    coords->i   = i   ;
    coords->j   = j   ;
    coords->k   = k   ;
    coords->pos = pos ;

    coord(i,j,k,pos,coords->xp);
    coord_of_xp( coords->xp, coords ) ;

#if( USE_MASK) 
    double *x = coords->x;
    if( x[1] < x1_min ) { x1_min = x[1] ; }
    if( x[1] > x1_max ) { x1_max = x[1] ; }
    if( x[2] < x2_min ) { x2_min = x[2] ; }
    if( x[2] > x2_max ) { x2_max = x[2] ; }
    if( x[3] < x3_min ) { x3_min = x[3] ; }
    if( x[3] > x3_max ) { x3_max = x[3] ; }
#endif

  }

  TRACE_END;
  return;
}

/***********************************************************************/
/***********************************************************************
  copy_coord():
  --------------------
   -- copies data from one of_coord (coords1) to another (coords2);
***********************************************************************/
void copy_coord(struct of_coord  *coords1, struct of_coord *coords2 )
{
  int i,j;

  coords2->i   = coords1->i  ;
  coords2->j   = coords1->j  ;
  coords2->k   = coords1->k  ;
  coords2->pos = coords1->pos;

  DLOOP1 { coords2->xp[i]         = coords1->xp[i]         ; }
  DLOOP1 { coords2->x[i]          = coords1->x[i]          ; }
  DLOOP1 { coords2->xcart[i]      = coords1->xcart[i]      ; }
  DLOOP2 { coords2->dx_dxp[i][j]  = coords1->dx_dxp[i][j]  ; }
  DLOOP2 { coords2->dxp_dx[i][j]  = coords1->dxp_dx[i][j]  ; }
  DLOOP2 { coords2->dxc_dxp[i][j] = coords1->dxc_dxp[i][j] ; }

  coords2->det_dx_dxp = coords1->det_dx_dxp; 

  return;
}

/***********************************************************************/
/***********************************************************************
   ROUTINES FOR   COORD_DIAGONAL2 
        -- requires th_length
***********************************************************************/
/***********************************************************************/

/* Step function from 0 to 1  with "strength" s at x=x0 */
#define STEP(x,x0,s) ( 0.5 + (atan((s)*((x)/(x0)-1.)))/M_PI )

/* Derivative of step() */
#define DSTEP(x,x0,s) ( ((s)*(x0))/(M_PI*((s)*(s)*((x) - (x0))*((x) - (x0)) + (x0)*(x0))) )

/*****************************************************************************************/
/*****************************************************************************************
  setup_diag2():
 -----------------------
  -- sets auxiliary constants needed for COORD_DIAGONAL2
  -- The coordinates x2 are defined by their discretization (dx2) as a function of xp2,
  -- dx2 is a set of C^\infinity piecewise-linear functions that are joined by a 
     continuous representation of the step function (to make the function differentiable)
  -- There can be n_diag2_lines, with slope ci_diag2[] and "y-intercepts" di_diag2[]; 
  -- ai_diag2[] and bi_diag2[] are derived quantities that are just the coefficients on the 
      linear and constant parts of each segment.  In other words, are desired function is 
 
      f(x) =   c_i x + d_i    if   x_i < x < x_{i+1} 

      and we derive a differentiable version of this this function by matching to 
 
      g(x) = \sum_i  (a_i x + b_i)*step(x,x_i,s_i)  

      where step(x,x_i,s_i) is our differentiable representation of the step function 
        going from 0 to 1  at x=x_i with "strength" s_i;  that is, the larger s_i 
        is the more like a step function step() is ;
    
   -- ai,bi are easily found to be linear combinations of ci,di respectively 
        (see calc_ai_bi_di_diag2())
  
*****************************************************************************************/
void setup_diag2(int type)
{
  int i; 
  double fe,fb;

  /* Calculate all auxiliary arrays (assume that they have already been allocated): */
  if( type == 0 ) { 
    calc_ai_bi_di_diag2();
  }
  else {
    calc_ai_bi_diag2();
  }    

  /* Need to renormalize coordinates to desired range : 
     -- move boundaries just within true boundaries since we remap xp2 to always lie within 
         x0 and xn 
  */
  fb = x2_of_xp2_diag2(xi_diag2[0            ]*(1.+SMALL)); 
  fe = x2_of_xp2_diag2(xi_diag2[n_diag2_lines]*(1.-SMALL)); 
  beg_diag2 = fb;
  end_diag2 = fe;
  len_diag2 = end_diag2 - beg_diag2;
  norm_diag2 = th_length/len_diag2;

  /* Because of how the step function is defined, we need to make sure that none of the 
     transition locations are at  x=0   */
  for(i=0; i<n_diag2_lines; i++) if(fabs(xi_diag2[i]) < SMALL) { xi_diag2[i] = SMALL ; }   

  fprintf(stdout,"ai = "); for(i=0;i<n_diag2_lines;i++) { fprintf(stdout," %28.18e ",ai_diag2[i]); }  fprintf(stdout,"\n");
  fprintf(stdout,"bi = "); for(i=0;i<n_diag2_lines;i++) { fprintf(stdout," %28.18e ",bi_diag2[i]); }  fprintf(stdout,"\n");
  fprintf(stdout,"ci = "); for(i=0;i<n_diag2_lines;i++) { fprintf(stdout," %28.18e ",ci_diag2[i]); }  fprintf(stdout,"\n");
  fprintf(stdout,"di = "); for(i=0;i<n_diag2_lines;i++) { fprintf(stdout," %28.18e ",di_diag2[i]); }  fprintf(stdout,"\n");
  fprintf(stdout,"si = "); for(i=0;i<n_diag2_lines;i++) { fprintf(stdout," %28.18e ",si_diag2[i]); }  fprintf(stdout,"\n");
  fprintf(stdout,"xi = "); for(i=0;i<n_diag2_lines+1;i++) { fprintf(stdout," %28.18e ",xi_diag2[i]); }  fprintf(stdout,"\n");

  fprintf(stdout,"beg_diag2  = %28.18e \n",beg_diag2);
  fprintf(stdout,"end_diag2  = %28.18e \n",end_diag2);
  fprintf(stdout,"len_diag2  = %28.18e \n",len_diag2);
  fprintf(stdout,"norm_diag2 = %28.18e \n",norm_diag2);
  fprintf(stdout,"theor norm = %28.18e \n",theoretical_norm_diag2());

  fflush(stdout);


  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  alloc_diag2(): 
 ------------------
  -- responsible for allocating global constant arrays that are used to 
       describe COORD_DIAGONAL2
*****************************************************************************************/
void alloc_diag2(void)
{ 
  if( (ai_diag2 = (double *) calloc(n_diag2_lines,sizeof(double))) == NULL ) { 
    fprintf(stderr,"alloc_diag2(): Cannot allocate ai_diag2 of size = %d \n", n_diag2_lines); fflush(stderr); fail( FAIL_BASIC, 0 ) ; 
  }
  if( (bi_diag2 = (double *) calloc(n_diag2_lines,sizeof(double))) == NULL ) { 
    fprintf(stderr,"alloc_diag2(): Cannot allocate bi_diag2 of size = %d \n", n_diag2_lines); fflush(stderr); fail( FAIL_BASIC, 0 ) ; 
  }
  if( (ci_diag2 = (double *) calloc(n_diag2_lines,sizeof(double))) == NULL ) { 
    fprintf(stderr,"alloc_diag2(): Cannot allocate ci_diag2 of size = %d \n", n_diag2_lines); fflush(stderr); fail( FAIL_BASIC, 0 ) ; 
  }
  if( (di_diag2 = (double *) calloc(n_diag2_lines,sizeof(double))) == NULL ) { 
    fprintf(stderr,"alloc_diag2(): Cannot allocate di_diag2 of size = %d \n", n_diag2_lines); fflush(stderr); fail( FAIL_BASIC, 0 ) ; 
  }
  if( (si_diag2 = (double *) calloc(n_diag2_lines,sizeof(double))) == NULL ) { 
    fprintf(stderr,"alloc_diag2(): Cannot allocate si_diag2 of size = %d \n", n_diag2_lines); fflush(stderr); fail( FAIL_BASIC, 0 ) ; 
  }
  if( (xi_diag2 = (double *) calloc(n_diag2_lines+1,sizeof(double))) == NULL ) { 
    fprintf(stderr,"alloc_diag2(): Cannot allocate xi_diag2 of size = %d \n", n_diag2_lines); fflush(stderr); fail( FAIL_BASIC, 0) ; 
  }
  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  dealloc_diag2():   -- unallocates all arrays used for COORD_DIAGONAL2
*****************************************************************************************/
void dealloc_diag2(void)
{ 
  free(ai_diag2);  free(bi_diag2);  free(ci_diag2);  free(di_diag2);  free(si_diag2);  free(xi_diag2);
  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  x2_of_xp2_diag2(): 
 ------------------
  -- Returns the non-normalized value of the coordinate function for COORD_DIAGONAL2
  -- makes sure that the discretization is periodic without making x2 periodic
  -- assumes equatorial symmetry;
*****************************************************************************************/
static double x2_of_xp2_diag2(double xp2)
{ 
  int j;
  double fout=0.;
  double intstep,intstepx,shift,xp2_old;
  xp2_old = xp2;
  xp2 = PERIODIC(xp2_old,xi_diag2[0],xi_diag2[n_diag2_lines]) ;
  for(j=0; j<n_diag2_lines; j++) { 
    int_step(xp2, xi_diag2[j], si_diag2[j], &intstep, &intstepx);
    fout +=   ai_diag2[j] * intstepx  +   bi_diag2[j] * intstep;
  }
  fout -=  len_diag2 * PERIODIC_SHIFT(xp2_old,xi_diag2[0],xi_diag2[n_diag2_lines]) ;
  return(fout);
}


/*****************************************************************************************/
/*****************************************************************************************
  dx2_dxp2_diag2(): 
 ------------------
  -- Returns the non-normalized value of the derivative of x2 w.r.t. xp2
      for COORD_DIAGONAL2
  -- makes sure that the discretization is periodic without making x2 periodic
  -- assumes equatorial symmetry;
*****************************************************************************************/
static double dx2_dxp2_diag2(double xp2)
{ 
  int j;
  double fout=0.;
  xp2 = PERIODIC(xp2,xi_diag2[0],xi_diag2[n_diag2_lines]) ;
  for(j=0; j<n_diag2_lines; j++) { 
    fout += (ai_diag2[j]*xp2 + bi_diag2[j]) * STEP(xp2,xi_diag2[j],si_diag2[j]);
  }
  return(fout);
}

/*****************************************************************************************/
/*****************************************************************************************
  dx2_dxp2_dxp2_diag2(): 
 ------------------
  -- Returns the non-normalized value of the 2nd derivative of x2 w.r.t. xp2
      for COORD_DIAGONAL2
  -- makes sure that the discretization is periodic without making x2 periodic
  -- assumes equatorial symmetry;
*****************************************************************************************/
static double dx2_dxp2_dxp2_diag2(double xp2)
{ 
  int j;
  double fout=0.;
  xp2 = PERIODIC(xp2,xi_diag2[0],xi_diag2[n_diag2_lines]) ;
  for(j=0; j<n_diag2_lines; j++) { 
    fout +=     ai_diag2[j]                    * STEP( xp2,xi_diag2[j],si_diag2[j]) 
             + (ai_diag2[j]*xp2 + bi_diag2[j]) * DSTEP(xp2,xi_diag2[j],si_diag2[j]) ; 
  }
  return(fout);
}

/*****************************************************************************************/
/*****************************************************************************************
  calc_ai_bi_di_diag2();
 ------------------
  --  Calculates the coefficients on the step functions that 
       result in a step-wise set of continuous linear curves. 
       (see comments in setup_diag2() for more details)

  --  We need only the slopes (ci), the locations of the 
       intersections (xi), and the value of the function 
       at x=xi[0]  (d0).   Calculates  ai,bi,di;  
   
  -- Assumes that all memory has been allocated;

*****************************************************************************************/
static void calc_ai_bi_di_diag2(void)
{
  int i;

  /* Adjust y-intercept to x0 */
  di_diag2[0] = d0_diag2 - ci_diag2[0]*xi_diag2[0]  ;

  /* Set di[] so that the lines are continuous at xi */
  for(i=1; i<n_diag2_lines; i++) { 
    di_diag2[i] = di_diag2[i-1] + xi_diag2[i]*(ci_diag2[i-1] - ci_diag2[i]);
  }

  calc_ai_bi_diag2();

  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  calc_ai_bi_diag2();
 ------------------
  --  Calculates the coefficients on the step functions that 
       result in a step-wise set of continuous linear curves. 
       (see comments in setup_diag2() for more details)

  --  We need only the slopes (ci), the locations of the 
       intersections (xi), and the value of the function 
       at x=xi[0]  (d0).   Calculates  ai,bi,di;  
   
  -- Assumes that all memory has been allocated;

*****************************************************************************************/
static void calc_ai_bi_diag2(void)
{
  int i,j;
  int e_ij;  /* transformation matrix from ci,di to ai,bi */

  for(i=0; i<n_diag2_lines; i++) { 
    ai_diag2[i] = bi_diag2[i] = 0.; 
    for(j=0; j<n_diag2_lines; j++) { 
      e_ij = (i == j) - (i == j+1) ;  /* linear transformation */
      ai_diag2[i] +=  e_ij * ci_diag2[j]; 
      bi_diag2[i] +=  e_ij * di_diag2[j]; 
    }
  }

  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  theoretical_norm_xp2():
 -----------------------
  -- returns the normalization of x2_of_xp2() that one should get in the limit that 
      step() approaches a true step function, or as   si -> infinity ; 

*****************************************************************************************/
static double theoretical_norm_diag2(void)
{
  int i; 
  double fout=0.; 
  for(i=0; i<n_diag2_lines; i++) { 
    fout += (xi_diag2[i+1] - xi_diag2[i]) * ( 0.5*(xi_diag2[i+1] + xi_diag2[i])*ci_diag2[i] 
					      + di_diag2[i]
					      )  ;
  }
  fout = M_PI/fout;
  return(fout);
}

/*****************************************************************************************/
/*****************************************************************************************
  int_step():
  ------------------
  -- calculates the integral of step() and x*step()   used for COORD_DIAGONAL2
  -- assumes that the step function is   atan(s*(x/x0-1))/M_PI + 0.5
*****************************************************************************************/
static void int_step(double x,double x0,double s,double *intstep, double *intstepx)
{
  double t1,t2,t3,t4,t5,t6,t7;
  x = PERIODIC(x,xi_diag2[0],xi_diag2[n_diag2_lines]) ;

  t1 = x - x0;
  t2 = x + x0;
  t3 = s*t1/x0;
  t4 = atan(t3);
  t5 = 0.5/(M_PI*s);
  t6 = x0*x0;
  t7 = s*s;

  *intstep = (s*(M_PI*x + 2*t1*t4) - x0*log(t7*t1*t1 + t6))*t5;
  
  *intstepx = (2*t4*(t6 + t7*t1*t2) 
	       + s*(t1*(-2*x0 + M_PI*s*t2) - 2*t6*log(1. + t3*t3))
	       )*t5*t5*M_PI;
  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  seesaw():
  ------------------
  -- intermediate function for COORD_DIAGONAL3 that resembles a seesaw pattern periodic 
      over  x0 to (xn+(xn-x0)) 
*****************************************************************************************/
static double seesaw(double x,double x0,double xn)
{
  double xn2,x1,retval;

  xn2 = 2*xn - x0;
  x1 = PERIODIC(x,x0,xn2);
  retval =  (x1 < xn)  ?  (x1-x0)  :  (xn2 - x1) ;
  return( retval );
}

/****************************************************************************************
 ===========================================
 set of functions for COORD_WARPED_SPHERICAL
 ===========================================

 -- stepfunc 
    --------
    \sigma in LaTeX doc

 -- square
    ------
    \tau in latex doc.
    Square function  with baseline = 0 , amplitude = 1
    also called "tophat" in latex document, "square" in latex document is something else;

 -- square_per
    ----------
    Periodic square function over  x \in [0,1]

 -- stepfunc_int
    ------------
    \Sigma in LaTeX doc 

 -- square_int
    ----------
    \Tau in LaTeX doc

 -- square_int_per
    --------------
    \tilde \Tau in LaTeX doc

 -- dsquare_per
    -----------
    \tilde \tau' in latex doc

 -- d2square_per
    -----------
    \tilde \tau''

 ****************************************************************************************/

static double stepfunc(double x, double x1, double hh) {
  return tanh(hh*(x - x1)) ;
}

static double square(double x, double x1, double h, double delta) {
  return 0.5*(stepfunc(x,x1-delta,h) - stepfunc(x,x1+delta,h)) ;
}

static double square_per(double x, double x1, double h, double delta) {
  return square(x,x1,h,delta) + square(x,x1-1.,h,delta) + square(x,x1+1.,h,delta) ;
}

static double stepfunc_int(double x, double x1, double h) {
  // double y = x-x1 ;
  // return ( log(cosh(h*(y-floor(y)-0.5)))/h ) ;
  // return ( log(cosh(h*(y)))/h  - log(cosh(h*(-x1)))/h ) ;
  return log(cosh(h*(x-x1)))/h ;
}

static double square_int(double x, double x1, double h, double delta) {
  return 0.5*(stepfunc_int(x,x1-delta,h) - stepfunc_int(x,x1+delta,h)) ;
}

static double square_int_per(double x, double x1, double h, double delta) {
  return square_int(x,x1,h,delta) + square_int(x,x1-1.,h,delta) + square_int(x,x1+1.,h,delta) ;
}


static double dsquare(double x, double x1, double h, double delta) {
  return 0.5*h*( sech2(h*(x-x1+delta)) - sech2(h*(x-x1-delta)) ) ;
}

static double dsquare_per(double x, double x1, double h, double delta) {
  return dsquare(x,x1,h,delta) + dsquare(x,x1-1.,h,delta) + dsquare(x,x1+1.,h,delta) ;
}


static double d2square(double x, double x1, double h, double delta) {
  return h*h * ( sech2(h*(x - x1 - delta)) * tanh(h*(x - x1 - delta)) 
                -sech2(h*(x - x1 + delta)) * tanh(h*(x - x1 + delta)) ) ;
}

static double d2square_per(double x, double x1, double h, double delta) {
  return d2square(x,x1,h,delta) + d2square(x,x1-1.,h,delta) + d2square(x,x1+1.,h,delta) ;
}



static void rfuncs(double *rf, double *drdy, double *d2rdy2, double *ar_i,
                   double y, double ybh_i, double Rout_i, double Rin_i, double br_i, double s_i) {

  *ar_i = ( Rout_i - Rin_i - br_i ) / ( 
                                        sinh(s_i * (1. - ybh_i))
                                      - sinh(-s_i*ybh_i)
                                      - s_i
                                       ) ;

//--orig    *rf = Rin_i + ( br_i - (*ar_i)*s_i ) * y
//--orig                + (*ar_i) * (   sinh(s_i * (y - ybh_i))
//--orig                           - sinh(-s_i * ybh_i) ) ;

  *rf = Rin_i + ( br_i ) * y
    + (*ar_i) * (  sinh_terms( s_i, y, ybh_i ) ) ;
 
  *drdy   = br_i + s_i * (*ar_i) * ( cosh( s_i*(y - ybh_i) ) - 1. ) ;
  *d2rdy2 = s_i*s_i * (*ar_i) * sinh( s_i*(y - ybh_i) ) ;

}

static double sinh_terms(double s, double y, double y0 )
{
  double x  = s * y; 
  double x0 = s * y0; 

  if( (fabs(x) < 2.e-2) && (fabs(x0) < 2.e-2) ) { 

    double t1 = x0 * x0;
    double t89 = (((7315660800 + (609638400 + (362880 * t1 + 20321280) * t1) * t1) * t1 + ((-7315660800 + (-1219276800 + (-60963840 + (-20160 * t1 - 1451520) * t1) * t1) * t1) * x0 + (2438553600 + (1219276800 + (101606400 + (60480 * t1 + 3386880) * t1) * t1) * t1 + ((-609638400 + (-101606400 + (-5080320 + (-1680 * t1 - 120960) * t1) * t1) * t1) * x0 + (121927680 + (60963840 + (5080320 + (3024 * t1 + 169344) * t1) * t1) * t1 + ((-20321280 + (-3386880 + (-169344 + (-56 * t1 - 4032) * t1) * t1) * t1) * x0 + (2903040 + (1451520 + (120960 + (72 * t1 + 4032) * t1) * t1) * t1 + ((-362880 + (-60480 + (-3024 + (-t1 - 72) * t1) * t1) * t1) * x0 + (40320 + (20160 + (1680 + (t1 + 56) * t1) * t1) * t1) * x) * x) * x) * x) * x) * x) * x) * x) * x) / 0.14631321600e11;

    return(t89); 
  }
  else { 
    return( sinh(s*(y-y0)) + sinh(x0) - x ) ; 
  }
}


/*****************************************************************************************/
/*****************************************************************************************
  update_coord_time_funcs():
  ------------------
  -- calculates the time-dependent variables for the coordinate system;
*****************************************************************************************/
void update_coord_time_funcs(double t_now)
{
  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;

  TRACE_BEG;

#if(!DYNAMIC_COORDINATES)
  TRACE_END;
  return;
#endif

#if( COORD_TYPE_CHOICE == COORD_DIAGONAL3_DYN_RAD ) 

# if( BBH_SPACETIME ) 
  /**********************************************************************************
     BBH Spacetime setup:  inner radial boundary follows the shrinking binary's orbit;
   **********************************************************************************/
  struct of_bbh_traj bbh_traj; 
  get_bbh_traj_data(&bbh_traj) ;

#if( METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW  ||\
     METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND     ||\
     METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND   )
  double r12_dot = bbh_traj.r12dot ;
#else
  double r12_dot  =  (  bbh_traj.v12x*(bbh_traj.xi1x-bbh_traj.xi2x) +
			bbh_traj.v12y*(bbh_traj.xi1y-bbh_traj.xi2y) +
			bbh_traj.v12z*(bbh_traj.xi1z-bbh_traj.xi2z))/bbh_traj.r12;
#endif

  /* Assumes that nz_params is set to the same time as t_now : */ 
  if( REL_DIFF_FUNC( (bbh_traj.tt), (t_now) ) > 1.e2*SMALL  ) { 
    fprintf(stderr,"update_coord_time_funcs(): BBH trajectory data seems out of sync with current time : %26.16e %26.16e \n", 
	    bbh_traj.tt, t_now ); 
    fflush(stderr); fail(FAIL_BASIC,0);
  }

  double reduction_factor = 10. / initial_bbh_separation;   /* using min. separation = 10M */
  coord_params->f0_r = 0.;
  coord_params->fa_r = 0.75;  
  Rin  = coord_params->fa_r * (initial_bbh_separation) ; 
  Rout = 13.*initial_bbh_separation / reduction_factor ;   
  coord_params->upsilon_r = Rout/Rin;
  coord_params->Rin_of_t = coord_params->fa_r * bbh_traj.r12; 
  coord_params->dln_Rin_of_t_dt = r12_dot  / bbh_traj.r12 ;  /* time-derivative of Rin_of_t  divided by Rin_of_t  */

# else 
  /**********************************************************************************
     Single BH Bondi setup:  -- a constant radial velocity set by dln_Rin_of_t_dt  
           at t=0 by init.*.c  and use that to update Rin : 
  **********************************************************************************/

  double reduction_factor = 10. / initial_bbh_separation;   /* using min. separation = 10M */
  coord_params->f0_r = 0.;
  coord_params->fa_r = 0.75;  
  Rin  = coord_params->fa_r * (initial_bbh_separation) ; 
  Rout = 13.*initial_bbh_separation / reduction_factor ;   
  coord_params->upsilon_r = Rout/Rin;

  double r12, r12_dot;
  const double factor = 64./5.;
  double f1_tmp     = initial_bbh_separation;
  double f1_dot_tmp = 0.;
  double tc         = pow(initial_bbh_separation, 4.) / factor  +  t_shrink_bbh;
  double delt       = fabs(tc - t_now); 
  double f2_tmp     = pow( factor*delt , 0.25 ); 
  double f2_dot_tmp =   -0.25*f2_tmp / delt ; 
  double t_match              = t_shrink_bbh; 
  const double step_strength  = 200.;
  double ftmp                 = step_strength*(t_now - t_match)/t_match;
  double tanh_term            = tanh(ftmp);
  double ftmp_dot             = step_strength/t_match;

  r12     = 0.5 * ( 
		   (f2_tmp - f1_tmp)*tanh_term + 
		   f1_tmp + f2_tmp 
		    );

  r12_dot = 0.5 * ( 
		   (f2_dot_tmp - f1_dot_tmp)*tanh_term + 
		   (f2_tmp - f1_tmp)*(1. - tanh_term*tanh_term)*ftmp_dot + 
		   f1_dot_tmp + f2_dot_tmp 
		    );

  coord_params->Rin_of_t = coord_params->fa_r * r12; 
  coord_params->dln_Rin_of_t_dt = r12_dot  / r12 ;  /* time-derivative of Rin_of_t  divided by Rin_of_t  */

# endif /* if( BBH_SPACETIME )  */


#elif( COORD_TYPE_CHOICE == COORD_WARPED_SPHERICAL )
  coord_params->t = t_now;
  // we first define the warping parameters:

#define delta0 0.1 
  delta_x1=delta0;
  delta_x2=delta0 ;
  delta_x3=0.18 ;
  delta_x4=0.18 ;
  delta_y1=0.18 ;
  delta_y2=0.18 ;
  delta_y3=delta0 ;
  delta_y4=delta0 ;
  delta_z1=0.4    ;
  delta_z2=0.4    ;

  a_x10=1.5 ;
  a_x20=1.5;
#if( EQUATORIAL_RUN ) 
  a_z10=0. ;
#else 
  a_z10=4.3 ;
#endif
  //a_z10=4.3 ;
  //a_z10=0. ;
  
  h_x1=20. ;
  h_x2=20 ;
  h_x3=20. ;
  h_x4=20 ;
  h_y1=20.   ;
  h_y2=20;
  h_z1=20. ;
  h_z2=20. ;
  h_y3=10. ;
  h_y4=10. ;

  //  we're keeping index '0' for entry '3' in the notes.
  s_[1]=1e-2 ;
  s_[2]=1e-2 ;
  /* s_[0]=4. ; */
  s_[0]=1e-2 ;
  br_[1]=6.5e0 ;
  br_[2]=6.5e0 ;
  /* br_[0]=40. ; */
  br_[0]=6.5e0 ;
  /***************************************/
  
#define Rin0  1.
#define Rout0 300.
  Rin_[1]=Rin0;
  Rin_[2]=Rin0;
  Rin_[0]=Rin0;
  Rout_[1]=Rout0;
  Rout_[2] = Rout0;
  Rout_[0] = Rout0;

  struct of_bbh_traj bbh_traj; 
  get_bbh_traj_data(&bbh_traj) ;

# if(1)
  double r12 = bbh_traj.r12 ;

#if( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW ||\
     METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND     ||\
     METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND  ||\
     METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP)
  double r12_dot = bbh_traj.r12dot ;
#else
  double r12_dot  =  (  bbh_traj.v12x*(bbh_traj.xi1x-bbh_traj.xi2x) +
			bbh_traj.v12y*(bbh_traj.xi1y-bbh_traj.xi2y) +
			bbh_traj.v12z*(bbh_traj.xi1z-bbh_traj.xi2z))/bbh_traj.r12;
#endif

  double phi1 = atan2( bbh_traj.xi1y, bbh_traj.xi1x ) ;
  if (phi1 < 0) phi1 += 2*M_PI ;
  double phi2 = atan2( bbh_traj.xi2y, bbh_traj.xi2x ) ;
  if (phi2 < 0) phi2 += 2*M_PI ;

  double omega1   = bbh_traj.omega  ;
  double omega2   = omega1 ;
# else
  bbh_traj.tt = t_now;
  double r12      = 20.;
  double r12_dot  = 0.;
  double omega1   = pow(20.,-1.5);
  double omega2   = omega1 ;

  double phi1 = fmod( omega1*t_now, 2.*M_PI );
  if (phi1 < 0) phi1 += 2*M_PI ;
  double phi2 = fmod( phi1 + M_PI , 2.*M_PI ); 
  if (phi2 < 0) phi2 += 2*M_PI ;

# endif

  /* Assumes that nz_params is set to the same time as t_now : */ 
  if( REL_DIFF_FUNC( (bbh_traj.tt), (t_now) ) > 1.e2*SMALL  ) { 
    fprintf(stderr,"update_coord_time_funcs(): BBH trajectory data seems out of sync with current time : %26.16e %26.16e \n", 
	    bbh_traj.tt, t_now ); 
    fflush(stderr); fail(FAIL_BASIC,0);
  }
  
  double x_1d[1];
  double dfdx_1d[1];
  int n, retval ;

  double dybhi_drbhi[3] ;


  // BHs' position in the xp coordinates
  xbh_[1] = phi1/(2.*M_PI) ;
  xbh_[2] = phi2/(2.*M_PI) ;

  zbh_[1]=0.5;
  zbh_[2]=zbh_[1];


  // for the 'y' coordinate we need a root-finding algorithm

  // rhs:
  rbh_[1] = m_bh2/m_bh_tot * r12 ;
  rbh_[2] = m_bh1/m_bh_tot * r12 ;
  rbh_[0] = MAX(rbh_[1], rbh_[2]) ;

  // initial guess
  ybh_[1] = 0.5 ;
  ybh_[2] = 0.5 ;
  ybh_[0] = 0.5 ;

  char valstr[100];

  if( FAST_AND_FURIOUS || (myid == special1_pid) ) { 
    for(n = 0; n<3; n++) {
      x_1d[0] = ybh_[n] ;
      retval = newt_raphs_warped_spherical(x_1d, dfdx_1d, 1, 
					   Rout_[n], Rin_[n], br_[n], s_[n], rbh_[n], 
					   newt_raphs_func);
      if( retval ) { 
	fprintf(stderr,"update_coord_time_funcs(): newt_raphs_warped_spherical() failed with error code=[%d], iter=[%d],  ybh_=[%26.16e,%26.16e,%26.16e]  rbh_=[%26.16e,%26.16e,%26.16e]  \n", 
		retval,n,ybh_[0],ybh_[1],ybh_[2],rbh_[0],rbh_[1],rbh_[2]); 
	fflush(stderr); fail(FAIL_BASIC,0);
      }
      ybh_[n] = x_1d[0];
      dybhi_drbhi[n] = 1./dfdx_1d[0] ; // we write the first derivative here, it will be needed
      //    fprintf(stdout,"y_bh[%d,%d,%d, %26.16e] iter=[%d] ybh=[%26.16e] rbh=[%26.16e]  retval = %d \n", myid,nstep,n_substep,t_now,n,ybh_[n],rbh_[n],retval);  fflush(stdout);
      //      sprintf(valstr,"ybh_[%d]",n);    mpi_global_compare_double(valstr,ybh_[n]);
    }
//  printf("\ny_bh = %.16lf %.16lf %.16lf \n", ybh_[1],  ybh_[2],  ybh_[0]) ;
//  printf("dybh_drbh = %.16lf %.16lf %.16lf  \n\n", dybhi_drbhi[1],  dybhi_drbhi[2],  dybhi_drbhi[0]) ;
//  fprintf(stdout,"y_bh retval = %d \n", retval);  fflush(stdout);
  }

#if( !FAST_AND_FURIOUS )
  sync_vect_from_rank(ybh_,3,special1_pid);
  sync_vect_from_rank(dybhi_drbhi,3,special1_pid);
#endif

  // time derivatives

  dxbh1_dt = omega1/(2.*M_PI) ;
  dxbh2_dt = omega2/(2.*M_PI) ;

  dybh1_dt =  m_bh2/m_bh_tot * r12_dot * dybhi_drbhi[1] ;
  dybh2_dt =  m_bh1/m_bh_tot * r12_dot * dybhi_drbhi[2] ;

  dybh3_dt = MAX( m_bh2/m_bh_tot, m_bh1/m_bh_tot ) * r12_dot * dybhi_drbhi[0] ;



  // second time derivatives

  // TODO! for now these are just placeholders
  double omega_dot = 0. ;
  double r12_dot_dot = 0. ;

  double d2ybhi_drbhi2[3] ;

  for(n=0; n<3; n++) {
    d2ybhi_drbhi2[n] =  (
                       ( -br_[n] - Rin_[n] + Rout_[n]) * s_[n]*s_[n] * (   
                                                                        sinh(s_[n]/2.) * ( s_[n]*ybh_[n] - sinh(s_[n]*ybh_[n]) )
                                                                            *( 
                                                                               - 2*s_[n] * cosh(s_[n]/2. - s_[n]*ybh_[n]) 
                                                                               + 6*sinh(s_[n]/2.)
                                                                               - sinh(s_[n]*(1.5 - 2*ybh_[n]))
                                                                               + sinh(s_[n]/2. - 2*s_[n]*ybh_[n])
                                                                               )
                                                                            - 2*(-1 + cosh(s_[n]*ybh_[n]))
                                                                               *(
                                                                                   cosh(s_[n]*ybh_[n])
                                                                                 - cosh(s_[n] - s_[n]*ybh_[n])
                                                                                 )*( - s_[n] + sinh(s_[n]*ybh_[n]) 
                                                                                     + sinh(s_[n] - s_[n]*ybh_[n]) )
                                                                            + sinh(s_[n]*ybh_[n])*pow(- s_[n] + sinh(s_[n]*ybh_[n]) + sinh(s_[n] - s_[n]*ybh_[n]),2)
                                                                           )
                        )
                         /pow(-s_[n] + sinh(s_[n]*ybh_[n]) + sinh(s_[n] - s_[n]*ybh_[n]),3) ;

    d2ybhi_drbhi2[n] *= -dybhi_drbhi[n] * dybhi_drbhi[n] * dybhi_drbhi[n] ;
  }


  d2xbh1_dt2 = omega_dot/(2.*M_PI) ;
  d2xbh2_dt2 = omega_dot/(2.*M_PI) ;

  d2ybh1_dt2 =  m_bh2/m_bh_tot * r12_dot_dot * dybhi_drbhi[1] + pow( m_bh2/m_bh_tot * r12_dot, 2) * d2ybhi_drbhi2[1] ;
  d2ybh2_dt2 =  m_bh1/m_bh_tot * r12_dot_dot * dybhi_drbhi[2] + pow( m_bh1/m_bh_tot * r12_dot, 2) * d2ybhi_drbhi2[2] ;

  d2ybh3_dt2 =  MAX( m_bh2/m_bh_tot, m_bh1/m_bh_tot) * r12_dot_dot * dybhi_drbhi[0] + pow( MAX( m_bh2/m_bh_tot, m_bh1/m_bh_tot) * r12_dot, 2) * d2ybhi_drbhi2[0] ;


  // MZ we fill the arrays now

  // first we fill all the time dependent (and grid-point independent) variables

  for(n=0; n<3; n++) {
    ar_[n] = ( Rout_[n] - Rin_[n] - br_[n] ) / (
                                                  sinh(s_[n] * (1. - ybh_[n]))
                                                - sinh(-s_[n]*ybh_[n])
                                                - s_[n]
                                                ) ;
  }


  coord_params->square_per_x11      = square_per(1,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_per_x12      = square_per(1,xbh_[2],h_x2,delta_x2) ;
  coord_params->square_per_x01      = square_per(0,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_per_x02      = square_per(0,xbh_[2],h_x2,delta_x2) ;
  
  coord_params->square_int_per_zbh1 = square_int_per(zbh_[1],zbh_[1],h_z1,delta_z1) ;
  coord_params->square_int_per_z0   = square_int_per(0,zbh_[1],h_z1,delta_z1) ;
  coord_params->square_int_per_z1   = square_int_per(1,zbh_[1],h_z1,delta_z1) ;

  coord_params->square_int_per_xbh1 = square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) ;
  coord_params->square_int_per_xbh2 = square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) ;

  coord_params->square_int_per_x11  = square_int_per(1,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_int_per_x12  = square_int_per(1,xbh_[2],h_x2,delta_x2) ;
  coord_params->square_int_per_x01  = square_int_per(0,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_int_per_x02  = square_int_per(0,xbh_[2],h_x2,delta_x2) ;

  int i,j,k,pos;
  double xp[NDIM] ;

  j = k = 0 ;
  N1ALL_LOOP POSLOOP {
    coord(i,j,k,pos,xp);
    coord_params->square_xp1_3[i][pos]  = square(xp[1],ybh_[1],h_y3,delta_y3) ;
    coord_params->square_xp1_4[i][pos]  = square(xp[1],ybh_[2],h_y4,delta_y4) ;

    coord_params->dsquare_xp1_3[i][pos] = dsquare(xp[1],ybh_[1],h_y3,delta_y3) ;
    coord_params->dsquare_xp1_4[i][pos] = dsquare(xp[1],ybh_[2],h_y4,delta_y4) ;

    for(n=0; n<3; n++) {
//--orig       coord_params->r_of_yt_[i][pos][n] =  Rin_[n] + ( br_[n] - ar_[n]*s_[n] ) * xp[1]
//--orig                                          + ar_[n] * (   sinh(s_[n] * (xp[1] - ybh_[n]))
//--orig                                                       - sinh(-s_[n] * ybh_[n]) ) ;

      coord_params->r_of_yt_[i][pos][n] =  Rin_[n] + ( br_[n] ) * xp[1]
      + ar_[n] * (  sinh_terms( s_[n], xp[1], ybh_[n] ) );

      coord_params->drdy_of_yt_[i][pos][n] =  br_[n] 
        + s_[n] * ar_[n] * ( cosh( s_[n]*(xp[1] - ybh_[n]) ) - 1. ) ;
    }
  }

  i = k = 0 ;
  N2ALL_LOOP POSLOOP {
    coord(i,j,k,pos,xp);
    coord_params->square_per_xp2_1[j][pos]     = square_per(xp[2],zbh_[1],h_z2,delta_z2) ;
    coord_params->square_per_xp2_2[j][pos]     = square_per(xp[2],zbh_[2],h_z2,delta_z2) ;

    coord_params->dsquare_per_xp2_1[j][pos]    = dsquare_per(xp[2],zbh_[1],h_z2,delta_z2) ;
    coord_params->dsquare_per_xp2_2[j][pos]    = dsquare_per(xp[2],zbh_[2],h_z2,delta_z2) ;

    coord_params->square_int_per_xp2_1[j][pos] = square_int_per(xp[2],zbh_[1],h_z2,delta_z2) ;
    coord_params->square_int_per_xp2_2[j][pos] = square_int_per(xp[2],zbh_[2],h_z2,delta_z2) ;
  }

  i = j = 0 ;
  N3ALL_LOOP POSLOOP {
    coord(i,j,k,pos,xp);
    coord_params->square_per_xp3_1[k][pos]     = square_per(xp[3],xbh_[1],h_x1,delta_x1) ;
    coord_params->square_per_xp3_2[k][pos]     = square_per(xp[3],xbh_[2],h_x2,delta_x2) ;
    coord_params->square_per_xp3_3[k][pos]     = square_per(xp[3],xbh_[1],h_x3,delta_x3) ;
    coord_params->square_per_xp3_4[k][pos]     = square_per(xp[3],xbh_[2],h_x4,delta_x4) ;

    coord_params->square_int_per_xp3_1[k][pos] = square_int_per(xp[3],xbh_[1],h_x1,delta_x1) ;
    coord_params->square_int_per_xp3_2[k][pos] = square_int_per(xp[3],xbh_[2],h_x2,delta_x2) ;
    coord_params->square_int_per_xp3_3[k][pos] = square_int_per(xp[3],xbh_[1],h_x3,delta_x3) ;
    coord_params->square_int_per_xp3_4[k][pos] = square_int_per(xp[3],xbh_[2],h_x4,delta_x4) ;

    coord_params->dsquare_per_xp3_3[k][pos]    = dsquare_per(xp[3],xbh_[1],h_x3,delta_x3) ;
    coord_params->dsquare_per_xp3_4[k][pos]    = dsquare_per(xp[3],xbh_[2],h_x4,delta_x4) ;
  }


  // we now check the errors we're committing
#if(0) // commenting this part out
  double err_xbh1, err_xbh2 ;
  double err_ybh1, err_ybh2 ;
  double x[NDIM], xploc[NDIM];

  xploc[0] = t;
  xploc[1] = ybh_[1]; 
  xploc[2] = zbh_[1]; 
  xploc[3] = xbh_[1]; 
  x_of_xp(x,xploc);

  err_xbh1 = phi1 - x[3];
  err_ybh1 = rbh_[1] - x[1];

  xploc[0] = t;
  xploc[1] = ybh_[2]; 
  xploc[2] = zbh_[2]; 
  xploc[3] = xbh_[2]; 
  x_of_xp(x,xploc);

  err_xbh2 = phi2 - x[3];
  err_ybh2 = rbh_[2] - x[1];

  err_xbh1  = (phi1 == 0) ?  fabs(err_xbh1) : fabs(err_xbh1/phi1);
  err_xbh2  = (phi2 == 0) ?  fabs(err_xbh2) : fabs(err_xbh2/phi2);
  fprintf(stdout,"error in xbhi = %26.16e, %26.16e  phis = %26.16e %26.16e %26.16e \n", err_xbh1, err_xbh2,phi1,phi2,x[3]);  
  fprintf(stdout,"error in ybhi = %26.16e, %26.16e  rs   = %26.16e %26.16e %26.16e \n", fabs(err_ybh1/rbh_[1]), fabs(err_ybh2/rbh_[2]),rbh_[1],rbh_[2],x[1]); 
  fflush(stdout);
#endif

#elif( COORD_TYPE_CHOICE == COORD_WARPED_CARTESIAN ) 
  coord_params->t = t_now;
  // we first define the warping parameters:

  /***************************************/
#define delta0 0.1 
  delta_x1=delta0;
  delta_x2=delta0;
  delta_x3=delta0;
  delta_x4=delta0;
  delta_y1=delta0;
  delta_y2=delta0;
  delta_y3=delta0;
  delta_y4=delta0;

  a_x10=1.;
  a_x20=1.;
  a_y10=1.;
  a_y20=1.;
  
  h_x1=20.;
  h_x2=20.;
  h_x3=20.;
  h_x4=20.;
  h_y1=20.;
  h_y2=20.;
  h_y3=20.;
  h_y4=20.;

  xmin=-1;   xmax=1;
  ymin=-0.5; ymax=0.5;

  /***************************************/
  
  double omega1 = 2.*M_PI/5.;
  double omega2 = omega1 ;

  double phi1 = fmod( omega1*t_now, 2.*M_PI );
  if (phi1 < 0) phi1 += 2*M_PI ;
  double phi2 = fmod( phi1 + M_PI , 2.*M_PI ); 
  if (phi2 < 0) phi2 += 2*M_PI ;

  // central warp location(s) in the xp coordinates
  /* xbh_[1] = 0.25 ; */
  /* xbh_[2] = 0.75 ; */

  /* ybh_[1] = 0.5 ; */
  /* ybh_[2] = 0.5 ; */

  xbh_[1] = 0.5 - 0.25*cos(phi1) ;
  xbh_[2] = 0.5 + 0.25*cos(phi1) ;

  ybh_[1] = 0.5 - 0.25*sin(phi1) ;
  ybh_[2] = 0.5 + 0.25*sin(phi1) ;

  /* // time derivatives */
  /* dxbh1_dt = 0. ; */
  /* dxbh2_dt = 0. ; */

  /* dybh1_dt = 0. ; */
  /* dybh2_dt = 0. ; */

  dxbh1_dt =  omega1*0.25*sin(phi1) ;
  dxbh2_dt = -omega1*0.25*sin(phi1) ;

  dybh1_dt = -omega1*0.25*cos(phi1) ;
  dybh2_dt =  omega1*0.25*cos(phi1) ;


  /* dxbh1_dt = omega1/(2.*M_PI) ; */
  /* dxbh2_dt = omega2/(2.*M_PI) ; */

  /* dybh1_dt =  m_bh2/m_bh_tot * r12_dot * dybhi_drbhi[1] ; */
  /* dybh2_dt =  m_bh1/m_bh_tot * r12_dot * dybhi_drbhi[2] ; */

  coord_params->square_per_x11      = square_per(1,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_per_x12      = square_per(1,xbh_[2],h_x2,delta_x2) ;
  coord_params->square_per_x01      = square_per(0,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_per_x02      = square_per(0,xbh_[2],h_x2,delta_x2) ;

  coord_params->square_per_y11      = square_per(1,ybh_[1],h_y1,delta_y1) ;
  coord_params->square_per_y12      = square_per(1,ybh_[2],h_y2,delta_y2) ;
  coord_params->square_per_y01      = square_per(0,ybh_[1],h_y1,delta_y1) ;
  coord_params->square_per_y02      = square_per(0,ybh_[2],h_y2,delta_y2) ;
  
  coord_params->square_int_per_xbh1 = square_int_per(xbh_[1],xbh_[1],h_x1,delta_x1) ;
  coord_params->square_int_per_xbh2 = square_int_per(xbh_[2],xbh_[2],h_x2,delta_x2) ;
  coord_params->square_int_per_x11  = square_int_per(1,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_int_per_x12  = square_int_per(1,xbh_[2],h_x2,delta_x2) ;
  coord_params->square_int_per_x01  = square_int_per(0,xbh_[1],h_x1,delta_x1) ;
  coord_params->square_int_per_x02  = square_int_per(0,xbh_[2],h_x2,delta_x2) ;

  coord_params->square_int_per_ybh1 = square_int_per(ybh_[1],ybh_[1],h_y1,delta_y1) ;
  coord_params->square_int_per_ybh2 = square_int_per(ybh_[2],ybh_[2],h_y2,delta_y2) ;
  coord_params->square_int_per_y11  = square_int_per(1,ybh_[1],h_y1,delta_y1) ;
  coord_params->square_int_per_y12  = square_int_per(1,ybh_[2],h_y2,delta_y2) ;
  coord_params->square_int_per_y01  = square_int_per(0,ybh_[1],h_y1,delta_y1) ;
  coord_params->square_int_per_y02  = square_int_per(0,ybh_[2],h_y2,delta_y2) ;

  int i,j,k,pos;
  double xp[NDIM] ;

  i = k = 0 ;
  N2ALL_LOOP POSLOOP {
    coord(i,j,k,pos,xp);

    coord_params->square_per_xp2_1[j][pos]  = square_per(xp[2],ybh_[1],h_y1,delta_y1) ;
    coord_params->square_per_xp2_2[j][pos]  = square_per(xp[2],ybh_[2],h_y2,delta_y2) ;
    coord_params->square_per_xp2_3[j][pos]  = square_per(xp[2],ybh_[1],h_y3,delta_y3) ;
    coord_params->square_per_xp2_4[j][pos]  = square_per(xp[2],ybh_[2],h_y4,delta_y4) ;

    coord_params->dsquare_per_xp2_3[j][pos] = dsquare_per(xp[2],ybh_[1],h_y3,delta_y3) ;
    coord_params->dsquare_per_xp2_4[j][pos] = dsquare_per(xp[2],ybh_[2],h_y4,delta_y4) ;

    coord_params->square_int_per_xp2_1[j][pos] = square_int_per(xp[2],ybh_[1],h_y1,delta_y1) ;
    coord_params->square_int_per_xp2_2[j][pos] = square_int_per(xp[2],ybh_[2],h_y2,delta_y2) ;
  }

  k = j = 0 ;
  N1ALL_LOOP POSLOOP {
    coord(i,j,k,pos,xp);

    coord_params->square_per_xp1_1[i][pos]     = square_per(xp[1],xbh_[1],h_x1,delta_x1) ;
    coord_params->square_per_xp1_2[i][pos]     = square_per(xp[1],xbh_[2],h_x2,delta_x2) ;

    coord_params->square_per_xp1_3[i][pos]     = square_per(xp[1],xbh_[1],h_x3,delta_x3) ;
    coord_params->square_per_xp1_4[i][pos]     = square_per(xp[1],xbh_[2],h_x4,delta_x4) ;

    coord_params->square_int_per_xp1_1[i][pos] = square_int_per(xp[1],xbh_[1],h_x1,delta_x1) ;
    coord_params->square_int_per_xp1_2[i][pos] = square_int_per(xp[1],xbh_[2],h_x2,delta_x2) ;

    coord_params->dsquare_per_xp1_3[i][pos]    = dsquare_per(xp[1],xbh_[1],h_x3,delta_x3) ;
    coord_params->dsquare_per_xp1_4[i][pos]    = dsquare_per(xp[1],xbh_[2],h_x4,delta_x4) ;
  }


#else 
  fprintf(stderr,"update_coord_time_funcs(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
  fflush(stderr); fail(FAIL_BASIC,0);
#endif

  TRACE_END;
 
  return;

}





static int newt_raphs_warped_spherical( FTYPE x[], FTYPE dfdx[], int n,
                                        double Rout_i, double Rin_i, double br_i, double s_i, double rbh_i, 
                                        void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
                                                       FTYPE [][NEWT_DIM], FTYPE *, 
                                                       FTYPE [], int,
                                                       double, double, double, double, double)
                                        )
{
  FTYPE f, df[NEWT_DIM], dx[NEWT_DIM], x_old[NEWT_DIM], x_older[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  int   keep_iterating, i_increase;


  // Initialize various parameters and variables:
  errx = 1. ; 
  /* df = f = 1.; */
  f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_older[id] = x_old[id] = x_orig[id] = x[id] ;

  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, df, n,
              Rout_i, Rin_i, br_i, s_i, rbh_i);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_older[id] = x_old[id];
      x_old[id] = x[id] ;
    }

    /* don't use line search : */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /* if we are cycling in a dimension, then take a half-step instead */
    for( id = 0; id < n ; id++) {
      if( x_older[id] == x[id] ) { 
	x[id] -= 0.5*dx[id];
      }
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);

    //    printf("errx = %lf \n", errx) ;


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  for( id = 0; id < n ; id++) {
    dfdx[id] = df[id]  ;
  }


  /*  Check for bad untrapped divergences : */
  //  if( (finite(f)==0) || (finite(df)==0)  ) {
  if( (finite(f)==0) ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}


/**********************************************************************/ 
/*********************************************************************************
   newt_raphs_func(): 

        -- calculates the residuals, and Newton step for newt_raphs_warped_spherical();
        -- for this method, x=ybh_i here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = df/dx;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void newt_raphs_func(FTYPE x[], FTYPE dx[], FTYPE resid[], 
                            FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE df[], int n, 
                            double Rout_i, double Rin_i, double br_i, double s_i, double rbh_i)
{

  double rf, drdy, d2rdy2, ar_i, ybh_i, dri_dybh_i ;

  ybh_i = x[0] ;

  // we compute r for the given guess of ybh_i
  rfuncs(&rf, &drdy, &d2rdy2, &ar_i, ybh_i, ybh_i, Rout_i, Rin_i, br_i, s_i) ;

  resid[0] = rf - rbh_i ;


  dri_dybh_i = br_i + ( 
                        (br_i + Rin_i - Rout_i)*s_i*( -s_i - s_i*(-1. + ybh_i)*cosh(s_i*ybh_i)
                                                           + s_i*ybh_i*cosh(s_i - s_i*ybh_i)
                                                           - sinh(s_i) + sinh(s_i*ybh_i) + sinh(s_i - s_i*ybh_i) )
                        )
                      / pow(-s_i + sinh(s_i*ybh_i) + sinh(s_i - s_i*ybh_i),2) ; 

  jac[0][0] = dri_dybh_i ;

  dx[0] = -resid[0]/dri_dybh_i;


  *f = 0.5*resid[0]*resid[0];
  /* *df = -2. * (*f); */

  df[0] = dri_dybh_i ;

}

//Dennis
/**********************************************************************/ 
/*********************************************************************************
   jac_xBL_to_xIZ(): 
        -- Generated by maple script BL2IZ.ms with mini-disk notes 
        -- Calculates the Jacobian to take a tensor from the x_{BL} coordinate
            basis to the x_{IZ} coordinate basis.
        -- Here the x_{IZ} are Cook-Shield Horizon penetrating coordinates
        -- x_{IZ}(x_{BL}) are presented in Gallouin et al 2012
        -- Let (T,X,Y,Z)   be x_{IZ} and a be the BH spin
        -- Let (t,r,th,ph) be x_{BL} 
        -- r_{\pm} = M \pm \sqrt{M^2 - a^2}

           --         T = t_{BL} + ( r^2_+ + a^2 ) / ( r_+ - r_- ) ln| (r_{BL} - r_+)/(r_{BL} - r_- ) |
           --    X + iY = ( r_{BL} - M + i a ) exp( i \phi_{IK} ) * sin(\theta_{BL})
           --         Z = (r_{BL} - M) * cos(\theta_{BL})
           -- \phi_{IK} = \phi_{BL} + ( a / ( r_+ - r_- ) ) * ln | ( r_{BL} - r_+ ) / ( r_{BL} - r_- ) |

     Arguments:
         Mass = BH mass ( either m_bh1 or m_bh2 );
         Spin = BH spin ( either a_bh1 or a_bh2 );
          xBL = BL coordinates in terms of Mass;
    dxIZ_dxBL = \frac{ \partial x^{\alpha}_{IZ} }{\partial x^{\beta}_{BL}} ( x_{BL} )
 *********************************************************************************/

void jac_xBL_to_xIZ( double Mass, double Spin, double xBL[NDIM], double dxIZ_dxBL[NDIM][NDIM])
{

  double t1;
  double t2;  
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double r, th, ph;
  double cg, cg1, cg3, cg5, cg7;

  
  /* Set the BL coordinates from the xBL array 
     -- Only need spacial coordinates for Jacobian */
  r  = xBL[RR];  
  th = xBL[TH];  
  ph = xBL[PH];
  
  /* Define the cg terms that BL2IZ.ms renamed as a 
     result of assuming the quantities are real */
  cg  = Mass;     
  cg1 = Spin;
  cg3 = ph;
  cg5 = r;
  cg7 = th;

  /* Calculate dxIZ_dxBL as a function of known BL
     coordinates */
  t1 = pow(cg1, 0.2e1);
  t2 = pow(cg, 0.2e1) - t1;
  t3 = pow(t2, -0.1e1 / 0.2e1);
  t2 = t2 * t3;
  t4 = -cg5 + cg + t2;
  t5 = -cg5 + cg - t2;
  t6 = 0.1e1 / t5;
  t7 = t4 * t6;
  t8 = fabs(t7) / t7;
  t4 = 0.1e1 / t4;
  t4 = fabs(t4 * t5);
  t5 = pow(t6, 0.2e1);
  t6 = sin(cg7);
  t3 = (cg3 * t2 + cg1 * log(fabs(t7)) / 0.2e1) * t3;
  t7 = cos(t3);
  t3 = sin(t3);
  t1 = pow(cg5, 0.2e1) + 0.2e1 * cg5 * t2 + 0.2e1 * (-cg5 + cg - t2) * cg - t1;
  t9 = -cg5 + cg;
  t10 = t9 * t3 - cg1 * t7;
  t11 = cos(cg7);
  t12 = t9 * t7;
  t13 = cg1 * t3;
  t14 = t12 + t13;
  dxIZ_dxBL[0][0] = 1;
  dxIZ_dxBL[0][1] = 0.2e1 * cg * (cg + t2) * t8 * t4 * t5;
  dxIZ_dxBL[0][2] = 0.;
  dxIZ_dxBL[0][3] = 0.;
  dxIZ_dxBL[1][0] = 0.;
  dxIZ_dxBL[1][1] = t6 * (cg1 * t10 * t4 * t8 + t7 * t1) * t5;
  dxIZ_dxBL[1][2] = -t11 * t14;
  dxIZ_dxBL[1][3] = t6 * t10;
  dxIZ_dxBL[2][0] = 0.;
  dxIZ_dxBL[2][1] = t6 * (cg1 * (-t12 - t13) * t4 * t8 + t3 * t1) * t5;
  dxIZ_dxBL[2][2] = -t11 * t10;
  dxIZ_dxBL[2][3] = -t6 * t14;
  dxIZ_dxBL[3][0] = 0.;
  dxIZ_dxBL[3][1] = t11;
  dxIZ_dxBL[3][2] = t9 * t6;
  dxIZ_dxBL[3][3] = 0.;
  
  return;
}


/**********************************************************************/ 
/*********************************************************************************
   jac_xIZ_to_xBL():
        -- Generated by maple script IZ2BL_approx.ms with mini-disk notes 
        -- Calculates the Jacobian to take a tensor from the x_{IZ} coordinate
            basis to the x_{BL} coordinate basis.
        -- Here the x_{IZ} are Cook-Shield Horizon penetrating coordinates
        -- x_{IZ}(x_{BL}) are presented in Gallouin et al 2012
        -- Let (T,X,Y,Z)   be x_{IZ} and a be the BH spin
        -- Let (t,r,th,ph) be x_{BL} 
        -- r_{\pm} = M \pm \sqrt{M^2 - a^2}

           --         T = t_{BL} + ( r^2_+ + a^2 ) / ( r_+ - r_- ) ln| (r_{BL} - r_+)/(r_{BL} - r_- ) |
           --    X + iY = ( r_{BL} - M + i a ) exp( i \phi_{IK} ) * sin(\theta_{BL})
           --         Z = (r_{BL} - M) * cos(\theta_{BL})
           -- \phi_{IK} = \phi_{BL} + ( a / ( r_+ - r_- ) ) * ln | ( r_{BL} - r_+ ) / ( r_{BL} - r_- ) |

        --Here we actually use the approximately inverted equations for (t,r,theta,phi) provided in
          Hergt and Schafer 2008 Phys Rev D 77 104001 equations 70-78

     Arguments:
         Mass = BH mass ( either m_bh1 or m_bh2 );
         Spin = BH spin ( either a_bh1 or a_bh2 );
          xIZ = IZ coordinates in terms of Mass;
    dxBL_dxIZ = \frac{ \partial x^{\alpha}_{BL} }{\partial x^{\beta}_{IZ}} ( x_{IZ} )
 *********************************************************************************/

void jac_xIZ_to_xBL( double Mass, double Spin, double xIZ[NDIM], double dxBL_dxIZ[NDIM][NDIM])
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t25;

  double X,Y,Z;
  double cg,cg1,cg3;
  /* Set the BL coordinates from the xBL array 
     -- Only need spacial coordinates for Jacobian */
  X = xIZ[XX];  
  Y = xIZ[YY];  
  Z = xIZ[ZZ];
  
  cg  = X;
  cg1 = Y;
  cg3 = Z;
    
  /* Calculate dxBL_dxIZ as a function of known IZ
     coordinates */
  t1 = Spin * Spin;
  t2 = Mass * Mass - t1;
  t3 = pow(t2, -0.1e1 / 0.2e1);
  t2 = t2 * t3;
  t4 = Mass + t2;
  t5 = sqrt(0.2e1);
  t6 = pow(cg3, 0.2e1);
  t7 = pow(cg, 0.2e1);
  t8 = pow(cg1, 0.2e1);
  t9 = pow(t7, 0.2e1) + pow(t1, 0.2e1) + pow(t6, 0.2e1) + pow(t8, 0.2e1) + 0.2e1 * (t7 - t1) * t8 - 0.2e1 * t7 * t1 + 0.2e1 * (t7 + t1 + t8) * t6;
  t10 = pow(t9, -0.1e1 / 0.2e1);
  t9 = t7 + t8 + t6 - t1 + t9 * t10;
  t11 = pow(t9, -0.3e1 / 0.2e1);
  t12 = t9 * t11;
  t13 = t7 + t8 + t6 - t1;
  t14 = 0.1e1 / 0.2e1;
  t10 = t14 * t10;
  t15 = 0.2e1 * cg + 0.4e1 * t10 * cg * t13;
  t16 = t14 * t5;
  t17 = t16 * pow(t9, 0.2e1) * t11;
  t18 = t17 + t2;
  t2 = t17 - t2;
  t17 = 0.1e1 / t18;
  t19 = 0.1e1 - t17 * t2;
  t20 = t5 * t12 / 0.4e1;
  t21 = t20 * t15;
  t22 = t21 * t17 * t19;
  t13 = 0.2e1 * cg1 + 0.4e1 * t10 * cg1 * t13;
  t23 = t20 * t13;
  t24 = t23 * t17 * t19;
  t7 = 0.2e1 * cg3 + 0.4e1 * t10 * cg3 * (t7 + t8 + t6 + t1);
  t10 = t20 * t7;
  t17 = t10 * t17 * t19;
  t2 = 0.1e1 / t2;
  t4 = t14 * (pow(t4, 0.2e1) + t1) * t3;
  t9 = 0.1e1 / t9;
  t6 = 0.1e1 - 0.2e1 * t6 * t9;
  t6 = pow(t6, -0.1e1 / 0.2e1);
  t16 = t16 * cg3 * t11;
  t19 = 0.1e1 / cg;
  t20 = pow(t19, 0.2e1);
  t8 = 0.1e1 + t8 * t20;
  t1 = 0.1e1 + 0.2e1 * t1 * t9;
  t8 = 0.1e1 / t8;
  t1 = 0.1e1 / t1;
  t9 = t5 * t11;
  t25 = t14 * Spin;
  dxBL_dxIZ[0][0] = 1.;
  dxBL_dxIZ[0][1] = -t4 * t22 * t2 * t18;
  dxBL_dxIZ[0][2] = -t4 * t24 * t2 * t18;
  dxBL_dxIZ[0][3] = -t4 * t17 * t2 * t18;
  dxBL_dxIZ[1][0] = 0.;
  dxBL_dxIZ[1][1] = t21;
  dxBL_dxIZ[1][2] = t23;
  dxBL_dxIZ[1][3] = t10;
  dxBL_dxIZ[2][0] = 0.;
  dxBL_dxIZ[2][1] = t16 * t15 * t6;
  dxBL_dxIZ[2][2] = t16 * t13 * t6;
  dxBL_dxIZ[2][3] = -t5 * (t12 - t14 * cg3 * t11 * t7) * t6;
  dxBL_dxIZ[3][0] = 0.;
  dxBL_dxIZ[3][1] = t25 * (t9 * t15 * t1 - t3 * t22 * t2 * t18) - cg1 * t20 * t8;
  dxBL_dxIZ[3][2] = t25 * (t9 * t13 * t1 - t3 * t24 * t2 * t18) + t19 * t8;
  dxBL_dxIZ[3][3] = t25 * (t9 * t7 * t1 - t3 * t17 * t2 * t18);

  
  return;
}


/**********************************************************************/ 
/*********************************************************************************
   jac_xNZ_to_xIZ(): 
        -- Generated by maple script CS2PN.mpl with mini-disk notes 
        -- Calculates the Jacobian to take a tensor from the x_{NZ} coordinate
            basis to the x_{IZ} coordinate basis.
        -- Here the x_{IZ} are Cook-Scheel Horizon penetrating coordinates

     Arguments:
          xNZ = NZ Cartesian coordinates in terms of Mass;
    dxIZ_dxBL = \frac{ \partial x^{\alpha}_{IZ} }{\partial x^{\beta}_{NZ}} ( x_{NZ} )
 *********************************************************************************/
 void jac_xNZ_to_xIZ( int BH, double xNZ[NDIM], double dxIZ_dxNZ[NDIM][NDIM])
{
  double m1,m2,phi,b,omega,r12dot;
  double t,x,y,z;
  double t1,t2,t3,t4,t5,t6,t7,t8,t9;
  double t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
  double t20,t21,t22,t23,t24,t25,t26,t27,t28,t29;
  double t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;
  double t40,t41,t42,t43,t44,t45,t46,t47,t48,t49;
  double t50,t51,t52,t53,t54,t55,t56,t57,t58,t59;
  double t60,t61,t62,t63,t64,t65,t66,t67,t68,t69;
  double t70,t71,t72,t73,t74,t75,t76,t77,t78,t79;
  double t80,t81,t82,t83,t84,t85,t86,t87,t88,t89;
  double t90,t91,t92,t93,t94,t95,t96,t97,t98,t99;
  double t101,t102,t121;
  double t104,t105,t106,t108,t109,t111,t117,t127;
  double t103, t122, t125, t129, t130, t134, t138;
  double t139, t143, t145, t159, t161, t164, t165;
  double t167, t169, t173, t177, t178, t180, t200;
  double t201, t202, t250, t267, t277, t110, t144;
  double t148, t100, t118, t141, t146, t150;

  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  struct of_bbh_traj bbh_traj;
  get_bbh_traj_data(&bbh_traj);

  if( BH == MBH1 ) {
    t = xNZ[TT]; x = xNZ[XX];
    y = xNZ[YY]; z = xNZ[ZZ];
    
    /* Set BH Traj Values */
    m1    = m_bh1;
    m2    = m_bh2;
    phi   = bbh_traj.phi;
    b     = bbh_traj.r12;
    r12dot= bbh_traj.r12dot;
    omega = bbh_traj.omega;
  }
  else if( BH == MBH2 ) {
    t = xNZ[TT] ; x = -xNZ[XX];
    y = -xNZ[YY]; z = xNZ[ZZ];
    
    /* Set BH Traj Values */
    m1    = m_bh2;
    m2    = m_bh1;
    phi   = bbh_traj.phi;
    b     = bbh_traj.r12;
    r12dot= bbh_traj.r12dot;
    omega = bbh_traj.omega;
  }


  /* Maple Generated Code */
  t1 = b * b;
  t2 = t1 * b;
  t3 = 0.1e1 / t2;
  t4 = m2 * t3;
  t5 = m2 * b;
  t6 = 0.1e1 / M;
  t7 = cos(phi);
  t8 = t7 * t7;
  t11 = sin(phi);
  t12 = t11 * t11;
  t13 = t6 * t12;
  t16 = pow(-t5 * t6 * t8 - t5 * t13, 0.2e1);
  t17 = m2 * m2;
  t18 = t17 * b;
  t19 = M * M;
  t20 = 0.1e1 / t19;
  t23 = t18 * t20 * t12 * r12dot;
  t24 = t20 * r12dot;
  t25 = t24 * t8;
  t41 = t17 * m2;
  t42 = t41 * t20;
  t43 = t42 * t12;
  t46 = 0.1e1 / t1;
  t47 = t46 * r12dot;
  t57 = t1 * t1;
  t59 = m2 / t57;
  t60 = t6 * t7;
  t63 = pow(x - t5 * t60, 0.2e1);
  t64 = t6 * t11;
  t67 = pow(y - t5 * t64, 0.2e1);
  t68 = z * z;
  t69 = t63 + t67 + t68;
  t73 = sqrt(t69);
  t78 = 0.1e1 / t73;
  t80 = t41 * t1;
  t82 = 0.1e1 / t19 / M;
  t83 = t12 * t12;
  t85 = t82 * t83 * t11;
  t89 = t41 * t2;
  t90 = t82 * t83;
  t94 = t12 * t11;
  t95 = t82 * t94;
  t99 = t82 * t12;
  t103 = t82 * r12dot;
  t122 = t17 * t6;
  t125 = -t122 - 0.4e1 * t122 * t12;
  t129 = m2 * t46;
  t130 = m2 * t6;
  t134 = -0.3e1 * t130 + 0.6e1 * t130 * t12;
  t138 = t17 * t3;
  t139 = r12dot * t6;
  t143 = t17 * t46;
  t145 = t6 * omega * t7;
  t159 = t11 * omega;
  t161 = m2 / b;
  t164 = sqrt(t161);
  t165 = 0.1e1 / t164;
  t167 = sqrt(t130);
  t169 = t20 * t83;
  t173 = t17 * t1;
  t177 = 0.6e1 * t173 * t20 * t94 * omega;
  t178 = 0.9e1 * t23;
  t180 = t159 * t7;
  t200 = m2 * r12dot;
  t201 = t200 * t64;
  t202 = t5 * t145;
  t250 = t64 * omega;
  t267 = t7 * omega;
  t277 = b * t20 * r12dot;
  dxIZ_dxNZ[0][0] = (-t4 * t16 * (0.2e1 * t23 + 0.2e1 * t18 * t25) / 0.2e1 + (0.1e1 + 0.5e1 / 0.128e3 * (0.2e1 * m1 + 0.3e1 * m2) * t20 / m1 * t1 * r12dot + ((t42 + t43) * t12 * t47 + ((t42 + 0.2e1 * t43) * t46 * r12dot + t41 * t46 * t25) * t8 - t59 * r12dot * t69) * t73) * t73) * t78 + ((t4 * (0.9e1 * t80 * t85 * r12dot + (0.3e1 * t89 * t90 * omega + (0.18e2 * t80 * t95 * r12dot + (0.6e1 * t89 * t99 * omega + (0.9e1 * t80 * t103 * t11 + 0.3e1 * t89 * t82 * omega * t7) * t7) * t7) * t7) * t7) / 0.3e1 + (t125 * t11 * t3 * r12dot + (t129 * t134 * omega / 0.3e1 + (-0.4e1 * t138 * t139 * t11 + 0.2e1 * t143 * t145) * t7) * t7) * t69) * t78 + (t7 * m2 * t47 / 0.2e1 + t159 * t161) * t165 * t167 + ((t4 * (-0.9e1 * t18 * t169 * r12dot + (-t177 + (-t178 - 0.6e1 * t173 * t20 * t180) * t7) * t7) / 0.3e1 + (-0.2e1 * t4 * t180 + 0.3e1 * t59 * t12 * r12dot) * t69) * t78 - t4 * t12 * t78 * (-0.2e1 * t201 - 0.2e1 * t202) * y / 0.2e1) * y) * y + ((t4 * (-0.3e1 * t89 * t85 * omega + (0.9e1 * t80 * t90 * r12dot + (-0.6e1 * t89 * t95 * omega + (0.18e2 * t80 * t99 * r12dot + (-0.3e1 * t89 * t82 * t11 * omega + 0.9e1 * t80 * t103 * t7) * t7) * t7) * t7) * t7) / 0.3e1 + (-t129 * t134 * t11 * omega / 0.3e1 + (t125 * t3 * r12dot + (-0.2e1 * t143 * t250 - 0.4e1 * t138 * t139 * t7) * t7) * t7) * t69) * t78 + (-t11 * m2 * t47 / 0.2e1 + t267 * t161) * t165 * t167 + ((t4 * (0.6e1 * t173 * t169 * omega + (-0.18e2 * t94 * t17 * t277 + (-0.18e2 * t11 * t17 * t277 - 0.6e1 * t173 * t20 * omega * t7) * t8) * t7) / 0.3e1 + (0.2e1 * t4 * t12 * omega + (0.6e1 * t59 * t11 * r12dot - 0.2e1 * t4 * t267) * t7) * t69) * t78 + t4 * (-0.3e1 * t5 * t6 * t94 * omega + (0.9e1 * t200 * t13 + 0.6e1 * t5 * t6 * t180) * t7) * t78 * y / 0.3e1) * y + ((t4 * (t177 + (-t178 + (0.6e1 * t173 * t20 * t11 * omega - 0.9e1 * t18 * t24 * t7) * t7) * t7) * t7 / 0.3e1 + (0.2e1 * t4 * t159 + 0.3e1 * t59 * r12dot * t7) * t7 * t69) * t78 + t4 * (-0.6e1 * t5 * t13 * omega + (0.9e1 * t201 + 0.3e1 * t202) * t7) * t7 * t78 * y / 0.3e1 - t4 * t8 * t78 * (-0.2e1 * t200 * t60 + 0.2e1 * t5 * t250) * x / 0.2e1) * x) * x;

  t1 = b * b;
  t4 = m2 / t1 / b;
  t5 = m2 * b;
  t6 = 0.1e1 / M;
  t7 = cos(phi);
  t8 = t7 * t7;
  t10 = t5 * t6 * t8;
  t11 = sin(phi);
  t12 = t11 * t11;
  t14 = t5 * t6 * t12;
  t15 = -t10 - t14;
  t16 = t15 * t15;
  t19 = b * t6 * t7;
  t22 = m2 * t6;
  t34 = pow(x - t5 * t6 * t7, 0.2e1);
  t38 = pow(y - t5 * t6 * t11, 0.2e1);
  t39 = z * z;
  t40 = t34 + t38 + t39;
  t43 = sqrt(t40);
  t44 = 0.1e1 / t43;
  t50 = sqrt(m2 / b);
  t51 = sqrt(t22);
  t57 = t11 * t7;
  t64 = m2 * m2;
  t74 = t64 * t1;
  t75 = M * M;
  t76 = 0.1e1 / t75;
  t77 = t12 * t12;
  dxIZ_dxNZ[0][1] = t4 * (0.3e1 * t16 * m2 * t19 + ((-0.3e1 * t22 + 0.6e1 * t22 * t12) * b + 0.6e1 * t10) * t7 * t40) * t44 / 0.3e1 + t50 * t51 * t11 + (t4 * (0.6e1 * t15 * t11 * m2 * t19 - 0.6e1 * t57 * t40) * t44 / 0.3e1 + t64 / t1 * t12 * t44 * t6 * t7 * y) * y + (t4 * ((-0.3e1 * t74 * t76 * t77 + (-0.12e2 * t74 * t76 * t12 - 0.9e1 * t74 * t76 * t8) * t8 + (0.3e1 - 0.6e1 * t8) * t40) * t44 + ((0.6e1 * t5 * t6 * t12 * t11 + 0.12e2 * t11 * t8 * t5 * t6) * t44 - 0.3e1 * t12 * t44 * y) * y) / 0.3e1 + (t4 * ((0.6e1 * t14 + 0.9e1 * t10) * t7 * t44 - 0.6e1 * t57 * t44 * y) / 0.3e1 - t4 * t8 * t44 * x) * x) * x;

  t1 = b * b;
  t4 = m2 / t1 / b;
  t5 = m2 * b;
  t6 = 0.1e1 / M;
  t7 = cos(phi);
  t8 = t7 * t7;
  t10 = t5 * t6 * t8;
  t11 = sin(phi);
  t12 = t11 * t11;
  t14 = t5 * t6 * t12;
  t15 = -t10 - t14;
  t16 = t15 * t15;
  t18 = b * t6;
  t22 = m2 * t6;
  t32 = 0.6e1 * t11 * t8 * t5 * t6;
  t37 = pow(x - t5 * t6 * t7, 0.2e1);
  t40 = y - t5 * t6 * t11;
  t41 = t40 * t40;
  t42 = z * z;
  t43 = t37 + t41 + t42;
  t46 = sqrt(t43);
  t47 = 0.1e1 / t46;
  t53 = sqrt(m2 / b);
  t54 = sqrt(t22);
  t57 = m2 * m2;
  t58 = t57 * t1;
  t59 = M * M;
  t60 = 0.1e1 / t59;
  t61 = t12 * t12;
  t99 = t11 * t7;
  dxIZ_dxNZ[0][2] = t4 * (0.3e1 * t16 * m2 * t18 * t11 + ((-0.3e1 * t22 + 0.6e1 * t22 * t12) * t11 * b + t32) * t43) * t47 / 0.3e1 - t53 * t54 * t7 + (t4 * (-0.9e1 * t58 * t60 * t61 + (-0.12e2 * t58 * t60 * t12 - 0.3e1 * t58 * t60 * t8) * t8 + (0.3e1 - 0.6e1 * t12) * t43) * t47 / 0.3e1 + (t4 * (0.9e1 * t5 * t6 * t12 * t11 + t32) * t47 / 0.3e1 - t4 * t12 * t47 * y) * y) * y + (t4 * ((0.6e1 * t15 * t11 * m2 * t18 * t7 - 0.6e1 * t99 * t43) * t47 + ((0.12e2 * t14 + 0.6e1 * t10) * t7 * t47 - 0.6e1 * t99 * t47 * y) * y) / 0.3e1 - t4 * t8 * t47 * t40 * x) * x;

  t1 = b * b;
  t5 = m2 * b;
  t6 = 0.1e1 / M;
  t7 = cos(phi);
  t8 = t7 * t7;
  t10 = t5 * t6 * t8;
  t11 = sin(phi);
  t12 = t11 * t11;
  t14 = t5 * t6 * t12;
  t15 = -t10 - t14;
  t16 = t15 * t15;
  t20 = pow(x - t5 * t6 * t7, 0.2e1);
  t24 = pow(y - t5 * t6 * t11, 0.2e1);
  t25 = z * z;
  t29 = sqrt(t20 + t24 + t25);
  t30 = 0.1e1 / t29;
  t34 = t30 * z;
  dxIZ_dxNZ[0][3] = m2 / t1 / b * ((-0.3e1 * t16 + 0.3e1 * t20 + 0.3e1 * t24 + 0.3e1 * t25) * t30 * z + (-0.6e1 * t15 * t11 * t34 - 0.3e1 * t12 * t30 * z * y) * y + (-0.6e1 * (-t10 - t14 + t11 * y) * t7 * t34 - 0.3e1 * t8 * t30 * z * x) * x) / 0.3e1;

  t1 = 0.1e1 / M;
  t2 = m2 * t1;
  t3 = m2 * m2;
  t4 = M * M;
  t5 = 0.1e1 / t4;
  t7 = sin(phi);
  t8 = t7 * t7;
  t9 = t3 * t5 * t8;
  t15 = t1 * t7;
  t21 = t3 * m2;
  t23 = t7 * omega;
  t24 = cos(phi);
  t25 = t23 * t24;
  t30 = b * b;
  t31 = 0.1e1 / t30;
  t32 = m2 * t31;
  t33 = t2 * t8;
  t34 = -0.1e1 / 0.2e1 - t33;
  t39 = t3 * t1 * t8;
  t41 = -m2 - 0.2e1 * t39;
  t43 = 0.1e1 / t30 / b;
  t46 = t3 * t31;
  t48 = t46 * t15 * omega;
  t53 = 0.2e1 * t3 * t43 * t1 * r12dot * t24;
  t59 = z * z;
  t61 = 0.1e1 / b;
  t62 = m2 * t61;
  t64 = t2 / 0.2e1 + t9;
  t80 = t21 * t61;
  t90 = m2 * t43;
  t92 = t90 * t8 * omega;
  t93 = t30 * t30;
  t95 = m2 / t93;
  t101 = 0.3e1 * t95 * t7 * r12dot - t90 * omega * t24;
  t103 = t92 + t101 * t24;
  t110 = 0.4e1 * t39;
  t144 = t90 * t23;
  t148 = 0.3e1 * t95 * r12dot * t24;
  dxIZ_dxNZ[1][0] = (m2 * (t2 + t9 / 0.2e1) * t7 + m2 * b * t15) * omega + (-m2 * r12dot * t1 + t21 * t5 * t25 / 0.2e1) * t24 + (t32 * t34 * t7 * omega + (t41 * t43 * r12dot + (-t48 - t53) * t24) * t24) * t59 + (t62 * t64 * t8 * omega + (t32 * t64 * t7 * r12dot + (-t3 * t61 * t1 * omega / 0.2e1 + (t21 * t31 * t5 * t7 * r12dot - t80 * t5 * omega * t24) * t24) * t24) * t24 + t103 * t59 + (t32 * (-0.1e1 / 0.2e1 - 0.2e1 * t33) * t7 * omega + ((-m2 - t110) * t43 * r12dot + (t48 - t53) * t24) * t24 + t103 * y) * y) * y + (-t32 * (0.1e1 + (0.3e1 / 0.2e1 * t2 + t9) * t8) * r12dot + (t62 * (t2 + 0.2e1 * t9) * t7 * omega + (-t32 * (t2 + t9) * r12dot + 0.2e1 * t80 * t5 * t25) * t24) * t24 + (0.2e1 * t144 + t148) * t24 * t59 + ((0.2e1 * m2 + t110) * t7 * t43 * r12dot + t32 * (-0.6e1 * t33 - 0.1e1) * omega * t24 + (-0.3e1 * t95 * t8 * r12dot + (0.4e1 * t144 + t148) * t24) * y) * y + (-t32 * t34 * t7 * omega + (-t41 * t43 * r12dot - 0.2e1 * t46 * t1 * t25) * t24 + (-t92 - t101 * t24) * y) * x) * x;

  t2 = m2 / M;
  t4 = m2 * m2;
  t5 = M * M;
  t8 = sin(phi);
  t9 = t8 * t8;
  t10 = t4 / t5 * t9;
  t16 = 0.1e1 / b;
  t18 = m2 * t16;
  t20 = cos(phi);
  t21 = t20 * t20;
  t24 = b * b;
  t27 = m2 / t24 / b;
  t28 = z * z;
  t31 = 0.1e1 / t24;
  t35 = -0.1e1 - 0.2e1 * t2 * t9;
  dxIZ_dxNZ[1][1] = (m2 * (0.1e1 + (0.3e1 / 0.2e1 * t2 + t10) * t9) + b) * t16 + t18 * (t2 + t10) * t21 - t27 * t21 * t28 + (m2 * t31 * t35 * t8 - t27 * (-t9 + t21) * y) * y + t18 * (t35 * t16 * t20 + 0.2e1 * t31 * t8 * t20 * y) * x;

  t1 = 0.1e1 / b;
  t2 = m2 * t1;
  t3 = 0.1e1 / M;
  t4 = m2 * t3;
  t6 = m2 * m2;
  t7 = M * M;
  t9 = t6 / t7;
  t10 = sin(phi);
  t11 = t10 * t10;
  t15 = cos(phi);
  t16 = t15 * t15;
  t21 = b * b;
  t22 = 0.1e1 / t21;
  t23 = t22 * t10;
  t24 = z * z;
  t27 = t4 * t11;
  dxIZ_dxNZ[1][2] = t2 * (((-t4 / 0.2e1 - t9 * t11) * t10 - t9 * t10 * t16) * t15 - t23 * t15 * t24 + (((0.1e1 + 0.4e1 * t27) * t1 + 0.2e1 * t2 * t3 * t16) * t15 - 0.3e1 * t23 * t15 * y) * y + ((-0.1e1 - 0.2e1 * t27) * t10 * t1 - t22 * (0.2e1 * t16 - 0.2e1 * t11) * y + t23 * t15 * x) * x);

  t1 = 0.1e1 / b;
  t2 = m2 * t1;
  t3 = 0.1e1 / M;
  t5 = sin(phi);
  t6 = t5 * t5;
  t11 = cos(phi);
  t12 = t11 * t11;
  t19 = b * b;
  t21 = 0.1e1 / t19 * z;
  dxIZ_dxNZ[1][3] = t2 * (((0.2e1 * m2 * t3 * t6 + 0.1e1) * t1 + 0.2e1 * t2 * t3 * t12) * t11 * z - 0.2e1 * t21 * t5 * t11 * y - 0.2e1 * t21 * t12 * x);

  t2 = 0.1e1 / M;
  t3 = sin(phi);
  t4 = t2 * t3;
  t6 = m2 * t2;
  t7 = m2 * m2;
  t8 = M * M;
  t9 = 0.1e1 / t8;
  t11 = t3 * t3;
  t12 = t7 * t9 * t11;
  t20 = t7 * m2;
  t22 = cos(phi);
  t23 = t22 * t22;
  t34 = b * b;
  t36 = 0.1e1 / t34 / b;
  t38 = (-m2 - 0.2e1 * t7 * t2 * t11) * t3 * t36 * r12dot;
  t39 = 0.1e1 / t34;
  t40 = m2 * t39;
  t41 = t6 * t11;
  t45 = t7 * t36;
  t46 = r12dot * t2;
  t48 = t45 * t46 * t3;
  t50 = t7 * t39;
  t51 = t2 * omega;
  t53 = t50 * t51 * t22;
  t54 = -0.2e1 * t48 + t53;
  t59 = z * z;
  t64 = 0.1e1 / b;
  t65 = m2 * t64;
  t75 = t20 * t64;
  t76 = t9 * t3;
  t80 = t20 * t39;
  t90 = m2 * t36;
  t95 = t34 * t34;
  t97 = m2 / t95;
  t100 = 0.3e1 * t97 * t11 * r12dot;
  t103 = t3 * r12dot;
  t118 = t6 / 0.2e1 + t12;
  t141 = t90 * t11 * omega;
  t146 = 0.3e1 * t97 * t103 - t90 * omega * t22;
  t148 = t141 + t146 * t22;
  t150 = t3 * omega;
  dxIZ_dxNZ[2][0] = -m2 * r12dot * t4 + ((m2 * (-t6 - t12 / 0.2e1) - m2 * b * t2) * omega - t20 * t9 * omega * t23 / 0.2e1) * t22 + (t38 + (t40 * (t41 + 0.1e1 / 0.2e1) * omega + t54 * t22) * t22) * t59 + (-t40 * (0.1e1 + t41) * r12dot + (t65 * (-t6 - 0.2e1 * t12) * t3 * omega + (-t40 * (0.3e1 / 0.2e1 * t6 + t12) * r12dot + (-0.2e1 * t75 * t76 * omega - t80 * t9 * r12dot * t22) * t22) * t22) * t22 + (-0.2e1 * t90 * t22 * t3 * omega + t100) * t59 + (t90 * t103 + (t40 * (0.2e1 * t41 - 0.1e1 / 0.2e1) * omega - t54 * t22) * t22) * y) * y + (t65 * t118 * t11 * omega + (t40 * t118 * t3 * r12dot + (-t7 * t64 * t51 / 0.2e1 + (t80 * t76 * r12dot - t75 * t9 * omega * t22) * t22) * t22) * t22 + t148 * t59 + (t40 * t150 + (0.2e1 * t90 * r12dot + (0.6e1 * t50 * t4 * omega + 0.4e1 * t45 * t46 * t22) * t22) * t22 + (-t141 - t146 * t22) * y) * y + (t38 + (t40 * (-t41 + 0.1e1 / 0.2e1) * omega + (-0.4e1 * t48 + 0.2e1 * t53) * t22) * t22 + (t100 + (-0.4e1 * t90 * t150 - 0.3e1 * t97 * r12dot * t22) * t22) * y + t148 * x) * x) * x;

  t1 = 0.1e1 / b;
  t2 = m2 * t1;
  t3 = 0.1e1 / M;
  t4 = m2 * t3;
  t6 = m2 * m2;
  t7 = M * M;
  t9 = t6 / t7;
  t10 = sin(phi);
  t11 = t10 * t10;
  t15 = cos(phi);
  t16 = t15 * t15;
  t21 = b * b;
  t22 = 0.1e1 / t21;
  t23 = t22 * t10;
  t24 = z * z;
  dxIZ_dxNZ[2][1] = t2 * (((-t4 / 0.2e1 - t9 * t11) * t10 - t9 * t10 * t16) * t15 - t23 * t15 * t24 + ((-t1 - 0.2e1 * t2 * t3 * t16) * t15 + t23 * t15 * y) * y + ((0.2e1 * t4 * t11 + 0.1e1) * t10 * t1 + 0.4e1 * t2 * t3 * t10 * t16 - t22 * (0.2e1 * t11 - 0.2e1 * t16) * y - 0.3e1 * t23 * t15 * x) * x);

  t1 = 0.1e1 / M;
  t2 = m2 * t1;
  t3 = sin(phi);
  t4 = t3 * t3;
  t9 = 0.1e1 / b;
  t11 = m2 * t9;
  t13 = m2 * m2;
  t14 = M * M;
  t15 = 0.1e1 / t14;
  t22 = cos(phi);
  t23 = t22 * t22;
  t28 = b * b;
  t31 = m2 / t28 / b;
  t32 = z * z;
  dxIZ_dxNZ[2][2] = (m2 * (0.1e1 + t2 * t4) + b) * t9 + (t11 * (0.3e1 / 0.2e1 * t2 + t13 * t15 * t4) + t13 * m2 * t9 * t15 * t23) * t23 - t31 * t4 * t32 + t11 * (-t3 * t9 - 0.2e1 * t11 * t1 * t3 * t23) * y + (t11 * ((-t9 - 0.2e1 * t11 * t1 * t23) * t22 + 0.2e1 / t28 * t3 * t22 * y) - t31 * (t4 - t23) * x) * x;

  t1 = 0.1e1 / b;
  t2 = m2 * t1;
  t3 = 0.1e1 / M;
  t5 = sin(phi);
  t6 = t5 * t5;
  t13 = cos(phi);
  t14 = t13 * t13;
  t20 = b * b;
  t22 = 0.1e1 / t20 * z;
  dxIZ_dxNZ[2][3] = t2 * (((0.2e1 * m2 * t3 * t6 + 0.1e1) * t5 * t1 + 0.2e1 * t2 * t3 * t5 * t14) * z - 0.2e1 * t22 * t6 * y - 0.2e1 * t22 * t5 * t13 * x);

  t1 = b * b;
  t2 = 0.1e1 / t1;
  t3 = m2 * t2;
  t4 = 0.1e1 / M;
  t5 = m2 * t4;
  t6 = m2 * m2;
  t7 = M * M;
  t8 = 0.1e1 / t7;
  t10 = sin(phi);
  t11 = t10 * t10;
  t12 = t6 * t8 * t11;
  t25 = cos(phi);
  t26 = t25 * t25;
  t37 = 0.2e1 * m2 + 0.4e1 * t6 * t4 * t11;
  t40 = 0.1e1 / t1 / b;
  t45 = -0.1e1 - 0.2e1 * t5 * t11;
  t48 = t6 * t40;
  t49 = r12dot * t4;
  t53 = t6 * t2;
  t64 = t1 * t1;
  t66 = m2 / t64;
  t70 = m2 * t40;
  dxIZ_dxNZ[3][0] = (-t3 * (0.1e1 + (t5 + t12) * t11) * r12dot + (-t3 * (t5 + 0.2e1 * t12) * r12dot - t6 * m2 * t2 * t8 * r12dot * t26) * t26) * z + ((t37 * t10 * t40 * r12dot + (t3 * t45 * omega + (0.4e1 * t48 * t49 * t10 - 0.2e1 * t53 * t4 * omega * t25) * t25) * t25) * z + (-0.3e1 * t66 * t11 * r12dot + 0.2e1 * t70 * t25 * t10 * omega) * z * y) * y + ((-t3 * t45 * t10 * omega + (t37 * t40 * r12dot + (0.2e1 * t53 * t4 * t10 * omega + 0.4e1 * t48 * t49 * t25) * t25) * t25) * z + (-0.2e1 * t70 * t11 * omega + (-0.6e1 * t66 * t10 * r12dot + 0.2e1 * t70 * omega * t25) * t25) * z * y + (-0.2e1 * t70 * t10 * omega - 0.3e1 * t66 * r12dot * t25) * t25 * z * x) * x;

  t1 = 0.1e1 / b;
  t2 = m2 * t1;
  t3 = 0.1e1 / M;
  t5 = sin(phi);
  t6 = t5 * t5;
  t11 = cos(phi);
  t12 = t11 * t11;
  t19 = b * b;
  t21 = 0.1e1 / t19 * z;
  dxIZ_dxNZ[3][1] = t2 * (((-0.1e1 - 0.2e1 * m2 * t3 * t6) * t1 - 0.2e1 * t2 * t3 * t12) * t11 * z + 0.2e1 * t21 * t5 * t11 * y + 0.2e1 * t21 * t12 * x);

  t1 = 0.1e1 / b;
  t2 = m2 * t1;
  t3 = 0.1e1 / M;
  t5 = sin(phi);
  t6 = t5 * t5;
  t13 = cos(phi);
  t14 = t13 * t13;
  t20 = b * b;
  t22 = 0.1e1 / t20 * z;
  dxIZ_dxNZ[3][2] = t2 * (((-0.1e1 - 0.2e1 * m2 * t3 * t6) * t5 * t1 - 0.2e1 * t2 * t3 * t5 * t14) * z + 0.2e1 * t22 * t6 * y + 0.2e1 * t22 * t5 * t13 * x);

  t1 = 0.1e1 / M;
  t2 = m2 * t1;
  t3 = m2 * m2;
  t4 = M * M;
  t5 = 0.1e1 / t4;
  t7 = sin(phi);
  t8 = t7 * t7;
  t9 = t3 * t5 * t8;
  t15 = 0.1e1 / b;
  t17 = m2 * t15;
  t23 = cos(phi);
  t24 = t23 * t23;
  t31 = -0.1e1 - 0.2e1 * t2 * t8;
  t40 = b * b;
  t43 = m2 / t40 / b;
  dxIZ_dxNZ[3][3] = (m2 * (0.1e1 + (t2 + t9) * t8) + b) * t15 + (t17 * (t2 + 0.2e1 * t9) + t3 * m2 * t15 * t5 * t24) * t24 + (t17 * (t31 * t7 * t15 - 0.2e1 * t17 * t1 * t7 * t24) + t43 * t8 * y) * y + (t17 * ((t31 * t15 - 0.2e1 * t17 * t1 * t24) * t23 + 0.2e1 / t40 * t7 * t23 * y) + t43 * t24 * x) * x;

  return;
}

/**********************************************************************/ 
/*********************************************************************************
   jac_xNZ_to_xIZ_r12const(): 
        -- Generated by maple script CS2PN.mpl with mini-disk notes 
        -- Assumes r12dot = 0
        -- Calculates the Jacobian to take a tensor from the x_{NZ} coordinate
            basis to the x_{IZ} coordinate basis.
        -- Here the x_{IZ} are Cook-Scheel Horizon penetrating coordinates

     Arguments:
          xNZ = NZ Cartesian coordinates in terms of Mass;
    dxIZ_dxBL = \frac{ \partial x^{\alpha}_{IZ} }{\partial x^{\beta}_{NZ}} ( x_{NZ} )
 *********************************************************************************/
 void jac_xNZ_to_xIZ_r12const( int BH, double xNZ[NDIM], double dxIZ_dxNZ[NDIM][NDIM])
{
  double m1,m2,phi,b,omega;
  double t,x,y,z;
  double t1,t2,t3,t4,t5,t6,t7,t8,t9;
  double t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
  double t20,t21,t22,t23,t24,t25,t26,t27,t28,t29;
  double t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;
  double t40,t41,t42,t43,t44,t45,t46,t47,t48,t49;
  double t50,t51,t52,t53,t54,t55,t56,t57,t58,t59;
  double t60,t61,t62,t63,t64,t65,t66,t67,t68,t69;
  double t70,t71,t72,t73,t74,t75,t76,t77,t78,t79;
  double t80,t81,t82,t83,t84,t85,t86,t87,t88,t89;
  double t90,t91,t92,t93,t94,t95,t96,t97,t98,t99;
  double t101,t102,t121;
  double t104,t105,t106,t108,t109,t111,t117,t127;

  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  struct of_bbh_traj bbh_traj;
  get_bbh_traj_data(&bbh_traj);

  if( BH == MBH1 ) {
    t = xNZ[TT]; x = xNZ[XX];
    y = xNZ[YY]; z = xNZ[ZZ];
    
    /* Set BH Traj Values */
    m1    = m_bh1;
    m2    = m_bh2;
    phi   = bbh_traj.phi;
    b     = bbh_traj.r12;
    //r12dot= bbh_traj->r12dot;
    omega = bbh_traj.omega;
  }
  else if( BH == MBH2 ) {
    t = xNZ[TT] ; x = -xNZ[XX];
    y = -xNZ[YY]; z = xNZ[ZZ];
    
    /* Set BH Traj Values */
    m1    = m_bh2;
    m2    = m_bh1;
    phi   = bbh_traj.phi;
    b     = bbh_traj.r12;
    //r12dot= bbh_traj->r12dot;
    omega = bbh_traj.omega;
  }

  /* Maple generated code */
  t1 = cos(phi);
  t3 = sin(phi);
  t6 = m2 * (x * t1 + y * t3);
  t7 = m2 * b;
  t8 = 0.1e1 / M;
  t9 = t8 * t1;
  t12 = pow(x - t7 * t9, 0.2e1);
  t13 = t8 * t3;
  t16 = pow(y - t7 * t13, 0.2e1);
  t17 = z * z;
  t19 = sqrt(t12 + t16 + t17);
  t28 = m2 * t8;
  t30 = t3 * t3;
  t31 = t28 * t30;
  t33 = -0.3e1 * t28 + 0.6e1 * t31;
  t35 = t1 * t1;
  t37 = t28 * t35 * omega;
  t42 = 0.1e1 / t19;
  t45 = y * y;
  t55 = t3 * omega;
  t58 = 0.6e1 * t28 * t55 * t35;
  t90 = -t28 * t35 - t31;
  t93 = t13 * omega;
  t96 = t9 * omega;
  t104 = sqrt(m2 / b);
  t105 = sqrt(t28);
  t106 = t104 * t105;
  t108 = m2 * m2;
  t109 = t90 * t90;
  t111 = t108 * t109 * t42;
  t127 = b * b;
  dxIZ_dxNZ[0][0] = (-0.2e1 * t6 * t19 * (-x * t3 * omega + y * t1 * omega) + (m2 * (((t33 * omega + 0.6e1 * t37) * t1 * t19 + 0.3e1 * t30 * t42 * m2 * t9 * omega * t45) * y + ((-t33 * t3 * omega - t58) * t19 + (-0.3e1 * t28 * t30 * t3 * omega + t58) * t42 * t45 + ((-0.6e1 * t28 * t30 * omega + 0.3e1 * t37) * t1 * t42 * y - 0.3e1 * t35 * t42 * m2 * t13 * omega * x) * x) * x) / 0.3e1 + (-t6 * t90 * t42 * (0.2e1 * x * m2 * t93 - 0.2e1 * y * m2 * t96) + (0.1e1 + (t106 * t55 + t111 * t96) * y + (t106 * t1 * omega - t111 * t93) * x) * b) * b) * b) / t127 / b;
  t1 = sin(phi);
  t2 = y * t1;
  t3 = m2 * b;
  t4 = 0.1e1 / M;
  t5 = cos(phi);
  t6 = t4 * t5;
  t9 = pow(x - t3 * t6, 0.2e1);
  t13 = pow(y - t3 * t4 * t1, 0.2e1);
  t14 = z * z;
  t16 = sqrt(t9 + t13 + t14);
  t20 = t5 * t5;
  t24 = y * y;
  t25 = t1 * t1;
  t27 = 0.1e1 / t16;
  t28 = t24 * t25 * t27;
  t43 = m2 * t4;
  t45 = t43 * t25;
  t46 = 0.6e1 * t45;
  t47 = t43 * t20;
  t74 = -t47 - t45;
  t80 = m2 * m2;
  t81 = M * M;
  t83 = t80 / t81;
  t84 = t25 * t25;
  t101 = sqrt(m2 / b);
  t102 = sqrt(t43);
  t105 = t74 * t74;
  t117 = b * b;
  dxIZ_dxNZ[0][1] = (m2 * (-0.6e1 * t2 * t16 * t5 + ((0.3e1 - 0.6e1 * t20) * t16 - 0.3e1 * t28 + (-0.6e1 * t2 * t5 * t27 - 0.3e1 * t20 * t27 * x) * x) * x) / 0.3e1 + (m2 * ((-0.3e1 * t43 + t46 + 0.6e1 * t47) * t5 * t16 + 0.3e1 * t28 * t43 * t5 + ((0.6e1 * t43 * t25 * t1 + 0.12e2 * t1 * t20 * t43) * t27 * y + (t46 + 0.9e1 * t47) * t5 * t27 * x) * x) / 0.3e1 + (m2 * (0.6e1 * t2 * t74 * t27 * m2 * t6 + (-0.3e1 * t83 * t84 + (-0.12e2 * t83 * t25 - 0.9e1 * t83 * t20) * t20) * t27 * x) / 0.3e1 + (t101 * t102 * t1 + t80 * t105 * t27 * t4 * t5) * b) * b) * b) / t117 / b;
  t1 = sin(phi);
  t2 = t1 * t1;
  t5 = m2 * b;
  t6 = 0.1e1 / M;
  t7 = cos(phi);
  t11 = pow(x - t5 * t6 * t7, 0.2e1);
  t15 = pow(y - t5 * t6 * t1, 0.2e1);
  t16 = z * z;
  t18 = sqrt(t11 + t15 + t16);
  t20 = y * y;
  t22 = 0.1e1 / t18;
  t34 = t7 * t7;
  t44 = m2 * t6;
  t46 = t44 * t2;
  t50 = t1 * t34;
  t52 = 0.6e1 * t50 * t44;
  t62 = t44 * t34;
  t77 = m2 * m2;
  t78 = M * M;
  t80 = t77 / t78;
  t81 = t2 * t2;
  t93 = -t62 - t46;
  t105 = sqrt(m2 / b);
  t106 = sqrt(t44);
  t109 = t93 * t93;
  t121 = b * b;
  dxIZ_dxNZ[0][2] = (m2 * (((0.3e1 - 0.6e1 * t2) * t18 - 0.3e1 * t20 * t2 * t22) * y + (-0.6e1 * t7 * t18 * t1 - 0.6e1 * t20 * t1 * t7 * t22 - 0.3e1 * t34 * t22 * y * x) * x) / 0.3e1 + (m2 * (((-0.3e1 * t44 + 0.6e1 * t46) * t1 + t52) * t18 + (0.9e1 * t44 * t2 * t1 + t52) * t22 * t20 + ((0.12e2 * t46 + 0.6e1 * t62) * t7 * t22 * y + 0.3e1 * t50 * t22 * t44 * x) * x) / 0.3e1 + (m2 * ((-0.9e1 * t80 * t81 + (-0.12e2 * t80 * t2 - 0.3e1 * t80 * t34) * t34) * t22 * y + 0.6e1 * t7 * t93 * t22 * t44 * t1 * x) / 0.3e1 + (-t105 * t106 * t7 + t77 * t109 * t22 * t6 * t1) * b) * b) * b) / t121 / b;
  t1 = b * b;
  t5 = m2 * b;
  t6 = 0.1e1 / M;
  t7 = cos(phi);
  t11 = pow(x - t5 * t6 * t7, 0.2e1);
  t12 = sin(phi);
  t16 = pow(y - t5 * t6 * t12, 0.2e1);
  t17 = z * z;
  t19 = sqrt(t11 + t16 + t17);
  t22 = y * y;
  t23 = t12 * t12;
  t25 = 0.1e1 / t19;
  t26 = t25 * z;
  t29 = y * t12;
  t34 = t7 * t7;
  t43 = m2 * t6;
  t46 = -t43 * t34 - t43 * t23;
  t50 = t46 * t46;
  dxIZ_dxNZ[0][3] = m2 / t1 / b * (0.3e1 * t19 * z - 0.3e1 * t22 * t23 * t26 + (-0.6e1 * t29 * t7 * t25 * z - 0.3e1 * t34 * t25 * z * x) * x + (-0.6e1 * (x * t7 + t29) * t46 * t26 - 0.3e1 * t50 * t25 * z * b) * b) / 0.3e1;
  t1 = sin(phi);
  t2 = t1 * t1;
  t4 = cos(phi);
  t5 = t4 * t4;
  t6 = t5 * omega;
  t7 = t2 * omega - t6;
  t8 = z * z;
  t10 = y * y;
  t15 = t1 * omega;
  t29 = 0.1e1 / M;
  t30 = m2 * t29;
  t31 = t30 * t2;
  t32 = -0.1e1 / 0.2e1 - t31;
  t35 = t15 * t5;
  t36 = t30 * t35;
  t61 = m2 * m2;
  t62 = M * M;
  t64 = t61 / t62;
  t65 = t64 * t2;
  t80 = t64 * t35;
  t105 = b * b;
  dxIZ_dxNZ[1][0] = (m2 * ((t7 * t8 + t7 * t10) * y + (0.2e1 * t8 * t4 * t15 + 0.4e1 * t10 * t1 * t4 * omega - t7 * y * x) * x) + (m2 * ((t32 * t1 * omega - t36) * t8 + ((-0.1e1 / 0.2e1 - 0.2e1 * t31) * t1 * omega + t36) * t10 + ((-0.6e1 * t31 - 0.1e1) * omega * t4 * y + (-t32 * t1 * omega - 0.2e1 * t36) * x) * x) + (m2 * (((t30 / 0.2e1 + t65) * t2 * omega + (-t30 * omega / 0.2e1 - t64 * t6) * t5) * y + ((t30 + 0.2e1 * t65) * t1 * omega + 0.2e1 * t80) * t4 * x) + (m2 * ((t30 + t65 / 0.2e1) * t1 * omega + t80 / 0.2e1) + m2 * b * t29 * t1 * omega) * b) * b) * b) / t105 / b;
  t1 = z * z;
  t2 = cos(phi);
  t3 = t2 * t2;
  t5 = sin(phi);
  t6 = t5 * t5;
  t8 = y * y;
  t17 = m2 / M;
  t20 = -0.1e1 - 0.2e1 * t17 * t6;
  t28 = m2 * m2;
  t29 = M * M;
  t32 = t28 / t29 * t6;
  t44 = b * b;
  dxIZ_dxNZ[1][1] = (m2 * (-t1 * t3 + (t6 - t3) * t8 + 0.2e1 * y * t5 * t2 * x) + (m2 * (t20 * t5 * y + t20 * t2 * x) + (m2 * (0.1e1 + (0.3e1 / 0.2e1 * t17 + t32) * t6 + (t32 + t17) * t3) + b) * b) * b) / t44 / b;
  t1 = b * b;
  t5 = z * z;
  t6 = cos(phi);
  t8 = sin(phi);
  t10 = y * y;
  t14 = t6 * t6;
  t15 = t8 * t8;
  t24 = m2 / M;
  t25 = t24 * t15;
  t37 = m2 * m2;
  t38 = M * M;
  t40 = t37 / t38;
  dxIZ_dxNZ[1][2] = m2 / t1 / b * (-t5 * t6 * t8 - 0.3e1 * t10 * t8 * t6 + ((-0.2e1 * t14 + 0.2e1 * t15) * y + t6 * t8 * x) * x + ((0.4e1 * t25 + 0.1e1 + 0.2e1 * t24 * t14) * t6 * y + (-0.1e1 - 0.2e1 * t25) * t8 * x + ((-t24 / 0.2e1 - t40 * t15) * t8 - t40 * t8 * t14) * t6 * b) * b);
  t1 = b * b;
  t5 = cos(phi);
  t7 = sin(phi);
  t14 = m2 / M;
  t15 = t7 * t7;
  t18 = t5 * t5;
  dxIZ_dxNZ[1][3] = m2 / t1 / b * (-0.2e1 * z * (x * t5 + y * t7) * t5 + (0.2e1 * t14 * t15 + 0.1e1 + 0.2e1 * t14 * t18) * t5 * z * b);
  t1 = z * z;
  t2 = cos(phi);
  t4 = sin(phi);
  t5 = t4 * omega;
  t9 = t4 * t4;
  t11 = t2 * t2;
  t12 = t11 * omega;
  t13 = t9 * omega - t12;
  t16 = y * y;
  t29 = 0.1e1 / M;
  t30 = m2 * t29;
  t31 = t30 * t9;
  t34 = t30 * t12;
  t44 = t5 * t11;
  t59 = m2 * m2;
  t60 = M * M;
  t62 = t59 / t60;
  t63 = t62 * t9;
  t79 = t62 * t12;
  t104 = b * b;
  dxIZ_dxNZ[2][0] = (m2 * (-0.2e1 * t1 * t2 * t5 * y + (t13 * t1 - t13 * t16 + (-0.4e1 * y * t4 * t2 * omega + t13 * x) * x) * x) + (m2 * (((0.1e1 / 0.2e1 + t31) * omega + t34) * t2 * t1 + ((0.2e1 * t31 - 0.1e1 / 0.2e1) * omega - t34) * t2 * t16 + ((t5 + 0.6e1 * t30 * t44) * y + ((-t31 + 0.1e1 / 0.2e1) * omega + 0.2e1 * t34) * t2 * x) * x) + (m2 * (((-t30 - 0.2e1 * t63) * t4 * omega - 0.2e1 * t62 * t44) * t2 * y + ((t30 / 0.2e1 + t63) * t9 * omega + (-t30 * omega / 0.2e1 - t79) * t11) * x) + (m2 * ((-t30 - t63 / 0.2e1) * omega - t79 / 0.2e1) * t2 - m2 * b * t29 * t2 * omega) * b) * b) * b) / t104 / b;
  t1 = b * b;
  t5 = z * z;
  t6 = cos(phi);
  t8 = sin(phi);
  t10 = y * y;
  t13 = t6 * t6;
  t14 = t8 * t8;
  t24 = m2 / M;
  t34 = t8 * t13;
  t40 = m2 * m2;
  t41 = M * M;
  t43 = t40 / t41;
  dxIZ_dxNZ[2][1] = m2 / t1 / b * (-t5 * t6 * t8 + t10 * t8 * t6 + ((0.2e1 * t13 - 0.2e1 * t14) * y - 0.3e1 * t6 * t8 * x) * x + ((-0.1e1 - 0.2e1 * t24 * t13) * t6 * y + ((0.2e1 * t24 * t14 + 0.1e1) * t8 + 0.4e1 * t34 * t24) * x + ((-t24 / 0.2e1 - t43 * t14) * t8 - t43 * t34) * t6 * b) * b);
  t1 = z * z;
  t2 = sin(phi);
  t3 = t2 * t2;
  t6 = cos(phi);
  t9 = t6 * t6;
  t18 = m2 / M;
  t32 = m2 * m2;
  t33 = M * M;
  t35 = t32 / t33;
  t47 = b * b;
  dxIZ_dxNZ[2][2] = (m2 * (-t1 * t3 + (0.2e1 * y * t2 * t6 + (t9 - t3) * x) * x) + (m2 * ((-t2 - 0.2e1 * t2 * t9 * t18) * y + (-0.1e1 - 0.2e1 * t18 * t9) * t6 * x) + (m2 * (0.1e1 + t18 * t3 + (0.3e1 / 0.2e1 * t18 + t35 * t3 + t35 * t9) * t9) + b) * b) * b) / t47 / b;
  t1 = b * b;
  t5 = cos(phi);
  t7 = sin(phi);
  t14 = m2 / M;
  t15 = t7 * t7;
  t20 = t5 * t5;
  dxIZ_dxNZ[2][3] = m2 / t1 / b * (-0.2e1 * z * (x * t5 + y * t7) * t7 + ((0.2e1 * t14 * t15 + 0.1e1) * t7 + 0.2e1 * t7 * t20 * t14) * z * b);
  t1 = b * b;
  t5 = cos(phi);
  t7 = sin(phi);
  t19 = m2 / M;
  t20 = t7 * t7;
  t23 = -0.1e1 - 0.2e1 * t19 * t20;
  t25 = t5 * t5;
  dxIZ_dxNZ[3][0] = m2 / t1 / b * (0.2e1 * z * (x * t5 + y * t7) * (-x * t7 * omega + y * t5 * omega) + ((t23 * omega - 0.2e1 * t19 * t25 * omega) * t5 * z * y + (-t23 * t7 * omega + 0.2e1 * t19 * t7 * omega * t25) * z * x) * b);
  t1 = b * b;
  t5 = cos(phi);
  t7 = sin(phi);
  t14 = m2 / M;
  t15 = t7 * t7;
  t18 = t5 * t5;
  dxIZ_dxNZ[3][1] = m2 / t1 / b * (0.2e1 * z * (x * t5 + y * t7) * t5 + (-0.1e1 - 0.2e1 * t14 * t15 - 0.2e1 * t14 * t18) * t5 * z * b);
  t1 = b * b;
  t5 = cos(phi);
  t7 = sin(phi);
  t14 = m2 / M;
  t15 = t7 * t7;
  t20 = t5 * t5;
  dxIZ_dxNZ[3][2] = m2 / t1 / b * (0.2e1 * z * (x * t5 + y * t7) * t7 + ((-0.1e1 - 0.2e1 * t14 * t15) * t7 - 0.2e1 * t7 * t20 * t14) * z * b);
  t1 = cos(phi);
  t3 = sin(phi);
  t6 = pow(x * t1 + y * t3, 0.2e1);
  t9 = m2 / M;
  t10 = t3 * t3;
  t12 = 0.2e1 * t9 * t10;
  t15 = t1 * t1;
  t28 = m2 * m2;
  t29 = M * M;
  t31 = t28 / t29;
  t32 = t31 * t10;
  t46 = b * b;
  dxIZ_dxNZ[3][3] = (t6 * m2 + (m2 * (((-0.1e1 - t12) * t3 - 0.2e1 * t3 * t15 * t9) * y + (-0.1e1 - t12 - 0.2e1 * t9 * t15) * t1 * x) + (m2 * (0.1e1 + (t32 + t9) * t10 + (t9 + 0.2e1 * t32 + t31 * t15) * t15) + b) * b) * b) / t46 / b;

  return;
}


/**********************************************************************/ 
/*********************************************************************************
   coordsBL_of_coordsIZ(): 
         -- Returns the approximate BL coordinates associated with the IZ
            Coordsinates (T,X,Y,Z)
         -- Uses approximate inversions given in  Hergt and Schafer 2008 
            Phys Rev D 77 104001 equations 70-78
         -- Naming convention variables with H are harmonic while BL
            are the BL coordinates variables

         -- Equations used:
             rp = M + sqrt(M^2 - a^2);
             rm = M - sqrt(M^2 - a^2);
             tH = tBL + ( rp^2 + a^2 ) / ( rp - rm ) log( abs( (rBL - rp) / ( rBL - rm ) ) );
            tBL = tH  - ( rp^2 + a^2 ) / ( rp - rm ) log( abs( (rBL - rp) / ( rBL - rm ) ) );
             rH = sqrt( xH*xH + yH*yH + zH*zH );
            rBL = 1/sqrt(2) * sqrt( rH^2 - a^2 + sqrt( (rH^2 - a^2)^2 + 4a^2 zH^2 ) ) + M;
            thBL = arccos( zH / ( rBL - M ) );
            phBL = phH - phA - ( a / ( rp - rm ) ) * log( abs( (rBL - rp)/(rBL - rm) ) );
     Arguments:
         Mass = BH mass ( either m_bh1 or m_bh2 );
         Spin = BH spin ( either a_bh1 or a_bh2 );
          xIZ = IZ coordinates in terms of Mass;
          xBL = BL coordinates in terms of Mass;
 *********************************************************************************/
void coordsBL_of_coordsIZ( double Mass, double Spin, double xIZ[NDIM], double xBL[NDIM])
{
  double rp, rm;
  double tH, xH, yH, zH;
  double rH, phH, phA;
  double tBL, rBL, thBL, phBL;
  double lnarg, lnterm;

  /* Define horizons */
  rp = Mass + sqrt(Mass*Mass - Spin*Spin);
  rm = Mass - sqrt(Mass*Mass - Spin*Spin);

  /* Define short hands for IZ coords */
  tH = xIZ[TT]; xH = xIZ[XX];
  yH = xIZ[YY]; zH = xIZ[ZZ];

  /* Build radial quantities from IZ coords */
  rH     = sqrt(xH*xH + yH*yH + zH*zH);
  rBL    = (1. / sqrt(2.)) * sqrt( rH*rH - Spin*Spin + sqrt( (rH*rH - Spin*Spin)*(rH*rH - Spin*Spin) + 4. * Spin*Spin * zH * zH ) ) + Mass;
  lnarg  = ( rBL - rp ) / ( rBL - rm );
  lnterm = log( fabs( lnarg ) );

  /* Build intermediate quantities */
  phH = atan2( yH, xH );
  phA = atan2( Spin, rBL - Mass );

  /* Calculate remaining coordinates */
  tBL  = tH - ( (rp*rp + Spin*Spin)/(rp - rm) ) * lnterm;
  thBL = acos( zH / ( rBL - Mass ) );
  phBL = phH - phA - ( Spin / ( rp - rm ) ) * lnterm;

  /* Return approximate BL coordinates */
  xBL[TT] = tBL;
  xBL[RR] = rBL;
  xBL[TH] = thBL;
  xBL[PH] = phBL;
  
  return;
}

/**********************************************************************/ 
/*********************************************************************************
   coordsIZ_of_coordsNZ(): --optimized slightly
         -- Returns the coordinates of the IZ ( Cook-Scheel )
            (T,X,Y,Z) given NZ PN harmonic coordinates (t,x,y,z)
         -- Taken from Brennan's maple
            init_data/spin/maple/All_in_one_IZ_spin_LTE.mpl
         Arguments:
          xNZ = NZ coordinates in cartesian coordinates. ( xcart[NDIM] )
          xIZ = IZ coordinates returned ;
          bbh_traj = structure containing bbh trajectory information
 *********************************************************************************/
 void coordsIZ_of_coordsNZ(int BH, double xNZ[NDIM], double xIZ[NDIM])
{
  double t1 ,t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9;
  double t10,t11,t12,t13,t14,t15,t16,t17,t18;
  double m1,m2,b,phi;
  double t,x,y,z;

  extern void get_bbh_traj_data(struct of_bbh_traj *bbh_traj) ;
  struct of_bbh_traj bbh_traj;
  get_bbh_traj_data(&bbh_traj);

  if( BH == MBH1 ) {
    /* Set BBH information */
    m1  = m_bh1;
    m2  = m_bh2;
    b   = bbh_traj.r12;
    phi = bbh_traj.phi;
    
    /* Rename NZ coordinates */
    t = xNZ[TT]; x = xNZ[XX];
    y = xNZ[YY]; z = xNZ[ZZ];
  }
  else if( BH == MBH2 ) {
    /* Set BBH information */
    m1  = m_bh2;
    m2  = m_bh1;
    b   = bbh_traj.r12;
    phi = bbh_traj.phi;
    
    /* Rename NZ coordinates */
    t = xNZ[TT]; x = -xNZ[XX];
    y = -xNZ[YY]; z = xNZ[ZZ];
  }
  
  t1 = 1./b;//0.1e1 / b;
  t2 = m2 * t1;
  t3 = m1 + m2;
  t3 = 1./t3;//0.1e1 / t3;
  t4 = m2 * t3;
  t5 = sin(phi);
  t6 = t4 * b;
  t7 = y - t6 * t5;
  t8 = cos(phi);
  t6 = x - t6 * t8;
  t9 = t7 * t8 - t6 * t5;
  t10 = t6*t6 + t7*t7 + z * z;//pow(t6, 0.2e1) + pow(t7, 0.2e1) + z * z;
  t11 = sqrt(t10);
  t12 = t6 * t8 + t7 * t5;
  t13 = t12*t12;//pow(t12, 0.2e1);
  t14 = t1*t1;//pow(t1, 0.2e1);
  t15 = t10 * t8;
  t16 = 1. - t12 * t1;//0.1e1 - t12 * t1;
  t17 = t4 * t9;
  t3 = b * t3;
  t18 = t10 * t5;
  xIZ[0] = t - sqrt(t2*t4) * t9 + m2 * t1 * t14 * (t10 * t11 - 3. * t13 * t11) / 3.;//t - sqrt(t2) * sqrt(t4) * t9 + m2 * t1 * t14 * (t10 * t11 - 0.3e1 * t13 * t11) / 0.3e1;
  xIZ[1] = (-t3 * t8 + t1 * (t15 * t1 * 0.5 - t17 * t5 * 0.5 - t14 * t12 * (t15 - t12 * t6) + t6 * t16)) * m2 + x;//(-t3 * t8 + t1 * (t15 * t1 / 0.2e1 - t17 * t5 / 0.2e1 - t14 * t12 * (t15 - t12 * t6) + t6 * t16)) * m2 + x;
  xIZ[2] = (-t3 * t5 + t1 * (t18 * t1 * 0.5 + t17 * t8 * 0.5 - t14 * t12 * (t18 - t12 * t7) + t7 * t16)) * m2 + y;//(-t3 * t5 + t1 * (t18 * t1 / 0.2e1 + t17 * t8 / 0.2e1 - t14 * t12 * (t18 - t12 * t7) + t7 * t16)) * m2 + y;
  xIZ[3] = z + t2 * z * (t14 * t13 + t16);
  return;
}


 /************************************************************************************
  Tranform rank1 contravariant vector from BL coordinates centered on a BH
  to PN Harmonic Coordinates (numeric)
  Arguments:
        BH: either 1 (BH1) or 2 (BH2)
    coords: local PN harmonic coordinates structure
      ucon: 4-vector to transform ( say 4-velocity )
 *************************************************************************************/
extern void transform_rank1con_BL2PN( int BH, struct of_coord *coords, double ucon_bl[NDIM] )
{
  double Mass, Spin;
  
  double xNZ[NDIM], xIZ[NDIM], xBL[NDIM];
  double dxIZ_dxBL[NDIM][NDIM], dxIZ_dxNZ[NDIM][NDIM], dxNZ_dxIZ[NDIM][NDIM];
  double dxs_dxc[NDIM][NDIM];
  double *det_tmp;
  ALLOC_ARRAY(det_tmp,1);


  /* Set various quantities based off which BH */
  if( BH == 1 ) {
    MASS_TYPE = MBH1;
    Mass = m_bh1;
    Spin = a_bh1;
  }
  else if( BH == 2 ) {
    MASS_TYPE = MBH2;
    Mass = m_bh2;
    Spin = a_bh2;
  }
  else{ 
    fprintf(stderr,"coord.c:transform_rank1con_BL2PN() shouldn't be here (BH == %d)\n",BH);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }
  /* Calculate local Cook-Scheel and BL coordinates */
  xNZ[TT] = coords->xcart[TT]; xNZ[XX] = coords->xcart[XX];
  xNZ[YY] = coords->xcart[YY]; xNZ[ZZ] = coords->xcart[ZZ];
  coordsIZ_of_coordsNZ(BH, xNZ, xIZ);
  coordsBL_of_coordsIZ(Mass, Spin, xIZ, xBL);

  /* Get Jacobians */
  jac_xBL_to_xIZ(Mass,Spin,xBL,dxIZ_dxBL);
#if( SET_TSHRINK )
  if( t < t_shrink_bbh ) { jac_xNZ_to_xIZ_r12const(MASS_TYPE, xNZ, dxIZ_dxNZ); } // Get Jacobian 
  else                   { jac_xNZ_to_xIZ(MASS_TYPE,xNZ,dxIZ_dxNZ)           ; } // Get Jacobian    
#else
  jac_xNZ_to_xIZ(MASS_TYPE,xNZ,dxIZ_dxNZ);
#endif
  invert_matrix2( dxIZ_dxNZ, dxNZ_dxIZ, det_tmp );
  xspher_of_xcart(coords->x, coords->xcart, dxs_dxc);

  /* Transform the 4-vector */
  transform_rank1con2(dxIZ_dxBL, ucon_bl); //BL --> CS
  transform_rank1con2(dxNZ_dxIZ, ucon_bl); //CS --> PN
  if( MASS_TYPE == MBH2 ) { ucon_bl[XX] *= -1.; ucon_bl[YY] *= -1.; } //correct Jacobian error for BH2
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
  transform_rank1con2(dxs_dxc, ucon_bl); //PN Cartesian --> PN Spherical
#endif

  transform_rank1con2(coords->dxp_dx, ucon_bl); //PN Physical --> PN Numeric

  MASS_TYPE = MTOT;
  return;
}
 /************************************************************************************
  Tranform rank1 covariant dualvector from BL coordinates centered on a BH
  to PN Harmonic Coordinates ( numeric )
  Arguments:
        BH: either 1 (BH1) or 2 (BH2)
    coords: local PN harmonic coordinates structure
      Acov: dual vector to transform ( say vector potential )
 *************************************************************************************/
extern void transform_rank1cov_BL2PN( int BH, struct of_coord *coords, double acov[NDIM] )
{
  double Mass, Spin;
  
  double xNZ[NDIM], xIZ[NDIM], xBL[NDIM];
  double dxIZ_dxBL[NDIM][NDIM], dxBL_dxIZ[NDIM][NDIM], dxIZ_dxNZ[NDIM][NDIM];
  double dxs_dxc[NDIM][NDIM], dxc_dxs[NDIM][NDIM];
  double *det_tmp;
  ALLOC_ARRAY(det_tmp,1);


  /* Set various quantities based off which BH */
  if( BH == 1 ) {
    MASS_TYPE = MBH1;
    Mass = m_bh1;
    Spin = a_bh1;
  }
  else if( BH == 2 ) {
    MASS_TYPE = MBH2;
    Mass = m_bh2;
    Spin = a_bh2;
  }
  else{ 
    fprintf(stderr,"coord.c:transform_rank1cov_BL2PN() shouldn't be here (BH == %d)\n",BH);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }
  /* Calculate local Cook-Scheel and BL coordinates */
  xNZ[TT] = coords->xcart[TT]; xNZ[XX] = coords->xcart[XX];
  xNZ[YY] = coords->xcart[YY]; xNZ[ZZ] = coords->xcart[ZZ];
  coordsIZ_of_coordsNZ(BH, xNZ, xIZ);
  coordsBL_of_coordsIZ(Mass, Spin, xIZ, xBL);

  /* Get Jacobians */
  jac_xBL_to_xIZ(Mass,Spin,xBL,dxIZ_dxBL);
  invert_matrix2(dxIZ_dxBL, dxBL_dxIZ, det_tmp);
#if( SET_TSHRINK )
  if( t < t_shrink_bbh ) { jac_xNZ_to_xIZ_r12const(MASS_TYPE, xNZ, dxIZ_dxNZ); } // Get Jacobian 
  else                   { jac_xNZ_to_xIZ(MASS_TYPE,xNZ,dxIZ_dxNZ)           ; } // Get Jacobian    
#else
  jac_xNZ_to_xIZ(MASS_TYPE,xNZ,dxIZ_dxNZ);
#endif
  xspher_of_xcart(coords->x, coords->xcart, dxs_dxc);
  invert_matrix2(dxs_dxc, dxc_dxs, det_tmp);

  /* Transform the covector */
  transform_rank1cov2(dxBL_dxIZ, acov); //BL --> CS
  transform_rank1cov2(dxIZ_dxNZ, acov); //CS --> PN
  if( MASS_TYPE == MBH2 ) { acov[XX] *= -1.; acov[YY] *= -1.; } //correct Jacobian error for BH2
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
  transform_rank1cov2(dxc_dxs, acov); //PN Cartesian --> PN Spherical
#endif

  transform_rank1cov2(coords->dx_dxp, acov); //PN Physical ---> PN Numeric

  MASS_TYPE = MTOT;  
  return;
}



#undef STEP
#undef DSTEP
#undef PERIODIC_SHIFT
#undef PERIODIC
