
#include "decs.h"
#include "recon.h"

/**************************************************************************************
  Directional arrays to handle logic of how to select array elements when reconstructing 
   in a give direction;
**************************************************************************************/
const int idir[] = {1,0,0};
const int jdir[] = {0,1,0};
const int kdir[] = {0,0,1};


static void rescale_prims_r( double ****p_h, int rescale_type ) ;
static void recon_regularize_prim( double ****p_h, int rescale_type ) ;

/**************************************************************************************/
/**************************************************************************************
  linear_recon_local():
  --------------

   -- uses various limiters to perform linear interpolation at a cell face;
   -- "local" version allows recon method to be selected dynamically;

**************************************************************************************/
static double linear_recon_local(int irecon, double y1,double y2,double y3) 
{
	double Dqm,Dqp,Dqc,aDqm,aDqp,aDqc,s ;

	switch( irecon ) { 

	/* woodward, or monotonized central, slope limiter */
	case  RECON_MC  :
	  Dqm = 2.*(y2 - y1) ;
	  Dqp = 2.*(y3 - y2) ;
	  Dqc = 0.5*(y3 - y1) ;
	  s = Dqm*Dqp ;
	  if(s <= 0.) return 0. ;
	  else {
	    aDqm = fabs(Dqm) ;
	    aDqp = fabs(Dqp) ;
	    aDqc = fabs(Dqc) ;
	    if( (aDqm <= aDqp) && (aDqm <= aDqc) ) { return( Dqm ); }
	    else if( aDqp <= aDqc )                { return( Dqp ); }
	    else                                   { return( Dqc ); }
	  }
	  break;

	  /* minmod slope limiter (crude but robust) */
	case  RECON_MINMOD  :
	  Dqm = (y2 - y1) ;
	  Dqp = (y3 - y2) ;
	  s = Dqm*Dqp ;
	  if(s <= 0.) return( 0. ) ;
	  else if(fabs(Dqm) < fabs(Dqp)) return( Dqm );
	  else return( Dqp ) ;
	  break;

	  /* van leer slope limiter */
	case  RECON_VANL  :
	  Dqm = (y2 - y1) ;
	  Dqp = (y3 - y2) ;
	  s = Dqm*Dqp ;
	  if(s <= 0.) return( 0. );
	  else return(2.*s/(Dqm+Dqp)) ;
	  break;

	default :
	  fprintf(stderr,"linear_recon(): Unknown limiter type : %d \n", irecon);
	  fflush(stderr); fail(FAIL_BASIC,0);
	  break;

	}

	return( 0. ) ;
}


/**************************************************************************************/
/**************************************************************************************
  linear_recon():
  --------------

   -- uses various limiters to perform linear interpolation at a cell face;
   -- static version of linear_recon()

**************************************************************************************/
double linear_recon(double y1,double y2,double y3) 
{
	double Dqm,Dqp,Dqc,s,aDqm,aDqp,aDqc ;


#if( RECON_TYPE_CHOICE == RECON_MC || RECON_TYPE_CHOICE == RECON_PPM_MC )
	/* woodward, or monotonized central, slope limiter */
	  Dqm = 2.*(y2 - y1) ;
	  Dqp = 2.*(y3 - y2) ;
	  Dqc = 0.5*(y3 - y1) ;
	  s = Dqm*Dqp ;
	  if(s <= 0.) return 0. ;
	  else {
	    aDqm = fabs(Dqm) ;
	    aDqp = fabs(Dqp) ;
	    aDqc = fabs(Dqc) ;
	    if( (aDqm <= aDqp) && (aDqm <= aDqc) ) { return( Dqm ); }
	    else if( aDqp <= aDqc )                { return( Dqp ); }
	    else                                   { return( Dqc ); }
	  }

#elif( RECON_TYPE_CHOICE == RECON_MINMOD || RECON_TYPE_CHOICE == RECON_PPM_MM )
	  /* minmod slope limiter (crude but robust) */
	  Dqm = (y2 - y1) ;
	  Dqp = (y3 - y2) ;
	  s = Dqm*Dqp ;
	  if(s <= 0.) return( 0. ) ;
	  else if(fabs(Dqm) < fabs(Dqp)) return( Dqm );
	  else return( Dqp ) ;

#elif( RECON_TYPE_CHOICE == RECON_VANL || RECON_TYPE_CHOICE == RECON_PPM_VL )
	  /* van leer slope limiter */
	  Dqm = (y2 - y1) ;
	  Dqp = (y3 - y2) ;
	  s = Dqm*Dqp ;
	  if(s <= 0.) return( 0. );
	  else return(2.*s/(Dqm+Dqp)) ;

#else
	  fprintf(stderr,"linear_recon(): Unknown limiter type : %d \n", RECON_TYPE_CHOICE);
	  fflush(stderr); fail(FAIL_BASIC,0);
#endif

	return( 0. ) ;
}


/**************************************************************************************/
/**************************************************************************************
  linear_recon_vect_[1-3]():
  --------------
   -- vector version of linear_recon();
   -- uses various limiters to perform linear interpolation at a cell face;
   -- static version of linear_recon()
   -- the differences between the _[1-3] versions is just the length of the loop;
      -- it is most efficient to do loops of constant bounds;
      -- horrendous repetition of code, but it's faster;

**************************************************************************************/
static void linear_recon_vect_1( double *y, double *Dq )
{
  register unsigned i;
  register double a, b, a_abs, b_abs, s, Dqi;
  register int i_end = N1TOT-1;

  /* Loop over the interior points, making sure to recycle values from last step */
  b = y[1] - y[0]; 
  for( i = 1; i < i_end; i++ ) { 
    a = b;   b = y[i+1] - y[i  ];  s = a*b;
    if( s <= 0. ) { Dqi =  0. ; }
    else {
#if( RECON_TYPE_CHOICE == RECON_MC || RECON_TYPE_CHOICE == RECON_PPM_MC )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if(      (3*a_abs) <= b_abs )  { Dqi =      a + a  ; }
      else if( (3*b_abs) <= a_abs )  { Dqi =      b + b  ; }
      else                           { Dqi = 0.5*(a + b) ; }
#elif( RECON_TYPE_CHOICE == RECON_MINMOD || RECON_TYPE_CHOICE == RECON_PPM_MM )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if( a_abs <= b_abs)  { Dqi = a ; }
      else                 { Dqi = b ; }
#elif( RECON_TYPE_CHOICE == RECON_VANL || RECON_TYPE_CHOICE == RECON_PPM_VL )
      Dqi = 2*s/(a+b);
#endif 
    }
    Dq[i] = Dqi; 
  }
  return ;
}

static void linear_recon_vect_2( double *y, double *Dq )
{
  register unsigned i;
  register double a, b, a_abs, b_abs, s, Dqi;
  register int i_end = N2TOT-1;

  /* Loop over the interior points, making sure to recycle values from last step */
  b = y[1] - y[0]; 
  for( i = 1; i < i_end; i++ ) { 
    a = b;   b = y[i+1] - y[i  ];  s = a*b;
    if( s <= 0. ) { Dqi =  0. ; }
    else {
#if( RECON_TYPE_CHOICE == RECON_MC || RECON_TYPE_CHOICE == RECON_PPM_MC )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if(      (3*a_abs) <= b_abs )  { Dqi =      a + a  ; }
      else if( (3*b_abs) <= a_abs )  { Dqi =      b + b  ; }
      else                           { Dqi = 0.5*(a + b) ; }
#elif( RECON_TYPE_CHOICE == RECON_MINMOD || RECON_TYPE_CHOICE == RECON_PPM_MM )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if( a_abs <= b_abs)  { Dqi = a ; }
      else                 { Dqi = b ; }
#elif( RECON_TYPE_CHOICE == RECON_VANL || RECON_TYPE_CHOICE == RECON_PPM_VL )
      Dqi = 2*s/(a+b);
#endif 
    }
    Dq[i] = Dqi; 
  }
  return ;
}

static void linear_recon_vect_3( double *y, double *Dq )
{
  register unsigned i;
  register double a, b, a_abs, b_abs, s, Dqi;
  register int i_end = N3TOT-1;

  /* Loop over the interior points, making sure to recycle values from last step */
  b = y[1] - y[0]; 
  for( i = 1; i < i_end; i++ ) { 
    a = b;   b = y[i+1] - y[i  ];  s = a*b;
    if( s <= 0. ) { Dqi =  0. ; }
    else {
#if( RECON_TYPE_CHOICE == RECON_MC || RECON_TYPE_CHOICE == RECON_PPM_MC )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if(      (3*a_abs) <= b_abs )  { Dqi =      a + a  ; }
      else if( (3*b_abs) <= a_abs )  { Dqi =      b + b  ; }
      else                           { Dqi = 0.5*(a + b) ; }
#elif( RECON_TYPE_CHOICE == RECON_MINMOD || RECON_TYPE_CHOICE == RECON_PPM_MM )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if( a_abs <= b_abs)  { Dqi = a ; }
      else                 { Dqi = b ; }
#elif( RECON_TYPE_CHOICE == RECON_VANL || RECON_TYPE_CHOICE == RECON_PPM_VL )
      Dqi = 2*s/(a+b);
#endif 
    }
    Dq[i] = Dqi; 
  }
  return ;
}

void linear_recon_vect( double *y, double *Dq, int ibeg, int iend)
{
  unsigned i;
  double a, b, a_abs, b_abs, s, Dqi,a3,b3;

  /* Loop over the interior points, making sure to recycle values from last step */
  b = y[1] - y[0]; 
  for( i = ibeg; i <= iend; i++ ) { 
    a = b;   b = y[i+1] - y[i  ];  s = a*b;
    if( s <= 0. ) { Dqi =  0. ; }
    else {
#if( RECON_TYPE_CHOICE == RECON_MC || RECON_TYPE_CHOICE == RECON_PPM_MC )
      a_abs = fabs(a);   b_abs = fabs(b); a3 = a_abs+a_abs+a_abs; b3 = b_abs+b_abs+b_abs;
      if(      (a3) <= b_abs )  { Dqi =      a + a  ; }
      else if( (b3) <= a_abs )  { Dqi =      b + b  ; }
      else                           { Dqi = 0.5*(a + b) ; }
#elif( RECON_TYPE_CHOICE == RECON_MINMOD || RECON_TYPE_CHOICE == RECON_PPM_MM )
      a_abs = fabs(a);   b_abs = fabs(b); 
      if( a_abs <= b_abs)  { Dqi = a ; }
      else                 { Dqi = b ; }
#elif( RECON_TYPE_CHOICE == RECON_VANL || RECON_TYPE_CHOICE == RECON_PPM_VL )
      Dqi = 2*s/(a+b);
#endif 
    }
    Dq[i] = Dqi; 
  }
  return ;
}


/**************************************************************************************/
/**************************************************************************************
  para_local():
  --------------

   -- implementation of Colella & Woodward's parabolic interpolation scheme;
   -- "local" version allow method to be chosen dynamically;

 * parabolic interpolation subroutin  
 * ref. Collella && Woodward's PPM paper
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *   
 * Written by Xiaoyue Guan 2003-5?  
 *    with a few trivial modifications by Scott Noble

**************************************************************************************/
static void para_local(int irecon, 
		 double x1, double x2, double x3, double x4, double x5, 
		 double *lout, double *rout)
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

         y[0]=x1;
         y[1]=x2;
         y[2]=x3;
         y[3]=x4;
         y[4]=x5;

         /*CW1.7 */
	 switch( irecon ) { 

	   /* PPM + MC */
	 case  RECON_PPM_MC  :
	   /* PPM + MC */
	   for(i=1 ; i<4 ; i++) {
	     Dqm = 2. *(y[i]-y[i-1]);
	     Dqp = 2. *(y[i+1]-y[i]);
	     //               Dqm = 2.5 *(y[i]-y[i-1]);
	     //               Dqp = 2.5 *(y[i+1]-y[i]);

	     s = Dqm*Dqp;
	     if( s <= 0. ) {  dq[i]=0.; }      //CW1.8  
	     else { 
	       Dqc = 0.5 *(y[i+1]-y[i-1]);
	       aDqm = fabs(Dqm) ;
	       aDqp = fabs(Dqp) ;
	       aDqc = fabs(Dqc) ;
	       if( (aDqm <= aDqp) && (aDqm <= aDqc) ) { dq[i] = Dqm; }
	       else if( aDqp <= aDqc )                { dq[i] = Dqp; }
	       else                                   { dq[i] = Dqc; }
	     }
	   }
	   break; 

	   /* PPM + Minmod */
	 case  RECON_PPM_MM  :
	   for(i=1 ; i<4 ; i++) {
	     Dqm = y[i]   - y[i-1];
	     Dqp = y[i+1] - y[i];
	     s = Dqm*Dqp;
	     if (s <=0.) dq[i]=0.;   
	     else if( fabs(Dqm) < fabs(Dqp) )   dq[i]= Dqm;
	     else dq[i] = Dqp;
	   }
	   break;

	   /* PPM + van leer */
	 case  RECON_PPM_VL  :
	   for(i=1 ; i<4 ; i++) {
	     Dqm = y[i]   - y[i-1];
	     Dqp = y[i+1] - y[i];
	     s = Dqm*Dqp;
	     if (s <=0.) dq[i]=0.;   
	     else dq[i] = 2.*s/(Dqm+Dqp);
	   }
	   break;

	 default : 
	   fprintf(stderr,"para(): Unknown limiter type : %d \n", irecon);
	   fflush(stderr);  fail(FAIL_BASIC,0);
	   break;
	 }	   

	 /* CW1.6 */
         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0; /* Finite volume form */
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0; /* Finite volume form */

//         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/8.0; /* Finite difference form */
//         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/8.0; /* Finite difference form */

         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));

         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

	 if( qd*(qd-qe) < 0.0 ) l=3.0*y[2]-2.0*r;

	 if( qd*(qd+qe) < 0.0 ) r=3.0*y[2]-2.0*l;

	 
         *lout=l;   //a_L,j
	 *rout=r;
         //*dw=r-l;                      //CW1.5
         //*w6=6.0*(y[2]-0.5*(l+r));
}


/**************************************************************************************/
/**************************************************************************************
  para():
  --------------

   -- implementation of Colella & Woodward's parabolic interpolation scheme;
   -- this is the global version, recon method is static;

 * parabolic interpolation subroutin  
 * ref. Collella && Woodward's PPM paper
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *   
 * Written by Xiaoyue Guan 2003-5?  
 *    with a few trivial modifications by Scott Noble

**************************************************************************************/
void para(double x1, double x2, double x3, double x4, double x5, 
		 double *lout, double *rout)
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

         y[0]=x1;
         y[1]=x2;
         y[2]=x3;
         y[3]=x4;
         y[4]=x5;

         /*CW1.7 */

#if( RECON_TYPE_CHOICE == RECON_PPM_MC )
	   /* PPM + MC */
	   for(i=1 ; i<4 ; i++) {
	     Dqm = 2. *(y[i  ] - y[i-1]);
	     Dqp = 2. *(y[i+1] - y[i  ]);
	     //               Dqm = 2.5 *(y[i]-y[i-1]);
	     //               Dqp = 2.5 *(y[i+1]-y[i]);

	     s = Dqm*Dqp;
	     if( s <= 0. ) {  dq[i]=0.; }      //CW1.8  
	     else { 
	       Dqc = 0.5 *(y[i+1]-y[i-1]);
	       aDqm = fabs(Dqm) ; aDqp = fabs(Dqp) ; aDqc = fabs(Dqc) ;
	       if( (aDqm <= aDqp) && (aDqm <= aDqc) ) { dq[i] = Dqm; }
	       else if( aDqp <= aDqc )                { dq[i] = Dqp; }
	       else                                   { dq[i] = Dqc; }
		   //	       dq[i]=MIN(aDqc,MIN(aDqm,aDqp))*SIGN(Dqc);
	     }
	   }

#elif( RECON_TYPE_CHOICE == RECON_PPM_MM )
	   /* PPM + Minmod */
	   for(i=1 ; i<4 ; i++) {
	     Dqm = y[i]   - y[i-1];
	     Dqp = y[i+1] - y[i];
	     s = Dqm*Dqp;
	     if( s <= 0. ) dq[i]=0.;   
	     else if( fabs(Dqm) < fabs(Dqp) )   dq[i]= Dqm;
	     else dq[i] = Dqp;
	   }

#elif( RECON_TYPE_CHOICE == RECON_PPM_VL )
	   /* PPM + van leer */
	   for(i=1 ; i<4 ; i++) {
	     Dqm = y[i]   - y[i-1];
	     Dqp = y[i+1] - y[i];
	     s = Dqm*Dqp;
	     if (s <=0.) dq[i]=0.;   
	     else dq[i] = 2.*s/(Dqm+Dqp);
	   }

#else 
	   fprintf(stderr,"para(): Should not be here!! \n");
	   fflush(stderr); fail(FAIL_BASIC,0); 
	   for(i=1 ; i<4 ; i++) { dq[i] = 0.;}
#endif


	 /* CW1.6 */
         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0; /* Finite volume form */
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0; /* Finite volume form */

//         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/8.0; /* Finite difference form */
//         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/8.0; /* Finite difference form */

         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));

         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

	 if( qd*(qd-qe) < 0.0 ) l=3.0*y[2]-2.0*r;

	 if( qd*(qd+qe) < 0.0 ) r=3.0*y[2]-2.0*l;

	 
         *lout=l;   //a_L,j
	 *rout=r;
         //*dw=r-l;                      //CW1.5
         //*w6=6.0*(y[2]-0.5*(l+r));
}

/**************************************************************************************/
/**************************************************************************************
  para_vect_[1-3]():
  --------------
   -- vector version of para(); 
   -- assumes dq (aka Dq in linear_recon_vect())  has been calculated 
       with linear_recon_vect(); 

   -- implementation of Colella & Woodward's parabolic interpolation scheme;
   -- this is the global version, recon method is static;
   -- the _[1-3] versions are for the different directions and use explicit
      array lengths in order to speed up the calculation (non-constant bounds in 
      a loop runs more slowly);

 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *   
 * Originally written by Xiaoyue Guan 2003-5?  
 *    with major changes by Scott Noble

**************************************************************************************/
static void para_vect_1(double *y, double *dq, double *yL, double *yR )
{
  register int i;
  register double l,r,qa, lnew, rnew,yi, r2,l2, y3;
  register double yip;
  register int i_end = N1TOT-2;

  r=0.5*(y[2] + y[1]) - PPM_FACTOR*(dq[2] - dq[1]); 
  yip = y[2]; 

  for( i = 2; i < i_end; i++ ) { 
    /* shift values */
    yi  = yip; yip = y[i+1];  l = r;
    r = 0.5*(yip + yi) - PPM_FACTOR*(dq[i+1] - dq[i]); /* CW1.6 */
    qa=(r-yi)*(yi-l);
    if( qa <=0. ) {  lnew = rnew = yi; }
    else { 
      r2 = r+r; l2 = l+l; y3 = 3*yi;
      lnew = l; rnew = r; 
      if( r > l ) { 
	if( y3 > (l +r2) ) { lnew = y3 - r2; }
	if( y3 < (l2+r ) ) { rnew = y3 - l2; }
      }
      else { 
	if( y3 < (l +r2) ) { lnew = y3 - r2; }
	if( y3 > (l2+r ) ) { rnew = y3 - l2; }
      }
    }
    yL[i] = lnew;   yR[i] = rnew;
  }
  return;
}

static void para_vect_2(double *y, double *dq, double *yL, double *yR )
{
  register int i;
  register double l,r,qa,lnew, rnew,yi,r2,l2, y3;
  register double yip;
  register int i_end = N2TOT-2;

  r=0.5*(y[2] + y[1]) - PPM_FACTOR*(dq[2] - dq[1]); 
  yip = y[2]; 

  for( i = 2; i < i_end; i++ ) { 
    /* shift values */
    yi  = yip;  yip = y[i+1];  l = r;
    r = 0.5*(yip + yi) - PPM_FACTOR*(dq[i+1] - dq[i]); /* CW1.6 */
    qa=(r-yi)*(yi-l);
    if( qa <=0. ) {  lnew = rnew = yi; }
    else { 
      r2 = r+r; l2 = l+l; y3 = 3*yi;
      lnew = l; rnew = r; 
      if( r > l ) { 
	if( y3 > (l +r2) ) { lnew = y3 - r2; }
	if( y3 < (l2+r ) ) { rnew = y3 - l2; }
      }
      else { 
	if( y3 < (l +r2) ) { lnew = y3 - r2; }
	if( y3 > (l2+r ) ) { rnew = y3 - l2; }
      }
    }
    yL[i] = lnew;   yR[i] = rnew;
  }
  return;
}

static void para_vect_3(double *y, double *dq, double *yL, double *yR )
{
  register int i;
  register double l,r,qa, lnew, rnew,yi, r2,l2, y3;
  register double yip;
  register int i_end = N3TOT-2;

  r=0.5*(y[2] + y[1]) - PPM_FACTOR*(dq[2] - dq[1]); 
  yip = y[2]; 

  for( i = 2; i < i_end; i++ ) { 
    /* shift values */
    yi  = yip;  yip = y[i+1];  l = r;
    r = 0.5*(yip + yi) - PPM_FACTOR*(dq[i+1] - dq[i]); /* CW1.6 */
    qa=(r-yi)*(yi-l);
    if( qa <=0. ) {  lnew = rnew = yi; }
    else { 
      r2 = r+r; l2 = l+l; y3 = 3*yi;
      lnew = l; rnew = r; 
      if( r > l ) { 
	if( y3 > (l +r2) ) { lnew = y3 - r2; }
	if( y3 < (l2+r ) ) { rnew = y3 - l2; }
      }
      else { 
	if( y3 < (l +r2) ) { lnew = y3 - r2; }
	if( y3 > (l2+r ) ) { rnew = y3 - l2; }
      }
    }
    yL[i] = lnew;   yR[i] = rnew;
  }
  return;
}

/* arbitrary bounds version (about just as fast---to within 2%) : */
void para_vect(double *y, double *dq, double *yL, double *yR , int ibeg, int iend)
{
  int i;
  double l,r,qa, lnew, rnew,yi, r2,l2, y3;
  double yip;

  r=0.5*(y[2] + y[1]) - PPM_FACTOR*(dq[2] - dq[1]); 
  yip = y[2]; 

  for( i = ibeg; i <= iend; i++ ) { 
    /* shift values */
    yi  = yip; yip = y[i+1];  l = r;
    r = 0.5*(yip + yi) - PPM_FACTOR*(dq[i+1] - dq[i]); /* CW1.6 */
    qa=(r-yi)*(yi-l);
    if( qa <=0. ) {  lnew = rnew = yi; }
    else { 
      r2 = r+r; l2 = l+l; y3 = yi+yi+yi;
      lnew = l; rnew = r; 
      if( r > l ) { 
	if( y3 > (l +r2) ) { lnew = y3 - r2; }
	if( y3 < (l2+r ) ) { rnew = y3 - l2; }
      }
      else { 
	if( y3 < (l +r2) ) { lnew = y3 - r2; }
	if( y3 > (l2+r ) ) { rnew = y3 - l2; }
      }
    }
    yL[i] = lnew;   yR[i] = rnew;
  }
  return;
}


/**************************************************************************************/
/**************************************************************************************
  reconstruct(): 
  --------------

   -- calculates the left and right states about each FACE of the every cell;
   -- use RECON_TYPE_CHOICE to handle choice of the reconstruction method ; 

**************************************************************************************/
void reconstruct( double ****p_h )
{
  int d,i,j,k,l,im,jm,km,im2,jm2,km2, irecon, iface,i_s,i_e,j_s,j_e,k_s,k_e;
  double dq,geomfactor;
  static int first_time = 1; 
  struct of_geom *geom;

  TRACE_BEG;

  /* Output exactly how the reconstruction will take place for posterity :  */
  if( first_time ) {
    if( myid == printer_pid ) { 
      fprintf(stdout,"\n#################################################################\n");
      fprintf(stdout,"reconstruct(): RESCALE_B         =    %d \n", RESCALE_B ); 
      fprintf(stdout,"reconstruct(): RESCALE_R         =    %d \n", RESCALE_R ); 
      fprintf(stdout,"reconstruct(): RECON_TYPE_CHOICE =    %d \n", RECON_TYPE_CHOICE ); 
      fprintf(stdout,"reconstruct(): DEGENERATE        =    %d \n", DEGENERATE); 
      fprintf(stdout,"reconstruct(): do_recon_dir      =  (%d,%d,%d) \n", 
	      do_recon_dir[0], do_recon_dir[1], do_recon_dir[2]); 
      fprintf(stdout,"#################################################################\n");
      fflush(stdout);
    }
    first_time = 0 ; 
  }

  /***********************************************************************************
   Scale B^i so we reconstruct gdet*B^i (i.e. the fields in the divergence constraint):
  *************************************************************************************/
#if( RESCALE_B )
  ALL_LOOP {
    get_geometry(i,j,k,CENT,ncurr,geom);
    geomfactor = geom->g;
    BLOOP { p_h[i][j][k][l] *= geomfactor; }
  }
#endif
 

  /***********************************************************************************
    Section using DYNAMIC  reconstruction method : 
  *************************************************************************************/
#if( USE_LOCAL_RECON_TYPE ) 

  FACE_LOOP if(do_recon_dir[d])   {

    im = idir[d];  
    jm = jdir[d];  
    km = kdir[d];  
    im2 = 2*im;
    jm2 = 2*jm;
    km2 = 2*km;

    i_s = N1_R_S + (BUF-1)*idir[d];  j_s = N2_R_S + (BUF-1)*jdir[d];  k_s = N3_R_S + (BUF-1)*kdir[d];
    i_e = N1_R_E - (BUF-1)*idir[d];  j_e = N2_R_E - (BUF-1)*jdir[d];  k_e = N3_R_E - (BUF-1)*kdir[d];


    {   /* Perform reconstruction in this direction */ 

      for(i=i_s;i<=i_e;i++) for(j=j_s;j<=j_e;j++) for(k=k_s;k<=k_e;k++)  { 
	irecon = recon_type[i][j][k][d] ;

	/* LINEAR RECONSTRUCTION : */
	if( irecon < RECON_PPM_MM ) { 
	  PLOOP {
	    dq = linear_recon_local( irecon,
				     p_h[i-im][j-jm][k-km][l],
				     p_h[i   ][j   ][k   ][l],
				     p_h[i+im][j+jm][k+km][l]
				     ) ;

	    p_L[i+im][j+jm][k+km][d][l] = p_h[i][j][k][l] + 0.5*dq ;
	    p_R[i   ][j   ][k   ][d][l] = p_h[i][j][k][l] - 0.5*dq ;
	  }
	}
	/* HIGHER-ORDER RECONSTRUCTION : */
	else {
	  /****************************************************************************
	    Note that para() returns, respectively, the states left and right w.r.t 
	     to the center of the cell, not the left/right states about its boundary 
	  ****************************************************************************/
	  PLOOP {
	    para_local( irecon, 
			p_h[i-2*im][j-2*jm][k-2*km][l],
			p_h[i-  im][j-  jm][k-  km][l],
			p_h[i     ][j     ][k     ][l],
			p_h[i+  im][j+  jm][k+  km][l],
			p_h[i+2*im][j+2*jm][k+2*km][l],
			&(p_R[i   ][j   ][k   ][d][l]), 
			&(p_L[i+im][j+jm][k+km][d][l]) );
	  }
	}
      }
    }
  }

  /***********************************************************************************
    Section using STATIC reconstruction method : 
  *************************************************************************************/
#else


      /* LINEAR RECONSTRUCTION : */
#if( RECON_TYPE_CHOICE < RECON_PPM_MM )

  FACE_LOOP   if(do_recon_dir[d])   {  
    im = idir[d];
    jm = jdir[d];
    km = kdir[d];

    i_s = N1_R_S + (BUF-1)*idir[d];  j_s = N2_R_S + (BUF-1)*jdir[d];  k_s = N3_R_S + (BUF-1)*kdir[d];
    i_e = N1_R_E - (BUF-1)*idir[d];  j_e = N2_R_E - (BUF-1)*jdir[d];  k_e = N3_R_E - (BUF-1)*kdir[d];

    for(i=i_s;i<=i_e;i++) for(j=j_s;j<=j_e;j++) for(k=k_s;k<=k_e;k++)  { 
      PLOOP {
	dq = linear_recon( p_h[i-im][j-jm][k-km][l],
			   p_h[i   ][j   ][k   ][l],
			   p_h[i+im][j+jm][k+km][l]
			   ) ;

	p_L[i+im][j+jm][k+km][d][l] = p_h[i][j][k][l] + 0.5*dq ;
	p_R[i   ][j   ][k   ][d][l] = p_h[i][j][k][l] - 0.5*dq ;
      }
    }
  }


  /* HIGHER-ORDER RECONSTRUCTION : */
#else 

  FACE_LOOP   if(do_recon_dir[d])   {  
    im = idir[d];
    jm = jdir[d];
    km = kdir[d];
    im2 = 2*im;
    jm2 = 2*jm;
    km2 = 2*km;

#if( RESCALE_R )
    if( d == 0 )  rescale_prims_r( p_h, 0 ) ;
#endif

    i_s = N1_R_S + (BUF-1)*idir[d];  j_s = N2_R_S + (BUF-1)*jdir[d];  k_s = N3_R_S + (BUF-1)*kdir[d];
    i_e = N1_R_E - (BUF-1)*idir[d];  j_e = N2_R_E - (BUF-1)*jdir[d];  k_e = N3_R_E - (BUF-1)*kdir[d];

    for(i=i_s;i<=i_e;i++) for(j=j_s;j<=j_e;j++) for(k=k_s;k<=k_e;k++)  { 
      /****************************************************************************
	    Note that para() returns, respectively, the states left and right w.r.t 
	     to the center of the cell, not the left/right states about its boundary 
      ****************************************************************************/
      PLOOP {
	para( 
	     p_h[i- im2][j- jm2][k- km2][l],
	     p_h[i-  im][j-  jm][k-  km][l],
	     p_h[i     ][j     ][k     ][l],
	     p_h[i+  im][j+  jm][k+  km][l],
	     p_h[i+ im2][j+ jm2][k+ km2][l],
	     &(p_R[i   ][j   ][k   ][d][l]), 
	     &(p_L[i+im][j+jm][k+km][d][l]) );
      }
    }

#if( RESCALE_R )
    if( d == 0 )  rescale_prims_r( p_h, 1 ) ;
#endif

  }

/* end of RECON_TYPE_CHOICE */
#endif 

/* end of USE_LOCAL_RECON_TYPE */
#endif


  /***********************************************************************************
    Copy L/R states into ghosts if one or more dimensions is degenerate (i.e. not used)
  *************************************************************************************/
  if( DEGENERATE ) {
    /* X1-dir */
    iface = 0 ;
    if( !do_recon_dir[iface] )   {  
      for(j=N2S-1;j<=N2E+1;j++) for(k=N3S-1;k<=N3E+1;k++)  { 
	  PLOOP p_L[N1S][j][k][iface][l] = p_R[N1S][j][k][iface][l] = p_h[N1S][j][k][l] ; 
	  FACE_LOOP PLOOP  p_L[N1S+1][j][k][d][l] = p_L[N1S-1][j][k][d][l] = p_L[N1S][j][k][d][l];
	  FACE_LOOP PLOOP  p_R[N1S+1][j][k][d][l] = p_R[N1S-1][j][k][d][l] = p_R[N1S][j][k][d][l];
	}
    }
    /* X2-dir */
    iface = 1 ;
    if( !do_recon_dir[iface] )   {  
      for(i=N1S-1;i<=N1E+1;i++) for(k=N3S-1;k<=N3E+1;k++)  { 
	  PLOOP p_L[i][N2S][k][iface][l] = p_R[i][N2S][k][iface][l] = p_h[i][N2S][k][l] ; 
	  FACE_LOOP PLOOP  p_L[i][N2S+1][k][d][l] = p_L[i][N2S-1][k][d][l] = p_L[i][N2S][k][d][l];
	  FACE_LOOP PLOOP  p_R[i][N2S+1][k][d][l] = p_R[i][N2S-1][k][d][l] = p_R[i][N2S][k][d][l];
	}
    }
    /* X3-dir */
    iface = 2 ;
    if( !do_recon_dir[iface] )   {  
      for(i=N1S-1;i<=N1E+1;i++)  for(j=N2S-1;j<=N2E+1;j++) {
	  PLOOP  p_L[i][j][N3S][iface][l] = p_R[i][j][N3S][iface][l] = p_h[i][j][N3S][l] ; 
	  FACE_LOOP PLOOP  p_L[i][j][N3S+1][d][l] = p_L[i][j][N3S-1][d][l] = p_L[i][j][N3S][d][l];
	  FACE_LOOP PLOOP  p_R[i][j][N3S+1][d][l] = p_R[i][j][N3S-1][d][l] = p_R[i][j][N3S][d][l];
	}
    }
  }


  /***********************************************************************************
    Need to scale back B^i at center and boundaries: 
  *************************************************************************************/
#if( RESCALE_B )
  ALL_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom);
    geomfactor = geom->g_inv;
    BLOOP { p_h[i][j][k][l] *= geomfactor; }
  }
  FACE_LOOP ALL_LOOP { 
    get_geometry(i,j,k,d,ncurr,geom);
    geomfactor = geom->g_inv;
    BLOOP { p_L[i][j][k][d][l] *= geomfactor; }
    BLOOP { p_R[i][j][k][d][l] *= geomfactor; }
  }
#endif

  TRACE_END;

  return;
}


/**************************************************************************************/
/**************************************************************************************
  reconstruct_fast(): 
  --------------
   -- optimized version of reconstruct();
   -- linear reconstruction is performed the same, but PPM is different;
   -- calculates the left and right states about each FACE of the every cell;
   -- use RECON_TYPE_CHOICE to handle choice of the reconstruction method ; 
   -- currently does not support the "USE_LOCAL_RECON_TYPE" option;

**************************************************************************************/
void reconstruct_fast( double ****p_h )
{
  register unsigned int d,i,j,k,l,im,jm,km;
  double geomfactor, dqtmp;
  static int first_time = 1; 
  struct of_geom *geom;

  TRACE_BEG;

  /* See if we have a valid reconstruction method set: */
#if( (RECON_TYPE_CHOICE < 0) || (RECON_TYPE_CHOICE > RECON_PPM_MC) )
  fprintf(stderr,"reconstruct_fast(): Invalid RECON_TYPE_CHOICE =  %d \n", RECON_TYPE_CHOICE);
  fflush(stderr);
  fail(FAIL_BASIC,0);
#endif
 

  /* Output exactly how the reconstruction will take place for posterity :  */
  if( first_time ) {
    if( myid == printer_pid ) { 
      fprintf(stdout,"\n#################################################################\n");
      fprintf(stdout,"reconstruct_fast(): RESCALE_B          =    %d \n", RESCALE_B ); 
      fprintf(stdout,"reconstruct_fast(): RESCALE_R          =    %d \n", RESCALE_R ); 
      fprintf(stdout,"reconstruct_fast(): RESCALE_REGULARIZE =    %d \n", RESCALE_REGULARIZE ); 
      fprintf(stdout,"reconstruct_fast(): RECON_TYPE_CHOICE  =    %d \n", RECON_TYPE_CHOICE ); 
      fprintf(stdout,"reconstruct_fast(): DEGENERATE         =    %d \n", DEGENERATE); 
      fprintf(stdout,"reconstruct_fast(): do_recon_dir       =  (%d,%d,%d) \n", 
	      do_recon_dir[0], do_recon_dir[1], do_recon_dir[2]); 
      fprintf(stdout,"#################################################################\n");
      fflush(stdout);
    }
    first_time = 0 ; 
  }

  /***********************************************************************************
   Scale B^i so we reconstruct gdet*B^i (i.e. the fields in the divergence constraint):
  *************************************************************************************/
#if( RESCALE_B )
  ALL_LOOP {
    get_geometry(i,j,k,CENT,ncurr,geom);
    geomfactor = geom->g;
    BLOOP { p_h[i][j][k][l] *= geomfactor; }
  }
#endif

#if( RESCALE_REGULARIZE )   
  recon_regularize_prim( p_h, 0 );
#endif

  /***********************************************************************************
    Section using DYNAMIC  reconstruction method : 
        -- currently unsupported for this fast version; 
  *************************************************************************************/
#if( USE_LOCAL_RECON_TYPE ) 
  fprintf(stderr,"reconstruct_fast(): Does not support USE_LOCAL_RECON_TYPE !!\n"); 
  fflush(stderr);
  fail(FAIL_BASIC,0);
#endif


  /************************************************************************************/
  /************************************************************************************
    Below, we reconstruct along lines in each direction.  This requires storing  
     the line of values in a work array that is passed to the reconstruction procedures. 
  *************************************************************************************/
  /************************************************************************************/
#if( !RECON_USE_PPM )
  FACE_LOOP   if(do_recon_dir[d])   {  
    im = idir[d];
    jm = jdir[d];
    km = kdir[d];

#if( RESCALE_R )
    if( d == 0 ) rescale_prims_r( p_h, 0 ) ;
#endif

    int i_s = N1_R_S + (BUF-1)*idir[d];  int j_s = N2_R_S + (BUF-1)*jdir[d];  int k_s = N3_R_S + (BUF-1)*kdir[d];
    int i_e = N1_R_E - (BUF-1)*idir[d];  int j_e = N2_R_E - (BUF-1)*jdir[d];  int k_e = N3_R_E - (BUF-1)*kdir[d];

    for(i=i_s;i<=i_e;i++) for(j=j_s;j<=j_e;j++) for(k=k_s;k<=k_e;k++)  { 
      PLOOP {
	dqtmp = 0.5*linear_recon( p_h[i-im][j-jm][k-km][l],
				  p_h[i   ][j   ][k   ][l],
				  p_h[i+im][j+jm][k+km][l]
				  ) ;
	p_L[i+im][j+jm][k+km][d][l] = p_h[i][j][k][l] + dqtmp ;
	p_R[i   ][j   ][k   ][d][l] = p_h[i][j][k][l] - dqtmp ;
      }
    }
#if( RESCALE_R )
    if( d == 0 )  rescale_prims_r( p_h, 1 ) ;
#endif
  }

#else 

  /************************************************************************************
    X1-dir reconstruction :  
      -- may need to rescale variables to eliminate non-polynomial profiles from functions;
  *************************************************************************************/
  if( do_recon_dir[0] ) { 
#if( RESCALE_R )
    rescale_prims_r( p_h, 0 ) ;
#endif
    for(j=N2_R_S;j<=N2_R_E;j++) for(k=N3_R_S;k<=N3_R_E;k++)  PLOOP { 
      N1ALL_LOOP { p_vect[i] = p_h[i][j][k][l]; }   /* Set values to 1d array    */
      linear_recon_vect(p_vect,dq_vect,                N1S-2,N1E+2);       /* For now always need linear recon. */
      para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N1S-1,N1E+1); 
      for(i=N1_R_S+BUF-1;i<=N1_R_E-BUF+1;i++) { 
	p_L[i+1][j  ][k  ][0][l] = pR_vect[i]; 
	p_R[i  ][j  ][k  ][0][l] = pL_vect[i]; 
      }
    }
#if( RESCALE_R )
    rescale_prims_r( p_h, 1 ) ;
#endif
  }


  /************************************************************************************
    X2-dir reconstruction : 
  *************************************************************************************/
  if( do_recon_dir[1] ) { 
    for(i=N1_R_S;i<=N1_R_E;i++) for(k=N3_R_S;k<=N3_R_E;k++)  PLOOP { 
      N2ALL_LOOP { p_vect[j] = p_h[i][j][k][l]; }   /* Set values to 1d array    */
      linear_recon_vect(p_vect,dq_vect,                N2S-2,N2E+2);       /* For now always need linear recon. */
      para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N2S-1,N2E+1); 

      for(j=N2_R_S+BUF-1;j<=N2_R_E-BUF+1;j++)  { 
	p_L[i  ][j+1][k  ][1][l] = pR_vect[j]; 
	p_R[i  ][j  ][k  ][1][l] = pL_vect[j]; 
      }
    }
  }


  /************************************************************************************
    X3-dir reconstruction : 
  *************************************************************************************/
  if( do_recon_dir[2] ) { 
    for(i=N1_R_S;i<=N1_R_E;i++) for(j=N2_R_S;j<=N2_R_E;j++) PLOOP { 
      N3ALL_LOOP { p_vect[k] = p_h[i][j][k][l]; }   /* Set values to 1d array    */
      linear_recon_vect(p_vect,dq_vect,                N3S-2,N3E+2);       /* For now always need linear recon. */
      para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N3S-1,N3E+1); 
      for(k=N3_R_S+BUF-1;k<=N3_R_E-BUF+1;k++)  { 
	p_L[i  ][j  ][k+1][2][l] = pR_vect[k]; 
	p_R[i  ][j  ][k  ][2][l] = pL_vect[k]; 
      }
    }
  }

#endif

  /***********************************************************************************
    Copy L/R states into ghosts if one or more dimensions is degenerate (i.e. not used)
  *************************************************************************************/
  if( DEGENERATE ) {
    /* X1-dir */
    if( !do_recon_dir[0] )   {  
      N2ALL_LOOP  N3ALL_LOOP  { 
	PLOOP p_L[N1S][j][k][0][l] = p_R[N1S][j][k][0][l] = p_h[N1S][j][k][l] ; 
	for(i=1;i<=BUF;i++) FACE_LOOP PLOOP  p_L[N1S-i][j][k][d][l] = p_L[N1S+i][j][k][d][l] = p_L[N1S][j][k][d][l];
	for(i=1;i<=BUF;i++) FACE_LOOP PLOOP  p_R[N1S-i][j][k][d][l] = p_R[N1S+i][j][k][d][l] = p_R[N1S][j][k][d][l];
      }
    }
    /* X2-dir */
    if( !do_recon_dir[1] )   {  
      N1ALL_LOOP  N3ALL_LOOP  { 
	PLOOP  p_L[i][N2S][k][1][l] = p_R[i][N2S][k][1][l] = p_h[i][N2S][k][l] ; 
	for(j=1;j<=BUF;j++) FACE_LOOP PLOOP  p_L[i][N2S-j][k][d][l] = p_L[i][N2S+j][k][d][l] = p_L[i][N2S][k][d][l];
	for(j=1;j<=BUF;j++) FACE_LOOP PLOOP  p_R[i][N2S-j][k][d][l] = p_R[i][N2S+j][k][d][l] = p_R[i][N2S][k][d][l];
      }
    }
    /* X3-dir */
    if( !do_recon_dir[2] )   {  
      N1ALL_LOOP  N2ALL_LOOP {
	PLOOP  p_L[i][j][N3S][2][l] = p_R[i][j][N3S][2][l] = p_h[i][j][N3S][l] ; 
	for(k=1;k<=BUF;k++) FACE_LOOP PLOOP  p_L[i][j][N3S-k][d][l] = p_L[i][j][N3S+k][d][l] = p_L[i][j][N3S][d][l];
	for(k=1;k<=BUF;k++) FACE_LOOP PLOOP  p_R[i][j][N3S-k][d][l] = p_R[i][j][N3S+k][d][l] = p_R[i][j][N3S][d][l];
      }
    }
  }


  /***********************************************************************************
    Need to scale back B^i at center and boundaries: 
  *************************************************************************************/
#if( RESCALE_B )
  ALL_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom);
    geomfactor = geom->g_inv;
    BLOOP { p_h[i][j][k][l] *= geomfactor; }
  }
  FACE_LOOP ALL_LOOP  { 
    get_geometry(i,j,k,d,ncurr,geom);
    geomfactor = geom->g_inv;
    BLOOP { p_L[i][j][k][d][l] *= geomfactor; }
    BLOOP { p_R[i][j][k][d][l] *= geomfactor; }
  }
#endif

#if( RESCALE_REGULARIZE )   
  recon_regularize_prim( p_h, 1 );
#endif


  TRACE_END;

  return;
}



/**************************************************************************************/
/**************************************************************************************
  rescale_prims_r();
  --------------
   -- rescales the primitive variables so that reconstruction is done on functions 
      tha should be polynomial (not inverse polynomial);

   -- rescale_type =  0   :  scale functions by appropriate power-law of r ; 
                   =  1   :  unscale functions ;

   -- should only be called when you want to reconstruct in r-direction;
   -- assumes that x1 is independent of x2
**************************************************************************************/
void rescale_prims_r( double ****p_h, int rescale_type ) 
{
  register unsigned int i,j,k,l;
  
  double rf0, rf1, rf2, rf3,rf4,rf5,rloc;
  struct of_coord *coords;

  static int first_time = 1; 

#if( DYNAMIC_COORDINATES )
  fprintf(stdout,"rescale_prims_r(): should-not-use-RESCALE_R-with-dynamic-coordinates-\n"); fflush(stdout); fail(FAIL_BASIC,0);                                                      
#endif

#if( \
  (COORD_TYPE_CHOICE != COORD_IDENTITY ) && \
  (COORD_TYPE_CHOICE != COORD_DIAGONAL ) && \
  (COORD_TYPE_CHOICE != COORD_DIAGONAL2) && \
  (COORD_TYPE_CHOICE != COORD_DIAGONAL3) )
  fprintf(stdout,"rescale_prims_r(): should-not-use-RESCALE_R-with-this coordinate system \n"); fflush(stdout); fail(FAIL_BASIC,0);                                                   
#endif

  if( first_time ) { 

    j = N2S;  k = N3S; 
    N1ALL_LOOP { 
      get_coord(i,j,k,CENT,ncurr,coords);
      rloc = coords->x[RR];
      rf2 = rloc*rloc;  rf0 = sqrt(rf2 * rloc);   rf1 = rf2*rloc;
      rscale[0][i] = rf0;
      rscale[1][i] = rf1;
      rscale[2][i] = rloc;

      get_coord(i,j,k,FACE1,ncurr,coords);
      rloc = coords->x[RR];
      rf2 = rloc*rloc;  rf0 = sqrt(rf2 * rloc);   rf1 = rf2*rloc;
      rscale[3+0][i] = rf0;
      rscale[3+1][i] = rf1;
      rscale[3+2][i] = rloc;

      inv_rscale[0][i] = 1/rscale[0][i];
      inv_rscale[1][i] = 1/rscale[1][i];
      inv_rscale[2][i] = 1/rscale[2][i];
      inv_rscale[3][i] = 1/rscale[3][i];
      inv_rscale[4][i] = 1/rscale[4][i];
      inv_rscale[5][i] = 1/rscale[5][i];
    }
    first_time = 0; 
  }
  

  /******************************************************************************
    Multiply the primitive variables by power-laws in r 
  ******************************************************************************/
  if( rescale_type == 0 ) { 
    N1ALL_LOOP {
      rf0 = rscale[0][i];
      rf1 = rscale[1][i];
      rf2 = rscale[2][i];
      N2ALL_LOOP N3ALL_LOOP { 
	p_h[i][j][k][RHO] *= rf0;
	p_h[i][j][k][UU ] *= rf2;
	p_h[i][j][k][U1 ] *= rf1;
	p_h[i][j][k][U2 ] *= rf1;
	p_h[i][j][k][U3 ] *= rf1;
	p_h[i][j][k][B1 ] *= rf1;
	p_h[i][j][k][B2 ] *= rf1;
	p_h[i][j][k][B3 ] *= rf1;
      }
    }
  }
  /******************************************************************************
    Unscale the primitive variables and the reconstructed values ; 
  ******************************************************************************/
  else { 
    N1ALL_LOOP {
      rf0 = inv_rscale[0][i];
      rf1 = inv_rscale[1][i];
      rf2 = inv_rscale[2][i];
      rf3 = inv_rscale[3][i];
      rf4 = inv_rscale[4][i];
      rf5 = inv_rscale[5][i];
      N2ALL_LOOP N3ALL_LOOP { 
	p_h[i][j][k][RHO] *= rf0;
	p_h[i][j][k][UU ] *= rf2;
	p_h[i][j][k][U1 ] *= rf1;
	p_h[i][j][k][U2 ] *= rf1;
	p_h[i][j][k][U3 ] *= rf1;
	p_h[i][j][k][B1 ] *= rf1;
	p_h[i][j][k][B2 ] *= rf1;
	p_h[i][j][k][B3 ] *= rf1;
      }

      N2ALL_LOOP N3ALL_LOOP { 
	p_L[i][j][k][0][RHO] *= rf3;
	p_L[i][j][k][0][UU ] *= rf5;
	p_L[i][j][k][0][U1 ] *= rf4;
	p_L[i][j][k][0][U2 ] *= rf4;
	p_L[i][j][k][0][U3 ] *= rf4;
	p_L[i][j][k][0][B1 ] *= rf4;
	p_L[i][j][k][0][B2 ] *= rf4;
	p_L[i][j][k][0][B3 ] *= rf4;
      }

      N2ALL_LOOP N3ALL_LOOP { 
	p_R[i][j][k][0][RHO] *= rf3;
	p_R[i][j][k][0][UU ] *= rf5;
	p_R[i][j][k][0][U1 ] *= rf4;
	p_R[i][j][k][0][U2 ] *= rf4; 
	p_R[i][j][k][0][U3 ] *= rf4;
	p_R[i][j][k][0][B1 ] *= rf4;
	p_R[i][j][k][0][B2 ] *= rf4; 
	p_R[i][j][k][0][B3 ] *= rf4;
      }
    }
  }

  return;
}

/**************************************************************************************/
/**************************************************************************************
  recon_regularize_prim();
  --------------
   -- first transforms vectors to "x" coordinates then 
      scales the vector components of the primitive variables 
      to eliminate any geometric factors;

   -- rescale_type =  0   :  performs regularization;
                   =  1   :  performs inverse regularization; 

   -- needs to be called before reconstruction in any direction;
**************************************************************************************/
static void recon_regularize_prim( double ****p_h, int rescale_type ) 
{
  register unsigned int i,j,k,pos;

#if( RESCALE_REGULARIZE )   
  --rescale_regularize-not-supported-anymore-------wont-compile
  if( rescale_type == 0 ) { 
    ALL_LOOP { 
      regularize_prim( p_h[i][j][k], regularize_gf[i][j][k][CENT], dx_dxp_gf[i][j][k][CENT] );
    }
  }
  else { 
    ALL_LOOP { 
      unregularize_prim( p_h[i][j][k], unregularize_gf[i][j][k][CENT], dxp_dx_gf[i][j][k][CENT] );
      for(pos=FACE1;pos<=FACE3;pos++) { 
	unregularize_prim_LR( p_L[i][j][k][pos], p_R[i][j][k][pos], unregularize_gf[i][j][k][pos], dxp_dx_gf[i][j][k][pos] );
      }
    }
  }
#endif

  return;
}


///**************************************************************************************/
///**************************************************************************************
//  recon_regularize_prim();
//  --------------
//   -- first transforms vectors to "x" coordinates then 
//      scales the vector components of the primitive variables 
//      to eliminate any geometric factors;
//
//   -- rescale_type =  0   :  performs regularization;
//                   =  1   :  performs inverse regularization; 
//
//   -- needs to be called before reconstruction in any direction;
//**************************************************************************************/
//static void recon_regularize_prim( double p_h[][N2TOT][N3TOT][NP], int rescale_type ) 
//{
//  register unsigned int i,j,k,pos;
//  register double rf0, rf1, rf2;
//  register double v1,v2,b1,b2,dxdxp11,dxdxp12,dxdxp21,dxdxp22;
//
//
//#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
//  /* Just regularize */
//  if( rescale_type == 0 ) { 
//    ALL_LOOP { 
//      rf0 = regularize_gf[i][j][k][CENT][0];
//      rf1 = regularize_gf[i][j][k][CENT][1];
//      rf2 = regularize_gf[i][j][k][CENT][2];
//      p_h[i][j][k][U1] *= rf0; 
//      p_h[i][j][k][U2] *= rf1; 
//      p_h[i][j][k][U3] *= rf2; 
//      p_h[i][j][k][B1] *= rf0; 
//      p_h[i][j][k][B2] *= rf1; 
//      p_h[i][j][k][B3] *= rf2; 
//    }
//  }
//  else { 
//    ALL_LOOP { 
//      rf0 = unregularize_gf[i][j][k][CENT][0];
//      rf1 = unregularize_gf[i][j][k][CENT][1];
//      rf2 = unregularize_gf[i][j][k][CENT][2];
//      p_h[i][j][k][U1] *= rf0; 
//      p_h[i][j][k][U2] *= rf1; 
//      p_h[i][j][k][U3] *= rf2; 
//      p_h[i][j][k][B1] *= rf0; 
//      p_h[i][j][k][B2] *= rf1; 
//      p_h[i][j][k][B3] *= rf2; 
//
//      for(pos=FACE1;pos<=FACE3;pos++) { 
//	rf0 = unregularize_gf[i][j][k][pos][0];
//	rf1 = unregularize_gf[i][j][k][pos][1];
//	rf2 = unregularize_gf[i][j][k][pos][2];
//	p_L[i][j][k][pos][U1] *= rf0;    
//	p_L[i][j][k][pos][U2] *= rf1;    
//	p_L[i][j][k][pos][U3] *= rf2;    
//	p_L[i][j][k][pos][B1] *= rf0;    
//	p_L[i][j][k][pos][B2] *= rf1;    
//	p_L[i][j][k][pos][B3] *= rf2;    
//	p_R[i][j][k][pos][U1] *= rf0; 
//	p_R[i][j][k][pos][U2] *= rf1; 
//	p_R[i][j][k][pos][U3] *= rf2; 
//	p_R[i][j][k][pos][B1] *= rf0; 
//	p_R[i][j][k][pos][B2] *= rf1; 
//	p_R[i][j][k][pos][B3] *= rf2; 
//      }
//    }
//  }
//
//#elif( COORD_TYPE_CHOICE == COORD_DIAGONAL || COORD_TYPE_CHOICE == COORD_DIAGONAL2 )
//  /* Transform to x coordinates and then regularize */
//  /* Only need to take into account diagonal factors in coordinate transformation: */
//  /* We assume also that  regularize_gf[0] = dx_dxp[3][3] = 1. :  */
//  if( rescale_type == 0 ) { 
//    ALL_LOOP { 
//      rf0 = dx_dxp_gf[i][j][k][CENT][RR][1];
//      rf1 = regularize_gf[i][j][k][CENT][1] * dx_dxp_gf[i][j][k][CENT][TH][2];
//      rf2 = regularize_gf[i][j][k][CENT][2] ;
//      p_h[i][j][k][U1] *= rf0; 
//      p_h[i][j][k][U2] *= rf1; 
//      p_h[i][j][k][U3] *= rf2; 
//      p_h[i][j][k][B1] *= rf0; 
//      p_h[i][j][k][B2] *= rf1; 
//      p_h[i][j][k][B3] *= rf2; 
//    }
//  }
//  else { 
//    ALL_LOOP { 
//      rf0 = dxp_dx_gf[i][j][k][CENT][1][RR];
//      rf1 = unregularize_gf[i][j][k][CENT][1] * dxp_dx_gf[i][j][k][CENT][2][TH];
//      rf2 = unregularize_gf[i][j][k][CENT][2] ;
//      p_h[i][j][k][U1] *= rf0; 
//      p_h[i][j][k][U2] *= rf1; 
//      p_h[i][j][k][U3] *= rf2; 
//      p_h[i][j][k][B1] *= rf0; 
//      p_h[i][j][k][B2] *= rf1; 
//      p_h[i][j][k][B3] *= rf2; 
//
//      for(pos=FACE1;pos<=FACE3;pos++) { 
//	rf0 = dxp_dx_gf[i][j][k][pos][1][RR];
//	rf1 = unregularize_gf[i][j][k][pos][1] * dxp_dx_gf[i][j][k][pos][2][TH];
//	rf2 = unregularize_gf[i][j][k][pos][2] ;
//	p_L[i][j][k][pos][U1] *= rf0;    
//	p_L[i][j][k][pos][U2] *= rf1;    
//	p_L[i][j][k][pos][U3] *= rf2;    
//	p_L[i][j][k][pos][B1] *= rf0;    
//	p_L[i][j][k][pos][B2] *= rf1;    
//	p_L[i][j][k][pos][B3] *= rf2;    
//	p_R[i][j][k][pos][U1] *= rf0; 
//	p_R[i][j][k][pos][U2] *= rf1; 
//	p_R[i][j][k][pos][U3] *= rf2; 
//	p_R[i][j][k][pos][B1] *= rf0; 
//	p_R[i][j][k][pos][B2] *= rf1; 
//	p_R[i][j][k][pos][B3] *= rf2; 
//      }
//    }
//  }
//
//
//#elif( COORD_TYPE_CHOICE == COORD_MIXED )
//  /* Transform to x coordinates and then regularize */
//  /* Need to consider off-diagonal elements between RR/1 and TH/2 */
//  /* We assume also that  regularize_gf[0] = dx_dxp[3][3] = 1. :  */
//  if( rescale_type == 0 ) { 
//    ALL_LOOP { 
//      v1 = p_h[i][j][k][U1];  v2 = p_h[i][j][k][U2];  
//      b1 = p_h[i][j][k][B1];  b2 = p_h[i][j][k][B2];
//      dxdxp11 = dx_dxp_gf[i][j][k][CENT][RR][1];   dxdxp12 = dx_dxp_gf[i][j][k][CENT][RR][2];
//      dxdxp21 = dx_dxp_gf[i][j][k][CENT][TH][1];   dxdxp22 = dx_dxp_gf[i][j][k][CENT][TH][2];
//      rf1 = regularize_gf[i][j][k][CENT][1];
//      rf2 = regularize_gf[i][j][k][CENT][2];
//      p_h[i][j][k][U1]  = dxdxp11*v1  +  dxdxp12*v2 ;
//      p_h[i][j][k][U2]  = rf1*(dxdxp21*v1  +  dxdxp22*v2); 
//      p_h[i][j][k][U3] *= rf2; 
//      p_h[i][j][k][B1]  = dxdxp11*b1  +  dxdxp12*b2 ; 
//      p_h[i][j][k][B2]  = rf1*(dxdxp21*b1  +  dxdxp22*b2); 
//      p_h[i][j][k][B3] *= rf2; 
//    }
//  }
//  else { 
//    ALL_LOOP { 
//      v1 = p_h[i][j][k][U1];  v2 = p_h[i][j][k][U2];  
//      b1 = p_h[i][j][k][B1];  b2 = p_h[i][j][k][B2];
//      dxdxp11 = dxp_dx_gf[i][j][k][CENT][1][RR];   dxdxp12 = dxp_dx_gf[i][j][k][CENT][1][TH];
//      dxdxp21 = dxp_dx_gf[i][j][k][CENT][2][RR];   dxdxp22 = dxp_dx_gf[i][j][k][CENT][2][TH];
//      rf1 = unregularize_gf[i][j][k][CENT][1];
//      rf2 = unregularize_gf[i][j][k][CENT][2];
//      p_h[i][j][k][U1]  = dxdxp11*v1  +  dxdxp12*v2 ;
//      p_h[i][j][k][U2]  = rf1*(dxdxp21*v1  +  dxdxp22*v2); 
//      p_h[i][j][k][U3] *= rf2; 
//      p_h[i][j][k][B1]  = dxdxp11*b1  +  dxdxp12*b2 ; 
//      p_h[i][j][k][B2]  = rf1*(dxdxp21*b1  +  dxdxp22*b2); 
//      p_h[i][j][k][B3] *= rf2; 
//
//      for(pos=FACE1;pos<=FACE3;pos++) { 
//	dxdxp11 = dxp_dx_gf[i][j][k][pos][1][RR];   dxdxp12 = dxp_dx_gf[i][j][k][pos][1][TH];
//	dxdxp21 = dxp_dx_gf[i][j][k][pos][2][RR];   dxdxp22 = dxp_dx_gf[i][j][k][pos][2][TH];
//	rf1 = unregularize_gf[i][j][k][pos][1];
//	rf2 = unregularize_gf[i][j][k][pos][2];
//
//	v1 = p_L[i][j][k][U1];  v2 = p_L[i][j][k][U2];  
//	b1 = p_L[i][j][k][B1];  b2 = p_L[i][j][k][B2];
//	p_L[i][j][k][pos][U1]  = dxdxp11*v1  +  dxdxp12*v2 ;
//	p_L[i][j][k][pos][U2]  = rf1*(dxdxp21*v1  +  dxdxp22*v2); 
//	p_L[i][j][k][pos][U3] *= rf2; 
//	p_L[i][j][k][pos][B1]  = dxdxp11*b1  +  dxdxp12*b2 ; 
//	p_L[i][j][k][pos][B2]  = rf1*(dxdxp21*b1  +  dxdxp22*b2); 
//	p_L[i][j][k][pos][B3] *= rf2; 
//
//	v1 = p_R[i][j][k][U1];  v2 = p_R[i][j][k][U2];  
//	b1 = p_R[i][j][k][B1];  b2 = p_R[i][j][k][B2];
//	p_R[i][j][k][pos][U1]  = dxdxp11*v1  +  dxdxp12*v2 ;
//	p_R[i][j][k][pos][U2]  = rf1*(dxdxp21*v1  +  dxdxp22*v2); 
//	p_R[i][j][k][pos][U3] *= rf2; 
//	p_R[i][j][k][pos][B1]  = dxdxp11*b1  +  dxdxp12*b2 ; 
//	p_R[i][j][k][pos][B2]  = rf1*(dxdxp21*b1  +  dxdxp22*b2); 
//	p_R[i][j][k][pos][B3] *= rf2; 
//      }
//    }
//  }
//
//
//#else 
//  fprintf(stderr,"recon_regularize_prim(): COORD_TYPE_CHOICE = %d  not supported! \n",COORD_TYPE_CHOICE);
//  fflush(stderr);  fail(FAIL_BASIC,0); 
//#endif
//
//
//  return;
//}




#undef N1_R_S
#undef N1_R_E
#undef N2_R_S
#undef N2_R_E
#undef N3_R_S
#undef N3_R_E
#undef BUF
#undef DEGENERATE
#undef ONE_SIXTH 
#undef RECON_USE_PPM
#undef PPM_FACTOR 
