
#include "decs.h"



/***********************************************************************************/
/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

kom5: alfven wave
***********************************************************************************/

void init(void)
{
  double ftmp;
  void init_data(void) ;

  /* Do some checks on the pre-compiler definitions : */ 
  ftmp = (double) N1TOT;  ftmp *= N2TOT  ;  ftmp *= N3TOT  ;  ftmp *= NP; 
  if( ftmp > ( (double) UINT_MAX ) ) { 
    fprintf(stderr,"init_base(): Max. no. of elements too large !!! N = %g \n", ftmp);
    fflush(stderr);
    fail(FAIL_BASIC,0);
  }

      
  init_base();  /* Set global grid parameters and static work arrays */

  init_data();  /* Set MHD grid functions */

  return;
}

/***********************************************************************************/
/***********************************************************************************
  init_base(): 
  ---------
   -- initializes grid structure and constant grid functions; 

Kom5
***********************************************************************************/
void init_base(void)
{
	int i,j,k ;
	double X[NDIM];

	void calc_all_geom(void) ;

	/* Set the parameters that define the grid and other constants : */ 
	gam = (4./3.) ;
	cour = 0.8 ;

	t = 0. ;                     /* Initial time */ 
	GridLength[0] = 2. ;        /* Length of X0 dimension (evolution period) */ 
	GridLength[1] = 2.5 ;        /* Length of X1 dimension */ 
	GridLength[2] = 2.5 ;        /* Length of X2 dimension */ 
	GridLength[3] = 2.5 ;        /* Length of X3 dimension */ 

	dx[0] = 1.e-5;               /* Discretization size in X0 direction */
	dx[1] = GridLength[1]/N1 ;   /* Discretization size in X1 direction */
	dx[2] = dx[1];   /* Discretization size in X2 direction */
	dx[3] = dx[1];   /* Discretization size in X3 direction */

	startx[0] =  0.;               /* Set the Physical Minimum boundary  */
	startx[1] = -1.;               /* Set the Physical Minimum boundary  */
	startx[2] = -1.;               /* Set the Physical Minimum boundary  */
	startx[3] = -1.;               /* Set the Physical Minimum boundary  */
	SDLOOP1 startx[i] -= NG*dx[i];  /* Adjust to numerical minimum boundary */

	dV = dx[1]*dx[2]*dx[3] ;     /* Coordinate volume */

	DT_out[OUT_ASCII]   = GridLength[0]/10. ;	/* dumping frequency */
	DT_out[OUT_IMAGE]   = GridLength[0]/40. ;	/* image frequency */
	DT_out[OUT_HISTORY] = 1. ;                	/* logfile frequency */
	DT_out[OUT_HDF5]    = GridLength[0]/200. ;      /* hdf frequency */
	DT_out[OUT_SDF]     = GridLength[0]/200. ;	/* sdf frequency */
	DT_out[OUT_STAT]    = 80.;                      /* statistics dumps frequency */
	DT_out[OUT_STAT2]   = DT_out[OUT_HDF5];         /* statistics dumps frequency */
	DT_out[OUT_RADFLUX] = DT_out[OUT_HDF5];         /* radiative flux dumps       */
	n_restart = 100;                /* number of time steps between restart dumps */
	

	calc_all_geom() ;  /* Calculate the static metric grid functions */

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
  int i, j, k ; 
  double X[NDIM], r;
  double fr(double R,double min,double max);
  double x1,x2,f1,f2;
  double quartic_interp(double x , double x1, double x2, double f1, double f2) ;
  double th, th_tot, rate;
  double ur[5], Br[5];
  double rho, pressure;
  double u1l, u2l, u3l, B1l, B2l, B3l;
  double u1r, u2r, u3r, B1r, B2r, B3r;
  void alfv( double ux, double Blx, double Bly, double p, double d, 
	     double teta, double ur[5], double Br[5], double Gamma ) ;


  x1 = -0.5;
  x2 = 0.;

  pressure = 1.;
  rho      = 1.;

  u1l = 0.;
  u3l = 0.;

  B1l = 3.;
  B2l = 3.;
  B3l = 0.;


  th_tot = M_PI ; 

  /* Get right states: */ 
  alfv( u1l, B1l, B2l, pressure, rho, th_tot, ur, Br, gam ) ;

  B1r = Br[1];
  B2r = Br[2];
  B3r = Br[3];
  
  u1r = ur[1];
  u2r = ur[2];
  u3r = ur[3];

  /* Loop over the grid  */
  LOOP {
    coord(i,j,k,CENT,X) ; 

    //r = sqrt(X[1]*X[1] + X[2]*X[2]) ;

    if( X[1] <= x1 ) { 
      p[i][j][k][U1]  = u1l;
      p[i][j][k][U2]  = u2l;
      p[i][j][k][U3]  = u3l;
      p[i][j][k][B1]  = B1l;
      p[i][j][k][B2]  = B2l;
      p[i][j][k][B3]  = B3l;
      p[i][j][k][RHO] = rho;
      p[i][j][k][UU]  = pressure/(gam - 1.) ;
    }
    else if( X[1] >= x2 ) { 
      p[i][j][k][U1]  = u1r;
      p[i][j][k][U2]  = u2r;
      p[i][j][k][U3]  = u3r;
      p[i][j][k][B1]  = B1r;
      p[i][j][k][B2]  = B2r;
      p[i][j][k][B3]  = B3r;
      p[i][j][k][RHO] = rho;
      p[i][j][k][UU]  = pressure/(gam - 1.) ;
    }
    else { 
      rate = (X[1] - x1) / (x2 - x1);
      th = th_tot * rate * rate * ( 3. - 2.*rate );
      alfv( u1l, B1l, B2l, pressure, rho, th, ur, Br, gam ) ;
      
      p[i][j][k][U1]  = ur[1];
      p[i][j][k][U2]  = ur[2];
      p[i][j][k][U3]  = ur[3];
      p[i][j][k][B1]  = Br[1];
      p[i][j][k][B2]  = Br[2];
      p[i][j][k][B3]  = Br[3];
      p[i][j][k][RHO] = rho;
      p[i][j][k][UU]  = pressure/(gam - 1.) ;

    }
  }

  fixup(p)  ;     /* Set floor, and correct unphysical states */
  bounds(p,0) ;	  /* enforce boundary conditions */

}

double fr(double R,double min_r,double max_r)
{
        double r_outer,r_inner,dR ;

        r_outer = GridLength[1]/12. ;
        r_inner = 0.8*r_outer ;
        dR = (r_outer - r_inner) ;

        if(R > r_outer) return(min_r) ;
        if(R < r_inner) return(max_r) ;
        else {
                return( exp(
                        log(max_r)*(r_outer - R)/dR +
                        log(min_r)*(R - r_inner)/dR)) ;
        }
}

double quartic_interp(double x , double x1, double x2, double f1, double f2) 
{
  double dx1,dx2,dx12;

  dx1  = x-x1;
  dx2  = x-x2;
  dx12 = x1-x2;
//  return( 
//	 (f2*dx1*dx1*(2*x + x1 - 3*x2) - f1*dx2*dx2*(2*x - 3*x1 + x2))/(dx12*dx12*dx12) 
//	 );
  return( (f1*x - f2*x + f2*x1 - f1*x2)/dx12   );
}



/*********************************************************************

  aflv(): 
 -------
 -- komissarov's Alfven wave solver  used in his 1997 paper and 
    the Aflven/Compound wave tests

 -- manually translated by Scott Noble into C ;

 -- so that we don't have to change all the array indices by -1, I am
    increasing the lengths of all arrays

-- note that the vectors are indexed ala:

   1 -> x
   2 -> y
   3 -> z
   4 -> t

*********************************************************************/


/*=============================================================
*  Analytical solution for relativistic Alfven wave (shock/continious)   
*  This gives the right state of the wave corresponding to rotation 
*  by angle TETA assuming that in the left state  Uy=Uz=0, Hz=0    
*============================================================= */
void alfv( double ux, double hlx, double hly, double p, double d, 
	   double teta, double ur[5], double hr[5], double Gamma ) 
{
  /*  input: */
  /*
    left state: ux -  x-component of 4-velocity
    hlx,hly -  x- and y-components of magnetic field (usual b)
    p,d - pressure and desity 
    del - fraction of electrons (not needed for polytropes) 

    right state: ur(4) - 4-velocity;
    hr(3) - usual magnetic field b     
    output:
  */

  /*  procedures:  */
  double  w(double p, double d, double Gamma ) ; //! enthalpy function

  /*  variabes:   */
  double bl[5],br[5],ul[5],ww,fi ;
  double tmp,ey,ez,ey2,ez2,kappa,dby,b2   ;
  double hl[4],lorx,tmp1,tmp2,va,vx  ; 
  double ult[5],urt[5],blt[5],brt[5] ;
  int    i   ;
  /*   left state: */
  ul[1] = ux ;
  ul[2] = 0. ;
  ul[3] = 0. ;
  ul[4] = sqrt(1.+ux*ux) ;
  hl[1] = hlx;
  hl[2] = hly;
  hl[3] = 0. ;
  tmp = 0. ;
  for( i = 1; i <= 3 ; i++ ) {
    tmp = tmp+ul[i]*hl[i] ;
  }
  bl[4] = tmp ;
  for( i = 1; i <= 3 ; i++ ) {
    bl[i] = hl[i]+ul[i]*bl[4]/ul[4]  ;
  }
  /*   alfwen speed:*/
  ww = w(p,d,Gamma)  ;
  b2 = -bl[4]*bl[4];
  for( i = 1; i <= 3 ; i++ ) {
    b2 = b2+bl[i]*bl[i];
  }
  fi = b2 + ww ;
  vx = ul[1]/ul[4] ;
  tmp = ul[4]*sqrt(fi);
  tmp1=(bl[1]-tmp*vx)/(bl[4]-tmp) ;
  tmp2=(bl[1]+tmp*vx)/(bl[4]+tmp) ;
  va = MAX(tmp1,tmp2) ;
  /* c       print *,'alfven speed -',va  */
  /*   transformation to the alfvenic wave frame:  */
  vx = va;
  lorx = 1./sqrt(1. - vx*vx)  ;
  for( i = 1; i <= 4 ; i++ ) {
    ult[i]=ul[i]  ;
    blt[i]=bl[i] ; 
  }
  ul[1] = lorx*(ult[1]-vx*ult[4]) ;
  ul[4] = lorx*(ult[4]-vx*ult[1]) ; 
  bl[1] = lorx*(blt[1]-vx*blt[4]) ;
  bl[4] = lorx*(blt[4]-vx*blt[1]) ; 
  /****** parameters of magnetic ellipse *********/
  ey2 = (bl[2]*fi/ww)*(bl[2]*fi/ww);
  ez2 = (fi/ww)*(bl[2]*ul[4])*(bl[2]*ul[4]);
  ey = sqrt(ey2) ;
  ez = sqrt(ez2)  ;
  dby = - bl[2]*b2/ww ;
  /**** wave frame solution ***************/
  br[1]=bl[1]   ;
  tmp = (cos(teta)/ey)*(cos(teta)/ey) + (sin(teta)/ez)*(sin(teta)/ez);
  tmp = sqrt(1./tmp)  ;
  br[2] = tmp*cos(teta) + dby  ;
  br[3] = tmp*sin(teta) ;
  kappa = ul[1]/bl[1]    ;
  ur[1] = ul[1]   ;
  ur[2] = kappa*(br[2]-bl[2]) ;
  ur[3] = kappa*(br[3]-bl[3])   ;
  ur[4] = sqrt(1.+ur[1]*ur[1]+ur[2]*ur[2]+ur[3]*ur[3]) ;
  tmp = ur[1]*br[1]+ur[2]*br[2]+ur[3]*br[3] ;
  br[4] = tmp/ur[4]  ;
  /*** transformation to the initial frame:  */
  for( i = 1; i <= 4 ; i++ ) {
    urt[i]=ur[i] ; 
    brt[i]=br[i]  ;
  }
  vx = -va;
  lorx = 1./sqrt(1. - vx*vx)  ;
  ur[1] = lorx*(urt[1]-vx*urt[4]) ;
  ur[4] = lorx*(urt[4]-vx*urt[1]) ; 
  br[1] = lorx*(brt[1]-vx*brt[4]) ;
  br[4] = lorx*(brt[4]-vx*brt[1]) ; 
  /**** h-field: */
  hr[1]= br[1]*ur[4] - ur[1]*br[4] ;
  hr[2]= br[2]*ur[4] - ur[2]*br[4] ;
  hr[3]= br[3]*ur[4] - ur[3]*br[4] ;

  return;
}

/* find the enthalpy  */
double w( double p, double d, double Gamma ) 
{
  return( d + Gamma*p/(Gamma-1.) ) ;
}
