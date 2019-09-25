/***********************************************************************************
***********************************************************************************/

#include "decs.h"

/* apply floors to density, internal energy */

void fixup(double ****pv)
{
   int i,j,k;
   void fixup1zone(double *pv, double *gamma, struct of_coord *coords)  ;
   double gamma;
   struct of_geom *geom ;
   struct of_coord *coords;

  TRACE_BEG;

   LOOP  { 
     /* Get geom now since icurr, etc.  are needed by fail() calls */ 
     get_geometry(i,j,k,CENT,ncurr,geom   ) ; 

     if( gamma_calc(pv[i][j][k],geom,&gamma) ) { 
       fail( FAIL_GAMMA_CALC,0 );
       gamma = 1.;
     }

     get_coord(   i,j,k,CENT,ncurr,coords) ; 
     fixup1zone(pv[i][j][k], &gamma,coords);
   }

  TRACE_END;

}

/*******************************************************************************************/
/*******************************************************************************************
  fixup1zone()
  ------------------
     -- routine that contains all the non-interpolation based fixup schemes;
     -- these currently include 
          1) imposing the floor on  rho and uu 
          2) setting the ceiling on gamma

*******************************************************************************************/
void fixup1zone(double *pv, double *gamma, struct of_coord *coords) 
{
  
  /**********************************************************************
    Impose the floor on the density and internal energy density:  
  **********************************************************************/

#if( ALLOW_NEGATIVE_DENSITIES )
  /* If we accept rho<0, we still need to set floor at t=0 */
  if( t < 0.4*dx[0] ) {
    fixup_floor(pv,coords); 
  }
#else 
  /* Use the floor if we refuse negative densities */
  fixup_floor(pv,coords); 
#endif


  /*****************************************************************************
   In regions where bsq >> p, large gradients in p can arise, and these 
   large gradients are more likely to to do thermal expansion when p >> rho. 
    p >> rho is only expected in the jet, which needs a floor anyway and 
    so this just sets another arbitrary prescription for the floor. 
  ******************************************************************************/

#if( USE_TEMPERATURE_MAX ) 
  if( check_Tmax(pv) ) {  fixup_Tmax(pv) ; }
#endif


  /**********************************************************************
    Impose ceiling on gamma :
  **********************************************************************/
  if(check_gamma(*gamma))  {  fixup_gamma(pv,gamma); }

  return;
}

/*******************************************************************************************/
/*******************************************************************************************
  check_floor()
  ------------------
     -- returns with non-zero status if the floor 
     -- see fixup_floor() for the routine that actually imposes the floor;
*******************************************************************************************/
int check_floor(double *pv, struct of_coord *coords)
{
  int ret; 

  ret = 0 ; 

#if( !ALLOW_NEGATIVE_DENSITIES ) 
  if( pv[RHO] < coords->rhoflr )   ret +=   1;
  if( pv[UU]  < coords->uuflr  )   ret +=   2;
#endif 

  return( ret ); 
}

/*******************************************************************************************/
/*******************************************************************************************
  fixup_floor()
  ------------------
     -- imposes the floor and temperature maximum on the primitive variables;
*******************************************************************************************/
void fixup_floor(double *pv, struct of_coord *coords)
{
  unsigned short int ret = 0;

#if( !ALLOW_NEGATIVE_DENSITIES ) 
  if( pv[RHO] < coords->rhoflr )  {  pv[RHO] = coords->rhoflr;  /* fprintf(stdout, "rho floor reached\n"); fflush(stdout); */ ret = 1;  }    
  if( pv[ UU] < coords->uuflr  )  {  pv[ UU] = coords->uuflr ;  /* fprintf(stdout, "uu floor reached\n");  fflush(stdout); */ ret = 1;  }    
  if( ret                      )  {  fail(FAIL_FLOOR,0);                  }
#endif

    return;
}

/*******************************************************************************************/
/*******************************************************************************************
  check_gamma()
  ------------------
     -- returns with non-zero status if gamma exceeds the either of the two 
        gamma thresholds;
     -- requires gamma;
     -- we assume that GAMMAMAX2 > GAMMAMAX 
     -- see fixup_gamma() for the routine that actually imposes the gamma ceiling;
*******************************************************************************************/
int check_gamma(double gamma)
{
  if( gamma <= GAMMAMAX )    return( 0 ); 
  else{ 
    if( gamma > GAMMAMAX2 )  return( 2 ); 
    else                     return( 1 );
  }
}

/*******************************************************************************************/
/*******************************************************************************************
  fixup_gamma()
  ------------------
     -- rescales velocities so that their gamma  matches GAMMAMAX ;
*******************************************************************************************/
void fixup_gamma(double *pv, double *gamma)
{
  double f; 
  f = sqrt( (GAMMAMAX-1.)*(GAMMAMAX+1.)/ ((*gamma-1.)*(*gamma+1.)) );
  pv[U1] *= f ;	
  pv[U2] *= f ;	
  pv[U3] *= f ;	
  *gamma = GAMMAMAX;
  fail( FAIL_GAMMA_MAX,0 ) ;
  return;
}

/*******************************************************************************************/
/*******************************************************************************************
  check_Tmax()
  ------------------
     -- checks to see if the temperature (really just u/rho) is greater than the ceiling;
*******************************************************************************************/
int check_Tmax(double *pv)
{
#if( USE_TEMPERATURE_MAX ) 
  return( pv[UU] > pv[RHO]*TEMPERATURE_MAX ); 
#else
  return(0);
#endif
}

/*******************************************************************************************/
/*******************************************************************************************
  fixup_Tmax()
  ------------------
     -- rescales u to match max. allowed temperature;
*******************************************************************************************/
void fixup_Tmax(double *pv)
{
  pv[UU] = pv[RHO]*TEMPERATURE_MAX ; 
  fail( FAIL_TMAX,0 ) ;
  return;
}



/*******************************************************************************************/
/*******************************************************************************************
  fixup_global()
  ------------------
     -- performs non-local or extended operations on the primitive variables 
        that must be done with the global arrays and must be done before the boundary 
        conditions ; 
     -- this routine currently can (if macros are set accordingly): 
         1) 
     -- these currently include 
          1) imposing the floor on  rho and uu 
          2) setting the ceiling on gamma

*******************************************************************************************/
void fixup_global(double ****prim ) 
{
  int i,j,k,l,g;
  double f,f2;
  double bsq, uu1, uu2,g1,g2,g3;
  struct of_geom *geom;

#if( (!FIX_CUTOUT_GAMMA)  && (!FIX_CUTOUT_PRESSURE) )
  return; 
#endif
  
  /*************************************************************************************
    Smooth the gamma,pressure in x2-direction at points adjacent to the cutout boundary when:
       (bsq > 10*p) AND (p[j=0] is very different from p[j=1])
     -- uses logarithmic (aka geometric) averaging ;
  *************************************************************************************/

#if( BC_TYPE_CHOICE == BC_SPHERICAL_CUTOUT ) 
  /* Lower face */
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    j = N2S;
    N1ALL_LOOP N3ALL_LOOP  {
      get_geometry(i,j,k,CENT,ncurr,geom); 
      uu1 = prim[i][N2S][k][UU];  uu2 = prim[i][N2S+1][k][UU];
      g1 = p_gamma[i][N2S][k]; g2 = p_gamma[i][N2S+1][k]; 
      bsq = -1. ;
#if( FIX_CUTOUT_PRESSURE ) 
      if( REL_DIFF_FUNC(uu1,uu2) > 0.1 ) { 
	bsq = bsq_calc(prim[i][N2S][k],geom);
	if( bsq > 10.*uu1 ) { 
	  prim[i][N2S][k][UU] = prim[i][N2S+1][k][UU] = sqrt(uu2*uu1); 
	  pflag[i][N2S][k] = pflag[i][N2S+1][k] = MAX(pflag[i][N2S][k],pflag[i][N2S+1][k]);
	  fail(FAIL_CUTOUT_PRESSURE,0);
	}
      }
#endif
#if( FIX_CUTOUT_GAMMA ) 
      if( REL_DIFF_FUNC(g1,g2) > 0.1 ) { 
	if( bsq < 0. ) { bsq = bsq_calc(prim[i][N2S][k],geom); } 
	if( bsq > 10.*uu1 ) { 
	  g3 = sqrt(g1*g2);  f  = g3/g1;  f2 = g3/g2;
	  p_gamma[i][N2S][k] = p_gamma[i][N2S+1][k] = g3;
	  prim[i][N2S  ][k][U1] *= f;
	  prim[i][N2S  ][k][U2] *= f;
	  prim[i][N2S  ][k][U3] *= f;
	  prim[i][N2S+1][k][U1] *= f2;
	  prim[i][N2S+1][k][U2] *= f2;
	  prim[i][N2S+1][k][U3] *= f2;
	  pflag[i][N2S][k] = pflag[i][N2S+1][k] = MAX(pflag[i][N2S][k],pflag[i][N2S+1][k]);
	  fail(FAIL_CUTOUT_GAMMA,0);
	}
      }
#endif
    }
  }

  /* Upper face */
  if( bc_pid[2][BCUP] == BC_PHYS ) { 
    j = N2E;
    N1ALL_LOOP  N3ALL_LOOP  {
      get_geometry(i,j,k,CENT,ncurr,geom);  /* Assuming axisymmetric spacetime */
      uu1 = prim[i][N2E][k][UU];  uu2 = prim[i][N2E-1][k][UU];
      g1 = p_gamma[i][N2E][k]; g2 = p_gamma[i][N2E-1][k]; 
      bsq = -1. ;
#if( FIX_CUTOUT_PRESSURE ) 
      if( REL_DIFF_FUNC(uu1,uu2) > 0.1 ) { 
	bsq = bsq_calc(prim[i][N2E][k],geom);
	if( bsq > 10.*uu1 ) { 
	  prim[i][N2E][k][UU] = prim[i][N2E-1][k][UU] = sqrt(uu2*uu1); 
	  pflag[i][N2E][k] = pflag[i][N2E-1][k] = MAX(pflag[i][N2E][k],pflag[i][N2E-1][k]);
	  fail(FAIL_CUTOUT_PRESSURE,0);
	}
      }
#endif
#if( FIX_CUTOUT_GAMMA ) 
      if( REL_DIFF_FUNC(g1,g2) > 0.1 ) { 
	if( bsq < 0. ) { bsq = bsq_calc(prim[i][N2E][k],geom); } 
	if( bsq > 10.*uu1 ) { 
	  g3 = sqrt(g1*g2);  f = g3/g1;  f2 = g3/g2;
	  p_gamma[i][N2E][k] = p_gamma[i][N2E-1][k] = g3;
	  prim[i][N2E  ][k][U1] *= f;
	  prim[i][N2E  ][k][U2] *= f;
	  prim[i][N2E  ][k][U3] *= f;
	  prim[i][N2E-1][k][U1] *= f2;
	  prim[i][N2E-1][k][U2] *= f2;
	  prim[i][N2E-1][k][U3] *= f2;
	  pflag[i][N2E][k] = pflag[i][N2E-1][k] = MAX(pflag[i][N2E][k],pflag[i][N2E-1][k]);
	  fail(FAIL_CUTOUT_GAMMA,0);
	}
      }
#endif
    }
  }

#endif
    
  return;
}

/******************************************************************************************/
/******************************************************************************************
  check_entropy_eq():
 ------------------
   -- integrates the entropy density and then calculates the primitive variables 
      from the new entropy 
   -- inversion ignores the total energy equation; 
   -- this routine should be called before p_L,p_R,F[] are written over;

******************************************************************************************/
int check_entropy_eq(double th, double uu, double bsq, double rho)
{
#if( USE_ENTROPY_EQ )
//return( (th <= M_PI/6.0) || (th >= (5.0/6.0)*M_PI) || (uu < BETA_MIN * bsq) || (rho < BETA_MIN * bsq) ); 
  return( (uu < BETA_MIN * bsq) );
#else
  return(0);
#endif
}

/******************************************************************************************/
/******************************************************************************************
  fixup_entropy_eq():
 ------------------
   -- integrates the entropy density and then calculates the primitive variables 
      from the new entropy ;
   -- inversion ignores the total energy equation; 
   -- this routine should be called before p_L,p_R,F[] are written over;
   -- prim_old[] are the primitive variables from the timestep we are integrating 
       i.e. the n^th step  ("p_old[]" array)
   -- performs the fix for just one cell indicated by the indices in the arguments;

******************************************************************************************/
int fixup_entropy_eq(int i, int j, int k, double *prim_old, double *prim, 
		      struct of_geom *geom, struct of_geom *geom_old, double *gamma_out)
{
  int d,l,ret=0;
  int ih,jh,kh;
  double s,S,S_old, Dt, gamma,cmin,cmax,ctop,prim_backup[NP];
  
  double Ftot_m[NDIM-1],Ftot_p[NDIM-1]; 
  double Fl,Fr,Ul,Ur;
  int Utoprim_1d_ee(double U[NP], double S, struct of_geom *geom, double prim[NP], double *gamma );
  int Utoprim_1d_ee2(double U[NP], double S, struct of_geom *geom, double prim[NP], double *gamma );

  static const int idir[] = {1,0,0};
  static const int jdir[] = {0,1,0};
  static const int kdir[] = {0,0,1};

  PLOOP prim_backup[l] = prim[l]; 

  /* Set the time step depending on what substep we are on */
  if( n_substep == (N0-1) ) {  Dt = dx[0];     }
  else                      {  Dt = 0.5*dx[0]; }

  /* Set local variables : */
  gamma_calc(prim_old,geom_old,&gamma);
  S_old = geom_old->g * (gam-1.)*prim_old[UU]*gamma/( geom_old->alpha * pow( prim_old[RHO], (gam-1.) ) );


  /* Lower face loop : */
  FACE_LOOP { 
    cmin = c_gf[i][j][k][d][0];  cmax = c_gf[i][j][k][d][1];  ctop = MAX(cmin,cmax); 
    s = (gam-1.) * p_L[i][j][k][d][UU] * pow(p_L[i][j][k][d][RHO],1.-gam);
    Ul = s * ucon_L[i][j][k][d][0];
    Fl = s * ucon_L[i][j][k][d][1];

    s = (gam-1.) * p_R[i][j][k][d][UU] * pow(p_R[i][j][k][d][RHO],1.-gam);
    Ur = s * ucon_R[i][j][k][d][0];
    Fr = s * ucon_R[i][j][k][d][1];

#if( FLUX_TYPE_CHOICE == FLUX_HLLE ) 	
    Ftot_m[d] = (cmax*Fl + cmin*Fr - cmax*cmin*(Ur - Ul))/(cmax+cmin+SMALL) ; /* HLLE flux: */
#elif( FLUX_TYPE_CHOICE == FLUX_LF ) 	
    Ftot_m[d] = 0.5*(Fl + Fr - ctop*(Ur - Ul)) ;	  /* Lax-Friedrichs flux:  */
#endif 
  }
    
  /* Upper face loop : */
  FACE_LOOP { 
    ih = i+idir[d]; jh = j+jdir[d]; kh = k+kdir[d]; 
    cmin = c_gf[ih][jh][kh][d][0];  cmax = c_gf[ih][jh][kh][d][1];  ctop = MAX(cmin,cmax); 

    s = (gam-1.) * p_L[ih][jh][kh][d][UU] * pow(p_L[ih][jh][kh][d][RHO],1.-gam);
    Ul = s * ucon_L[ih][jh][kh][d][0];
    Fl = s * ucon_L[ih][jh][kh][d][1];

    s = (gam-1.) * p_R[ih][jh][kh][d][UU] * pow(p_R[ih][jh][kh][d][RHO],1.-gam);
    Ur = s * ucon_R[ih][jh][kh][d][0];
    Fr = s * ucon_R[ih][jh][kh][d][1];

#if( FLUX_TYPE_CHOICE == FLUX_HLLE ) 	
    Ftot_p[d] = (cmax*Fl + cmin*Fr - cmax*cmin*(Ur - Ul))/(cmax+cmin+SMALL) ; /* HLLE flux: */
#elif( FLUX_TYPE_CHOICE == FLUX_LF ) 	
    Ftot_p[d] = 0.5*(Fl + Fr - ctop*(Ur - Ul)) ;	  /* Lax-Friedrichs flux:  */
#endif 
  }

  /* Update entropy : */
  S = S_old + Dt*(
		  - ( Ftot_p[0] - Ftot_m[0] ) * invdx[1] 
		  - ( Ftot_p[1] - Ftot_m[1] ) * invdx[2] 
		  - ( Ftot_p[2] - Ftot_m[2] ) * invdx[3] 
		  );

  /* Perform the inversion: */
  ret = Utoprim_1d_ee(U_gf[n_U][i][j][k],S,geom,prim,&gamma);

  if( (!ret) && (check_gamma(gamma)==2) ) { /* check here so we can try for another solution: */
    PLOOP { prim[l] = prim_backup[l]; }   
    ret = 1; 
    gamma = 1.;
  }   
  if( ret ) { 
    ret = Utoprim_1d_ee2(U_gf[n_U][i][j][k],S,geom,prim,&gamma);

    if( (!ret) && (check_gamma(gamma)==2) ) { 
      PLOOP { prim[l] = prim_backup[l]; }   
      ret = 1; 
    } 
    if( ret ) { 
      gamma = 1.; /* Don't use bad gamma */
    }

  }
  
//  if( ret ) { 
//    fprintf(stdout,"ee0 %d %d %d %d : %d \n",nstep,(globalpos[1]+icurr-N1S),
//	    (globalpos[2]+jcurr-N2S),(globalpos[3]+kcurr-N3S),ret);
//    fprintf(stdout,"ee1 ");
//    PLOOP fprintf(stdout," %28.18e ",prim_old[l]);
//    fprintf(stdout,"\n");
//    fprintf(stdout,"ee2 ");
//    PLOOP fprintf(stdout," %28.18e ",prim[l]);
//    fprintf(stdout,"\n");
//    fprintf(stdout,"ee3 ");
//    PLOOP fprintf(stdout," %28.18e ",U_gf[n_U][i][j][k][l]);
//    fprintf(stdout,"\n");
//    fprintf(stdout,"ee4 %28.18e  %28.18e \n",S_old,S);
//  }
//  else { 
//    fprintf(stdout,"eetest(%d,%d,%d,%d) Pn : ",
//	    nstep,(globalpos[1]+i-N1S),(globalpos[2]+j-N2S),(globalpos[3]+k-N3S));
//    PLOOP { fprintf(stdout," %26.16e",prim[l]); }
//    fprintf(stdout,"\n");
//  }

  *gamma_out = gamma;

  if( !ret )  fail(FAIL_USE_EE,0) ; 

  return(ret);
}


/*******************************************************************************************/
/*******************************************************************************************
  fixup_interp_prim(): 

  -- figures out (w/ pflag[]) which stencil to use to interpolate bad point from neighbors;

  -- interpolation is done using a 3x3x3 cube of parent (interpolating) cells 
     centered about the child (interpolated) cell; 

  -- In order to describe the stencils, let us use the following index notation 
     that is defined relative to the child cell.  A give cell has index 
      (d_i + 3*d_j + 9*d_k)   where "d_[i,j,k]"  is the change in [i,j,k] 
     indices from the center cell.  The figure below shows slices  in the x1,x2
     plane at different x3.  The  "0" marks the child cell. 


  ^ x2           d_k=-1            d_k=0             d_k=1
  |            -7   -6   -5       2   3   4       11   12   13
  |           -10   -9   -8      -1   0   1        8    9   10       
  |           -13  -12  -11      -4  -3  -2        5    6    7
  ----> x1   


  -- by numbering the cells in this way, we can exploit the parity symmetry of this 
     interpolation scheme;

  -- note that this routine uses pflag[] in the ghost zones, so be sure to set 
     pflag[]'s boundary conditions

  -- ngood_min[] is set so that this routine performs an interpolation when there are 
     an appropriate number of good neighboring cells;
      -- 1D runs: we require there to be at least two good neighbors ;
      -- 2D runs: we require there to be 2 unique (nondegenerate) lines (i.e. criss-cross)
      -- 3D runs: we are more strict in order to ensure that interpolation "samples" 
                  all dimensions and require that there be at least 4 lines;

  -- Note that ngood_min[] represents each dimensions number of  lines including degeneracies;
         for example, in 2D each "diagonal" line is thrice degenerate plus there is a line 
         connect two ghosts of the center cell.

 *******************************************************************************************/

void fixup_interp_prim( double ****pv )
{
  int i, j, k, l, ind1, i_ret;
  int ii, jj, kk, ngood, ngood2;
  int pf[27];
  int o=13;
  double pold[NP];
  static const int ngood_min[] = { 13, 9, 6, 4 };
  double gamma, f, gamma_old;
  struct of_geom *geom, *geom2 ;


  TRACE_BEG;

  if( !failure_exists ) { 
    TRACE_END;
    return;
  }


  /* Fix the interior points first */
  LOOP {
    if( pflag[i][j][k] == PFLAG_INTERP_PRIM ) { 
      icurr = i; jcurr = j; kcurr = k;
 
      /* Collect the local flags and make pf true if it is a good point */
      FLOOP pold[l] = pv[i][j][k][l]; 

      for(ii=-1; ii<=1; ii++) for(jj=-1; jj<=1; jj++) for(kk=-1; kk<=1; kk++) { 
	ind1 = ii + 3*(jj + 3*kk);
	pf[o+ind1] = !pflag[i+ii][j+jj][k+kk];
      }

      /* Cancel mirror cells, and count the number of good pairs */
      ngood = 0;
      for(ii=1; ii<=13 ; ii++) { 
	if( pf[o+ii] && pf[o-ii]  ) { ngood++ ;                }
	else                        { pf[o+ii] = pf[o-ii] = 0; }
      }
      
      /* Leave unchanged from previous step's value if there are too few good neighbors */
      get_geometry(i,j,k,CENT,ncurr,geom) ; 
      if( ngood < ngood_min[n_spatial_dims] ) { 
	gamma_calc(pv[i][j][k],geom,&(p_gamma[i][j][k]));
	fail( FAIL_INTERP_PRIM,0 ) ; 
      }
      /* Only fix via interpolation if it is sufficiently covered */
      else{ 

	/* Average the good points : */ 
	FLOOP pv[i][j][k][l] = 0.;
	for(ii=-1; ii<=1; ii++) for(jj=-1; jj<=1; jj++) for(kk=-1; kk<=1; kk++) { 
	  if( pf[o + ii + 3*(jj + 3*kk)] ) { 
	    FLOOP  pv[i][j][k][l] += pv[i+ii][j+jj][k+kk][l];
	  }
	}
	ngood *= 2;
	FLOOP  pv[i][j][k][l]  /=   ngood; 

	/* Make sure that we can calculate a valid lorentz factor : */
	i_ret = gamma_calc(pv[i][j][k],geom,&gamma);

	/* Interpolate over gamma if fixed gamma is bad and if there is a valid stencil: */
	/* note that we cannot rescale v^i if (i_ret!=0)  since v^i may be bad */
	if( check_gamma(gamma) ) { 
	  gamma_old = gamma;
	  gamma = 0.;
	  ngood2 = 0;
	  for(ii=-1; ii<=1; ii++) for(jj=-1; jj<=1; jj++) for(kk=-1; kk<=1; kk++) { 
	    if( pf[o + ii + 3*(jj + 3*kk)] ) { 
	      get_geometry(i+ii,j+jj,k+kk,CENT,ncurr,geom2) ; 
	      i_ret = gamma_calc(pv[i+ii][j+jj][k+kk],geom2,&f);  
	      /* only use good gamma's */
	      if( !i_ret ) {  
		gamma += f ;  
		ngood2++;  
	      }
	    }
	  }
	  gamma /= ngood2;
	  f = sqrt( (gamma-1.)*(gamma+1.)/ ((gamma_old-1.)*(gamma_old+1.)) );
	  pv[i][j][k][U1] *= f ;	
	  pv[i][j][k][U2] *= f ;	
	  pv[i][j][k][U3] *= f ;	
	  i_ret = gamma_calc(pv[i][j][k],geom,&gamma);  /* Verify the new gamma just in case */
	}

	/* If interpolation cannot be done or the interpolation is bad, then set to old values */ 
	if( i_ret ) { 
	  VLOOP  pv[i][j][k][l] = pold[l]; 
	  gamma_calc(pv[i][j][k],geom,&gamma);
	  fail( FAIL_INTERP_PRIM,0 ) ; 
	}

	fail( FAIL_USE_INTERP_PRIM, 0 ); /* need to do this after we are sure gamma can be calculated */

	if( (!i_ret) && check_gamma(gamma) ) {
	  fixup_gamma(pv[i][j][k],&gamma);
	}

	p_gamma[i][j][k] = gamma;

	if( check_Tmax(pv[i][j][k]) )  fixup_Tmax(pv[i][j][k]) ; 

      }
    }

  }

  TRACE_END;
  return;
}


/*******************************************************************************************/
/*******************************************************************************************
  fixup_interp_v(): 

  -- just like fixup_interp_prim() but only performs the interpolation on the velocity 
          components; 

  -- 
  -- figures out (w/ pflag[]) which stencil to use to interpolate bad point from neighbors;

  -- interpolation is done using a 3x3x3 cube of parent (interpolating) cells 
     centered about the child (interpolated) cell; 

  -- In order to describe the stencils, let us use the following index notation 
     that is defined relative to the child cell.  A give cell has index 
      (d_i + 3*d_j + 9*d_k)   where "d_[i,j,k]"  is the change in [i,j,k] 
     indices from the center cell.  The figure below shows slices  in the x1,x2
     plane at different x3.  The  "0" marks the child cell. 


  ^ x2           d_k=-1            d_k=0             d_k=1
  |            -7   -6   -5       2   3   4       11   12   13
  |           -10   -9   -8      -1   0   1        8    9   10       
  |           -13  -12  -11      -4  -3  -2        5    6    7
  ----> x1   


  -- by numbering the cells in this way, we can exploit the parity symmetry of this 
     interpolation scheme;

  -- note that this routine uses pflag[] in the ghost zones, so be sure to set 
     pflag[]'s boundary conditions

  -- ngood_min[] is set so that this routine performs an interpolation when there are 
     an appropriate number of good neighboring cells;
      -- 1D runs: we require there to be at least two good neighbors ;
      -- 2D runs: we require there to be 2 unique (nondegenerate) lines (i.e. criss-cross)
      -- 3D runs: we are more strict in order to ensure that interpolation "samples" 
                  all dimensions and require that there be at least 4 lines;

  -- Note that ngood_min[] represents each dimensions number of  lines including degeneracies;
         for example, in 2D each "diagonal" line is thrice degenerate plus there is a line 
         connect two ghosts of the center cell.

 *******************************************************************************************/

void fixup_interp_v( double ****pv )
{
  int i, j, k, l, ind1, i_ret;
  int ii, jj, kk, ngood, ngood2;
  int pf[27];
  int o=13;
  double pold[NP];
  static const int ngood_min[] = { 13, 9, 6, 4 };
  double gamma, f, gamma_old;
  struct of_geom *geom, *geom2 ;

  TRACE_BEG;

  if( !failure_exists ) {
    TRACE_END;
    return;
  }


  /* Fix the interior points first */
  LOOP {
    if( pflag[i][j][k] == PFLAG_INTERP_V ) { 
      icurr = i; jcurr = j; kcurr = k;
 
      /* Collect the local flags and make pf true if it is a good point */
      VLOOP pold[l] = pv[i][j][k][l]; 

      for(ii=-1; ii<=1; ii++) for(jj=-1; jj<=1; jj++) for(kk=-1; kk<=1; kk++) { 
	ind1 = ii + 3*(jj + 3*kk);
	pf[o+ind1] = !pflag[i+ii][j+jj][k+kk];
      }

      /* Cancel mirror cells, and count the number of good pairs */
      ngood = 0;
      for(ii=1; ii<=13 ; ii++) { 
	if( pf[o+ii] && pf[o-ii]  ) { ngood++ ;                }
	else                        { pf[o+ii] = pf[o-ii] = 0; }
      }
      
      /* Just limit gamma if there are too few good neighbors : */
      if( ngood < ngood_min[n_spatial_dims] ) { 
	fixup_gamma(pv[i][j][k],&(p_gamma[i][j][k]));
	fail( FAIL_INTERP_V,0 ) ; 
      }
      /* Only fix via interpolation if it is sufficiently covered */
      else{ 

	/* Average the good points : */ 
	VLOOP pv[i][j][k][l] = 0.;
	for(ii=-1; ii<=1; ii++) for(jj=-1; jj<=1; jj++) for(kk=-1; kk<=1; kk++) { 
	  if( pf[o + ii + 3*(jj + 3*kk)] ) { 
	    VLOOP  pv[i][j][k][l] += pv[i+ii][j+jj][k+kk][l];
	  }
	}
	ngood *= 2;
	VLOOP  pv[i][j][k][l]  /=   ngood; 

	/* Make sure that we can calculate a valid lorentz factor : */
	get_geometry(i,j,k,CENT,ncurr,geom) ; 
	i_ret = gamma_calc(pv[i][j][k],geom,&gamma);

	/* Interpolate over gamma if fixed gamma is bad and if there is a valid stencil: */
	/* note that we cannot rescale v^i if (i_ret!=0)  since v^i may be bad */
	if( check_gamma(gamma) ) { 
	  gamma_old = gamma;
	  gamma = 0.;
	  ngood2 = 0;
	  for(ii=-1; ii<=1; ii++) for(jj=-1; jj<=1; jj++) for(kk=-1; kk<=1; kk++) { 
	    if( pf[o + ii + 3*(jj + 3*kk)] ) { 
	      get_geometry(i+ii,j+jj,k+kk,CENT,ncurr,geom2) ; 
	      i_ret = gamma_calc(pv[i+ii][j+jj][k+kk],geom2,&f);  
	      /* only use good gamma's */
	      if( !i_ret ) {  
		gamma += f ;  
		ngood2++;  
	      }
	    }
	  }
	  gamma /= ngood2;
	  f = sqrt( (gamma-1.)*(gamma+1.)/ ((gamma_old-1.)*(gamma_old+1.)) );
	  pv[i][j][k][U1] *= f ;	
	  pv[i][j][k][U2] *= f ;	
	  pv[i][j][k][U3] *= f ;	
	  i_ret = gamma_calc(pv[i][j][k],geom,&gamma);  /* Verify the new gamma just in case */
	}
      
	/* If interpolation cannot be done or the interpolation is bad, then try other ways: */ 
	if( i_ret || check_gamma(gamma) ) { 
	  VLOOP  pv[i][j][k][l] = pold[l];
	  fixup_gamma(pv[i][j][k],&(p_gamma[i][j][k]));
	  fail( FAIL_INTERP_V,0 ) ; 
	}
	else { 
	  fail( FAIL_USE_INTERP_V, 0 );
	  p_gamma[i][j][k] = gamma;
	}
      }
    }
  }

  TRACE_END;
  return;
}


#undef VLOOP 
#undef FLOOP 

#if( USE_SIMPLE_EOS )
/***********************************************************************
   set_Katm():

       -- sets the EOS constant used for Font's fix. 

       -- see utoprim_1dfix1.c and utoprim_1dvsq2fix1.c  for more
           information. 

       -- uses the initial floor values of rho/u determined by fixup1zone()

       -- we assume here that Constant X1,r is independent of theta,X2

***********************************************************************/
void set_Katm( void )
{
  int i, j, k, l, G_type ;
  double prim[NP], G_tmp, gamma;
  int get_G_ATM( double *G_tmp );
  struct of_coord *coords;

  double x[NDIM], xp[NDIM], r, c_s_sq, Omega_K;

  G_type = get_G_ATM( &G_tmp );

  fflush(stdout);
  fprintf(stdout,"G_tmp = %26.20e \n", G_tmp );
  fflush(stdout);

  /* Isentropic : */
  if( G_type ) { 
    j = N2S;  k = N3S;
    for( i = N1S ; i < N1E; i++ ) { 

      PLOOP prim[l] = 0.;
      prim[RHO] = prim[UU] = -1.;

      get_coord(i,j,k,CENT,ncurr,coords);      
      fixup1zone(prim, &gamma,coords);
      Katm[i] = (gam - 1.) * prim[UU] / pow( prim[RHO], G_tmp ) ;
    
      fflush(stdout);
      fprintf(stdout,"Katm[%d] = %26.20e \n", i, Katm[i] );
      fflush(stdout);
    }
  } 
  /* Isothermal : */
  else { 

    /* Extra factor is to compensate for the fact that MM08 use H = Gaussian Std. Dev. instead of scaleheight;
       See Noble, Krolik & Hawley 2010 */
    const double h_o_r = 0.1 * sqrt(2./M_PI);  

    fprintf(stdout,"set_Katm():  using  h_o_r = %26.16e \n", h_o_r); fflush(stdout);

    for( i = N1S ; i < N1E; i++ ) { 
      get_coord(i,N2S,N3S,CENT,ncurr,coords);      
      Katm[i]      = target_temperature( coords->x[RR], h_o_r );
      nu_visc[i]   = alpha_visc * Katm[i] / omega_circular_equatorial(coords->x[RR]);
    }
  }



  return;
}
#endif




