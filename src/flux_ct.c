
#include "decs.h"

/********************************************************************************************************/
/********************************************************************************************************
  flux_ct(): 
  ----------

   -- the CT flux at (i,j,k)  needs the fluxes from the following cells: 
         i-1,j  ,k  
         i-1,j  ,k+1  
         i-1,j+1,k
         i  ,j-1,k  
         i  ,j-1,k+1  
         i  ,j  ,k-1  
         i  ,j  ,k  
         i  ,j  ,k+1  
         i  ,j+1,k-1  
         i  ,j+1,k
         i+1,j-1,k  
         i+1,j  ,k-1  
         i+1,j  ,k  
********************************************************************************************************/
/********************************************************************************************************/

void flux_ct( double *****F ) 
{
	int i,j,k,d ;

  TRACE_BEG;

	/* sign convention on EMF: d_t B = curl (EMF) */
	/* Toth approach: average */
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  emf[0][i][j][k] = 0.25*(F[i][j][k][1][B3] + F[i][j  ][k-1][1][B3]
			        - F[i][j][k][2][B2] - F[i][j-1][k  ][2][B2]) ;
	  emf[1][i][j][k] = 0.25*(F[i][j][k][2][B1] + F[i-1][j][k  ][2][B1]
			        - F[i][j][k][0][B3] - F[i  ][j][k-1][0][B3]) ;
	  emf[2][i][j][k] = 0.25*(F[i][j][k][0][B2] + F[i  ][j-1][k][0][B2]
			        - F[i][j][k][1][B1] - F[i-1][j  ][k][1][B1]) ;
	}

	/* zero fluxes */
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++)	for(k=N3S;k<=N3E+1;k++) {
	  F[i][j][k][0][B1] = 0. ;
	  F[i][j][k][0][B2] = 0. ;
	  F[i][j][k][0][B3] = 0. ;
	  F[i][j][k][1][B1] = 0. ;
	  F[i][j][k][1][B2] = 0. ;
	  F[i][j][k][1][B3] = 0. ;
	  F[i][j][k][2][B1] = 0. ;
	  F[i][j][k][2][B2] = 0. ;
	  F[i][j][k][2][B3] = 0. ;
	}

#if(0 && USE_MASK) 
    /* Set the fluxes of the components of the induction equations in the excision region to zero so that we can 
       update the magnetic field there without having any problems:  
        -- Need to probe ghost zones for excision state in order to correctly set the fluxes on the physical boundary.
    */
        for(i=N1S-1;i<=N1E+1;i++)
        for(j=N2S-1;j<=N2E+1;j++) 
        for(k=N3S-1;k<=N3E+1;k++) 
	  if( evol_mask[ncurr][i][j][k] == MASK_EXCISED ) {
	    /* Lower xp1 face: */
	    emf[2][i  ][j  ][k  ] = 0.;
	    emf[2][i  ][j+1][k  ] = 0.;
	    emf[1][i  ][j  ][k  ] = 0.;
	    emf[1][i  ][j  ][k+1] = 0.;
	    /* Upper xp1 face: */
	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[1][i+1][j  ][k  ] = 0.;
	    emf[1][i+1][j  ][k+1] = 0.;
	    /* Lower xp2 face: */
//	    emf[2][i  ][j  ][k  ] = 0.;
//	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[0][i  ][j  ][k  ] = 0.; 
	    emf[0][i  ][j  ][k+1] = 0.;
	    /* Upper xp2 face: */
//	    emf[2][i  ][j+1][k  ] = 0.;
//	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[0][i  ][j+1][k  ] = 0.; 
	    emf[0][i  ][j+1][k+1] = 0.;
	    /* Lower xp3 face: */
//	    emf[1][i  ][j  ][k  ] = 0.;
//	    emf[1][i+1][j  ][k  ] = 0.;
//	    emf[0][i  ][j  ][k  ] = 0.;
//	    emf[0][i  ][j+1][k  ] = 0.;
	    /* Upper xp3 face: */
//	    emf[1][i  ][j  ][k+1] = 0.;
//	    emf[1][i+1][j  ][k+1] = 0.;
//	    emf[0][i  ][j  ][k+1] = 0.;
//	    emf[0][i  ][j+1][k+1] = 0.;
	  }
#endif

	/* signs check on EMFs; just compare w/ defns above */

        for(i=N1S;i<=N1E+1;i++)
        for(j=N2S;j<=N2E  ;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][0][B2] +=  0.5*(emf[2][i][j][k] + emf[2][i  ][j+1][k  ]) ;
		F[i][j][k][0][B3] += -0.5*(emf[1][i][j][k] + emf[1][i  ][j  ][k+1]) ;
        }

        for(i=N1S;i<=N1E  ;i++)
        for(j=N2S;j<=N2E+1;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][1][B1] += -0.5*(emf[2][i][j][k] + emf[2][i+1][j  ][k  ]) ;
		F[i][j][k][1][B3] +=  0.5*(emf[0][i][j][k] + emf[0][i  ][j  ][k+1]) ;
        }

	for(i=N1S;i<=N1E  ;i++)
	for(j=N2S;j<=N2E  ;j++)
        for(k=N3S;k<=N3E+1;k++) 
	{
		F[i][j][k][2][B1] +=  0.5*(emf[1][i][j][k] + emf[1][i+1][j  ][k  ]) ;
		F[i][j][k][2][B2] += -0.5*(emf[0][i][j][k] + emf[0][i  ][j+1][k  ]) ;
	}

  TRACE_END;

  return;

}

/********************************************************************************************************/
/********************************************************************************************************
  flux_ct_para(): 
  ----------

   -- the parabolic  CT flux at (i,j,k)  needs the fluxes from the following cells: 
            i-3,j  ,k
            i-3,j+1,k
            i-3,j  ,k+1
            i-2,j  ,k
            i-2,j+1,k
            i-2,j  ,k+1
            i-1,j  ,k
            i-1,j+1,k
            i-1,j  ,k+1
            i  ,j-3,k
            i  ,j-3,k+1
            i  ,j-2,k
            i  ,j-2,k+1
            i  ,j-1,k
            i  ,j-1,k+1
            i  ,j  ,k  
            i  ,j  ,k-3
            i  ,j  ,k-2
            i  ,j  ,k-1
            i  ,j  ,k+1
            i  ,j  ,k+2
            i  ,j  ,k+3
            i  ,j+1,k-3
            i  ,j+1,k-2
            i  ,j+1,k-1
            i  ,j+1,k
            i  ,j+1,k+1
            i  ,j+1,k+2
            i  ,j+2,k
            i  ,j+2,k+1
            i  ,j+3,k
            i+1,j-3,k
            i+1,j-2,k
            i+1,j-1,k
            i+1,j  ,k-3
            i+1,j  ,k-2
            i+1,j  ,k-1
            i+1,j  ,k
            i+1,j  ,k+1
            i+1,j  ,k+2
            i+1,j+1,k
            i+1,j+2,k
            i+2,j  ,k
            i+2,j  ,k+1
            i+2,j+1,k
            i+3,j  ,k

********************************************************************************************************/
/********************************************************************************************************/
void flux_ct_para( double *****F )
{
	int i,j,k,d ;
	double fl1,fr1,fl2,fr2,dq;

	void para(double x1, double x2, double x3, double x4, double x5, 
		  double *lout, double *rout);
	
	double linear_recon(double y1,double y2,double y3) ;

	void comp_emf( double ***emf[NDIM-1], double *****F  ) ;

  TRACE_BEG;


	/* sign convention on EMF: d_t B = curl (EMF) */
	/* fluxes only defined  one cell beyond physical domain */
	/* emfs calculated on all physical cells plus one past the end (e.g. N1S to N1E+1) */

  /*****************************************************************************************
    We use 3rd-order reconstruction using PPM for the EMFs on the interior points. 
    Averaging (the old way) is used for the boundary points since neighboring subdomains (even
    the phi boundary condition) need to have a consistent way of calculating the shared
    EMFs since they are not communicated.  Linear (slope-limiting) interpolation is 
    used for the L-state part of the average of the first interior cells.  An interior cell
    is one that does NOT abut a boundary. 
  *****************************************************************************************/


  /* zero-out emf's since we average different reconstructions disjointly below: */
  FACE_LOOP ALL_LOOP {  emf[d][i][j][k] = 0.; }
	   
	/******************************************************************************
            EMF 0 
	*******************************************************************************/
	/* first term : */
	for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++) {
	  k = N3S; 
	  dq = linear_recon( F[i][j][k-1][1][B3], 
			     F[i][j][k  ][1][B3], 
			     F[i][j][k+1][1][B3]);
	  emf[0][i][j][k  ] += 0.25 *(F[i][j][k][1][B3] + F[i][j][k-1][1][B3] );  
	  emf[0][i][j][k+1] += 0.25*(F[i][j][k][1][B3] + 0.5*dq );

	  k = N3E; 
	  dq = linear_recon( F[i][j][k-1][1][B3], 
			     F[i][j][k  ][1][B3], 
			     F[i][j][k+1][1][B3]);
	  emf[0][i][j][k  ] += 0.25*(F[i][j][k][1][B3] - 0.5*dq );  
	  emf[0][i][j][k+1] += 0.25 *(F[i][j][k+1][1][B3] + F[i][j][k][1][B3] ); 
	}
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++) for(k=N3S+1;k<=N3E-1;k++) {
	  para(F[i][j][k-2][1][B3],
	       F[i][j][k-1][1][B3],
	       F[i][j][k  ][1][B3],
	       F[i][j][k+1][1][B3],
	       F[i][j][k+2][1][B3],
	       &fl1,&fr1);
	  emf[0][i][j][k  ] += 0.25*fl1;
	  emf[0][i][j][k+1] += 0.25*fr1;
	}
	/* second term : */
	for(i=N1S;i<=N1E+1;i++) for(k=N3S;k<=N3E+1;k++)  {
	  j = N2S; 
	  dq = linear_recon( F[i][j-1][k][2][B2], 
			     F[i][j  ][k][2][B2], 
			     F[i][j+1][k][2][B2]);
	  emf[0][i][j  ][k] -=  0.25*(F[i][j][k][2][B2] + F[i][j-1][k][2][B2]);
	  emf[0][i][j+1][k] -=  0.25*(F[i][j][k][2][B2] + 0.5*dq );

	  j = N2E; 
	  dq = linear_recon( F[i][j-1][k][2][B2], 
			     F[i][j  ][k][2][B2], 
			     F[i][j+1][k][2][B2]);
	  emf[0][i][j  ][k] -= 0.25*(F[i][j][k][2][B2] - 0.5*dq );  
	  emf[0][i][j+1][k] -= 0.25*(F[i][j+1][k][2][B2] + F[i][j][k][2][B2]);
	}
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S+1;j<=N2E-1;j++) for(k=N3S;k<=N3E+1;k++) {
	  para(F[i][j-2][k][2][B2],
	       F[i][j-1][k][2][B2],
	       F[i][j  ][k][2][B2],
	       F[i][j+1][k][2][B2],
	       F[i][j+2][k][2][B2],
	       &fl2,&fr2);
	  emf[0][i][j  ][k] -= 0.25*fl2;
	  emf[0][i][j+1][k] -= 0.25*fr2;
	}

	/******************************************************************************
            EMF 1 
	*******************************************************************************/
	/* first term : */
	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  i = N1S; 
	  dq = linear_recon( F[i-1][j][k][2][B1], 
			     F[i  ][j][k][2][B1], 
			     F[i+1][j][k][2][B1]);
	  emf[1][i  ][j][k] += 0.25*(F[i][j][k][2][B1] + F[i-1][j][k][2][B1]);
	  emf[1][i+1][j][k] += 0.25*(F[i][j][k][2][B1] + 0.5*dq );
	  
	  i = N1E; 
	  dq = linear_recon( F[i-1][j][k][2][B1], 
			     F[i  ][j][k][2][B1], 
			     F[i+1][j][k][2][B1]);
	  emf[1][i  ][j][k] += 0.25*(F[i][j][k][2][B1] - 0.5*dq );  
	  emf[1][i+1][j][k] += 0.25*(F[i+1][j][k][2][B1] + F[i][j][k][2][B1]);
	}
	for(i=N1S+1;i<=N1E-1;i++)  for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  para(F[i-2][j][k][2][B1],
	       F[i-1][j][k][2][B1],
	       F[i  ][j][k][2][B1],
	       F[i+1][j][k][2][B1],
	       F[i+2][j][k][2][B1],
	       &fl1,&fr1);
	  emf[1][i  ][j][k] += 0.25*fl1;
	  emf[1][i+1][j][k] += 0.25*fr1;
	}
	/* second term : */
	for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++) {
	  k = N3S; 
	  dq = linear_recon( F[i][j][k-1][0][B3], 
			     F[i][j][k  ][0][B3], 
			     F[i][j][k+1][0][B3]);
	  emf[1][i][j][k  ] -= 0.25*(F[i][j][k][0][B3] + F[i][j][k-1][0][B3]);
	  emf[1][i][j][k+1] -= 0.25*(F[i][j][k][0][B3] + 0.5*dq );

	  k = N3E; 
	  dq = linear_recon( F[i][j][k-1][0][B3], 
			     F[i][j][k  ][0][B3], 
			     F[i][j][k+1][0][B3]);
	  emf[1][i][j][k  ] -= 0.25*(F[i][j][k][0][B3] - 0.5*dq );  
	  emf[1][i][j][k+1] -= 0.25*(F[i][j][k+1][0][B3] + F[i][j][k][0][B3]);
	}
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++) for(k=N3S+1;k<=N3E-1;k++) {
	  para(F[i][j][k-2][0][B3],
	       F[i][j][k-1][0][B3],
	       F[i][j][k  ][0][B3],
	       F[i][j][k+1][0][B3],
	       F[i][j][k+2][0][B3],
	       &fl2,&fr2);
	  emf[1][i][j][k  ] -= 0.25*fl2;
	  emf[1][i][j][k+1] -= 0.25*fr2;
	}

	/******************************************************************************
            EMF 2
	*******************************************************************************/
	/* first term : */
	for(i=N1S;i<=N1E+1;i++) for(k=N3S;k<=N3E+1;k++)  {
	  j = N2S; 
	  dq = linear_recon( F[i][j-1][k][0][B2], 
			     F[i][j  ][k][0][B2], 
			     F[i][j+1][k][0][B2]);
	  emf[2][i][j  ][k] += 0.25*(F[i][j][k][0][B2] + F[i][j-1][k][0][B2]);
	  emf[2][i][j+1][k] += 0.25*(F[i][j][k][0][B2] + 0.5*dq );
	  
	  j = N2E; 
	  dq = linear_recon( F[i][j-1][k][0][B2], 
			     F[i][j  ][k][0][B2], 
			     F[i][j+1][k][0][B2]);
	  emf[2][i][j  ][k] += 0.25*(F[i][j][k][0][B2] - 0.5*dq );  
	  emf[2][i][j+1][k] += 0.25*(F[i][j+1][k][0][B2] + F[i][j][k][0][B2]);
	}
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S+1;j<=N2E-1;j++) for(k=N3S;k<=N3E+1;k++) {
	  para(F[i][j-2][k][0][B2],
	       F[i][j-1][k][0][B2],
	       F[i][j  ][k][0][B2],
	       F[i][j+1][k][0][B2],
	       F[i][j+2][k][0][B2],
	       &fl1,&fr1);
	  emf[2][i][j  ][k] += 0.25*fl1;
	  emf[2][i][j+1][k] += 0.25*fr1;
	}
	/* second term : */
	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  i = N1S; 
	  dq = linear_recon( F[i-1][j][k][1][B1], 
			     F[i  ][j][k][1][B1], 
			     F[i+1][j][k][1][B1]);
	  emf[2][i  ][j][k] -= 0.25*(F[i][j][k][1][B1] + F[i-1][j][k][1][B1]);
	  emf[2][i+1][j][k] -= 0.25*(F[i][j][k][1][B1] + 0.5*dq );

	  i = N1E; 
	  dq = linear_recon( F[i-1][j][k][1][B1], 
			     F[i  ][j][k][1][B1], 
			     F[i+1][j][k][1][B1]);
	  emf[2][i  ][j][k] -= 0.25*(F[i][j][k][1][B1] - 0.5*dq );  
	  emf[2][i+1][j][k] -= 0.25*(F[i+1][j][k][1][B1] + F[i][j][k][1][B1]);
	}
	for(i=N1S+1;i<=N1E-1;i++)  for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  para(F[i-2][j][k][1][B1],
	       F[i-1][j][k][1][B1],
	       F[i  ][j][k][1][B1],
	       F[i+1][j][k][1][B1],
	       F[i+2][j][k][1][B1],
	       &fl2,&fr2);
	  emf[2][i  ][j][k] -= 0.25*fl2;
	  emf[2][i+1][j][k] -= 0.25*fr2;
	}

	//	comp_emf(  emf , F  ) ;

	/******************************************************************************
	   zero fluxes
	*******************************************************************************/
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++)	for(k=N3S;k<=N3E+1;k++) {
	  F[i][j][k][0][B1] = 0. ;
	  F[i][j][k][0][B2] = 0. ;
	  F[i][j][k][0][B3] = 0. ;
	  F[i][j][k][1][B1] = 0. ;
	  F[i][j][k][1][B2] = 0. ;
	  F[i][j][k][1][B3] = 0. ;
	  F[i][j][k][2][B1] = 0. ;
	  F[i][j][k][2][B2] = 0. ;
	  F[i][j][k][2][B3] = 0. ;
	}

#if(0 && USE_MASK) 
    /* Set the fluxes of the components of the induction equations in the excision region to zero so that we can 
       update the magnetic field there without having any problems:  
        -- Need to probe ghost zones for excision state in order to correctly set the fluxes on the physical boundary.
    */
        for(i=N1S-1;i<=N1E+1;i++)
        for(j=N2S-1;j<=N2E+1;j++) 
        for(k=N3S-1;k<=N3E+1;k++) 
	  if( evol_mask[ncurr][i][j][k] == MASK_EXCISED ) {
	    /* Lower xp1 face: */
	    emf[2][i  ][j  ][k  ] = 0.;
	    emf[2][i  ][j+1][k  ] = 0.;
	    emf[1][i  ][j  ][k  ] = 0.;
	    emf[1][i  ][j  ][k+1] = 0.;
	    /* Upper xp1 face: */
	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[1][i+1][j  ][k  ] = 0.;
	    emf[1][i+1][j  ][k+1] = 0.;
	    /* Lower xp2 face: */
//	    emf[2][i  ][j  ][k  ] = 0.;
//	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[0][i  ][j  ][k  ] = 0.; 
	    emf[0][i  ][j  ][k+1] = 0.;
	    /* Upper xp2 face: */
//	    emf[2][i  ][j+1][k  ] = 0.;
//	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[0][i  ][j+1][k  ] = 0.; 
	    emf[0][i  ][j+1][k+1] = 0.;
	    /* Lower xp3 face: */
//	    emf[1][i  ][j  ][k  ] = 0.;
//	    emf[1][i+1][j  ][k  ] = 0.;
//	    emf[0][i  ][j  ][k  ] = 0.;
//	    emf[0][i  ][j+1][k  ] = 0.;
	    /* Upper xp3 face: */
//	    emf[1][i  ][j  ][k+1] = 0.;
//	    emf[1][i+1][j  ][k+1] = 0.;
//	    emf[0][i  ][j  ][k+1] = 0.;
//	    emf[0][i  ][j+1][k+1] = 0.;
	  }
#endif

	/******************************************************************************
	   signs check on EMFs; just compare w/ defns above 
	*******************************************************************************/

        for(i=N1S;i<=N1E+1;i++)
        for(j=N2S;j<=N2E  ;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][0][B2] +=  0.5*(emf[2][i][j][k] + emf[2][i  ][j+1][k  ]) ;
		F[i][j][k][0][B3] += -0.5*(emf[1][i][j][k] + emf[1][i  ][j  ][k+1]) ;
        }

        for(i=N1S;i<=N1E  ;i++)
        for(j=N2S;j<=N2E+1;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][1][B1] += -0.5*(emf[2][i][j][k] + emf[2][i+1][j  ][k  ]) ;
		F[i][j][k][1][B3] +=  0.5*(emf[0][i][j][k] + emf[0][i  ][j  ][k+1]) ;
        }

	for(i=N1S;i<=N1E  ;i++)
	for(j=N2S;j<=N2E  ;j++)
        for(k=N3S;k<=N3E+1;k++) 
	{
		F[i][j][k][2][B1] +=  0.5*(emf[1][i][j][k] + emf[1][i+1][j  ][k  ]) ;
		F[i][j][k][2][B2] += -0.5*(emf[0][i][j][k] + emf[0][i  ][j+1][k  ]) ;
	}

  TRACE_END;

  return;

}



/*****************************************************************************************/
/*****************************************************************************************
  flux_ct_para_fast():
  --------------------
      -- faster version of flux_ct_para();
      -- using vectorized reconstruction methods for EMFs
      -- moved the factor of 0.25 to the flux calculations at the end 
         (saves many multiplications);

      -- We use 3rd-order reconstruction using PPM for the EMFs on the interior points. 
         Averaging (the old way) is used for the boundary points since neighboring 
         subdomains (even the phi boundary condition) need to have a consistent way 
         of calculating the shared EMFs since they are not communicated.  Linear 
         (slope-limiting) interpolation is used for the L-state part of the average 
         of the first interior cells.  An interior cell is one that does NOT abut a boundary. 

      --  sign convention on EMF: d_t B = curl (EMF) 
      --  fluxes only defined  one cell beyond physical domain 
      --  emfs calculated on all physical cells plus one past the end (e.g. N1S to N1E+1) 

*****************************************************************************************/
void flux_ct_para_fast( double *****F )
{
	int i,j,k,d ;
	double fl1,fr1,fl2,fr2,dq;

	void linear_recon_vect( double *y, double *Dq, int ibeg, int iend);
	void para_vect(double *y, double *dq, double *yL, double *yR , int ibeg, int iend);

	//	void comp_emf( double ***emf[NDIM-1] , double *****F  ) ;

  TRACE_BEG;


  /* zero-out emf's since we average different reconstructions disjointly below: */
  FACE_LOOP ALL_LOOP {  emf[d][i][j][k] = 0.; }
	   
	/******************************************************************************
            EMF 0   N3S-1  N3E+1
	*******************************************************************************/
	/* first term : */
	for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++) {
	  for(k=N3S-1;k<=N3E+1;k++) { p_vect[k] = F[i][j][k][1][B3]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N3S  ,N3E  );  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N3S+1,N3E-1); 
	    
	  k = N3S; 
	  emf[0][i][j][k  ] += (F[i][j][k][1][B3] + F[i][j][k-1][1][B3] );  
	  emf[0][i][j][k+1] += (F[i][j][k][1][B3] + 0.5*dq_vect[k] );

	  k = N3E; 
	  emf[0][i][j][k  ] += (F[i][j][k][1][B3] - 0.5*dq_vect[k] );  
	  emf[0][i][j][k+1] += (F[i][j][k+1][1][B3] + F[i][j][k][1][B3] ); 

	  for(k=N3S+1;k<=N3E-1;k++) { 
	    emf[0][i][j][k  ] += pL_vect[k];  
	    emf[0][i][j][k+1] += pR_vect[k]; 
	  }
	}
	/* second term : */
	for(i=N1S;i<=N1E+1;i++) for(k=N3S;k<=N3E+1;k++)  {
	  for(j=N2S-1;j<=N2E+1;j++) { p_vect[j] = F[i][j][k][2][B2]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N2S  ,N2E  );  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N2S+1,N2E-1); 
	    
	  j = N2S; 
	  emf[0][i][j  ][k] -=  (F[i][j][k][2][B2] + F[i][j-1][k][2][B2]);
	  emf[0][i][j+1][k] -=  (F[i][j][k][2][B2] + 0.5*dq_vect[j] );

	  j = N2E; 
	  emf[0][i][j  ][k] -= (F[i][j][k][2][B2] - 0.5*dq_vect[j] );  
	  emf[0][i][j+1][k] -= (F[i][j+1][k][2][B2] + F[i][j][k][2][B2]);

	  for(j=N2S+1;j<=N2E-1;j++) {
	    emf[0][i][j  ][k] -= pL_vect[j];
	    emf[0][i][j+1][k] -= pR_vect[j];
	  }
	}

	/******************************************************************************
            EMF 1 
	*******************************************************************************/
	/* first term : */
	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  for(i=N1S-1;i<=N1E+1;i++) { p_vect[i] = F[i][j][k][2][B1]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N1S  ,N1E  );  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N1S+1,N1E-1); 
	    
	  i = N1S; 
	  emf[1][i  ][j][k] += (F[i][j][k][2][B1] + F[i-1][j][k][2][B1]);
	  emf[1][i+1][j][k] += (F[i][j][k][2][B1] + 0.5*dq_vect[i] );

	  i = N1E; 
	  emf[1][i  ][j][k] += (F[i][j][k][2][B1] - 0.5*dq_vect[i] );  
	  emf[1][i+1][j][k] += (F[i+1][j][k][2][B1] + F[i][j][k][2][B1]);

	  for(i=N1S+1;i<=N1E-1;i++) {
	    emf[1][i  ][j][k] += pL_vect[i];
	    emf[1][i+1][j][k] += pR_vect[i];
	  }
	}
	/* second term : */
	for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++)  {
	  for(k=N3S-1;k<=N3E+1;k++)  { p_vect[k] = F[i][j][k][0][B3]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N3S  ,N3E  );  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N3S+1,N3E-1); 
	    
	  k = N3S; 
	  emf[1][i][j][k  ] -= (F[i][j][k][0][B3] + F[i][j][k-1][0][B3]);
	  emf[1][i][j][k+1] -= (F[i][j][k][0][B3] + 0.5*dq_vect[k] );

	  k = N3E; 
	  emf[1][i][j][k  ] -= (F[i][j][k][0][B3] - 0.5*dq_vect[k] );  
	  emf[1][i][j][k+1] -= (F[i][j][k+1][0][B3] + F[i][j][k][0][B3]);

	  for(k=N3S+1;k<=N3E-1;k++) {
	    emf[1][i][j][k  ] -= pL_vect[k];
	    emf[1][i][j][k+1] -= pR_vect[k];
	  }
	}

	/******************************************************************************
            EMF 2
	*******************************************************************************/
	/* first term : */
	for(i=N1S;i<=N1E+1;i++) for(k=N3S;k<=N3E+1;k++) {
	  for(j=N2S-1;j<=N2E+1;j++) { p_vect[j] = F[i][j][k][0][B2]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N2S  ,N2E  );  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N2S+1,N2E-1); 
	    
	  j = N2S; 
	  emf[2][i][j  ][k] += (F[i][j][k][0][B2] + F[i][j-1][k][0][B2]);
	  emf[2][i][j+1][k] += (F[i][j][k][0][B2] + 0.5*dq_vect[j] );
	  
	  j = N2E; 
	  emf[2][i][j  ][k] += (F[i][j][k][0][B2] - 0.5*dq_vect[j] );  
	  emf[2][i][j+1][k] += (F[i][j+1][k][0][B2] + F[i][j][k][0][B2]);

	  for(j=N2S+1;j<=N2E-1;j++) {
	    emf[2][i][j  ][k] += pL_vect[j];
	    emf[2][i][j+1][k] += pR_vect[j];
	  }
	}
	/* second term : */
	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  for(i=N1S-1;i<=N1E+1;i++) { p_vect[i] = F[i][j][k][1][B1]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N1S  ,N1E  );  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N1S+1,N1E-1); 
	    
	  i = N1S; 
	  emf[2][i  ][j][k] -= (F[i][j][k][1][B1] + F[i-1][j][k][1][B1]);
	  emf[2][i+1][j][k] -= (F[i][j][k][1][B1] + 0.5*dq_vect[i] );

	  i = N1E; 
	  emf[2][i  ][j][k] -= (F[i][j][k][1][B1] - 0.5*dq_vect[i] );  
	  emf[2][i+1][j][k] -= (F[i+1][j][k][1][B1] + F[i][j][k][1][B1]);

	  for(i=N1S+1;i<=N1E-1;i++) {
	    emf[2][i  ][j][k] -= pL_vect[i];
	    emf[2][i+1][j][k] -= pR_vect[i];
	  }
	}

	//	comp_emf(  emf , F  ) ;

	/******************************************************************************
	   zero fluxes
	*******************************************************************************/
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++)	for(k=N3S;k<=N3E+1;k++) {
	  F[i][j][k][0][B1] = 0. ;
	  F[i][j][k][0][B2] = 0. ;
	  F[i][j][k][0][B3] = 0. ;
	  F[i][j][k][1][B1] = 0. ;
	  F[i][j][k][1][B2] = 0. ;
	  F[i][j][k][1][B3] = 0. ;
	  F[i][j][k][2][B1] = 0. ;
	  F[i][j][k][2][B2] = 0. ;
	  F[i][j][k][2][B3] = 0. ;
	}

#if(0 && USE_MASK) 
    /* Set the fluxes of the components of the induction equations in the excision region to zero so that we can 
       update the magnetic field there without having any problems:  
        -- Need to probe ghost zones for excision state in order to correctly set the fluxes on the physical boundary.
    */
        for(i=N1S-1;i<=N1E+1;i++)
        for(j=N2S-1;j<=N2E+1;j++) 
        for(k=N3S-1;k<=N3E+1;k++) 
	  if( evol_mask[ncurr][i][j][k] == MASK_EXCISED ) {
	    /* Lower xp1 face: */
	    emf[2][i  ][j  ][k  ] = 0.;
	    emf[2][i  ][j+1][k  ] = 0.;
	    emf[1][i  ][j  ][k  ] = 0.;
	    emf[1][i  ][j  ][k+1] = 0.;
	    /* Upper xp1 face: */
	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[1][i+1][j  ][k  ] = 0.;
	    emf[1][i+1][j  ][k+1] = 0.;
	    /* Lower xp2 face: */
//	    emf[2][i  ][j  ][k  ] = 0.;
//	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[0][i  ][j  ][k  ] = 0.; 
	    emf[0][i  ][j  ][k+1] = 0.;
	    /* Upper xp2 face: */
//	    emf[2][i  ][j+1][k  ] = 0.;
//	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[0][i  ][j+1][k  ] = 0.; 
	    emf[0][i  ][j+1][k+1] = 0.;
	    /* Lower xp3 face: */
//	    emf[1][i  ][j  ][k  ] = 0.;
//	    emf[1][i+1][j  ][k  ] = 0.;
//	    emf[0][i  ][j  ][k  ] = 0.;
//	    emf[0][i  ][j+1][k  ] = 0.;
	    /* Upper xp3 face: */
//	    emf[1][i  ][j  ][k+1] = 0.;
//	    emf[1][i+1][j  ][k+1] = 0.;
//	    emf[0][i  ][j  ][k+1] = 0.;
//	    emf[0][i  ][j+1][k+1] = 0.;
	  }
#endif

	/******************************************************************************
	   signs check on EMFs; just compare w/ defns above 
	*******************************************************************************/

        for(i=N1S;i<=N1E+1;i++)
        for(j=N2S;j<=N2E  ;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][0][B2] +=  0.125*(emf[2][i][j][k] + emf[2][i  ][j+1][k  ]) ;
		F[i][j][k][0][B3] += -0.125*(emf[1][i][j][k] + emf[1][i  ][j  ][k+1]) ;
        }

        for(i=N1S;i<=N1E  ;i++)
        for(j=N2S;j<=N2E+1;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][1][B1] += -0.125*(emf[2][i][j][k] + emf[2][i+1][j  ][k  ]) ;
		F[i][j][k][1][B3] +=  0.125*(emf[0][i][j][k] + emf[0][i  ][j  ][k+1]) ;
        }

	for(i=N1S;i<=N1E  ;i++)
	for(j=N2S;j<=N2E  ;j++)
        for(k=N3S;k<=N3E+1;k++) 
	{
		F[i][j][k][2][B1] +=  0.125*(emf[1][i][j][k] + emf[1][i+1][j  ][k  ]) ;
		F[i][j][k][2][B2] += -0.125*(emf[0][i][j][k] + emf[0][i  ][j+1][k  ]) ;
	}

  TRACE_END;

  return;

}

/*****************************************************************************************/
/*****************************************************************************************
  flux_ct_para_fast2():
  --------------------
      -- faster version of flux_ct_para();
      -- uses PPM everywhere for EMFS, unlike  flux_ct_para_fast();

      -- using vectorized reconstruction methods for EMFs
      -- moved the factor of 0.25 to the flux calculations at the end 
         (saves many multiplications);

      -- We use 3rd-order reconstruction using PPM for the EMFs on the interior points. 
         Averaging (the old way) is used for the boundary points since neighboring 
         subdomains (even the phi boundary condition) need to have a consistent way 
         of calculating the shared EMFs since they are not communicated.  Linear 
         (slope-limiting) interpolation is used for the L-state part of the average 
         of the first interior cells.  An interior cell is one that does NOT abut a boundary. 

      --  sign convention on EMF: d_t B = curl (EMF) 
      --  fluxes only defined  one cell beyond physical domain 
      --  emfs calculated on all physical cells plus one past the end (e.g. N1S to N1E+1) 

*****************************************************************************************/
void flux_ct_para_fast2( double *****F )
{
	int i,j,k,d ;
	double fl1,fr1,fl2,fr2,dq;

	void linear_recon_vect( double *y, double *Dq, int ibeg, int iend);
	void para_vect(double *y, double *dq, double *yL, double *yR , int ibeg, int iend);

	//	void comp_emf( double ***emf[NDIM-1] , double *****F  ) ;

  TRACE_BEG;


  /* zero-out emf's since we average different reconstructions disjointly below: */
  FACE_LOOP ALL_LOOP {  emf[d][i][j][k] = 0.; }

	   
	/******************************************************************************
            EMF 0   N3S-1  N3E+1
	*******************************************************************************/
	/* first term : */
	for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++) {
	  N3ALL_LOOP { p_vect[k] = F[i][j][k][1][B3]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N3S-2,N3E+2);  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N3S-1,N3E+1); 
	  for(k=N3S-1;k<=N3E+1;k++) { 
	    emf[0][i][j][k  ] += pL_vect[k];  
	    emf[0][i][j][k+1] += pR_vect[k]; 
	  }
	}
	/* second term : */
	for(i=N1S;i<=N1E+1;i++) for(k=N3S;k<=N3E+1;k++)  {
	  N2ALL_LOOP { p_vect[j] = F[i][j][k][2][B2]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N2S-2,N2E+2);  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N2S-1,N2E+1); 
	  for(j=N2S-1;j<=N2E+1;j++) {
	    emf[0][i][j  ][k] -= pL_vect[j];
	    emf[0][i][j+1][k] -= pR_vect[j];
	  }
	}

	/******************************************************************************
            EMF 1 
	*******************************************************************************/
	/* first term : */
	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  N1ALL_LOOP { p_vect[i] = F[i][j][k][2][B1]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N1S-2,N1E+2);  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N1S-1,N1E+1); 
	  for(i=N1S-1;i<=N1E+1;i++) {
	    emf[1][i  ][j][k] += pL_vect[i];
	    emf[1][i+1][j][k] += pR_vect[i];
	  }
	}
	/* second term : */
	for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++)  {
	  N3ALL_LOOP  { p_vect[k] = F[i][j][k][0][B3]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N3S-2,N3E+2);  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N3S-1,N3E+1); 
	  for(k=N3S-1;k<=N3E+1;k++) {
	    emf[1][i][j][k  ] -= pL_vect[k];
	    emf[1][i][j][k+1] -= pR_vect[k];
	  }
	}

	/******************************************************************************
            EMF 2
	*******************************************************************************/
	/* first term : */
	for(i=N1S;i<=N1E+1;i++) for(k=N3S;k<=N3E+1;k++) {
	  N2ALL_LOOP { p_vect[j] = F[i][j][k][0][B2]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N2S-2,N2E+2);  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N2S-1,N2E+1); 
	  for(j=N2S-1;j<=N2E+1;j++) {
	    emf[2][i][j  ][k] += pL_vect[j];
	    emf[2][i][j+1][k] += pR_vect[j];
	  }
	}
	/* second term : */
	for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
	  N1ALL_LOOP { p_vect[i] = F[i][j][k][1][B1]; }  /* load recon array */
	  linear_recon_vect(p_vect,dq_vect,                N1S-2,N1E+2);  
	  para_vect(        p_vect,dq_vect,pL_vect,pR_vect,N1S-1,N1E+1); 
	  for(i=N1S-1;i<=N1E+1;i++) {
	    emf[2][i  ][j][k] -= pL_vect[i];
	    emf[2][i+1][j][k] -= pR_vect[i];
	  }
	}

	//	comp_emf(  emf , F  ) ;

	/******************************************************************************
	   zero fluxes
	*******************************************************************************/
	for(i=N1S;i<=N1E+1;i++)	for(j=N2S;j<=N2E+1;j++)	for(k=N3S;k<=N3E+1;k++) {
	  F[i][j][k][0][B1] = 0. ;
	  F[i][j][k][0][B2] = 0. ;
	  F[i][j][k][0][B3] = 0. ;
	  F[i][j][k][1][B1] = 0. ;
	  F[i][j][k][1][B2] = 0. ;
	  F[i][j][k][1][B3] = 0. ;
	  F[i][j][k][2][B1] = 0. ;
	  F[i][j][k][2][B2] = 0. ;
	  F[i][j][k][2][B3] = 0. ;
	}

#if(0 && USE_MASK) 
    /* Set the fluxes of the components of the induction equations in the excision region to zero so that we can 
       update the magnetic field there without having any problems:  
        -- Need to probe ghost zones for excision state in order to correctly set the fluxes on the physical boundary.
    */
        for(i=N1S-1;i<=N1E+1;i++)
        for(j=N2S-1;j<=N2E+1;j++) 
        for(k=N3S-1;k<=N3E+1;k++) 
	  if( evol_mask[ncurr][i][j][k] == MASK_EXCISED ) {
	    /* Lower xp1 face: */
	    emf[2][i  ][j  ][k  ] = 0.;
	    emf[2][i  ][j+1][k  ] = 0.;
	    emf[1][i  ][j  ][k  ] = 0.;
	    emf[1][i  ][j  ][k+1] = 0.;
	    /* Upper xp1 face: */
	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[1][i+1][j  ][k  ] = 0.;
	    emf[1][i+1][j  ][k+1] = 0.;
	    /* Lower xp2 face: */
//	    emf[2][i  ][j  ][k  ] = 0.;
//	    emf[2][i+1][j  ][k  ] = 0.;
	    emf[0][i  ][j  ][k  ] = 0.; 
	    emf[0][i  ][j  ][k+1] = 0.;
	    /* Upper xp2 face: */
//	    emf[2][i  ][j+1][k  ] = 0.;
//	    emf[2][i+1][j+1][k  ] = 0.;
	    emf[0][i  ][j+1][k  ] = 0.; 
	    emf[0][i  ][j+1][k+1] = 0.;
	    /* Lower xp3 face: */
//	    emf[1][i  ][j  ][k  ] = 0.;
//	    emf[1][i+1][j  ][k  ] = 0.;
//	    emf[0][i  ][j  ][k  ] = 0.;
//	    emf[0][i  ][j+1][k  ] = 0.;
	    /* Upper xp3 face: */
//	    emf[1][i  ][j  ][k+1] = 0.;
//	    emf[1][i+1][j  ][k+1] = 0.;
//	    emf[0][i  ][j  ][k+1] = 0.;
//	    emf[0][i  ][j+1][k+1] = 0.;
	  }
#endif

	/******************************************************************************
	   signs check on EMFs; just compare w/ defns above 
	*******************************************************************************/

        for(i=N1S;i<=N1E+1;i++)
        for(j=N2S;j<=N2E  ;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][0][B2] +=  0.125*(emf[2][i][j][k] + emf[2][i  ][j+1][k  ]) ;
		F[i][j][k][0][B3] += -0.125*(emf[1][i][j][k] + emf[1][i  ][j  ][k+1]) ;
        }

        for(i=N1S;i<=N1E  ;i++)
        for(j=N2S;j<=N2E+1;j++) 
        for(k=N3S;k<=N3E  ;k++) 
	{
                F[i][j][k][1][B1] += -0.125*(emf[2][i][j][k] + emf[2][i+1][j  ][k  ]) ;
		F[i][j][k][1][B3] +=  0.125*(emf[0][i][j][k] + emf[0][i  ][j  ][k+1]) ;
        }

	for(i=N1S;i<=N1E  ;i++)
	for(j=N2S;j<=N2E  ;j++)
        for(k=N3S;k<=N3E+1;k++) 
	{
		F[i][j][k][2][B1] +=  0.125*(emf[1][i][j][k] + emf[1][i+1][j  ][k  ]) ;
		F[i][j][k][2][B2] += -0.125*(emf[0][i][j][k] + emf[0][i  ][j+1][k  ]) ;
	}

  TRACE_END;

  return;

}

void comp_emf( double ***emf[NDIM-1], double *****F )
{
  int i,j,k,l;
  double f1,f2,rd,rderr;

  TRACE_BEG;

  rderr = 0.5;
  for(i=N1S;i<=N1E+1;i++) for(j=N2S;j<=N2E+1;j++) for(k=N3S;k<=N3E+1;k++) {
    l = 0;  f1 = emf[l][i][j][k]; 
    f2 =  0.25*(F[i][j][k][1][B3] + F[i][j  ][k-1][1][B3]
		- F[i][j][k][2][B2] - F[i][j-1][k  ][2][B2]) ;
    if( fabs(f1) > 1.e-30 && fabs(f2) > 1.e-30 ) { 
      rd = fabs( 2. * (f2 - f1) / (fabs(f1)+fabs(f2)) );
      if( rd > rderr ) { 
	fprintf(stdout,"emfrd %d %d %d %d %26.16e %26.16e %26.16e \n",l,i,j,k,f1,f2,rd); fflush(stdout);
      }
    }

    l = 1;  f1 = emf[l][i][j][k]; 
    f2 = 0.25*(F[i][j][k][2][B1] + F[i-1][j][k  ][2][B1]
	       - F[i][j][k][0][B3] - F[i  ][j][k-1][0][B3]) ;
    if( fabs(f1) > 1.e-30 && fabs(f2) > 1.e-30 ) { 
      rd = fabs( 2. * (f2 - f1) / (fabs(f1)+fabs(f2)) );
      if( rd > rderr ) { 
	fprintf(stdout,"emfrd %d %d %d %d %26.16e %26.16e %26.16e \n",l,i,j,k,f1,f2,rd); fflush(stdout);
      }
    }

    l = 2;  f1 = emf[l][i][j][k]; 
    f2 = 0.25*(F[i][j][k][0][B2] + F[i  ][j-1][k][0][B2]
	       - F[i][j][k][1][B1] - F[i-1][j  ][k][1][B1]) ;
    if( fabs(f1) > 1.e-30 && fabs(f2) > 1.e-30 ) { 
      rd = fabs( 2. * (f2 - f1) / (fabs(f1)+fabs(f2)) );
      if( rd > rderr ) { 
	fprintf(stdout,"emfrd %d %d %d %d %26.16e %26.16e %26.16e \n",l,i,j,k,f1,f2,rd); fflush(stdout);
      }
    }
  }

  TRACE_END;

  return;
}

/*****************************************************************************************/
/*****************************************************************************************
  curl_of_A():
  --------------------
      -- calculates the curl of a vector potential to get the magnetic field components;
      -- uses the finite difference stencil consistent with FluxCT's divergence operator, 
         so that the resultant magnetic field will have zero divergence (w.r.t that 
         divergence stencil); 

*****************************************************************************************/
void curl_of_A( double ****A , double ****prim ) 
{
  int i,j,k,l,p;
  double f0,b1,b2,b3;
  struct of_geom *geom;

  double f1 = 0.25*invdx[1];
  double f2 = 0.25*invdx[2];
  double f3 = 0.25*invdx[3];
  

  N1_LOOP  N2_LOOP   N3_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom); 


    /* B^1 = \partial_2 A_3 - \partial_3 A_2  (averaged over 4 edges) */
    b1 =  ( 
	   f2*(
	       A[i  ][j+1][k  ][B3] - A[i  ][j  ][k  ][B3] + 
	       A[i+1][j+1][k  ][B3] - A[i+1][j  ][k  ][B3] + 
	       A[i  ][j+1][k+1][B3] - A[i  ][j  ][k+1][B3] + 
	       A[i+1][j+1][k+1][B3] - A[i+1][j  ][k+1][B3] 
	       ) - 
	   f3*(
	       A[i  ][j  ][k+1][B2] - A[i  ][j  ][k  ][B2] + 
	       A[i+1][j  ][k+1][B2] - A[i+1][j  ][k  ][B2] + 
	       A[i  ][j+1][k+1][B2] - A[i  ][j+1][k  ][B2] + 
	       A[i+1][j+1][k+1][B2] - A[i+1][j+1][k  ][B2] 
	       )
	    );

    /* B^2 = \partial_3 A_1 - \partial_1 A_3  (averaged over 4 edges) */
    b2 =  ( 
	   f3*(
	       A[i  ][j  ][k+1][B1] - A[i  ][j  ][k  ][B1] + 
	       A[i+1][j  ][k+1][B1] - A[i+1][j  ][k  ][B1] + 
	       A[i  ][j+1][k+1][B1] - A[i  ][j+1][k  ][B1] + 
	       A[i+1][j+1][k+1][B1] - A[i+1][j+1][k  ][B1] 
	       ) - 
	   f1*(
	       A[i+1][j  ][k  ][B3] - A[i  ][j  ][k  ][B3] + 
	       A[i+1][j+1][k  ][B3] - A[i  ][j+1][k  ][B3] + 
	       A[i+1][j  ][k+1][B3] - A[i  ][j  ][k+1][B3] + 
	       A[i+1][j+1][k+1][B3] - A[i  ][j+1][k+1][B3] 
	       )
	    );
    
    /* B^3 = \partial_1 A_2 - \partial_2 A_1  (averaged over 4 edges) */
    b3 =  ( 
	   f1*(
	       A[i+1][j  ][k  ][B2] - A[i  ][j  ][k  ][B2] + 
	       A[i+1][j  ][k+1][B2] - A[i  ][j  ][k+1][B2] + 
	       A[i+1][j+1][k  ][B2] - A[i  ][j+1][k  ][B2] + 
	       A[i+1][j+1][k+1][B2] - A[i  ][j+1][k+1][B2] 
	       ) - 
	   f2*(
	       A[i  ][j+1][k  ][B1] - A[i  ][j  ][k  ][B1] + 
	       A[i+1][j+1][k  ][B1] - A[i+1][j  ][k  ][B1] + 
	       A[i  ][j+1][k+1][B1] - A[i  ][j  ][k+1][B1] + 
	       A[i+1][j+1][k+1][B1] - A[i+1][j  ][k+1][B1] 
	       )
	    );

    f0 = 1./geom->g; 
    prim[i][j][k][B1] = b1 * f0;
    prim[i][j][k][B2] = b2 * f0;
    prim[i][j][k][B3] = b3 * f0;
    
  }

  return;
}
