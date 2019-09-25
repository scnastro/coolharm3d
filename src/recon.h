
/*************************************************************************************/
/******  HEADER  FOR  RECONSTRUCTION ROUTINES IN recon.c              ****************/
/*************************************************************************************/


/*********************************************************************************
  Choose the factor to be used in the PPM routine  
*********************************************************************************/
#define ONE_SIXTH  (0.16666666666666666667)
#define PPM_FACTOR (ONE_SIXTH)  /* finite volume form */
//#define PPM_FACTOR (0.125)    /* finite difference form */


/**************************************************************************************
  In general, one reconstructs in the nearest ghost zones for the CT method.  But, 
   if we are doing a 2D problem, we do not have to reconstruct in the degenerate 
   direction nor do we need to reconstruct in the other directions in that 
   degenerate direction's ghost zones.  
    -- do_recon_dir[idir] = 1  if we are to reconstruct  in direction "idir" 
**************************************************************************************/

#if( USE_FLUX_CT_PARA ) 
#  define BUF (NG)
#else
#  define BUF (1)
#endif

#define N1_R_S     (N1_R_S_glob)
#define N1_R_E     (N1_R_E_glob)
#define N2_R_S     (N2_R_S_glob)
#define N2_R_E     (N2_R_E_glob)
#define N3_R_S     (N3_R_S_glob)
#define N3_R_E     (N3_R_E_glob)
#define DEGENERATE (DEGENERATE_glob)
