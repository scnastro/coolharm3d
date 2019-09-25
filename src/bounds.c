
#include "decs.h"

#if( USEMPI )
extern void bounds_mpi( double ****prim_arr );
extern void bounds_pflag_mpi( int ***pflag_arr );
extern void bounds_mpi_fast( double ****prim_arr );
extern void bounds_pflag_mpi_fast( int ***pflag_arr );
#endif

static void inflow_check(double *pr, int ii, int jj, int kk, int dir, int type );
static void inflow_check2(double *pr, int ii, int jj, int kk, int dir, int type );
static  void limit_gamma_in_horizon( double ****prim  );
static  void set_cutout_boundary( double ****prim );
static  void set_cutout_boundary2(double ****prim );
static  void set_cutout_boundary3(double ****prim );
static  void set_cutout_boundary4(double ****prim );
static  void cutout_cons_interp(  double ****prim , int type );

/*********************************************************************************/
/*********************************************************************************
   bounds():
   ---------

     -- boundary condition loops done in order of slowest index to fastest, 
        that's why "GLOOP" is in different places depending on the face. 

*********************************************************************************/
/*********************************************************************************/

void bounds(double ****prim_arr, int skip_mpi )
{
  unsigned int i,j,k,l,g,j2,k2;
  double geomfactor;
  struct of_geom *geom;


  TRACE_BEG;

  /********************************************************************************
    Scale B^i by gdet in real cells so that ghosts are set with geometrized B-field.
     (could optimize this more so that we only geometrize the boundary values, but 
      it is complicated for all grid geometries). 
  ********************************************************************************/
#if( RESCALE_B )
  LOOP  { 
    get_geometry(i,j,k,CENT,ncurr,geom);
    geomfactor = geom->g;
    BLOOP { prim_arr[i][j][k][l] *= geomfactor; }
  }
#endif

  if( (!skip_mpi) || boundary_phys_pflag ) { 

    /********************************************************************************
     CARTESIAN + OUTFLOW ON ALL SIDES
    ********************************************************************************/
#if( BC_TYPE_CHOICE  == BC_CARTESIAN_OUTFLOW )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP  {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP  {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP  {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP  {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face  */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP    {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3S][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP    {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3S] ;
	}
      }    
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP    {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3E][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP    {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3E] ;
	}
      }    
    }

#endif


    /********************************************************************************
     CARTESIAN + PERIODIC ON ALL SIDES
    ********************************************************************************/
#if( BC_TYPE_CHOICE  == BC_CARTESIAN_PERIODIC )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1E-g][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP {
	  pflag[N1S-1-g][j][k] = pflag[N1E-g][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1S+g][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP {
	  pflag[N1E+1+g][j][k] = pflag[N1S+g][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2E-g][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP {
	  pflag[i][N2S-1-g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2S+g][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP {
	  pflag[i][N2E+1+g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }


#endif

    /********************************************************************************
     CARTESIAN + REFLECTION BOUNDARIES ON ALL SIDES
    ********************************************************************************/
#if( BC_TYPE_CHOICE  == BC_CARTESIAN_REFLECT )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S+g][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	prim_arr[N1S-1-g][j][k][U1] *= -1 ;
	prim_arr[N1S-1-g][j][k][B1] *= -1 ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S+g][j][k] ;
	}
      }
    }

    /* 1-lower face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E-g][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	prim_arr[N1E+1+g][j][k][U1] *= -1 ;
	prim_arr[N1E+1+g][j][k][B1] *= -1 ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E-g][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP  {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
      }
      N1ALL_LOOP      GLOOP   N3_LOOP  {
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP   N3_LOOP  {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP   N3_LOOP  {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
      }
      N1ALL_LOOP      GLOOP   N3_LOOP  {
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP   N3_LOOP  {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	prim_arr[i][j][N3S-1-g][U3] *= -1 ;
	prim_arr[i][j][N3S-1-g][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3S+g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	prim_arr[i][j][N3E+1+g][U3] *= -1 ;
	prim_arr[i][j][N3E+1+g][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3E-g] ;
	}
      }
    }


#endif


    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, constant IN x2, PERIODIC IN x3
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, AXI-SYM IN x2, PERIODIC IN x3
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW2 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, REFLECTIVE AXI-SYM IN x2, PERIODIC IN x3
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW3 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, REFLECTIVE AXI-SYM IN x2, PERIODIC IN x3
        -- no inflow_check();
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW4 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif


    /**************************************************************************************
     SPHERICAL + OUTFLOW (without inflow_check() ) IN RADIAL x1 DIRECTION, AXI-SYM IN x2, PERIODIC IN x3
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW5 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, REFLECTIVE AXI-SYM IN x2, PERIODIC IN x3
        -- no inflow_check() at inner radial boundary, but inflow_check() at outer radial boundary;
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW6 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif


    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, REFLECTIVE AXI-SYM IN x2, PERIODIC IN x3
        -- no inflow_check() at inner radial boundary, but inflow_check2() at outer radial boundary;
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW7 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check2(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif


    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, COPIES x2, PERIODIC IN x3
        -- no inflow_check() at inner radial boundary, but inflow_check() at outer radial boundary;
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW8 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

#if(!EQUATORIAL_RUN) 
    why-are-we-using-this?
#endif

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif


    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, COPIES x2, PERIODIC IN x3
        -- no inflow_check() at inner radial boundary, but inflow_check() at outer radial boundary;
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW9 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

#if(!EQUATORIAL_RUN) 
    why-are-we-using-this?
#endif

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif


    /*******************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, REFLECTION ALONG CUTOUT IN x2, PERIODIC IN x3
    ********************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_CUTOUT )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      //    limit_gamma_in_horizon( prim_arr ); 
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face (set pflag[] in set_cutout_boundary() call below ) */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
    }

    /* 2-upper face (set pflag[] in set_cutout_boundary() call below ) */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /*******************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, REFLECTION ALONG CUTOUT IN x2, PERIODIC IN x3 w/ regularization
    ********************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_CUTOUT2 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	j = N2S;   j2 = N2S-1-g; 
#if(RESCALE_REGULARIZE==2)	
	regularize_prim( prim_arr[i][j][k], regularize_gf[i][j][k][CENT], dx_dxp_gf[i][j][k][CENT] );
#endif
	PLOOP prim_arr[i][j2][k][l] = prim_arr[i][j][k][l] ;
	prim_arr[i][j2][k][U2] *= -1 ;	prim_arr[i][j2][k][B2] *= -1 ;
#if(RESCALE_REGULARIZE==2)	
	unregularize_prim( prim_arr[i][j ][k], unregularize_gf[i][j ][k][CENT], dxp_dx_gf[i][j ][k][CENT] );
	unregularize_prim( prim_arr[i][j2][k], unregularize_gf[i][j2][k][CENT], dxp_dx_gf[i][j2][k][CENT] );
#endif
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	j = N2E;   j2 = N2E+1+g; 
#if(RESCALE_REGULARIZE==2)	
	regularize_prim( prim_arr[i][j][k], regularize_gf[i][j][k][CENT], dx_dxp_gf[i][j][k][CENT] );
#endif
	PLOOP prim_arr[i][j2][k][l] = prim_arr[i][j][k][l] ;
	prim_arr[i][j2][k][U2] *= -1 ;	prim_arr[i][j2][k][B2] *= -1 ;
#if(RESCALE_REGULARIZE==2)	
	unregularize_prim( prim_arr[i][j ][k], unregularize_gf[i][j ][k][CENT], dxp_dx_gf[i][j ][k][CENT] );
	unregularize_prim( prim_arr[i][j2][k], unregularize_gf[i][j2][k][CENT], dxp_dx_gf[i][j2][k][CENT] );
#endif
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /*******************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, NO REFLECTION ALONG CUTOUT IN x2, PERIODIC IN x3 w/ regularization
    ********************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_CUTOUT3 )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	j = N2S;   j2 = N2S-1-g; 
#if(RESCALE_REGULARIZE==2)	
	regularize_prim( prim_arr[i][j][k], regularize_gf[i][j][k][CENT], dx_dxp_gf[i][j][k][CENT] );
#endif
	PLOOP prim_arr[i][j2][k][l] = prim_arr[i][j][k][l] ;
#if(RESCALE_REGULARIZE==2)	
	unregularize_prim( prim_arr[i][j ][k], unregularize_gf[i][j ][k][CENT], dxp_dx_gf[i][j ][k][CENT] );
	unregularize_prim( prim_arr[i][j2][k], unregularize_gf[i][j2][k][CENT], dxp_dx_gf[i][j2][k][CENT] );
#endif
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	j = N2E;   j2 = N2E+1+g; 
#if(RESCALE_REGULARIZE==2)	
	regularize_prim( prim_arr[i][j][k], regularize_gf[i][j][k][CENT], dx_dxp_gf[i][j][k][CENT] );
#endif
	PLOOP prim_arr[i][j2][k][l] = prim_arr[i][j][k][l] ;
#if(RESCALE_REGULARIZE==2)	
	unregularize_prim( prim_arr[i][j ][k], unregularize_gf[i][j ][k][CENT], dxp_dx_gf[i][j ][k][CENT] );
	unregularize_prim( prim_arr[i][j2][k], unregularize_gf[i][j2][k][CENT], dxp_dx_gf[i][j2][k][CENT] );
#endif
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL + OUTFLOW IN RADIAL x1 DIRECTION, AXI-SYM IN x2 (interp. U[]), PERIODIC IN x3
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_CONSINTERP )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      cutout_cons_interp( prim_arr , BCDN );    
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      cutout_cons_interp( prim_arr , BCUP );    
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif


    /********************************************************************************
     CYLINDRICAL + OUTFLOW IN x1 and x2 DIRECTIONS, PERIODIC IN x3
    ********************************************************************************/
#if( BC_TYPE_CHOICE  == BC_CYLINDRICAL_OUTFLOW )
    /* 1-lower face */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1S-1-g][j][k], (N1S-1-g), j, k, RR, BCDN );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	}
      }
    }

    /* 1-upper face */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
    }
    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL + STATIC IN RADIAL x1 DIRECTION, REFLECTIVE AXI-SYM IN x2, PERIODIC IN x3
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_STATIC )
    /* 1-lower face */
    /*  nothing     */

    /* 1-upper face */
    /*  nothing     */

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /********************************************************************************
     Keep the ghost zones the same as the initial solution (useful, e.g., for Bondi evolution)
    ********************************************************************************/
#if( (BC_TYPE_CHOICE  == BC_STATIC_ALL) )
    /* Do not set any boundary conditions  */
#endif


    /********************************************************************************
     USER SPECIFIED BOUNDARY CONDITION
    ********************************************************************************/
#if( BC_TYPE_CHOICE  == BC_USER_SUPPLIED )
    fprintf(stderr,"\n\nbounds(): Please supply special boundary condition!! \n\n");
    fflush(stderr);
    fail(FAIL_BASIC,0);
#endif

    /**************************************************************************************
     SPHERICAL with special origin cell, matching at origin and axis, outflow IN RADIAL x1 DIRECTION
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_FULLSPHERE )

    /* Note that if the subdomain covers an entire dimension, then its bc_pid is set to BC_PHYS 
       and it must perform the symmetry operation using its own data of course;   */
    
       
    /* 1-lower face :  this should call the routine that calculates the fluxes along the edge, updates
       the origin cell, and then fills the ghost zone space with the correct values.  */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
    }

    /* 1-upper face :  typical outflow boundary conditions */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k2][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][U3] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k2] ;
	}
      }
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k2][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][U3] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k2] ;
	}
      }
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL with full interdomain matching at origin and axis, outflow IN RADIAL x1 DIRECTION
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_FULLSPHERE2 )

    /* Note that if the subdomain covers an entire dimension, then its bc_pid is set to BC_PHYS 
       and it must perform the symmetry operation using its own data of course;   */

    /* 1-lower face :   */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      if( totalsize[2] != N2 ) { 
	fprintf(stdout,"bounds(): 1 totalsize[2] and N2 should match :  %d %d \n", totalsize[2], N2); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      GLOOP     N2_LOOP   { 
	j2 = totalsize[2] - 1 - j + 2*NG; 
	N3_LOOP  {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S+g][j2][k2][l] ;
	  prim_arr[N1S-1-g][j][k][U1] *= -1 ;
	  prim_arr[N1S-1-g][j][k][U2] *= -1 ;
	  prim_arr[N1S-1-g][j][k][B1] *= -1 ;
	  prim_arr[N1S-1-g][j][k][B2] *= -1 ;
	}
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   { 
	  j2 = totalsize[2] - 1 - j + 2*NG; 
	  N3_LOOP  {
	    k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	    pflag[N1S-1-g][j][k] = pflag[N1S+g][j2][k2] ;
	  }
	}
      }
    }

    /* 1-upper face :  typical outflow boundary conditions */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
#if(EQUATORIAL_RUN)
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
#else
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 1 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k2][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][U3] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k2] ;
	}
      }
#endif
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
#if(EQUATORIAL_RUN)
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
#else
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 2 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k2][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][U3] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k2] ;
	}
      }
#endif
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 3 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 4 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL with full interdomain matching at origin and axis, setting v=0 at origin, outflow IN RADIAL x1 DIRECTION
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_FULLSPHERE3 )

    /* Note that if the subdomain covers an entire dimension, then its bc_pid is set to BC_PHYS 
       and it must perform the symmetry operation using its own data of course;   */
    
       
    /* 1-lower face :   */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      if( totalsize[2] != N2 ) { 
	fprintf(stdout,"bounds(): 1 totalsize[2] and N2 should match :  %d %d \n", totalsize[2], N2); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      GLOOP     N2_LOOP   { 
	j2 = totalsize[2] - 1 - j + 2*NG; 
	N3_LOOP  {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  prim_arr[N1S-1-g][j][k][RHO] = prim_arr[N1S+g][j2][k2][RHO] ;
	  prim_arr[N1S-1-g][j][k][UU ] = prim_arr[N1S+g][j2][k2][UU ] ;
	  prim_arr[N1S-1-g][j][k][U1 ] = 0. ;
	  prim_arr[N1S-1-g][j][k][U2 ] = 0. ;
	  prim_arr[N1S-1-g][j][k][U3 ] = 0. ;
	  prim_arr[N1S-1-g][j][k][B1 ] = 0. ;
	  prim_arr[N1S-1-g][j][k][B2 ] = 0. ;
	  prim_arr[N1S-1-g][j][k][B3 ] = 0. ;
	}
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   { 
	  j2 = totalsize[2] - 1 - j + 2*NG; 
	  N3_LOOP  {
	    k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	    pflag[N1S-1-g][j][k] = pflag[N1S+g][j2][k2] ;
	  }
	}
      }
    }

    /* 1-upper face :  typical outflow boundary conditions */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
#if(EQUATORIAL_RUN)
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
#else
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k2][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][U3] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k2] ;
	}
      }
    }
#endif

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
#if(EQUATORIAL_RUN)
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
#else
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k2][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][U3] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k2] ;
	}
      }
    }
#endif

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif

    /**************************************************************************************
     SPHERICAL with outflow in inner and outer radial boundaries 
    ***************************************************************************************/
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_FULLSPHERE4 )

    /* Note that if the subdomain covers an entire dimension, then its bc_pid is set to BC_PHYS 
       and it must perform the symmetry operation using its own data of course;   */

    /* 1-lower face :   */
    if( bc_pid[1][BCDN] == BC_PHYS ) {
      if( totalsize[2] != N2 ) { 
	fprintf(stdout,"bounds(): 1 totalsize[2] and N2 should match :  %d %d \n", totalsize[2], N2); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      GLOOP     N2_LOOP   { 
	N3_LOOP  {
	  PLOOP prim_arr[N1S-1-g][j][k][l] = prim_arr[N1S][j][k][l] ;
//	  prim_arr[N1S-1-g][j][k][U1] *= -1 ;
//	  prim_arr[N1S-1-g][j][k][U3] *= -1 ;
//	  prim_arr[N1S-1-g][j][k][B1] *= -1 ;
//	  prim_arr[N1S-1-g][j][k][B3] *= -1 ;
	}
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   { 
	  N3_LOOP  {
	    pflag[N1S-1-g][j][k] = pflag[N1S][j][k] ;
	  }
	}
      }
    }

    /* 1-upper face :  typical outflow boundary conditions */
    if( bc_pid[1][BCUP] == BC_PHYS ) {
      GLOOP     N2_LOOP   N3_LOOP  {
	PLOOP prim_arr[N1E+1+g][j][k][l] = prim_arr[N1E][j][k][l] ;
      }
      GLOOP     N2_LOOP   N3_LOOP  {
	inflow_check(prim_arr[N1E+1+g][j][k], (N1E+1+g), j, k, RR, BCUP );
      }
      if( failure_exists ) { 
	GLOOP     N2_LOOP   N3_LOOP  {
	  pflag[N1E+1+g][j][k] = pflag[N1E][j][k] ;
	}
      }
    }

    /* 2-lower face */
    if( bc_pid[2][BCDN] == BC_PHYS ) {
#if(EQUATORIAL_RUN)
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2S-1-g][k] = pflag[i][N2S][k] ;
	}
      }
#else
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 1 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2S-1-g][k][l] = prim_arr[i][N2S+g][k2][l] ;
	prim_arr[i][N2S-1-g][k][U2] *= -1 ;
	prim_arr[i][N2S-1-g][k][U3] *= -1 ;
	prim_arr[i][N2S-1-g][k][B2] *= -1 ;
	prim_arr[i][N2S-1-g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2S-1-g][k] = pflag[i][N2S+g][k2] ;
	}
      }
#endif
    }

    /* 2-upper face */
    if( bc_pid[2][BCUP] == BC_PHYS ) {
#if(EQUATORIAL_RUN)
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E][k][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  pflag[i][N2E+1+g][k] = pflag[i][N2E][k] ;
	}
      }
#else
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 2 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      GLOOP    N3_LOOP   {
	k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	PLOOP prim_arr[i][N2E+1+g][k][l] = prim_arr[i][N2E-g][k2][l] ;
	prim_arr[i][N2E+1+g][k][U2] *= -1 ;
	prim_arr[i][N2E+1+g][k][U3] *= -1 ;
	prim_arr[i][N2E+1+g][k][B2] *= -1 ;
	prim_arr[i][N2E+1+g][k][B3] *= -1 ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      GLOOP    N3_LOOP   {
	  k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
	  pflag[i][N2E+1+g][k] = pflag[i][N2E-g][k2] ;
	}
      }
#endif
    }

    /* 3-lower face */
    if( bc_pid[3][BCDN] == BC_PHYS ) {
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 3 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3S-1-g][l] = prim_arr[i][j][N3E-g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3S-1-g] = pflag[i][j][N3E-g] ;
	}
      }
    }

    /* 3-upper face */
    if( bc_pid[3][BCUP] == BC_PHYS ) {
      if( totalsize[3] != N3 ) { 
	fprintf(stdout,"bounds(): 4 totalsize[3] and N3 should match :  %d %d \n", totalsize[3], N3); 
	fflush(stdout);     fail(FAIL_BASIC,0) ; 
      }
      N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	PLOOP prim_arr[i][j][N3E+1+g][l] = prim_arr[i][j][N3S+g][l] ;
      }
      if( failure_exists ) { 
	N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
	  pflag[i][j][N3E+1+g] = pflag[i][j][N3S+g] ;
	}
      }
    }

#endif




  }  // end of   if( (!skip_mpi) !! boundary_phys_pflag) ) { 


    /********************************************************************************
     MPI BOUNDARY COMMUNICATION ROUTINES:
  ********************************************************************************/
#if( USEMPI ) 

  if( (!skip_mpi) || boundary_mpi_pflag ) { 
    //    fprintf(stdout,"Pre bounds_mpi() pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);
    bounds_mpi(prim_arr);
    //    bounds_mpi_fast(prim_arr);
    //    fprintf(stdout,"Post bounds_mpi() pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);
    //    fprintf(stdout,"Pre bounds_pflag_mpi() pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);
    bounds_pflag_mpi(pflag);
    //    bounds_pflag_mpi_fast(pflag);
    //    fprintf(stdout,"Post bounds_pflag_mpi() pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);
  }

  if( !skip_mpi ) { 
    mpi_global_imax2(&boundary_mpi_pflag,&boundary_phys_pflag);
  }

  //  fprintf(stdout,"Post bounds_pflag_mpi() %d \n", myid); fflush(stdout);
#endif
  

  /********************************************************************************
     ADDITIONAL BOUNDARY CONDITIONS THAT REQUIRE MPI INFORMATION : 
  ********************************************************************************/
  /* We need B1 and B3 set in all ghost zones to impose divergence constraint */
#if( BC_TYPE_CHOICE  == BC_SPHERICAL_CUTOUT )
//  if( (!skip_mpi) || boundary_phys_pflag) { 
//    if( (bc_pid[2][BCUP] == BC_PHYS) || (bc_pid[2][BCDN] == BC_PHYS) ) {
//      //    set_cutout_boundary( prim_arr );
//      //    set_cutout_boundary2( prim_arr );
//      set_cutout_boundary3( prim_arr );
//      //    set_cutout_boundary4( prim_arr );
//    }
//  }
#endif  


  /********************************************************************************
     Remove gdet scale factor (need to do this for all cells including ghosts): 
  ********************************************************************************/
#if( RESCALE_B )
  ALL_LOOP { 
    get_geometry(i,j,k,CENT,ncurr,geom);
    geomfactor = geom->g_inv;
    BLOOP { prim_arr[i][j][k][l] *= geomfactor; }
  }
#endif

  /*******************************************************************************
    Special boundary conditions for excision outside of horizons ( mask )
  *********************************************************************************/
#if( USE_MASK && ( EXCISION_TYPE_CHOICE == EXCISE_IZ || EXCISION_TYPE_CHOICE == EXCISE_ISCO ) )
  struct of_coord *coords;
  LOOP {
    if( evol_mask[ncurr][i][j][k] != MASK_NORMAL ) {
      get_coord(i,j,k,CENT,ncurr,coords);
      /* Ramp down density and Ramp down uu */
      prim_arr[i][j][k][RHO] /= 100.;
      prim_arr[i][j][k][UU] /= 10000.;
      if( prim_arr[i][j][k][RHO] < coords->rhoflr ) { prim_arr[i][j][k][RHO] = coords->rhoflr; }
      if( prim_arr[i][j][k][UU ] < coords->uuflr  ) { prim_arr[i][j][k][UU ] = coords->uuflr;  }
    }
  }
#endif


  TRACE_END;

  return; 

}



/**********************************************************************************/
/**********************************************************************************
  fix_flux(): 
  -------------
    -- impose the coordinate system's symmetries in the flux functions; 
    -- uses BC_TYPE_CHOICE  do determine which symmetry to impose: 

**********************************************************************************/
void fix_flux( void ) 
{
  int i,j,k,l ;

  TRACE_BEG;

  /*******************************************************************************
    Special axisymmetry conditions along the axis (in those w/ a symmetry axis) 
  *********************************************************************************/
#if(BC_TYPE_CHOICE==BC_CYLINDRICAL_OUTFLOW||BC_TYPE_CHOICE==BC_SPHERICAL_OUTFLOW||BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT||BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT2||BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT3||BC_TYPE_CHOICE==BC_SPHERICAL_OUTFLOW2||BC_TYPE_CHOICE==BC_SPHERICAL_OUTFLOW3)
//  if( bc_pid[2][BCDN] == BC_PHYS ) { 
//    N1ALL_LOOP N3ALL_LOOP PLOOP  { 
//      F[i][N2S  ][k][1][l] = 0. ;   /* No flux thru coord. singularity at the axis : */
//      //    F[i][N2S-1][k][0][B2] = -F[i][N2S][k][0][B2] ;  /* EMF Axi-symmetry Condition ?? */ 
//    }
//  }
//  if( bc_pid[2][BCUP] == BC_PHYS ) { 
//    N1ALL_LOOP N3ALL_LOOP PLOOP  { 
//      F[i][N2E+1][k][1][l] = 0. ;  /* No flux thru coord. singularity at the axis : */
//    //    F[i][N2E+1][k][0][B2] = -F[i][N2E][k][0][B2] ;  /* EMF Axi-symmetry Condition ?? */ 
//    }
//  }
#endif
    
  /*******************************************************************************
    Radial flux condition to make sure that there is no mass flux along the 
     boundaries at the radial extremes: 
  *********************************************************************************/
#if(BC_TYPE_CHOICE==BC_CYLINDRICAL_OUTFLOW||BC_TYPE_CHOICE==BC_SPHERICAL_OUTFLOW||BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT||BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT2||BC_TYPE_CHOICE==BC_SPHERICAL_CUTOUT3||BC_TYPE_CHOICE==BC_SPHERICAL_OUTFLOW2||BC_TYPE_CHOICE==BC_SPHERICAL_OUTFLOW3||BC_TYPE_CHOICE==BC_SPHERICAL_CONSINTERP)
  if( bc_pid[1][BCDN] == BC_PHYS ) { 
    N2ALL_LOOP  N3ALL_LOOP { 
      if(F[N1S  ][j][k][0][RHO]>0) F[N1S  ][j][k][0][RHO]=0.;
    }
  }
  if( bc_pid[1][BCUP] == BC_PHYS ) { 
    N2ALL_LOOP  N3ALL_LOOP { 
      if(F[N1E+1][j][k][0][RHO]<0) F[N1E+1][j][k][0][RHO]=0.;
    }
  }
#endif 

#if(BC_TYPE_CHOICE==BC_SPHERICAL_FULLSPHERE||BC_TYPE_CHOICE==BC_SPHERICAL_FULLSPHERE2||BC_TYPE_CHOICE==BC_SPHERICAL_FULLSPHERE3||BC_TYPE_CHOICE==BC_SPHERICAL_FULLSPHERE4)
//   if( bc_pid[1][BCDN] == BC_PHYS ) { 
//     N2ALL_LOOP  N3ALL_LOOP  PLOOP {  F[N1S][j][k][0][l]=0.;   }
//   }
#endif


  TRACE_END;

  return;

}

/**********************************************************************************/
/**********************************************************************************
  inflow_check():
  -------------
    -- makes sure the dir-component of the fluid's 4-velocity is pointing 
        outwards at the specified boundary; "dir" must be a spatial component
        in "x" (or "non-numerical") coordinates; 

    -- will eventually need to transform 4-vel. here to "x" coordinates, 
       but we are leaving it as is for now in order to compare to harm2d;

    -- impose the coordinate system's symmetries in the flux functions; 

       type = 0  : DN/lower (aka "left", "inner", "horizon", etc.) boundary
            = 1  : UP/upper (aka "right", "outer", etc.) boundary

**********************************************************************************/
#define GAMMA_MAX_X1DN (5.) 

void inflow_check(double *pr, int ii, int jj, int kk, int dir, int type )
{
  int j,k ;
  double ucon[NDIM] ;
  double gamma,vsq,factor ;
  struct of_geom *geom ;
	
  if( dir < 1 || dir > 3 ) { 
    fprintf(stderr,"inflow_check(): invalid value of dir = %d \n", dir); 
    fflush(stderr);
    fail(FAIL_BASIC,0) ; 
  }

  get_geometry(ii,jj,kk,CENT,ncurr,geom) ;
  ucon_calc(pr, geom, ucon) ;

  if( ((ucon[dir] > 0.) && (type==0)) || ((ucon[dir] < 0.) && (type==1)) ) { 
    /* find gamma and remove it from primitives */
    gamma = geom->alpha * ucon[TT];
    pr[U1] /= gamma ;
    pr[U2] /= gamma ;
    pr[U3] /= gamma ;

    /* reset radial velocity so radial 4-velocity
     * is zero */
    pr[U1+dir-1] = geom->beta[dir-1]*geom->ncon[0] ;

    /* now find new gamma and put it back in */
    vsq = 0. ;
    for(j=1; j<NDIM; j++) for(k=1; k<NDIM; k++) { 
      vsq += geom->gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
    }
    //    if( fabs(vsq) < 1.e-13 )  vsq = 1.e-13;
    if( fabs(vsq) < 1.e-13 )  vsq = fabs(vsq);
    if( vsq < 0. )  vsq = 0. ; 
    if( vsq >= 1. ) { 
      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX) ;
    }
    gamma = 1./sqrt(1. - vsq) ;

#if( USE_GAMMA_CEILING_X1DN_BC ) 
    if( (type==0) &&  (gamma > GAMMA_MAX_X1DN) ) { gamma = GAMMA_MAX_X1DN; } 
#endif

    pr[U1] *= gamma ;
    pr[U2] *= gamma ;
    pr[U3] *= gamma ;

    /* done */
  }
  else if( type==0 ) { 
#if( USE_GAMMA_CEILING_X1DN_BC ) 
    gamma = geom->alpha * ucon[TT];
    if( gamma > GAMMA_MAX_X1DN ) {  
      factor = sqrt( (GAMMA_MAX_X1DN-1.)*(GAMMA_MAX_X1DN+1.) / ((gamma-1.)*(gamma+1.)) );
      pr[U1] *= factor ;
      pr[U2] *= factor ;
      pr[U3] *= factor ;
    }
#endif
  }

  return;

}


/**********************************************************************************/
/**********************************************************************************
  inflow_check2():
  -------------
    -- makes sure the dir-component of the fluid's 4-velocity is pointing 
        outwards at the specified boundary; "dir" must be a spatial component
        in "x" (or "non-numerical") coordinates; 

    -- only different than inflow_check() in that it sets pr[U1] to zero, not ucon[1]; 

    -- will eventually need to transform 4-vel. here to "x" coordinates, 
       but we are leaving it as is for now in order to compare to harm2d;

    -- impose the coordinate system's symmetries in the flux functions; 

       type = 0  : DN/lower (aka "left", "inner", "horizon", etc.) boundary
            = 1  : UP/upper (aka "right", "outer", etc.) boundary

**********************************************************************************/
#define GAMMA_MAX_X1DN (5.) 

void inflow_check2(double *pr, int ii, int jj, int kk, int dir, int type )
{
  int j,k ;
  double ucon[NDIM] ;
  double gamma,vsq,factor ;
  struct of_geom *geom ;
	
  if( dir < 1 || dir > 3 ) { 
    fprintf(stderr,"inflow_check(): invalid value of dir = %d \n", dir); 
    fflush(stderr);
    fail(FAIL_BASIC,0) ; 
  }

  get_geometry(ii,jj,kk,CENT,ncurr,geom) ;
  ucon_calc(pr, geom, ucon) ;

  if( ((ucon[dir] > 0.) && (type==0)) || ((ucon[dir] < 0.) && (type==1)) ) { 
    /* find gamma and remove it from primitives */
    gamma = geom->alpha * ucon[TT];
    pr[U1] /= gamma ;
    pr[U2] /= gamma ;
    pr[U3] /= gamma ;

    /* reset radial velocity so radial 4-velocity
     * is zero */
    pr[U1+dir-1] = 0.;

    /* now find new gamma and put it back in */
    vsq = 0. ;
    for(j=1; j<NDIM; j++) for(k=1; k<NDIM; k++) { 
      vsq += geom->gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
    }
    //    if( fabs(vsq) < 1.e-13 )  vsq = 1.e-13;
    if( fabs(vsq) < 1.e-13 )  vsq = fabs(vsq);
    if( vsq < 0. )  vsq = 0. ; 
    if( vsq >= 1. ) { 
      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX) ;
    }
    gamma = 1./sqrt(1. - vsq) ;

#if( USE_GAMMA_CEILING_X1DN_BC ) 
    if( (type==0) &&  (gamma > GAMMA_MAX_X1DN) ) { gamma = GAMMA_MAX_X1DN; } 
#endif

    pr[U1] *= gamma ;
    pr[U2] *= gamma ;
    pr[U3] *= gamma ;

    /* done */
  }
  else if( type==0 ) { 
#if( USE_GAMMA_CEILING_X1DN_BC ) 
    gamma = geom->alpha * ucon[TT];
    if( gamma > GAMMA_MAX_X1DN ) {  
      factor = sqrt( (GAMMA_MAX_X1DN-1.)*(GAMMA_MAX_X1DN+1.) / ((gamma-1.)*(gamma+1.)) );
      pr[U1] *= factor ;
      pr[U2] *= factor ;
      pr[U3] *= factor ;
    }
#endif
  }

  return;

}


/**********************************************************************************/
/**********************************************************************************
  limit_gamma_in_horizon():
  ---------------
    -- Rescales the velocity so that gamma never surpasses GAMMA_MAX_X1DN ;
    -- uses the same rescaling procedure as used in fixup1zone(); 
    -- all cells but the ones abutting the horizon are rescaled;
    -- this routine only rescales the physical cells, so it must be called
        before the boundary conditions are set; 
    -- if there is a failure of any kind in calculating/rescaling the 
       velocities, then they are not rescaled;  (there should not be a problem 
       at this stage of the update since all cells should have been fixed to 
       valid values by now); 
 
**********************************************************************************/
void limit_gamma_in_horizon( double ****prim  ) 
{
  int i,j,k,n1end;
  double gamma, ftmp, f;
  struct of_geom *geom ;

  n1end = N1S + n_within_horizon - 1;

  ftmp = (GAMMA_MAX_X1DN-1.)*(GAMMA_MAX_X1DN+1.);

  for(i=N1S; i<n1end; i++) N2_LOOP  N3_LOOP { 
      get_geometry(i,j,k,CENT,ncurr,geom) ; 
      if( !gamma_calc(prim[i][j][k], geom, &gamma) ) { 
	if( gamma > GAMMA_MAX_X1DN ) {
	  f = sqrt( ftmp / ((gamma-1.)*(gamma+1.)) );
	  prim[i][j][k][U1] *= f ;	
	  prim[i][j][k][U2] *= f ;	
	  prim[i][j][k][U3] *= f ;
	}
      }
    }

  return;
}


#undef GAMMA_MAX_X1DN


/**********************************************************************************/
/**********************************************************************************
  set_cutout_boundary(): 
  ---------------
    -- Reflection boundary conditions are imposed at the cutout surface 
       that surrounds the axis.  The specific conditions used in this routine 
       are zero-derivatives in all primitive variables except for U2, which 
       is enforced to be zero at the boundary. 

    -- The boundary state is the following : 
           1) (gdet*B^i) or (B^i) , depending on RESCALE_B, are constant 
           2) rho,u, u^1,u^3 are constant through boundary
           3)  u^2 = 0   (from mass flux condition) and ghost zones set to reconstruct same

    -- The ghost zones are then set so that minmod reconstruction will 
        result in no poloidal flux along the boundary.  This involves
        writing over the first real cell's value if it is too different 
        from the desired value at the boundary.  Linear extrapolation 
        is then used with these boundary values to set the ghost values.
        Thus, this routine assumes that we are using the minmod limiter 
        in the X2 direction along the X2 face.

    -- Assumes that (rho,u,B1,B3,u1,u3) have already been set in all ghost zones; 

    -- This routine is also responsible for setting pflag[] in X2-dir appropriately.

    -- Julian calls this "Option 4"
 
**********************************************************************************/
void set_cutout_boundary( double ****prim )
{
  int i,j,k,l,g,ii,d, sig_g;
  int jm2,jm1,jp1,jp2,pflag1,pflag2,pflag_out;
  double f0,f1,f2;
  double df10,df21;
  double p_0[NP];

  static int first_call=1;


  /* Implicitly uses the fact that there are at least 3 real cells in x^2 direction */
  if( N2 < 3 )  { 
    fprintf(stderr,"set_cutout_boundary(): N2 needs to be larger than 2:  %d \n", N2);
    fflush(stderr); 
    fail(FAIL_BASIC,0);
  }

  /*************************************************************************************
    Set the reconstruction mask function that determines where to use which slope limiter
      -- for now use MINMOD limiter along the first 2 cells (number of cells here  
         should be equal to the number of ghost cells required for the limiter used 
         to set the boundary values). 
      -- only need to use special limiter in the X2-direction to be consisten;
  *************************************************************************************/
  /* Limiter used only along X2-faces in X2-direction */
  if( first_call ) { 
#if( USE_LOCAL_RECON_TYPE )
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP GLOOP N3ALL_LOOP { recon_type[i][N2S+g][k][FACE2] = RECON_MINMOD ; }
    }
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP GLOOP N3ALL_LOOP { recon_type[i][N2E-g][k][FACE2] = RECON_MINMOD ; }
    }
#endif
    first_call = 0;
  }

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-LOWER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
      icurr = i; jcurr = jp1; kcurr = k; 
      pflag1 = pflag2 = 0;

      PLOOP p_0[l] = prim[i][jp1][k][l];
      p_0[U2] = 0.; 

      l = U2;
      f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1); }

      if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
      pflag1 *= pflag[i][jp1][k];
      pflag2 *= pflag[i][jp2][k];
      pflag_out = MAX(pflag1,pflag2);
      pflag[i][jp1][k] = pflag_out;
      GLOOP { pflag[i][jm1-g][k] = pflag_out; } 
    }
  }

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-UPPER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCUP] == BC_PHYS ) {
    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
      icurr = i; jcurr = jm1; kcurr = k; 
      pflag1 = pflag2 = 0;

      PLOOP p_0[l] = prim[i][jm1][k][l];
      p_0[U2] = 0.; 

      l = U2;
      f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

      if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
      pflag1 *= pflag[i][jm1][k];
      pflag2 *= pflag[i][jm2][k];
      pflag_out = MAX(pflag1,pflag2);
      pflag[i][jm1][k] = pflag_out; 
      GLOOP { pflag[i][jp1+g][k] = pflag_out; } 
    }
  }

  /*************************************************************************************
     Debugging output: 
  *************************************************************************************/
//  for(i=1;i<N1TOT-1;i++) for(k=1;k<N3TOT-1;k++) {
//    fprintf(stdout,"############################################################# \n");
//    fprintf(stdout,"########### PRIM  TEST                ####################### \n");
//    fprintf(stdout,"############################################################# \n");
//    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
//    PLOOP{ 
//      fprintf(stdout,"prim(%d,%d,%d,%d):  %28.18e %28.18e %28.18e %28.18e \n",i,jm2,k,l,
//	      prim[i][jm2][k][l],prim[i][jm1][k][l],prim[i][jp1][k][l],prim[i][jp2][k][l]); 
//    }
//    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
//    PLOOP{ 
//      fprintf(stdout,"prim(%d,%d,%d,%d):  %28.18e %28.18e %28.18e %28.18e \n",i,jm2,k,l,
//	      prim[i][jm2][k][l],prim[i][jm1][k][l],prim[i][jp1][k][l],prim[i][jp2][k][l]); 
//    }
//  }
//  fflush(stdout);

  return; 
}

/**********************************************************************************/
/**********************************************************************************
  set_cutout_boundary2(): 
  ---------------
    -- Reflection boundary conditions are imposed at the cutout surface 
       that surrounds the axis.  These conditions are primarily on the 
       mass, energy and momentum fluxes (i.e. their poloidal components 
       must be 0).  In order to be consistent with the ghost zone values, 
       we set the ghosts zones values such that the flux conditions are 
       satisfied.  Since the numerical flux involves a sum of the left and 
       right reconstructed fluxes, then we must ensure that the left and right
       ones are the same, as well as the left/right conserved variables and 
       primitive variables.  We must therefore make sure that the limiter 
       reconstructs in the way we want, so we need to interpolate over 
       the first physical cell nearest the cutout and force the code to 
       use one type of limiter near there (since we do not want to figure 
       this out for every limiter and the minmod limiter will probably 
       work best there since it is so diffusive).  

    -- The boundary state is the following : 
           1) B^1 = B^3 = 0 at boundary
           2) rho, v^{1,3}, B^2 are constant through boundary;
           3)  u^2 = u_2 = 0   (from mass flux condition)
           4)  p = bsq/2   (momentum flux)

    -- The ghost zones are then set so that minmod reconstruction will 
        result in no poloidal flux along the boundary.  This involves
        writing over the first real cell's value if it is too different 
        from the desired value at the boundary.  Linear extrapolation 
        is then used with these boundary values to set the ghost values.
        Thus, this routine assumes that we are using the minmod limiter 
        in the X2 direction along the X2 face.

    -- Assumes that (rho,v1,v3,B2) have already been set in all ghost zones; B1,B3 are 
        needed to impose divergence constraint for B2. Hence, only need to do 
         steps 2 - 7. 

    -- Assumes that x^2 is orthogonal to x^1 and x^3 in that g_{2 a} = g^{2 a} = 0 
          for all a!=2 ;  

    -- This routine is also responsible for setting pflag[] in X2-dir appropriately.

    -- Julian calls this "Option 3"
 
**********************************************************************************/
void set_cutout_boundary2( double ****prim )
{
  int i,j,k,l,g,ii,d,sig_g;
  int jm2,jm1,jp1,jp2,pflag1,pflag2,pflag_out;
  double f0,f1,f2;
  double df10,df21;
  double p_0[NP], bsq_bd, gamma, u_flr;
  struct of_geom *geom;
  struct of_coord *coords;


  /* Implicitly uses the fact that there are at least 3 real cells in x^2 direction */
  if( N2 < 3 )  { 
    fprintf(stderr,"set_cutout_boundary(): N2 needs to be larger than 2:  %d \n", N2);
    fflush(stderr); 
    fail(FAIL_BASIC,0);
  }

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-LOWER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
	get_coord(   i,jp1,k,FACE2,ncurr,coords);
	get_geometry(i,jp1,k,FACE2,ncurr,geom);

	p_0[RHO] = prim[i][jp1][k][RHO];
	p_0[U1 ] = prim[i][jp1][k][U1 ];
	p_0[U2 ] = 0.;
	p_0[U3 ] = prim[i][jp1][k][U3 ];
	p_0[B1 ] = 0.;
	p_0[B2 ] = prim[i][jp1][k][B2 ];
	p_0[B3 ] = 0.;

	if( gamma_calc(p_0,geom,&gamma) ) { 
	  fprintf(stderr,"set_cutout_boundary2(): bad gamma: %28.18e : %28.18e %28.18e %28.18e \n",
		  gamma,p_0[U1],p_0[U2],p_0[U3]); 
	  fflush(stderr);
	}

	bsq_bd  =  geom->alpha*p_0[B2]/gamma;
	bsq_bd  *=  bsq_bd * geom->gcov[2][2];
	p_0[UU] = 0.5*bsq_bd / ( gam - 1.);

	u_flr = UUMIN*pow(coords->x[RR],UUPOWER);
	p_0[UU] = MAX( p_0[UU] , u_flr );

	//      fprintf(stdout,"p_0(%d,%d,%d): %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e \n",
	//	      i,jp1,k,p_0[RHO],p_0[UU],p_0[U1],p_0[U2],p_0[U3],p_0[B1],p_0[B2],p_0[B3]);
	//      fflush(stdout);


	/*************************************************************************************
       Now extrapolate into the ghost zones (extrap. v^i) : 
	*************************************************************************************/
	icurr = i; jcurr = jp1; kcurr = k; 
	pflag1 = pflag2 = 0; 

	l = UU; 
	f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  pflag2 = 1;
	}
	else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = MAX( f0 + sig_g*(f0-f1) , u_flr ); }

	l = U2;
	f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  pflag2 = 1;
	}
	else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1); }

	l = B1;
	f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  //	pflag2 = 1;
	}
	//      else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1); }

	l = B3;
	f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  //	pflag2 = 1;
	}
	//      else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1); }

	if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
	pflag1 *= pflag[i][jp1][k];
	pflag2 *= pflag[i][jp2][k];
	pflag_out = MAX(pflag1,pflag2);
	pflag[i][jp1][k] = pflag_out;
	GLOOP { pflag[i][jm1-g][k] = pflag_out; } 
      }
  }

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-UPPER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCUP] == BC_PHYS ) {
    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
	get_coord(   i,jp1,k,FACE2,ncurr,coords);
	get_geometry(i,jp1,k,FACE2,ncurr,geom);

	p_0[RHO] = prim[i][jp1][k][RHO];
	p_0[U1 ] = prim[i][jp1][k][U1 ];
	p_0[U2 ] = 0.;
	p_0[U3 ] = prim[i][jp1][k][U3 ];
	p_0[B1 ] = 0.;
	p_0[B2 ] = prim[i][jp1][k][B2 ];
	p_0[B3 ] = 0.;

	if( gamma_calc(p_0,geom,&gamma) ) { 
	  fprintf(stderr,"set_cutout_boundary2(): bad gamma: %28.18e : %28.18e %28.18e %28.18e \n",
		  gamma,p_0[U1],p_0[U2],p_0[U3]); 
	  fflush(stderr);
	}

	bsq_bd  =  geom->alpha*p_0[B2]/gamma;
	bsq_bd  *=  bsq_bd * geom->gcov[2][2];
	p_0[UU] = 0.5*bsq_bd / ( gam - 1.);

	u_flr = UUMIN*pow(coords->x[RR],UUPOWER);
	p_0[UU] = MAX( p_0[UU] , u_flr );

	/*************************************************************************************
       Now extrapolate into the ghost zones (extrap. v^i) : 
	*************************************************************************************/
	icurr = i; jcurr = jm1; kcurr = k; 
	pflag1 = pflag2 = 0; 

	l = UU; 
	f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  pflag2 = 1;
	}
	else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = MAX( f0 + sig_g*(f0-f1) , u_flr );  }

	l = U2;
	f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  pflag2 = 1;
	}
	else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

	l = B1;
	f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  //	pflag2 = 1;
	}
	//      else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

	l = B3;
	f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
	df10 = f1-f0;   df21 = f2-f1;  
	/* Change first real cell if it does not lead to desired recon. at boundary: */
	if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	  f1 = (2*f0 + f2)/3.; 
	  //	pflag2 = 1;
	}
	//      else{ pflag1 = 1; }
	for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

	if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
	pflag1 *= pflag[i][jm1][k];
	pflag2 *= pflag[i][jm2][k];
	pflag_out = MAX(pflag1,pflag2);
	pflag[i][jm1][k] = pflag_out; 
	GLOOP { pflag[i][jp1+g][k] = pflag_out; } 
      }
  }

  /*************************************************************************************
     Debugging output: 
  *************************************************************************************/
//  for(i=1;i<N1TOT-1;i++) for(k=1;k<N3TOT-1;k++) {
//    fprintf(stdout,"############################################################# \n");
//    fprintf(stdout,"########### PRIM  TEST                ####################### \n");
//    fprintf(stdout,"############################################################# \n");
//    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
//    PLOOP{ 
//      fprintf(stdout,"prim(%d,%d,%d,%d):  %28.18e %28.18e %28.18e %28.18e \n",i,jm2,k,l,
//	      prim[i][jm2][k][l],prim[i][jm1][k][l],prim[i][jp1][k][l],prim[i][jp2][k][l]); 
//    }
//    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
//    PLOOP{ 
//      fprintf(stdout,"prim(%d,%d,%d,%d):  %28.18e %28.18e %28.18e %28.18e \n",i,jm2,k,l,
//	      prim[i][jm2][k][l],prim[i][jm1][k][l],prim[i][jp1][k][l],prim[i][jp2][k][l]); 
//    }
//  }
//  fflush(stdout);

  return; 
}

/**********************************************************************************/
/**********************************************************************************
  set_cutout_boundary3(): 
  ---------------
    -- Reflection boundary conditions are imposed at the cutout surface 
       that surrounds the axis.  The specific conditions used in this routine 
       are zero-derivatives in all primitive variables except for U2 & B2, which 
       are enforced to be zero at the boundary. 

    -- The boundary state is the following : 
           1) (gdet*B^{1,3}) or (B^{1,3}) , depending on RESCALE_B, are constant 
           2) b^2 = B^2 = 0 at boundary ; 
           3) rho,u, u^1,u^3 are constant through boundary
           4)  u^2 = 0   (from mass flux condition) and ghost zones set to reconstruct same

    -- The ghost zones are then set so that minmod reconstruction will 
        result in no poloidal flux along the boundary.  This involves
        writing over the first real cell's value if it is too different 
        from the desired value at the boundary.  Linear extrapolation 
        is then used with these boundary values to set the ghost values.
        Thus, this routine assumes that we are using the minmod limiter 
        in the X2 direction along the X2 face.

    -- Assumes that (rho,u,B1,B3,u1,u3) have already been set in all ghost zones; 

    -- This routine is also responsible for setting pflag[] in X2-dir appropriately.

    -- Julian calls this "Option 1"

**********************************************************************************/
void set_cutout_boundary3( double ****prim )
{
  int i,j,k,l,g,ii,d,sig_g;
  int jm2,jm1,jp1,jp2,pflag1,pflag2,pflag_out;
  double f0,f1,f2;
  double df10,df21;
  double p_0[NP];

  static int first_call=1;

  //  fprintf(stdout,"set_cutout_boundary3 beg pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);

  /* Implicitly uses the fact that there are at least 3 real cells in x^2 direction */
  if( N2 < 3 )  { 
    fprintf(stderr,"set_cutout_boundary(): N2 needs to be larger than 2:  %d \n", N2);
    fflush(stderr); 
    fail(FAIL_BASIC,0);
  }

  /*************************************************************************************
    Set the reconstruction mask function that determines where to use which slope limiter
      -- for now use MINMOD limiter along the first 2 cells (number of cells here  
         should be equal to the number of ghost cells required for the limiter used 
         to set the boundary values). 
      -- only need to use special limiter in the X2-direction to be consisten;
  *************************************************************************************/
  /* Limiter used only along X2-faces in X2-direction */
  if( first_call ) { 
#if( USE_LOCAL_RECON_TYPE )
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP GLOOP N3ALL_LOOP { recon_type[i][N2S+g][k][FACE2] = RECON_MINMOD ; }
    }
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP GLOOP N3ALL_LOOP { recon_type[i][N2E-g][k][FACE2] = RECON_MINMOD ; }
    }
#endif
    first_call = 0;
  }

  //  fprintf(stdout,"set_cutout_boundary3 here1 pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-LOWER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
      icurr = i; jcurr = jp1; kcurr = k; 
      pflag1 = pflag2 = 0;

      PLOOP p_0[l] = prim[i][jp1][k][l];
      p_0[U2] = 0.; 
      p_0[B2] = 0.; 

      l = U2;
      f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1);}

      l = B2;
      f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      /* we need not set pflag for B2 since it is never interpolated over */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	//	pflag2 = 1;
      }
      //      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1); }

      if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
      pflag1 *= pflag[i][jp1][k];
      pflag2 *= pflag[i][jp2][k];
      pflag_out = MAX(pflag1,pflag2);
      pflag[i][jp1][k] = pflag_out;
      GLOOP { pflag[i][jm1-g][k] = pflag_out; } 
    }
  }

  //  fprintf(stdout,"set_cutout_boundary3 here2 pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-UPPER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCUP] == BC_PHYS ) {
    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
      icurr = i; jcurr = jm1; kcurr = k; 
      pflag1 = pflag2 = 0;

      PLOOP p_0[l] = prim[i][jm1][k][l];
      p_0[U2] = 0.; 
      p_0[B2] = 0.; 

      l = U2;
      f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

      l = B2;
      f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	//	pflag2 = 1;
      }
      //      else{ pflag1 = 1; }
      for(g=-1; g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

      if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
      pflag1 *= pflag[i][jm1][k];
      pflag2 *= pflag[i][jm2][k];
      pflag_out = MAX(pflag1,pflag2);
      pflag[i][jm1][k] = pflag_out; 
      GLOOP { pflag[i][jp1+g][k] = pflag_out; } 
    }
  }

  /*************************************************************************************
     Debugging output: 
  *************************************************************************************/
//  for(i=1;i<N1TOT-1;i++) for(k=1;k<N3TOT-1;k++) {
//    fprintf(stdout,"############################################################# \n");
//    fprintf(stdout,"########### PRIM  TEST                ####################### \n");
//    fprintf(stdout,"############################################################# \n");
//    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
//    PLOOP{ 
//      fprintf(stdout,"prim(%d,%d,%d,%d):  %28.18e %28.18e %28.18e %28.18e \n",i,jm2,k,l,
//	      prim[i][jm2][k][l],prim[i][jm1][k][l],prim[i][jp1][k][l],prim[i][jp2][k][l]); 
//    }
//    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
//    PLOOP{ 
//      fprintf(stdout,"prim(%d,%d,%d,%d):  %28.18e %28.18e %28.18e %28.18e \n",i,jm2,k,l,
//	      prim[i][jm2][k][l],prim[i][jm1][k][l],prim[i][jp1][k][l],prim[i][jp2][k][l]); 
//    }
//  }
//  fflush(stdout);

//  fprintf(stdout,"set_cutout_boundary3 end pid=%d  nstep=%d \n", myid,nstep); fflush(stdout);

  return; 
}

/**********************************************************************************/
/**********************************************************************************
  set_cutout_boundary4(): 
  ---------------
    -- just like set_cutout_boundary3() except uses Krolik's idea of diminishing
        pressure so that the total energy does not increase when changing 
        the first real cell's value of U2 and B2; 

    -- The boundary state is the following : 
           1) (gdet*B^{1,3}) or (B^{1,3}) , depending on RESCALE_B, are constant 
           2) b^2 = B^2 = 0 at boundary ; 
           3) rho,u, u^1,u^3 are constant through boundary
           4)  u^2 = 0   (from mass flux condition) and ghost zones set to reconstruct same

    -- The ghost zones are then set so that minmod reconstruction will 
        result in no poloidal flux along the boundary.  This involves
        writing over the first real cell's value if it is too different 
        from the desired value at the boundary.  Linear extrapolation 
        is then used with these boundary values to set the ghost values.
        Thus, this routine assumes that we are using the minmod limiter 
        in the X2 direction along the X2 face.

    -- Assumes that (rho,u,B1,B3,u1,u3) have already been set in all ghost zones; 

    -- This routine is also responsible for setting pflag[] in X2-dir appropriately.

    -- Julian calls this "Option 2"
 
**********************************************************************************/
void set_cutout_boundary4( double ****prim )
{
  int i,j,k,l,g,ii,d,sig_g;
  int jm2,jm1,jp1,jp2,pflag1,pflag2,pflag_out;
  double f0,f1,f2;
  double df10,df21;
  double p_0[NP];

  struct of_state q;
  struct of_geom *geom; 
  double U_1[NP], U_2[NP], p_1[NP], p_2[NP], E_tot1, E_tot2, E_test, u_try, A, B, u_flr ;
  struct of_coord *coords;


  static int first_call=1;


  /* Implicitly uses the fact that there are at least 3 real cells in x^2 direction */
  if( N2 < 3 )  { 
    fprintf(stderr,"set_cutout_boundary(): N2 needs to be larger than 2:  %d \n", N2);
    fflush(stderr); 
    fail(FAIL_BASIC,0);
  }

  /*************************************************************************************
    Set the reconstruction mask function that determines where to use which slope limiter
      -- for now use MINMOD limiter along the first 2 cells (number of cells here  
         should be equal to the number of ghost cells required for the limiter used 
         to set the boundary values). 
      -- only need to use special limiter in the X2-direction to be consisten;
  *************************************************************************************/
  /* Limiter used only along X2-faces in X2-direction */
  if( first_call ) { 
#if( USE_LOCAL_RECON_TYPE )
    if( bc_pid[2][BCDN] == BC_PHYS ) {
      N1ALL_LOOP GLOOP N3ALL_LOOP { recon_type[i][N2S+g][k][FACE2] = RECON_MINMOD ; }
    }
    if( bc_pid[2][BCUP] == BC_PHYS ) {
      N1ALL_LOOP GLOOP N3ALL_LOOP { recon_type[i][N2E-g][k][FACE2] = RECON_MINMOD ; }
    }
#endif
    first_call = 0;
  }

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-LOWER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    jm2 = N2S - 2;  jm1 = N2S-1;  jp1 = N2S;  jp2 = N2S+1;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
      icurr = i; jcurr = jp1; kcurr = k; 
      pflag1 = pflag2 = 0;

      PLOOP p_1[l] = p_0[l] = prim[i][jp1][k][l];
      p_0[U2] = 0.; 
      p_0[B2] = 0.; 

      l = U2;
      f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1);}

      l = B2;
      f0 = p_0[l];  f1 = prim[i][jp1][k][l];  f2 = prim[i][jp2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      //      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jm1-g][k][l] = f0 + sig_g*(f0-f1); }

      
      /* If there was any change, make sure that the total energy did not increase: */
      if( pflag2 ) { 
	PLOOP p_2[l] = prim[i][jp1][k][l];
	get_geometry(i,jp1,k,CENT,ncurr,geom); 

	get_state(p_1, geom, &q );  primtoflux( p_1, &q, 0, geom, U_1);
	get_state(p_2, geom, &q );  primtoflux( p_2, &q, 0, geom, U_2);

	E_tot1 = -(U_1[UU] - U_1[RHO])/geom->g; 
	E_tot2 = -(U_2[UU] - U_2[RHO])/geom->g; 

	/* Can express   E =  A + B*u , solve for new  u  ,  note that E = -T^t_t */
	/* Assumes  P = (gam-1)*u  */
	if( E_tot2 > E_tot1 ) { 
	  A = -(p_2[RHO] + q.bsq)*q.ucon[TT]*q.ucov[TT] - 0.5*q.bsq + q.bcon[TT]*q.bcov[TT];
	  B = -( gam*q.ucon[TT]*q.ucov[TT] + gam - 1. );
	  u_try = (E_tot1 - A) / B; 

	  E_test = A + B * u_try;

	  get_coord(i,jp1,k,FACE2,ncurr,coords);
	  u_flr = UUMIN*pow(coords->x[RR],UUPOWER);

	  if( u_try < u_flr )  u_try = u_flr; 
	  
	  prim[i][jp1][k][UU] = u_try; 
//	  fprintf(stdout,"DN E_tot1 = %20.10e , E_tot2 = %20.10e , E_test = %20.10e, u_old = %20.10e , u_new = %20.10e \n",
//		  E_tot1, E_tot2, E_test, p_2[UU], u_try);
//	  fprintf(stdout,"DN u_t = %20.10e ,  u^t = %20.10e , A = %20.10e , B = %20.10e \n",
//		  q.ucov[TT], q.ucon[TT], A, B); 
//	  fflush(stdout);
	}
      }
      if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
      pflag1 *= pflag[i][jp1][k];
      pflag2 *= pflag[i][jp2][k];
      pflag_out = MAX(pflag1,pflag2);
      pflag[i][jp1][k] = pflag_out;
      GLOOP { pflag[i][jm1-g][k] = pflag_out; } 
    }
  }

  /*************************************************************************************/
  /*************************************************************************************
    Calculate boundary values and then extrapolate into ghost zones: 
      -- 2-UPPER FACE:
  *************************************************************************************/
  if( bc_pid[2][BCUP] == BC_PHYS ) {
    jm2 = N2E-1;  jm1 = N2E;  jp1 = N2E+1;  jp2 = N2E+2;
    for(i=1;i<(N1TOT-1);i++) for(k=1;k<(N3TOT-1);k++)  {   
      icurr = i; jcurr = jm1; kcurr = k; 
      pflag1 = pflag2 = 0;

      PLOOP p_1[l] = p_0[l] = prim[i][jm1][k][l];
      p_0[U2] = 0.; 
      p_0[B2] = 0.; 

      l = U2;
      f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      else{ pflag1 = 1; }
      for(g=-1 ;g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

      l = B2;
      f0 = p_0[l];  f1 = prim[i][jm1][k][l];  f2 = prim[i][jm2][k][l]; 
      df10 = f1-f0;   df21 = f2-f1;  
      /* Change first real cell if it does not lead to desired recon. at boundary: */
      if( ((2*fabs(df10)) > fabs(df21)) || ((df10*df21) < 0.) ) { 
	f1 = (2*f0 + f2)/3.; 
	pflag2 = 1;
      }
      //      else{ pflag1 = 1; }
      for(g=-1; g<NG; g++) { sig_g = 2*g+1;  prim[i][jp1+g][k][l] = f0 + sig_g*(f0-f1); }

      /* If there was any change, make sure that the total energy did not increase: */
      if( pflag2 ) { 
	PLOOP p_2[l] = prim[i][jm1][k][l];
	get_geometry(i,jm1,k,CENT,ncurr,geom); 

	get_state(p_1, geom, &q );  primtoflux( p_1, &q, 0, geom, U_1);
	get_state(p_2, geom, &q );  primtoflux( p_2, &q, 0, geom, U_2);

	E_tot1 = -(U_1[UU] - U_1[RHO])/geom->g; 
	E_tot2 = -(U_2[UU] - U_2[RHO])/geom->g; 

	/* Can express   E =  A + B*u , solve for new  u  ,  note that E = -T^t_t */
	/* Assumes  P = (gam-1)*u  */
	if( E_tot2 > E_tot1 ) { 
	  A = -(p_2[RHO] + q.bsq)*q.ucon[TT]*q.ucov[TT] - 0.5*q.bsq + q.bcon[TT]*q.bcov[TT];
	  B = -( gam*q.ucon[TT]*q.ucov[TT] + gam - 1. );
	  u_try = (E_tot1 - A) / B; 

	  E_test = A + B * u_try;

	  get_coord(i,jm1,k,FACE2,ncurr,coords);
	  u_flr = UUMIN*pow(coords->x[RR],UUPOWER);

	  if( u_try < u_flr )  u_try = u_flr; 
	  
	  prim[i][jm1][k][UU] = u_try; 
//	  fprintf(stdout,"UP E_tot1 = %20.10e , E_tot2 = %20.10e , E_test = %20.10e, u_old = %20.10e , u_new = %20.10e \n",
//		  E_tot1, E_tot2, E_test, p_2[UU], u_try);
//	  fprintf(stdout,"UP u_t = %20.10e ,  u^t = %20.10e , A = %20.10e , B = %20.10e \n",
//		  q.ucov[TT], q.ucon[TT], A, B); 
//	  fflush(stdout);
	}
      }
      if( pflag2 )  fail(FAIL_CUTOUT_BC,0);
      pflag1 *= pflag[i][jm1][k];
      pflag2 *= pflag[i][jm2][k];
      pflag_out = MAX(pflag1,pflag2);
      pflag[i][jm1][k] = pflag_out; 
      GLOOP { pflag[i][jp1+g][k] = pflag_out; } 
    }
  }

  return; 
}


/**********************************************************************************/
/**********************************************************************************
  cutout_cons_interp(): 
  ---------------
    -- performs 0th-order extrapolation of the densitized conserved variables 
       into the ghost zones along the X2 boundaries; 
    -- uses 0th-order extrapolation of the primitives as guesses;
    -- uses the primitive magnetic fields and not the "conserved" one;
    -- if there is a problem inverting for the primitive variables in the ghost 
        zones, then an error message is printed and the extrapolated values 
        of the primitives are used instead. 
    --        

    -- This routine is also responsible for setting pflag[] in X2-dir appropriately.
 
**********************************************************************************/
void cutout_cons_interp( double ****prim, int type )
{
  int i,j,k,l,g,ret, j_add,j_phys,ibeg,iend;

  struct of_state q;
  struct of_geom *geom, *geom_g; 
  double pg[NG][NP], Ug[NG][NP],gtmp;

  int Utoprim_2d_fast(double *U, struct of_geom *geom, double *prim);
  int    Utoprim_1d( double  U[NP], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
		     double gdet, double prim[NP]);

  /*************************************************************************************
   LOWER X2 BOUNDARY parameters: 
  *************************************************************************************/
  if( type == BCDN ) { 
    j_phys = N2S;    /* coordinate of closest physical cell */
    j_add = 0;       /* amount to add to ghost cell's coordinate */
  }
  /*************************************************************************************
   UPPER X2 BOUNDARY parameters: 
  *************************************************************************************/
  else if( type == BCUP ) { 
    j_phys = N2E;    /* coordinate of closest physical cell */
    j_add  = N2E+1;  /* amount to add to ghost cell's coordinate */
  }
  else { 
    fprintf(stderr,"cutout_cons_interp(): Incorrect boundary type = %d \n", type ) ; 
    fflush(stderr);
    fail(FAIL_BASIC,0); 
  }

  ibeg = N1S;  iend = N1E; 
  if( bc_pid[1][BCDN] == BC_PHYS )  ibeg = 0; 
  if( bc_pid[1][BCUP] == BC_PHYS )  iend = N1TOT-1;


  /*************************************************************************************
    This is the generic boundary condition that is tailored to the UP/DN boundary by 
     the above parameters.  
  *************************************************************************************/
  for( i=ibeg; i<=iend; i++)  N3_LOOP {  /* Note that if we do not */
      get_geometry(i,j_phys,k,CENT,ncurr,geom); 
      /* Get  states of P[] and U[] of first physical cell : */
      PLOOP pg[0][l] = prim[i][j_phys][k][l]; 
      get_state(  pg[0], geom, &q );
      primtoflux( pg[0], &q, 0, geom, Ug[0]) ;
	
      /* 0th-order extrapolation to the ghosts : */
      for( g=1; g<NG; g++ )  PLOOP  { 
	  pg[g][l] = pg[0][l]; 
	  Ug[g][l] = Ug[0][l]; 
	}

      /* Invert ghost conserved var's to get final ghost prim's : */
      GLOOP { 
	j=j_add+g;
	get_geometry(i,j,k,CENT,ncurr,geom_g); 
	
	/* Interpolate non-densitized conserved variables: */
	gtmp = geom_g->g / geom->g;
	PLOOP  { 
	  Ug[g][l] *= gtmp;
	}

	ret = Utoprim_2d_fast(Ug[g], geom_g, pg[g]);
	if( ret ) { 
	  ret = Utoprim_1d(Ug[g], geom_g->gcov, geom_g->gcon, geom_g->g, pg[g]);
	}
	if( ret ) { 
	  fprintf(stderr,
		  "cutout_cons_interp(): inversion failure at (%d,%d,%d,%d,%d) : ret=%d\n",
		  i,j,k,nstep,myid,ret); 
	  fprintf(stderr, "cutout_cons_interp(): U = ");
	  PLOOP fprintf(stderr," %28.18e ", Ug[g][l]);
	  fprintf(stderr,"\n");
	  fprintf(stderr, "cutout_cons_interp(): P = ");
	  PLOOP fprintf(stderr," %28.18e ", pg[g][l]);
	  fprintf(stderr,"\n");
	}
	PLOOP { prim[i][j][k][l] = pg[g][l] ; } 
      }
    }

  return; 

}  /* End of  cutout_cons_interp():  */


