
static const double mink[NDIM][NDIM] = {{-1.,0.,0.,0.},
                                        { 0.,1.,0.,0.},
					{ 0.,0.,1.,0.},
					{ 0.,0.,0.,1.} };

static const int    idel[NDIM][NDIM] = {{ 1,0,0,0},
                                        { 0,1,0,0},
					{ 0,0,1,0},
					{ 0,0,0,1} };

/*Coefficient distribution array (coeff_d2[i][j]):

  -1/144   1/18   -1/12   -1/18   1/144
   1/18   -4/9     4/3     4/9   -1/18
  -1/12    4/3    -5/2     4/3   -1/12
  -1/18    4/9     4/3    -4/9    1/18
   1/144  -1/18   -1/12    1/18  -1/144

Factorized:

/                          \
|  -1    8   -12   -8   1  |
|   8  -64   192   64  -8  |
| -12  192  -360  192 -12  | x 1/144
|  -8   64   192  -64   8  |
|   1   -8   -12    8  -1  |
\                          /

*/

static const double   coeff_d2[5][5] = {
                                        {  -1.0,   8.0, -12.0,  -8.0,   1.0},
                                        {   8.0, -64.0, 192.0,  64.0,  -8.0},
                                        { -12.0, 192.0,-360.0, 192.0, -12.0},
                                        {  -8.0,  64.0, 192.0, -64.0,   8.0},
                                        {   1.0,  -8.0, -12.0,   8.0,  -1.0}
                                       };


extern int invert_matrix( double Am[][NDIM], double Aminv[][NDIM], double *f )  ;
extern int invert_matrix2( double Am[][NDIM], double Aminv[][NDIM], double *f )  ;
extern void  calc_coord( int i, int j, int k, int pos, struct of_coord *coords);
extern void transform_all(struct of_coord *coords, struct of_geom *geom);


static double w1_v[NDIM];
static double w2_v[NDIM];

static double w_2lo[5], w_2hi[5]; 

/* BBH trajectory data structure in order to resolve differences between the various definitions of "nz_params" : */
struct of_bbh_traj { 
  double xi1x;
  double xi1y;
  double xi2x;
  double xi2y;
  double v1x;
  double v1y;
  double v2x;
  double v2y;
  double v1;
  double v2;
  double v12x;
  double v12y;
  double v21x;
  double v21y;
  double v12;
  double v21;
  double v1v2; // dot product of v1 and v2
  double v2v1; // dot product of v2 and v1 (different from v1v2 due to PN approximation)
  double t_c;  // time to merger from t=0
  double phi;  // orbital phase change from t=0
  double omega;// current orbital phase rate of change
  double r12;  // current separation
  double r21;  // current separation (from BH2's perspective)
  double r12dot; // for now, this element is only set when using
                 // METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW 
  double xi1z;   //  
  double xi2z;   //  
  double v1z ;   //  
  double v2z ;   //  
  double v12z;   //  
  double v21z;   //  
  double r1T0;   // Inner1-Near radius trans. func. parameter 
  double w1T0;   // Inner1-Near width trans. func. parameter 
  double r2T0;   // Inner2-Near radius trans. func. parameter 
  double w2T0;   // Inner2-Near width trans. func. parameter 
  double xNT0;   // Near-Inner radius trans. func. parameter 
  double wNT0;   // Near-Inner width trans. func. parameter 
  double lambda; // Near-Far trans. func. parameter 
  double tt;
  double n12x; 
  double n12y; 
  double n12z; 
  double n12v12; 
  double n12v1; 
  double n12v2; 
};

/******************************************************************************
  general_gcov_func():    -- special gcov routine:
******************************************************************************/
extern void null_time_funcs_setup(double t);
extern void dyn_fullpn_nz_gcov_func(double *xx, double gcov[][NDIM]);


//#define gcon_func invert_matrix 
#define gcon_func invert_matrix2

#if(   METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_35PN_NZ       )
  extern void dyn_35pn_nz_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   null_time_funcs_setup
  #define general_gcov_func  dyn_35pn_nz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ     )
  #define time_funcs_setup   null_time_funcs_setup
  #define general_gcov_func  dyn_fullpn_nz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_SIMPLE_NEWTONIAN   )
  extern  void dyn_simple_newtonian_gcov_func_setup(double t);
  extern  void dyn_simple_newtonian_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_simple_newtonian_gcov_func_setup
  #define general_gcov_func  dyn_simple_newtonian_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_FAST)
  extern void dyn_fullpn_opt_nz_gcov_func_setup(double t);
  extern void dyn_fullpn_opt_nz_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_fullpn_opt_nz_gcov_func_setup 
  #define general_gcov_func  dyn_fullpn_opt_nz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ)
  extern  void dyn_fullpn_nz_iz_gcov_func_setup(double t);
  extern  void dyn_fullpn_nz_iz_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_fullpn_nz_iz_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nz_iz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ || \
       METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW )
  extern  void dyn_fullpn_nz_iz_fz_gcov_func_setup(double t);
  extern  void dyn_fullpn_nz_iz_fz_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_fullpn_nz_iz_fz_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nz_iz_fz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND)
  extern  void dyn_fullpn_nz_iz_2nd_gcov_func_setup(double t);
  extern  void dyn_fullpn_nz_iz_2nd_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_fullpn_nz_iz_2nd_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nz_iz_2nd_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND)
  extern  void dyn_fullpn_nz_iz_fz_2nd_gcov_func_setup(double t);
  extern  void dyn_fullpn_nz_iz_fz_2nd_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_fullpn_nz_iz_fz_2nd_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nz_iz_fz_2nd_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_DROP)
  extern  void dyn_fullpn_nz_gcov_func_setup(double t);
  #define time_funcs_setup   dyn_fullpn_nz_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NR_FZ)
  extern  void dyn_fullpn_nr_fz_gcov_func_setup(double t);
  extern  void dyn_fullpn_nr_fz_gcov_func(double *xx, double gcov[][NDIM]);
  #define time_funcs_setup   dyn_fullpn_nr_fz_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nr_fz_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_KS_SPHERICAL  )
  #define time_funcs_setup   null_time_funcs_setup
  #define general_gcov_func  ks_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_MINK_CARTESIAN )
  #define time_funcs_setup   null_time_funcs_setup
  #define general_gcov_func  mink_cartesian_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_MINK_SPHERICAL )
  #define time_funcs_setup   null_time_funcs_setup
  #define general_gcov_func  mink_spherical_gcov_func
#elif( METRIC_DYNAMIC_TYPE_CHOICE  == METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP)
  extern  void dyn_fullpn_nz_iz_fz_aspin_gcov_func_setup(double t);
//#undef  time_funcs_setup
  #define time_funcs_setup   dyn_fullpn_nz_iz_fz_aspin_gcov_func_setup
  #define general_gcov_func  dyn_fullpn_nz_iz_fz_aspin_gcov_func
#else
   invalid-value-of-METRIC_DYNAMIC_TYPE_CHOICE
#endif

#if((COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL)|| (COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD))
# define ks_conn_func ks_conn_func_2
#else
# define ks_conn_func ks_conn_func_1
#endif


#define mink_cartesian_gcon_func  mink_cartesian_gcov_func


#if( METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND \
  || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND )

  #define SET_IZ_PARAMS(_iparams) {			\
     iz1_params = iz1_params_all[(_iparams)];		\
     iz2_params = iz2_params_all[(_iparams)];		\
     iz1_A2_Et = iz1_A2_Et_all[(_iparams)];		\
     iz2_A2_Et = iz2_A2_Et_all[(_iparams)];		\
     iz1_A2_Etdot = iz1_A2_Etdot_all[(_iparams)];	\
     iz2_A2_Etdot = iz2_A2_Etdot_all[(_iparams)];	\
     iz1_A2_Etdot0 = iz1_A2_Etdot0_all[(_iparams)];	\
     iz2_A2_Etdot0 = iz2_A2_Etdot0_all[(_iparams)];	\
     iz1_A3_Et = iz1_A3_Et_all[(_iparams)];		\
     iz2_A3_Et = iz2_A3_Et_all[(_iparams)];		\
     iz1_A3_Ct = iz1_A3_Ct_all[(_iparams)];		\
     iz2_A3_Ct = iz2_A3_Ct_all[(_iparams)];		\
     iz1_A3_Ctdot = iz1_A3_Ctdot_all[(_iparams)];	\
     iz2_A3_Ctdot = iz2_A3_Ctdot_all[(_iparams)];	\
     iz1_A3_Ctdot0 = iz1_A3_Ctdot0_all[(_iparams)];	\
     iz2_A3_Ctdot0 = iz2_A3_Ctdot0_all[(_iparams)];	\
     iz1_A4_Ct = iz1_A4_Ct_all[(_iparams)];		\
     iz2_A4_Ct = iz2_A4_Ct_all[(_iparams)];		\
   }

#else 
     #define SET_IZ_PARAMS(_iparams) {}
#endif


   /* For your new dynamic metric, tell us how to set the local pointer to the static memory: */
#if( METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )
  #define SPECIAL_TIME_PARAM_POINTER(_iparams)  { \
     bbh_params = &(bbh_params_all[(_iparams)]);  \
   }
#else 
  #define SPECIAL_TIME_PARAM_POINTER(_iparams)  { \
     nz_params  = nz_params_all[(_iparams)];	  \
   }
#endif


   /* The following define statements represent what is done to set
      the time-dependent-only functions for the metric and
      coordinates.  It performs all the operations for setting up a
      new in time of coordinates and metric.  Note that all these
      commands should be permissible even if we are not running with
      DYNAMIC_COORDINATES or DYNAMIC_SPACETIME. */
#define SET_GEOM_ONLY_TIME_FUNCS(_tloc,_iparams) {      \
     SET_IZ_PARAMS((_iparams));				\
     SPECIAL_TIME_PARAM_POINTER((_iparams));            \
     time_funcs_setup((_tloc));				\
   }

#define SET_GEOM_COORD_TIME_POINTERS(_iparams) {	\
     SET_IZ_PARAMS((_iparams));				\
     SPECIAL_TIME_PARAM_POINTER((_iparams));            \
     coord_params = &(coord_params_all[(_iparams)]);	\
   }

#define SET_GEOM_COORD_TIME_FUNCS(_tloc,_iparams) {	\
     SET_GEOM_COORD_TIME_POINTERS((_iparams));		\
     time_funcs_setup((_tloc));				\
     update_coord_time_funcs((_tloc));			\
   }

