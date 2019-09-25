
/****************************************************************************************/
/****************************************************************************************
 Macros for caclulating the connection or Christoffel coeficients: 
*****************************************************************************************/
/*****************************************************************************************/

/* 2nd-order time derivative of the metric for the connection : */ 
#if(DYNAMIC_SPACETIME) 
# define CONN_2ND_ORDER_TIME_DERIVATIVE                                                                        \
   {                                                                                                           \
         k=0 ;                                                                                                 \
         struct of_coord coords_loc;                                                                           \
	 SET_GEOM_COORD_TIME_POINTERS(2);                                                                      \
         calc_coord(ii,jj,kk,coords->pos,&coords_loc);                                                         \
         gcov_func(&coords_loc,gh) ;                                                                           \
	 SET_GEOM_COORD_TIME_POINTERS(0);                                                                      \
         calc_coord(ii,jj,kk,coords->pos,&coords_loc);                                                         \
         gcov_func(&coords_loc,gl) ;                                                                           \
         for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  { conn[i][j][k] = (gh[i][j] - gl[i][j])*half_inv_dt_conn; }  \
	 SET_GEOM_COORD_TIME_POINTERS(1);                                                                      \
   }
#else
# define CONN_2ND_ORDER_TIME_DERIVATIVE	 { k=0;  for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  { conn[i][j][k] = 0.; } }
#endif



/* Often used array operations to transform metric derivatives to connection coefficients:  */
#define DGCOV_TO_CONN                                                         \
  {                                                                           \
  /* now rearrange to find \Gamma_{ijk} */                                    \
  for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++) {	              \
	tmp[i][j][k] = 0.5*(conn[j][i][k] + conn[k][i][j] - conn[k][j][i]) ;} \
    /* finally, raise index */                                                \
    for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++) {	      \
	  conn[i][j][k] = 0. ;                                                \
	  for(l=0;l<NDIM;l++) conn[i][j][k] += geom->gcon[i][l]*tmp[l][j][k]; \
	}                                                                     \
  }


/* 4th-order(or third?)  spatial derivative of the metric w.r.t.  xp1 */
#define FOURTH_ORDER_METRIC_DERIVATIVE1 \
  {                                                                                                                     \
	k = 1;                                                                                                          \
	w1 = w1_v[k];                                                                                                   \
	w2 = w2_v[k];                                                                                                   \
	get_geometry(ii-1,jj  ,kk  ,CENT ,ncurr,geom_ll);                                                               \
	get_geometry(ii  ,jj  ,kk  ,FACE1,ncurr,geom_l );                                                               \
	get_geometry(ii+1,jj  ,kk  ,FACE1,ncurr,geom_h );                                                               \
	get_geometry(ii+1,jj  ,kk  ,CENT ,ncurr,geom_hh);                                                               \
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {                                                                      \
	    conn[i][j][k] = w1*(geom_ll->gcov[i][j]-geom_hh->gcov[i][j]) + w2*(geom_h->gcov[i][j] - geom_l->gcov[i][j]);\
	}                                                                                                               \
  }

/* 4th-order(or third?)  spatial derivative of the metric w.r.t.  xp2 */
#define FOURTH_ORDER_METRIC_DERIVATIVE2 \
  {                                                                                                                     \
	k = 2;                                                                                                          \
	w1 = w1_v[k];                                                                                                   \
	w2 = w2_v[k];                                                                                                   \
	get_geometry(ii  ,jj-1,kk  ,CENT ,ncurr,geom_ll);                                                               \
	get_geometry(ii  ,jj  ,kk  ,FACE2,ncurr,geom_l );                                                               \
	get_geometry(ii  ,jj+1,kk  ,FACE2,ncurr,geom_h );                                                               \
	get_geometry(ii  ,jj+1,kk  ,CENT ,ncurr,geom_hh);                                                               \
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {                                                                      \
	    conn[i][j][k] = w1*(geom_ll->gcov[i][j]-geom_hh->gcov[i][j]) + w2*(geom_h->gcov[i][j] - geom_l->gcov[i][j]);\
	}                                                                                                               \
  }

/* 4th-order(or third?)  spatial derivative of the metric w.r.t.  xp3 */
#define FOURTH_ORDER_METRIC_DERIVATIVE3 \
  {                                                                                                                     \
	k = 3;                                                                                                          \
	w1 = w1_v[k];                                                                                                   \
	w2 = w2_v[k];                                                                                                   \
	get_geometry(ii  ,jj  ,kk-1,CENT ,ncurr,geom_ll);                                                               \
	get_geometry(ii  ,jj  ,kk  ,FACE3,ncurr,geom_l );                                                               \
	get_geometry(ii  ,jj  ,kk+1,FACE3,ncurr,geom_h );                                                               \
	get_geometry(ii  ,jj  ,kk+1,CENT ,ncurr,geom_hh);                                                               \
	for(i=0;i<NDIM;i++) for(j=0;j<NDIM;j++)  {                                                                      \
	    conn[i][j][k] = w1*(geom_ll->gcov[i][j]-geom_hh->gcov[i][j]) + w2*(geom_h->gcov[i][j] - geom_l->gcov[i][j]);\
	}                                                                                                               \
  }

