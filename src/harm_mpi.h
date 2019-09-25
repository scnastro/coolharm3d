/* 

   Note: cannot have any cpu write to same file or pipe to same file
   without parallel I/O or blocking through each CPU

*/


#if(USEMPI==1)
#include "mpi.h"
#endif

///////////////////
//
// MNEMOICS  (NM = N for MPI)
//
//////////////////////

/* Total number of local physical cells along each direction */
#define NM1 (NM1_glob)
#define NM2 (NM2_glob)
#define NM3 (NM3_glob)

/* Number of boundary points to transfer for boundary conditions */
#define NM1_BND (NG)
#define NM2_BND (NG)
#define NM3_BND (NG)

#define NM1S (NM1_BND)
#define NM2S (NM2_BND)
#define NM3S (NM3_BND)

#define NM1E (NM1E_glob)
#define NM2E (NM2E_glob)
#define NM3E (NM3E_glob)

#define NM1_TOT (NM1_TOT_glob)
#define NM2_TOT (NM2_TOT_glob)
#define NM3_TOT (NM3_TOT_glob)

#define NM1_LOOP for(i=NM1S;i<=NM1E;i++) 
#define NM2_LOOP for(j=NM2S;j<=NM2E;j++) 
#define NM3_LOOP for(k=NM3S;k<=NM3E;k++) 

#define NM1_DN_LOOP for(i=NM1S;i<NM1S+NM1_BND;i++) 
#define NM2_DN_LOOP for(j=NM2S;j<NM2S+NM2_BND;j++) 
#define NM3_DN_LOOP for(k=NM3S;k<NM3S+NM3_BND;k++) 

#define NM1_UP_LOOP for(i=NM1_UP_BEG_glob;i<=NM1E;i++) 
#define NM2_UP_LOOP for(j=NM2_UP_BEG_glob;j<=NM2E;j++) 
#define NM3_UP_LOOP for(k=NM3_UP_BEG_glob;k<=NM3E;k++) 

/* Boundary mnemonics */
#define IBEG    (0)
#define IEND    (1)


/* max. allowed cpus in our group: */
#define MAXCPUS     (10000)
#define MAXFILENAME (1000)

/* Starting values for tags used to indicate output and boundary condition transfers  */
#define BASE_TAG_DUMP   (10*MAXCPUS)
#define BASE_TAG_BOUNDS (20*MAXCPUS)
#define BASE_TAG_HIST   (21*MAXCPUS)
#define BASE_TAG_HIST2  (22*MAXCPUS)
#define BASE_TAG_SURF   (23*MAXCPUS)

/* Filename suffix to indicate cpu ID number */
#define CPUTXT ".%04d"
