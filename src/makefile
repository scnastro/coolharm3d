.PHONY : failure all clean newrun vclean runit runit2 \
	emacs emn tags newtest archive redo testinterp

#***********************************************************************************
#### Set USEMPI to 1 if you want use MPI compiler options. 
###### These options will undoubtedly be system dependent, so beware. 
#***********************************************************************************
USEMPI = 1

#***********************************************************************************
#***********************************************************************************
#### If mpi, we need to set cluster-specific compiler options:
###### USECLUSTER = 0 (Default settings),  1 (TACC/Stampede), 2 (RIT workstation) , 3 (RIT/newhorizons) , 4 (Gauss/JHU), 5 (BlueWaters/NCSA), 6 (OpenSuse 13.1 Laptop), 7 (RIT/bluesky - sandy), 8 (RIT/bluesky - nehalem), 9 (Datura/AEI), 10 (gordon), 11 (Quasar), 12 (tandy/Tulsa), 13 (Mac laptop)
#***********************************************************************************
USECLUSTER = 5


#***********************************************************************************
####  COMPILER :  0 (gcc)   1  (icc)    2  (pgcc)  
#***********************************************************************************
COMPILER = 0

#####  set MAKESDFS to 0 if you do not want to use the sdf libraries:
MAKESDFS = 0

#####  set MAKEHDFS to 0 if you do not want to use the hdf5 libraries:
MAKEHDFS = 1

#####  set USE_GSL to 0 if you do not want to use the GSL libraries:
USE_GSL = 0


#***********************************************************************************
#***********************************************************************************
# Basic compiler options and settings:
#
#  OPTIONS: Set the makefile variables that are to be sent to the
#            pre-compiler: If you add an option above, please add
#            it here so that you can use the flag as a
#            pre-compiler macro in the source code.
#
#***********************************************************************************
#***********************************************************************************
OPTIONS = 
OTHER_CFLAGS = 
OTHER_LFLAGS = 
EXTRA_LIBS = -lm

#***********************************************************************************
# Intel compiler options
#***********************************************************************************
ifeq ($(COMPILER),1)
COMP       = icc
#CCFLAGS  = -g 
CCFLAGS = -O3 -DINTEL_HEADER 
OTHER_CFLAGS +=  -I$(INTEL_INC_DIR)
OTHER_LFLAGS +=  -L$(INTEL_LIB_DIR)
endif

#***********************************************************************************
# PGI  compiler options
#***********************************************************************************
ifeq ($(COMPILER),2)
COMP       = pgcc
endif

#***********************************************************************************
# gcc  compiler options
#***********************************************************************************
ifeq ($(COMPILER),0)
COMP       = gcc
CCFLAGS  = -O3
#CCFLAGS  = -pg
#CCFLAGS = -O0 -g
endif

#***********************************************************************************
#***********************************************************************************
# Specialized compiler options, flags, libraries, headers:
#***********************************************************************************
#***********************************************************************************

#***********************************************************************************
# MPI libraries: MPI_HOME needs to be set by environment or below in
# cluster-specific sections:
#***********************************************************************************
ifeq ($(USEMPI),1)

ifndef MPI_HOME
ifdef MPI_DIR
MPI_HOME=$(MPI_DIR)
else 
ifeq ($(USE_MODULE_ENVIRONMENT),0)
failure:
	@echo "MPI_HOME needs to be set in user's environment or in cluster-specific environment in the makefile"
endif
endif
endif

OPTIONS += -DUSEMPI
MPI_INCLUDE   = $(MPI_HOME)/include
MPI_LIB       = $(MPI_HOME)/lib
MPI_BIN       = $(MPI_HOME)/bin
OTHER_CFLAGS += -I$(MPI_INCLUDE)
OTHER_LFLAGS += -L$(MPI_LIB)
endif

#***********************************************************************************
# sdf libraries
#***********************************************************************************
ifeq ($(MAKESDFS),1)
OTHER_CFLAGS +=  -I${HOME}/local/include
OTHER_LFLAGS +=  -L${HOME}/local/lib 
EXTRA_LIBS += -lrnpl -lxvs -lsvs
endif

#***********************************************************************************
# Hdf5 libraries
#***********************************************************************************
ifeq ($(MAKEHDFS),1)

ifndef HDF_HOME
ifdef HDF_DIR
HDF_HOME=$(HDF_DIR)
else 
ifeq ($(USE_MODULE_ENVIRONMENT),0)
failure:
	@echo "HDF_HOME needs to be set in user's environment or in cluster-specific environment in the makefile"
endif
endif
endif

OPTIONS      += -DH5_USE_16_API  
EXTRA_LIBS   += -lz -lhdf5_hl -lhdf5
HDF_INCLUDE   = $(HDF_HOME)/include
HDF_LIB       = $(HDF_HOME)/lib
HDF_BIN       = $(HDF_HOME)/bin
OTHER_CFLAGS +=  -I$(HDF_INCLUDE)
OTHER_LFLAGS +=  -L$(HDF_LIB)

endif

#***********************************************************************************
# GSL libraries:  GSL_HOME needs to be set by environment
#***********************************************************************************
ifeq ($(USE_GSL),1)

ifndef GSL_HOME
ifdef GSL_DIR
GSL_HOME=$(GSL_DIR)
else 
ifeq ($(USE_MODULE_ENVIRONMENT),0)
failure:
	@echo "GSL_HOME needs to be set in user's environment or in cluster-specific environment in the makefile"
endif
endif
endif

OPTIONS      += -DUSE_GSL
EXTRA_LIBS   += -lgsl -lgslcblas -lm
GSL_INCLUDE   = $(GSL_HOME)/include
GSL_LIB       = $(GSL_HOME)/lib
OTHER_CFLAGS += -I$(GSL_INCLUDE)
OTHER_LFLAGS += -L$(GSL_HOME)
endif


#***********************************************************************************
#***********************************************************************************
# Cluster/system-specific changes:  
#   -- Override previously set options/variables here):
#   -- Only specify overrides here:  if the default is fine, then keep it as is. 
#   -- If the defaults are fine for a given cluster, please still 
#      
#  -- set USE_MODULE_ENVIRONMENT=1 if your system relies on "modules"
#       to setup the environment, otherwise set 
#	USE_MODULE_ENVIRONMENT=0  if you have to manually specify
#       header and library paths. 
# 
#***********************************************************************************
CC = $(COMP) 

USE_MODULE_ENVIRONMENT = 0 

# Default options -- you do not have to set anything if the above options work
ifeq ($(USECLUSTER),0) 
endif

# Stampede cluster at TACC
ifeq ($(USECLUSTER),1) 
USE_MODULE_ENVIRONMENT = 1
ifeq ($(USEMPI),1)
ifeq ($(MAKEHDFS),1)
CC = h5pcc
else 
CC = mpicc
endif
endif
CCFLAGS  =   -O2 -xhost 
HDF_INC = $(TACC_HDF5_INC)
HDF_LIB = $(TACC_HDF5_LIB)
MPI_HOME = $(MPICH_HOME)
endif

# RIT  Workstation
ifeq ($(USECLUSTER),2) 
ifeq ($(USEMPI),1)
MPI_HOME=/opt/mpich-1.2.7p1
COMP=$(MPI_BIN)/mpicc
HDF5_CC=$(COMP)
HDF5_CLINKER=$(COMP)
CC=$(HDF_BIN)/h5cc  
else
HDF5_CC=$(COMP)
HDF5_CLINKER=$(HDF5_CC)
CC=$(HDF_BIN)/h5cc  
endif
endif

# NewHorizons cluster at RIT 
ifeq ($(USECLUSTER),3) 
ifeq ($(USEMPI),1)
MPI_HOME=/usr/mpi/qlogic
MPI_LIB=$(MPI_HOME)/lib64
CC=mpicc -cc=icc
CCFLAGS  =  -O2 -xhost 
endif
endif

# Gauss at JHU
ifeq ($(USECLUSTER),4) 
ifeq ($(USEMPI),1)
HDF5_CC = mpicc     
HDF5_CLINKER = $(HDF5_CC)
MCC=h5pcc
CC=$(MCC) 
endif
endif

# BlueWaters at NCSA
ifeq ($(USECLUSTER),5)
USE_MODULE_ENVIRONMENT = 1
ifeq ($(USEMPI),1)
#Cray version:  (max. optimizations for static code that does not seg fault):
#CCFLAGS  =  -h cache2,scalar0,vector0,ipa3,fp2,thread0 -D H5_USE_16_API
#PGI version:   
#CCFLAGS  = -fast -D H5_USE_16_API
CC=cc
endif
endif

# OpenSuse 13.1 laptop
ifeq ($(USECLUSTER),6) 
ifeq ($(USEMPI),1)
CC=mpicc
endif
MPI_HOME=/usr/lib64/gcc/openmpi
MPI_LIB=$(MPI_HOME)/lib64
HDF_HOME=/usr
HDF_LIB=$(HDF_HOME)/lib64
CCFLAGS= -O2 
OTHER_LFLAGS +=  -lmpi
endif

# bluesky cluster at RIT - sandy
ifeq ($(USECLUSTER),7) 
ifeq ($(USEMPI),1)
MPI_HOME=/usr/mpi/gcc/openmpi-1.4.3-qlc
MPI_LIB=$(MPI_HOME)/lib64
CC=icc
#CCFLAGS=-O2 -D H5_USE_16_API  
CCFLAGS = -xAVX -align  
OTHER_LFLAGS +=  -Wl,--export-dynamic,--allow-shlib-undefined  -lmpi
endif
endif

# bluesky cluster at RIT - nehalem
ifeq ($(USECLUSTER),8) 
ifeq ($(USEMPI),1)
MPI_HOME=/usr/mpi/gcc/openmpi-1.4.3-qlc
MPI_LIB=$(MPI_HOME)/lib64
CC=icc
CCFLAGS= -xSSE4.2 -align  
OTHER_LFLAGS += -Wl,--export-dynamic,--allow-shlib-undefined  -lmpi
#CC=mpicc -cc=icc
#CCFLAGS=-g -traceback -debug all -xSSE4.2 -align  
endif
endif

# datura at AEI
ifeq ($(USECLUSTER),9)
ifeq ($(USEMPI),1)
COMP=$(MPI_BIN)/mpicc
endif
HDF5_CC=$(COMP)
HDF5_CLINKER=$(COMP)
CC=$(HDF_BIN)/h5cc   
endif


# gordon
ifeq ($(USECLUSTER),10) 
ifeq ($(USEMPI),1)
MPI_HOME=/opt/mvapich2/intel/ib
HDF_HOME=/opt/hdf5/intel/mvapich2/ib
CC=mpicc
CCFLAGS=-g -traceback -debug all -xSSE4.2 -align 
endif
endif

# quasar
ifeq ($(USECLUSTER),11)
ifeq ($(USEMPI),1)
COMP=$(MPI_BIN)/mpicc
OTHER_CFLAGS +=  -Wl,-rpath -Wl,$(MPI_LIB) -Wl,--enable-new-dtags -lmpi -lrt -lpthread
OTHER_LFLAGS +=  -Wl,-rpath -Wl,$(MPI_LIB) -Wl,--enable-new-dtags -lmpi -lrt -lpthread
else
OTHER_CFLAGS +=  -Wl,-rpath -Wl,$(MPI_LIB) -Wl,--enable-new-dtags -lmpi -lrt -lpthread
OTHER_LFLAGS +=  -Wl,-rpath -Wl,$(MPI_LIB) -Wl,--enable-new-dtags -lmpi -lrt -lpthread
endif
HDF5_CC=$(COMP)
HDF5_CLINKER=$(COMP)
CC=$(HDF_BIN)/h5cc  
endif

# Tandy supercomputer cluster
ifeq ($(USECLUSTER),12)
ifeq ($(USEMPI),1)
COMP=mpicc
HDF5_CC=$(COMP)
HDF5_CLINKER=$(HDF5_CC)
CC=$(COMP)   
OTHER_CFLAGS += -lmpi 
OTHER_LFLAGS += -lmpi 
endif
endif

# Mac laptop
ifeq ($(USECLUSTER),13)
ifeq ($(USEMPI),1)
COMP=/opt/local/bin/mpicc-openmpi-mp
HDF5_CC=$(COMP)
HDF5_CLINKER=$(HDF5_CC)
CC=$(COMP)
OTHER_CFLAGS += -D H5_USE_16_API -I/opt/local/include -I/usr/include 
OTHER_LFLAGS += -L/opt/local/lib -L/usr/lib  -lhdf5_hl -lhdf5 -D H5_USE_16_API
endif
endif

#***********************************************************************************
#***********************************************************************************
# Begin of compiler and option independent section
#***********************************************************************************
#***********************************************************************************
CCCFLAGS = $(OPTIONS) $(OTHER_CFLAGS)  -c
CCLFLAGS = $(OPTIONS) $(OTHER_LFLAGS)  $(EXTRA_LIBS)

CC_COMPILE  = $(CC) $(CCFLAGS) $(CCCFLAGS) 
CC_LOAD     = $(CC) $(CCFLAGS) $(CCLFLAGS) 

EXE = harm3d
MACH = mach.list 

SRCS = \
bounds.c coord.c diag.c dump.c fixup.c flux_ct.c image.c \
init.gen_axi_disk.c lu.c main.c metric.c phys.c \
ranc.c restart.c recon.c step_ch.c \
u2p_util.c utoprim_1d_ee.c utoprim_1d_ee2.c \
utoprim_1dfix1.c utoprim_1dvsq2fix1.c \
utoprim_1d.c utoprim_2d_fast.c harm_mpi.c \
dump_hdf.c dump_history.c dump_stat.c dump_stat2.c \
dump_radflux.c dump_surface.c dump_timers.c decomp.c dump_photosphere.c \
interp_data.c  corona.c  corona2.c

INCS = decs.h  defs.h  metric.h u2p_defs.h u2p_util.h harm_mpi.h recon.h

BASE = $(basename $(SRCS) )
OBJS = $(addsuffix .o, $(BASE) )


#***********************************************************************************
# Dependencies and Make Rules
#***********************************************************************************
.c.o: 
	$(CC_COMPILE) $*.c

all: $(EXE) 

defs.h  : globals.h 
	egrep "^extern|^#if|^#elif|^#else|^#endif" globals.h | sed s/^extern//  > defs.h


$(OBJS) : $(INCS) makefile

$(EXE): $(OBJS) 
	$(CC_LOAD) $(OBJS) $(EXTRA_LIBS) -o $(EXE)


clean:
	/bin/rm -f *.o 
	/bin/rm -f *.oo 
	/bin/rm -f *.il 
	/bin/rm -f defs.h 
	/bin/rm -f TAGS
	/bin/rm -f $(EXE) 


newrun:
	/bin/rm -rf dumps
	/bin/rm -rf images
	/bin/rm -f  ener.out
	/bin/rm -f  *.sdf
	/bin/rm -rf history
	/bin/rm -f *.h5
	/bin/rm -f *.out 

vclean: 
	make newrun
	make clean

runit: 
	make newrun
	make clean
	make
	./harm3d > output_log.dat 2>&1 

runit2: 
	nohup mpirun -np 8 ./harm3d 2 4 1 256 128 1 > output_log.dat 2>&1 &


emacs: 
	emacs *.c *.h README CHANGES makefile  &

emn: 
	emacs -nw *.c *.h README CHANGES makefile

tags: $(INCS) $(SRCS)
	etags $(INCS) $(SRCS)

newtest:
	make newrun
	make clean
	make > make.out 2>&1
	./harm3d 1 1 1 128 128  1  > output_log.dat 2>&1  & 
	ln -s . dumps

archive:
	tar zcvf code.tgz *.c *.h makefile README* CHANGES 

redo:
	make newrun
	make clean
	make > make.out 2>&1
	cp ../cart/KDHARM0.000000.h5 rdump_start_other.h5 
	./harm3d  1 1 1 256 1 256 > output_log.dat  2>&1  & 
	ln -s . dumps


##  This is used as an example for setting up environment variable for interpolating calculation:
testinterp:
	cp ../../diag3/KDHARM0.000000.h5 ../diag3
	make newrun
	make clean
	make > make.out 2>&1
	export DEST_INTERP_FILENAME="../diag3/KDHARM0.000000.h5" ; mpirun -n 4 ./harm3d  2 1 2 64 1 64  > output_log.dat  2>&1  
	ln -s . dumps


#cleanmpi: 

