#==========================================================================================#
#==========================================================================================#
#    Makefile include.intel.mk                                                             #
#                                                                                          #
#    Compilation controls optimised for Gent Univeristy Cluster                            #
#------------------------------------------------------------------------------------------#

#----- Define make (gnu make works best). -------------------------------------------------#
MAKE=/usr/bin/make
#------------------------------------------------------------------------------------------#


#----- Main path for compilation. ---------------------------------------------------------#
BASE=$(ED_ROOT)/build/
#------------------------------------------------------------------------------------------#


#------ Detect current system. ------------------------------------------------------------#
UNAME_S := $(shell uname -s)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    HDF 5 libraries.                                                                      #
#                                                                                          #
#    Since ED-2.1, this is no longer optional for real simulations.  You must have the     #
# HDF5 libraries compiled with the same compiler you set for F_COMP and C_COMP.  You may   #
# still be able to compile without HDF5 but the code is not going to run.                  #
#------------------------------------------------------------------------------------------#
USE_HDF5=1
ifeq ($(UNAME_S),Linux)
	HDF5_LIBS=-lz -lhdf5 -lhdf5_fortran -lhdf5_hl
endif
ifeq ($(UNAME_S),Darwin)
	HDF5_INCS=-I/usr/local/hdf5_mio/include
	HDF5_LIBS=-lz -L/usr/local/hdf5_mio/lib -lhdf5 -lhdf5_fortran -lhdf5_hl
endif
#------------------------------------------------------------------------------------------#



#################################### COMPILER SETTINGS #####################################
CMACH=PC_LINUX1
F_COMP=ifort
C_COMP=icc
LOADER=ifort
MOD_EXT=mod
#------------------------------------------------------------------------------------------#



##################################### COMPILER OPTIONS #####################################
#------------------------------------------------------------------------------------------#
# A. Pickiest - Use this whenever you change arguments on functions and subroutines.       #
#             This will do all possible checks. Interface is no longer an issue in this    #
#             version of the code. double compilation is no longer needed).                #
# B. Profiler - This is a debug build with flags for profiling. To produce the coverage    #
#             run the application, then in the bin folder execute profmerge (this will     #
#             merge the *.dyn files. Then run                                              #
#             codecov -prj <build/project_name> -spi <bin/file.spi> -dpi <bin/file.dpi>    #
# C. Fast     - This is all about performance, use only when you are sure that the model   #
#             has no code problem, and you want results asap. This will not check for any  #
#             problems, which means that this is an option suitable for end users, not de- #
#             velopers.                                                                    #
#------------------------------------------------------------------------------------------#
ifeq ($(KIND_COMP),)
	KIND_COMP=C
endif
#------------------------------------------------------------------------------------------#
#################################    DEBUG BUILD   #########################################
ifeq ($(KIND_COMP),A)
	F_OPTS= -check -g -debug extended -debug-parameters -traceback -u -fp-stack-check -warn
	C_OPTS= -g -traceback -debug extended -debug-parameters -warn
	#---------------------------------------------------------------------------------------#
endif
########################   DEBUG BUILD WITH INTEL PROFILING   ###############################
ifeq ($(KIND_COMP),B)
#	F_OPTS= -check -g -debug extended -debug-parameters -traceback -u -fp-stack-check -warn \
#			-prof-gen=srcpos
#	C_OPTS= -g -traceback -debug extended -debug-parameters -warn -prof-gen=srcpos
	#---------------------------------------------------------------------------------------#
endif
######################################   OPTIMIZED BUILD   ##################################
ifeq ($(KIND_COMP),C)
#	F_OPTS= -O3 -xHost -g -u -qopenmp -guide -qopt-report -parallel
#	C_OPTS= -O3 -xHost -g -qopenmp -guide -qopt-report -parallel
#	F_LOWO_OPTS= -O1 -xHost -g -u -qopenmp -guide -qopt-report -parallel
	F_OPTS= -O3 -xHost -g -u
	C_OPTS= -O3 -xHost -g
	F_LOWO_OPTS= -O2 -xHost -g -u
	#---------------------------------------------------------------------------------------#
endif
#-------------------------------------------------------------------------------------------#
############################################################################################



#------------------------------------------------------------------------------------------#
#     If you have a version of hdf5 compiled in parallel, then you may benefit from        #
# collective I/O, then use this flag = 1.  Otherwise, set it to zero.                      #
#------------------------------------------------------------------------------------------#
USE_COLLECTIVE_MPIO=0
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Archive options.                                                                     #
#------------------------------------------------------------------------------------------#
ifeq ($(UNAME_S),Linux)
	ARCHIVE=xiar crs
endif
ifeq ($(UNAME_S),Darwin)
	ARCHIVE=xilibtool -c -static -o
endif
#------------------------------------------------------------------------------------------#


