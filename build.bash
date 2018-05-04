#!/bin/bash
#
# svn $Id: build.bash 752 2015-01-07 23:01:08Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build.bash [options]                                             :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

parallel=1
clean=1

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
      ;;

    -noclean )
      shift
      clean=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
      ;;
  esac
done


export   ROMS_APPLICATION=FIEBERLING
export        MY_ROOT_DIR=/gscratch/stf/bperfect/ROMS
export     MY_PROJECT_DIR=/gscratch/stf/bperfect
export       MY_ROMS_SRC=${MY_ROOT_DIR}
export         COMPILERS=${MY_ROMS_SRC}/Compilers

export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAVERAGES"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DOUT_DOUBLE"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DPOSITIVE_ZERO"

#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDIAGNOSTICS_TS"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDIAGNOSTICS_UV"

# Set deprecated lateral boundary conditions CPP flags for backward
# compatibility with older versions of the code.


# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

export           USE_MPI=on            # distributed-memory parallelism
export        USE_MPIF90=on            # compile with mpif90 script
#export         which_MPI=mpich         # compile with MPICH library
export         which_MPI=mpich2        # compile with MPICH2 library
#export         which_MPI=openmpi       # compile with OpenMPI library
#export        USE_OpenMP=on            # shared-memory parallelism
export              FORT=ifort
#export         USE_DEBUG=on            # use Fortran debugging flags
export         USE_LARGE=on            # activate 64-bit compilation
export       USE_NETCDF4=on            # compile with NetCDF-4 library
#export   USE_PARALLEL_IO=on            # Parallel I/O with Netcdf-4/HDF5
#export       USE_MY_LIBS=on            # use my library paths below

export I_MPI_F90=ifort

export NC_CONFIG=/sw/netcdf-fortran+c-4.4.1.1_icc-17/bin/nc-config
export NETCDF_INCDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/include

if [ -n "${USE_MPIF90:+1}" ]; then
      if [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/gfortransoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/usr/include/openmpi/:$PATH
      fi
fi


if [ -n "${USE_MY_LIBS:+1}" ]; then
      export             ESMF_OS=Linux
      export       ESMF_COMPILER=gfortran
      export           ESMF_BOPT=O
      export            ESMF_ABI=64
      export           ESMF_COMM=mpich
      export           ESMF_SITE=default

      export       ARPACK_LIBDIR=/opt/gfortransoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        if [ "${which_MPI}" = "mpich2" ]; then
          export        ESMF_DIR=/opt/gfortransoft/mpich2/esmf
          export      MCT_INCDIR=/opt/gfortransoft/mpich2/mct/include
          export      MCT_LIBDIR=/opt/gfortransoft/mpich2/mct/lib
          export  PARPACK_LIBDIR=/opt/gfortransoft/mpich2/PARPACK
        elif [ "${which_MPI}" = "openmpi" ]; then
          export        ESMF_DIR=/opt/gfortransoft/openmpi/esmf
          export      MCT_INCDIR=/opt/gfortransoft/openmpi/mct/include
          export      MCT_LIBDIR=/opt/gfortransoft/openmpi/mct/lib
          export  PARPACK_LIBDIR=/opt/gfortransoft/openmpi/PARPACK
        fi
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
          if [ "${which_MPI}" = "mpich2" ]; then
            export     NC_CONFIG=/opt/gfortransoft/mpich2/netcdf4/bin/nc-config
            export NETCDF_INCDIR=/opt/gfortransoft/mpich2/netcdf4/include
          elif [ "${which_MPI}" = "openmpi" ]; then
            export     NC_CONFIG=/opt/gfortransoft/openmpi/netcdf4/bin/nc-config
            export NETCDF_INCDIR=/opt/gfortransoft/openmpi/netcdf4/include
          fi
        else
          export       NC_CONFIG=/opt/gfortransoft/serial/netcdf4/bin/nc-config
          export   NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf4/include
        fi
      else
        export     NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf3/include
        export     NETCDF_LIBDIR=/opt/gfortransoft/serial/netcdf3/lib
      fi
fi

# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

 export     MY_HEADER_DIR=${MY_PROJECT_DIR}

 export MY_ANALYTICAL_DIR=${MY_ROOT_DIR}/Functionals

# Put the binary to execute in the following directory.

 export            BINDIR=${MY_PROJECT_DIR}

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

 export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Remove build directory.

if [ $clean -eq 1 ]; then
  make clean
fi

# Compile (the binary will go to BINDIR set above).

if [ $parallel -eq 1 ]; then
  make $NCPUS
else
  make
fi
