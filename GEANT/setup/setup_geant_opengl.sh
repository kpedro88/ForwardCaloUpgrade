#/bin/csh -v
#  Proto-script to setup environment to use CERN /afs binary releases
#
#  First version: 12 April 2002, J.A.
#  Modified:      24 May   2002, J.A.
#  Modified:      16 June  2003, I.M. 
#  Modified:      18 March 2004, I.M. 
#  Version 0.4.0
#                                       Eventual CVS ver no: $Id: setup_geant_temp.csh,v 1.2 2011/12/16 16:01:40 eno Exp $
# 

echo "Configuration script for Geant4 binary releases at CERN, Version 1.2 December 2005"
#*********************************************************************
#* start customization area
#*********************************************************************
# Configuration parameters, for release
#
#g4version="geant4.9.4.p02"
#g4releases=/data/users/eno/geant
#clhep=/data/users/eno/clhep
#XercesC=/data/users/eno/xerces
g4version="geant4.9.4.p02-openGL"
g4releases=/data/users/jtemple/SarahGeant/geant
clhep=/data/users/jtemple/SarahGeant/clhep
XercesC=/data/users/jtemple/SarahGeant/xercesC
g4data=/data/users/jtemple/SarahGeant/geant/data/

# Configuration parameters, for system and compiler
#
gccversion="4.3.4"
os="x86_64-slc5-gcc43"
EXTRALIBS=-L/usr/X11R6/lib64
export EXTRALIBS
clhepversion="x86_64-slc5-gcc43-opt"
XercesCversion="xerces-c-3.1.1-x86_64-linux-gcc-3.4"


#*********************************************************************
#* end customization area
#*********************************************************************

g++ --version | grep $gccversion > /dev/null
echo `g++ --version` | grep $gccversion
if [ $? != 0 ]
then
  echo "It looks like your compiler settings are not suitable"
  echo "The Operating system is expected to be $os"
  echo    "The compiler version should be g++ (GCC) $gccversion"
  echo -n "The system reports that it is  "; g++ --version
  echo "Please your PATH and LD_LIBRARY_PATH environment variables"
  echo "You may use the setup script "
  echo " source /afs/cern.ch/sw/lcg/contrib/gcc/4.3.2/x86_64-slc5-gcc43-opt/setup.sh"
  echo "  to set your environment for this compiler"

else

  echo "Setting up the environment for $g4version"

  export G4SYSTEM=Linux-g++
  export G4INSTALL=$g4releases/$g4version
  export G4LIB=$g4releases/$g4version/$os/lib
#  export CLHEP_BASE_DIR=$clhep/$clhepversion/${os}-opt
  export CLHEP_BASE_DIR=$clhep/$clhepversion
#  export XERCESCROOT=$XercesC/$XercesCversion/${os}-opt
  export XERCESCROOT=$XercesC/$XercesCversion

  export G4ABLADATA=$g4data/G4ABLA3.0
  export G4LEDATA=$g4data/G4EMLOW6.19
  export G4NEUTRONHPDATA=$g4data/G4NDL3.14
  export G4LEVELGAMMADATA=$g4data/PhotonEvaporation2.1
  export G4RADIOACTIVEDATA=$g4data/RadioactiveDecay3.3
  export G4ELASTICDATA=$g4data/G4ELASTIC
  export G4REALSURFACEDATA=$g4data/RealSurface1.0
  export G4NEUTRONXSDATA=$g4data/G4NEUTRONXS1.0
  export G4PIIDATA=$g4data/G4PII1.2

  # Geant 4 interface, visualisation and other variables
  export G4UI_USE_TERMINAL=1
  export G4UI_USE_TCSH=1
  export G4UI_USE_GAG=1
  export G4UI_USE_XAW=1
#  export G4UI_USE_XM=1
  #
  export G4VIS_USE_DAWN=1
  export G4VIS_USE_DAWNFILE=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_USE_OPENGLXM=1
  export G4VIS_USE_RAYTRACER=1
  export G4VIS_USE_RAYTRACERX=1
  export G4VIS_USE_VRML=1
  export G4VIS_USE_VRMLFILE=1
  #
  # Geant 4 build variables
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4UI_BUILD_XAW_SESSION=1
#  export G4UI_BUILD_XM_SESSION=1
  export G4LIB_BUILD_G3TOG4=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_BUILD_RAYTRACERX_DRIVER=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4UI_BUILD_TERMINAL_SESSION=1
  export G4UI_BUILD_GAG_SESSION=1
  export G4VIS_BUILD_RAYTRACER_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
  export G4LIB_BUILD_GDML=1

  export G4LIB_USE_G3TOG4=1
  # The following is used to store your executables 
  #            (in subdirectories, one per system)

  if [ $LD_LIBRARY_PATH ]
  then
     LD_LIBRARY_PATH=${G4LIB}:${CLHEP_BASE_DIR}/lib:${XERCESCROOT}/lib:${LD_LIBRARY_PATH}
  else
     LD_LIBRARY_PATH=${G4LIB}:${CLHEP_BASE_DIR}/lib:${XERCESCROOT}/lib
  fi
  export LD_LIBRARY_PATH

  if [ $G4WORKDIR ]
  then
    echo "G4WORKDIR already set"   
  else
    echo "Setting G4WORKDIR"
    export G4WORKDIR=$HOME/geant4/$g4version/$os
  fi

  if [ ! -d $G4WORKDIR ]
  then
     echo "Creating the G4WORKDIR directory $G4WORKDIR"
     mkdir -p $G4WORKDIR
  fi

  # The following is used to store your executables 
  #            (in subdirectories, one per system)
  if [ $G4BIN ]
  then
    echo "G4BIN already set"
  else
    export G4BIN=$G4WORKDIR/bin
  fi

fi
