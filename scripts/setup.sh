#!/bin/sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh

BASEDIR=`pwd`
echo "base directory: "$BASEDIR

export LD_LIBRARY_PATH=$BASEDIR/lib:$BASEDIR/DynamicTTree/lib/:$BASEDIR/CfgManager/lib/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$BASEDIR/lib:$BASEDIR/DynamicTTree/lib/:$BASEDIR/CfgManager/lib/:$DYLD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$BASEDIR/interface:$BASEDIR/DynamicTTree/interface/:$BASEDIR/CfgManager/interface/:$ROOT_INCLUDE_PATH
