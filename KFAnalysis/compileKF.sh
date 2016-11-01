#!/bin/bash

p=$PWD
cd kfparticle
test -d build || mkdir build
cd build
source /cvmfs/alice.gsi.de/kfparticle/SetVc.sh
cmake ../
make -j2
cp libKFParticle.so $p
cd $p
