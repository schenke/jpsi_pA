#/usr/bin/env bash

flag=$1

mkdir -p build
cd build
rm -fr *
if [ "$flag" == "KNL" ]; then
    CC=mpiicc CXX=mpiicpc cmake .. -DKNL=ON
else
    cmake ..
fi

make -j4
make install
cd ..
rm -fr build
