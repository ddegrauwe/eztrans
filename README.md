EZTRANS

A mimimal LAM transformation package

installation (on kili):

module load foss CMake
export FC=mpif90
mkdir -p build/eztrans
cd build/eztrans
rm -rf ../../build/eztrans/*
cmake -DCMAKE_INSTALL_PREFIX=../../install/eztrans ../../sources/eztrans
make
make install
