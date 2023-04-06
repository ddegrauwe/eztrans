# EZTRANS

A mimimal LAM spectral transform package

## Introduction

GPUs are powerful but stupid. The underlying idea of this repository is that by reducing the spectral transforms to their core (transpositions + FFTs), they may be easier to port to GPUs and give better performance. Some features that are missing w.r.t. the full code:
* splitting of latitudes
* splitting of longitudes
* different types of fields (wind, vorticity/divergence, scalar fields)
* integrated calculation of derivatives
* different MPI distribution in gridpoint and spectral space


## Installation (on kili):

```
module load foss CMake
export FC=mpif90
mkdir -p build/eztrans
cd build/eztrans
rm -rf ../../build/eztrans/*
cmake -DCMAKE_INSTALL_PREFIX=../../install/eztrans ../../sources/eztrans
make
make install
```

## TODO:
* introduce FFTs
* improve performance on OpenMP (and possibly OpenACC) by making sure that different threads (i.e. outer loops in filling send/recv buffers) *write* to different memory regions, rather than reading from different memory regions
* introduce profiling with drhook
* introduce unit tests for MPI_alltoallv, FFT, individual transpositions
* check performance
* if performance is promising, introduce derivatives, uvtovd. This will require to dynamically change nfld. Note that not all fields need to go to S space: if only derivatives are needed, L-space is sufficient.