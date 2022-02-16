# CMakeList.txt : CMake project for L1D2, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("L1D2")

# Add source to this project's executable.
add_executable (l1d2 "l1d2.c")
configure_file (lorenz.l1d2 lorenz.l1d2 COPYONLY)
configure_file (lorenz.dat lorenz.dat COPYONLY)

