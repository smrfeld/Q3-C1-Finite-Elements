# Q3 C1 Finite Elements

Simple implementation of Q3 C1 finite elements in 1D/2D/3D (line/square/cube).

## Installation

Use the makefile: `make` and `make install` to build the shared library. The default location is `/usr/local/lib` for the library and `/usr/local/include` for the headers. 

There is a convenient header `/include/q3c1` such that we can just use `#include <q3c1>` after installation.

## Linking

It is a shared library - use `g++ -std=c++14 -lq3c1 my_program.cpp -o my_program.o`.

## Namespace

The namespace is: `q3c1`.