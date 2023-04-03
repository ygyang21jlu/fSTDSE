Last updated: 2015 July 11
------------

IMPORTANT: Quadruple-precision support is not enabled by default. You do not
IMPORTANT: need to follow instructions in this file unless you are building
IMPORTANT: the quadruple-precision version of the code.

Since LAPACK does not support Fortran quadruple precision by default, building 
with quadruple-precision support requires a special build of the LAPACK and BLAS 
libraries. The special build of the LAPACK exports quadruple-precision entry
points using the usual double-precision function and subroutine names, prefixed
by "quad_". The single-precision entry points are of undefined semantics and 
should not be used; to avoid clashes with the standard LAPACK build, these are 
also prefixed by "quad_".

The quad-precision LAPACK builds depend on the compiler's type-promotion features; 
since these features are not standardized, some experimentation may be necessary 
for the compiler you are using.

Once the build is complete, it is essential to verify that *all* LAPACK
and BLAS entry points in the quad-precision LAPACK build have the
quad_ prefix. If the command:

nm libquadlapack_*.a | egrep ' [T|t] ' | grep -v 'quad_'

produces a non-empty output for the combined library, the extra entry
points must be renamed using objcopy. Failure to do so will lead to
unpredictable execution results.

Examples, tested with LAPACK 3.5.0 (http://www.netlib.org/lapack/):

  gfortran 4.8.2
  --------------

Build and test LAPACK following instruction in LAPACK distribution. Use build
settings (make.inc) supplied in ./extras/make_gfortran_quad.inc as the starting
point. The "-freal-8-real-16" flag instructs gfortran to promote double to
quadruple precision. If you experience problems building quad-precision LAPACK
library using settings in make_gfortran_quad.inc, you may also try the build
settings in ./extras/make_gfortran_vanilla_quad.inc, which avoids advanced
optimization settings.

Once the build is complete, combine LAPACK and reference BLAS libraries,
and rename the entry points to include the quad_ prefix.

Assuming that liblapack.a, librefblas.a and rename-lapack-syms_gfortran.lst
(found in the extras/ subdirectory of the distribution) are in the current 
directory:

mkdir xxx
cd xxx
ar x ../librefblas.a
ar x ../liblapack.a
ar r libquadlapack_gfortran.a *.o
ar s libquadlapack_gfortran.a
objcopy --redefine-syms ../rename-lapack-syms_gfortran.lst libquadlapack_gfortran.a ../libquadlapack_gfortran.a
cd ..
rm -rf xxx/

  Intel Fortran 15.0.1
  --------------------

Follow instructions for gfortran above, but use build options in make_ifort_quad.inc 
and rename list in rename-lapack-syms_intel.lst, both found in the ./extras/
subdirectory. The Intel Fortran flag "-r16" instructs the compiler to promote
default REAL kind (single precision) to quadruple precision.
