                                                                                                                          CPC Librarian's Note
____________________________________________________________________________________________________________________________________________

I had some problems installing the quad version of Lapack on my ubuntu-14.04 machine using gfortran, 4.8.2

I wrote to the author and he sent the following information to help me find the problem. I found it useful and thought it might help others as well.


The SCID-TDSE "extras/make_gfortran_quad.inc" file assumes that a 64-bit AVX static-build environment is available.

Below are some checks you can make to ascertain if your environment conforms to this.

0. Make a list of all compiler configuration files you'll be using. This includes:

    a. LAPACK make.inc (only if you decide to continue with building the quad-precision library)
    b. scid-tdse build configuration file. This will be one of the files in the configs/ sub-directory; to
        find (and possibly change) which one, open the main Makefile, and look for the first uncommented
        line after "# System-specific overrides"

1. Does your kernel support 64-bit binaries?
     To check, run the command:

     uname -m

     If it reports "x86_64" or "amd64", your CPU and kernel support 64-bit binaries; proceed to #2 below.

     If it reports anything else, your system is not 64-bit capable; in all build configuration files, delete all
     instances of "-m64" command option.

2. Does your CPU and kernel have Intel AVX support?
    To check, run the command:

     grep avx /proc/cpuinfo

     If this command returns non-empty output, you have AVX support; proceed to #3 below.

     If the output is empty or you get an error message, your system is 64-bit, but does not support AVX; in all
     build configuration files, delete all instances of the "-mavx" command option.

3. Does your gfortran installation support building programs?
    To check, create a file "hello.f90", containing:
     --- cut here ---
 program hello
   print *, 'Hello, world!'
 end program hello
     --- cut here ---

    Compile if using the command (delete -m64 and/or -mavx as necessary):

    gfortran -m64 -mavx hello.f90

    If this command completes without errors, and you are able to run the resulting executable, proceed to #4 below.

    Otherwise, you are missing Fortran support libraries. Install the appropriate libraries using your favorite package
    manager.

4. Does your gfortran installation support building programs statically?

    To check, compile the hello.f90 test using the command (delete -m64 and/or -mavx as necessary):

    gfortran -m64 -mavx -static hello.f90

    If this command completes without errors, and you are able to run the resulting executable, proceed to #5 below.

    Otherwise, either install the necessary static-build support files, or edit all build control files you use, deleting all
    instances of the option "-static".
