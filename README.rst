An example of using f2py to call Fortran module from Python
=============================================================================================

f2py is a command line that is included with Numpy. So, we can use f2py without any installing.

Please note that before running the python file, we have to compile the Fortran file to the Python using the command line:
f2py -c -m fortran_lib 10fortran_lib.f90

Or we can run the bash script 1compile_fortran_sub.script in this repository.





