# MSci Project
Theoretical Physics MSci project repo for the LHAPDF toolkit for CERN. Contains: 1D/2D cubic spline interpolation, Chebyshev polynomials and future reaction data

## Compilation and Execution Guide
There is no CMAKE for this, or any other MakeFile scripts, as there are not too many dependancies. Compiling the code is simple, for example, to compile the code to test the differing effect of boundary conditions on the PDF-like function just run
```
 g++ main_pdflike.cpp bicubicspline.cpp cubicspline.cpp boundaryderivatives.cpp extrapolate.cpp funcs.cpp -o ../executables/pdflike_test.out
```
Then to execute just enter
```
../executables/pdflike_test.out
```
This will save the output data in a csv within the outputs subdirectory, with whatever name was given in the code. The automatic name is "pdfLike_bicubic.csv". 

Of course, to run any other piece of code such as the derivatives test, all you need to do is replace the first and last command line arguments with the program of choice and then whatever you want the executable to be called.

To run the plots just run the python scripts within the Python directory. They should be set up to automatically work, but if you update any of the grid parameters or filenames, this will have to then be updated in the scripts.