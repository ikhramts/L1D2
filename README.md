
# L1D2 - Largest Lyapunov exponent and correlation dimension

This is a GitHubized clone of the 
[original program](https://www.physionet.org/content/lyapunov/1.0.0/) 
posted on physionet.org. 

A recompiled Windows x64 executable is available in the Relases
section in the sidebar.

## Introduction

This document describes how to use the program L1D2 for estimating the
largest Lyapunov exponent and correlation dimension of a reconstructed
attractor.  The input to the program is a scalar time series along
with several parameters that control the operation of the algorithm.
More information can be found in the following papers as well as many
other publications on nonlinear time series analysis:

> M.T. Rosenstein, J.J. Collins, and C.J. De Luca. A practical
> method for calculating largest Lyapunov exponents from small
> data sets. Physica D 65:117-134, 1993.
> 
> M.T. Rosenstein, J.J. Collins, and C.J. De Luca. Reconstruction
> expansion as a geometry-based framework for choosing proper
> delay times. Physica D 73:82-98, 1994.

## WARRANTY and DISCLAIMER

This software package is provided AS IS without guarantees or
warranty, expressed or implied.  THE AUTHOR SHALL NOT BE LIABLE IN ANY
WAY, EVEN FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES, IN CONNECTION WITH,
OR ARISING OUT OF THE FURNISHING, THE PERFORMANCE, OR THE USE OF THIS
SOFTWARE PACKAGE.  The author is under no obligation to provide
technical support, corrections, or upgrades to this software package.

## About the source code

The core algorithm implemented in L1D2 has been around for a while and
has been tested, although the source code itself is in a (messy)
preliminary form.  Modifying the code is NOT recommended.  Use it only
to compile the program AS IS or to get some ideas for developing your
own implementation.

## Compiling the program

On Posix systems, the easiest way to compile the program is to use gcc
to create the executable file L1D2:

```
gcc -lm -o L1D2 l1d2.c
```

There is also a CMake project definition for use with IDEs.

The program is written in generic C code, so it should compile on other
systems as well.  For instance, Metrowerks CodeWarrior provides the
SIOUX console library that makes it easy to implement text-based
programs for graphics-based operating systems like MS Windows and MacOS.


## Sample program run

The following is a sample run of the program when executed from the
Unix command-line prompt (details below).

```
L1D2> nice +15 L1D2
*** L1D2: generates scaling information for L1 and D2 ***
          v1.0 Copyright (c) 1999 M.T. Rosenstein
                                  (mtr@cs.umass.edu)
          reference: M.T. Rosenstein, J.J. Collins, C.J. De Luca,
                     A practical method for calculating largest
                     Lyapunov exponents from small data sets,
                     Physica D 65:117-134, 1993.
Enter test file name: lorenz.l1d2
Reading Test File...
detected 2 tests
Enter output file root (no extension): lorenz
Processing test 1 of 2: 
0....1....2....3....4....5....6....7....8....9....10
Processing test 2 of 2: 
0....1....2....3....4....5....6....7....8....9....10
Saving data for largest Lyapunov exponent...
Saving data for correlation dimension...
Success!
```
## Program notes

The following list is meant to act as a preliminary user's guide to
the program:

1. Type L1D2 at the command prompt to run the program.  Some runs can
be computationally intensive, so if you share your machine with others
you may want to 'nice' the program as shown above.

2. After displaying some information, the program prompts you to
enter the name of a test file.  (This example uses the test file
`lorenz.l1d2` which is included in the zip archive along with a sample
data file called `lorenz.dat`.)  The test file contains a list of all
the information needed to run a set of tests.  More on this below.

3. After reading the test file and reporting the number of tests
detected, the program prompts you to enter the file root for the
output.  The example above uses `lorenz` and so the program will save
its results in two files: `lorenz.d2` for the correlation dimension,
and `lorenz.l1` for the largest Lyapunov exponent.

4. Next the program shows the progress on each test.  Each dot
represents roughly 2% of the computation and the numbers mark the 10%
complete point, the 20% complete point, and so on.

5. The output files are ASCII text that can be read into a spread
sheet program (or a custom program) for further analysis.  The first
column shows the independent variable for the scaling plots.  For D2,
the independent variable is `ln(R)`, where `R` refers to the length scale
or radius for the nearest neighbor test.  For L1, the independent
variable is divergence time; the units are "samples" so you need to
multiply by your data sampling period to get "real" time.

6. The next group of columns in the data file contain the dependent
variable for each test, in the order given by the test file.  For D2,
the dependent variable is the (normalized) natural log of the
correlation sum, and for L1 the dependent variable is the natural log
of the divergence.

7. The next block of columns in the output files is essentially a
repeat of the first group of columns.  This time, however, the
dependent variable is the slope of the line computed by a 7-point
least-squares fit to the "raw" data from the first group of columns.
(The 7 points were centered about the corresponding independent
variable.)  Plots showing the slope are sometimes helpful with the
analysis, although the choice of a 7-point window may not be best for
some applications.  The slope data are included only as a convenience.

8. The input data file should be ASCII text with one time series per
column and with white space separating columns (tabs are okay).  Each
line of the input file contains one data sample vector.  You can
include header information only at the top of the file by making the
first character of a line by an asterisk `*` as in the test file.

9. Each time you run the program it creates a test file called
`sample.l1d2`.  You can edit `sample.l1d2` to make a new test file,
but be sure to change the file name because the program will overwrite
the old sample file.  Here's the contents of the sample test file:

```
* Header info starts with an asterisk.
*
* Each line of the test file contains the name of a data file
* followed by the parameters for the test:
*   file_name series# startIndex stopIndex m J W divergeT
*
*      file_name  = name of the data file
*      series#    = time series number to use for delay
*                      reconstruction
*      startIndex = index of first data point to read (usually 1)
*      stopIndex  = index of last data point to read
*                      (enter 0 for maximum)
*      m          = embedding dimension
*      J          = delay in samples
*      W          = window size for skipping temporally close
*                      nearest neighbors
*      divergT    = total divergence time in samples
*   example: lorenz.dat 1 1 0 5 7 100 300
```

10. The sample test file contains header information that describes
the file format.  You MUST be familiar with all of the parameters,
like the embedding dimension and the window size from Theiler's work
(not to be confused with the embedding window size).  Here's the
contents of the test file (`lorenz.l1d2`) that's included in the zip
archive:

```
*   file_name series# startIndex stopIndex m J W divergeT
lorenz.dat 1 1 0 5 11 300 500
lorenz.dat 2 1 0 7 7 300 500
```

11. The above test file contains two tests.  The first is for the
x-coordinate time series of the Lorenz system (time series #1 in the
first column of the data file) and the second is for the y-coordinate
series (series #2 in the second column of the data file).  The zeros
indicate that the program should read all the available data, although
you can change this to limit the number of data points.  You can also
change the start index if you need the program to extract data from
the middle of the data file.

12. The divergence time is usually a guess at how long it will take
the nearest neighbors to reach maximum divergence (on average).
Guessing a value that is too big will slow down the computation, but
not by a whole lot.  Guessing a value that is too small may keep you
from seeing the scaling region in the L1 results.

## For more information
Please direct questions, comments, and bug reports to Mike Rosenstein
via email to mtr@cs.umass.edu

(**NOTE:** the above line was copied from the [original code's 00README.txt](https://www.physionet.org/content/lyapunov/1.0.0/).
Mike Rosenstein is not the owner of this repository, and might not be aware that it exists.)