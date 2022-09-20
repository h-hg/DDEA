+----------------------------------------------------------------------------+
|        Benchmark Functions for the CEC'2013 Special Session and            |
|             Competition on Large-Scale Global Optimization                 |
|    Xiaodong Li, Ke Tang, Mohammad N. Omidvar, Zhenyu Yang, and Kai Qin     |
|                                                                            |
| Version    : 1.1                                                           |
| Developers : Mohammad N. Omidvar                                           |
+----------------------------------------------------------------------------+


This is the Matlab/Octave version of test suite for the Special Session
on Large Scale Global Optimization for The IEEE Congress on Evolutionary
Computation 2013.


+----------------------------------------------------------------------------+
| Code Structure                                                             |
+----------------------------------------------------------------------------+
README.txt:           This file.
demoEA.m:             A demo program (random walk) to show how to use the
                      test suite.
benchmark_func.m:     Main program to compute the fitness values for all
functions
computeRotation.m:    A helper function for generating rotation matrices.
computeShiftVector.m: A helper function for generating shift vectors.
datafiles:            Shift and rotation data used by benchmark_func.m
LICENSE.txt:          The license file.


+----------------------------------------------------------------------------+
| Testing                                                                    |
+----------------------------------------------------------------------------+
The test suite has been tested on Windows/Linux OS with Matlab 7.x.


+----------------------------------------------------------------------------+
| Contact                                                                    |
+----------------------------------------------------------------------------+
Questions and bug reports should be sent to:
Mohammad Nabi Omidvar: mohammad.omidvar@rmit.edu.au
or
Xiaodong Li: xiaodong.li@rmit.edu.au



+----------------------------------------------------------------------------+
| Change Log                                                                 |
+----------------------------------------------------------------------------+
2013-03-23: The global optimum for f5 and f12 was outside the feasible region.
            The shift vectors for f5 and f12 have been updated to be within
            the feasible region.

2013-12-24: The inconsistency between the transformations on Ackley's funciton
            between C++, MATLAB code is resolved.
