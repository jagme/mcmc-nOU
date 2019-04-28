
MCMCnOU
==========

This is the Matlab/C++-Mex implementation of the MCMC algorithms presented in

“Bayesian calibration and number of jump components in electricity spot price models” by
Jhonny Gonzalez, John Moriarty & Jan Palczewski (https://www.sciencedirect.com/science/article/pii/S0140988317301299)

The software performs Bayesian inference for sum of Gaussian and non-Gaussian Ornstein-Uhlenbeck processes, and demonstrates their use for modelling the UK and EEX electricity markets.

Please start with file ‘mainMCMC0.m’, which contains the full series of commands to reproduce the numerical results and plots in the article.

The commands in ‘mainTest.m’ test the algorithms on simulated data. 

Requirements
============

You will need Statistics and Machine Learning Toolbox to run the algorithms, and optionally Econometrics Toolbox to plot sample ACFs and the Parallel Computing Toolbox to speed up the code in repeatAnalysis.m. 

The price data under study is in the folder ‘data’.

To improve computational performance, the code in file getY2.m has been compiled 
into a C++ Mex file and is available for Unix and Linux platforms with 
Matlab 2014b or later. Windows platforms will use the pure Matlab implementation automatically. If calling getY2 function under Unix/Linux produces an error, 
probably your platform/version is not compatible.  A workaround is to rename 
the files getY2.mexmaci64 (Unix) and/or getY2.mexa64 (Linux), e.g.
getY2.mexmaci64 ->>> getY2_.mexmaci64 so the mex files are not called.

The latest version of this software can be found at https://github.com/jagme/mcmc-nOU

Email me with comments at: jhonny.gonzalez at manchester dot ac dot uk 
                           

License
============
See ‘license.txt` for licensing information.
