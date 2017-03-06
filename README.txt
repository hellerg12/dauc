This directory contains the fortran source code (dauc.f) for conducting inference on the difference in AUCs from nested binary regression models. It can handle a maximum of 1000 observations and 20 variables. A simulated data example is provided (sim_ipmn.dat) along with the results obtained from the analysis (sim_ipmn.res).

Compile the fortran code to obtain the executable. For example "gfortran -o dauc dauc.f" will generate the executable named dauc

The program requires the data to be in a space delimited text file with the binary response variable as the first column and the covariates in the remaining columns.

When the program is executed it will prompt for the following inputs through STDIN (standard input).

1. sample size
2. number of existing covariates
3. number of new covariates
4. column indicators of existing covariates one at a time
5. column indicators of new covariates one at a time
6. seed
7. File where data is kept

NOTE: The column indicators for the covariates are labeled from 1 to k (think of response as column zero). The anchor variable is typically the variable with the largest logistic regression coefficient in order to achieve numerical stability of the MRC estimates. seed is a negative integer. File name is given as character in single quote.

See sim_ipmn.res for an example of the session transcript of running the program on the example data.

