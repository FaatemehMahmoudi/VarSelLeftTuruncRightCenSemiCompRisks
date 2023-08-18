# VarSelLeftTuruncRightCenSemiCompRisks
Variable selection using Broken Adaptive Ridge regression for potentially left-truncated right-censored data in semi-competing risks data. 

BAseline hazard functions are modeled in two ways: 
Parametric and semiparametric. 

files with "weibull" word in them are for the parametric model, and the files with "bernstein" are related to the semiparametric modeling as Bernstein polynomials have been used. 

functions related to data generation and unpenalized estimation are in "sim-....R" files, and functions for variable selection using Lasso, Adaptive Lasso and BAR are in the files starting with "varsel....R". 

