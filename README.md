# CT_meanfield-analysis-and-spatial-mutual-information
# Set up instructions
This is a repository for the paper "Critical transitions in spatial systems induced by Ornstein-Uhlenbeck noise: Spatial mutual information as a precursor" by  Smita Deb and Partha Sharathi Dutta.
# User instructions
Compilation of the given codes require Matlab version R2020b and Fortran.

# Codes and results
spatial_bif_analytic.m can solve the self consistent equation and also yield in the bifurcation diagram by solving the system for continuous evolution of the bifurcation parameter with other set of parameter values fixed. 

spdf.m calculates the stationary probability density function corresponding to the Fokker Planck equation of the spatially extended system.

numerical_data.f90 can simulate the sequence of spatial patterns along the bifurcation gradient for all other parameters kept fixed. This will generate pat_tau0.1.dat which can be further used for calculating spatial mutual information using mutualinfo.m and other spatial early warning signals using R package (spatial early warning master, Kefi et al, 2014). An analog of the fortran program is possible by repeating the algorithm in matlab. We use fortran for its high computation efficiency in generating large sequence of snapshots.

mutualinfo.m computes the normalized spatial mutual information using sequence of spatial pattern generated using numerical_data.f90

