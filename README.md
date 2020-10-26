# Spetra_Estimators

The codes available here can be used to estimate different statistics of the matter density field in an N-body simulation using different window functions (NGP, CIC, spherical, or exponential) with different smoothing scales. The codes compute the spectra for the overdensity field delta, the overdensity field with the window function deconvolved, the non-linear field A=ln(1+delta), and the marked field m.

This folder contains the following codes:
-Density_Grid.c: This code computes the matter density grid for a given window function and smoothing scale;
-PowerSpectrum.c: This code computes the power-spectrum for the density grid computed with the code above;
-BiSpectrum.c: This code computes the bi-spectrum for the density grid computed with the code above;
-TriSpectrum.c: This code computes the tri-spectrum for the density grid computed with the code above;
-Compile.sh: This is a very simple bash file used to compute the codes above.

The complete list of input options for each code is given by the code when started without any input.

The Density_Grid.c assumes a file with particles where the first line has the total number of particles and the other lines have the 3D position and velocity for each particle. The expected format for this file is C binary.

If you use any of the codes present here, please cite the paper: -----.
