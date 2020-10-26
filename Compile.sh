gcc P_perturbation.c -o P_perturbation -L$HOME/local/lib -I$HOME/local/include -lm -lgsl -lgslcblas -fopenmp

gcc Density_Grid.c -o Density_Grid -L$HOME/local/lib -I$HOME/local/include -lm
gcc PowerSpectrum.c -o PowerSpectrum -L$HOME/local/lib -I$HOME/local/include -lm -lfftw3
gcc Propagator.c -o Propagator -L$HOME/local/lib -I$HOME/local/include -lm -lfftw3
gcc BiSpectrum.c -o BiSpectrum -L$HOME/local/lib -I$HOME/local/include -lm -lfftw3
gcc BiSpectrum_old.c -o BiSpectrum_old -L$HOME/local/lib -I$HOME/local/include -lm -lfftw3
gcc TriSpectrum.c -o TriSpectrum -L$HOME/local/lib -I$HOME/local/include -lm -lfftw3
