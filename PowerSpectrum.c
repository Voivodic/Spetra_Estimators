#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#define PI 3.141592653589793

/*Evaluate the ciclic sum of x and y*/
int mod(int x, int y, int nd){
	int resp;

	resp = x + y;

	if(resp<0) resp += nd;
	else if(resp>=nd) resp -= nd;

	return resp;
}

/*Define the cyclic sum for floats*/
double cysumf(double x, double y, double L){
	double resp;

	resp = x + y;
	if(resp>=L)	resp -= L;
	if(resp<0)	resp += L;

	return resp;
}

/*Give a indice for the partricle*/
void ind(double x[], int xt[], double Ld, int nd){
	xt[0] = floor(x[0]/Ld);
	xt[1] = floor(x[1]/Ld);
	xt[2] = floor(x[2]/Ld);
	if(xt[0]==nd)	xt[0] -=1;
	if(xt[1]==nd)	xt[1] -=1;
	if(xt[2]==nd)	xt[2] -=1;
}

/*Define de sinc function*/
double sinc(double x){
	double resp;
	if(x!=0.0) resp = sin(x)/x;
	else resp = 1.0;	

	return resp;
}

/*Define window function for NGP and CIC*/
double W(double k1, double k2, double k3, double Lb, double R, int window){
	double resp, kmod;

	if(window == 0)
		resp = sinc(k1*Lb/2.0)*sinc(k2*Lb/2.0)*sinc(k3*Lb/2.0);
	if(window == 1)
		resp = pow(sinc(k1*Lb/2.0)*sinc(k2*Lb/2.0)*sinc(k3*Lb/2.0), 2.0);
	if(window == 2){
		kmod = sqrt(k1*k1 + k2*k2 + k3*k3);	
		resp = 3.0/kmod/kmod/R/R*(sinc(kmod*R) - cos(kmod*R));
	}
	if(window == 3){
		kmod = k1*k1 + k2*k2 + k3*k3;	
		resp = exp(-kmod*R*R/2.0);
	}

	return resp;
}

/*Define the bin for the mode*/
int Indice(double k, double kmin, double dk){
	int resp;
	double tmp;

	tmp = (k - kmin)/dk;

	resp = floor(tmp);

	return resp;
}

/*Define the marked transformation*/
double Marked(double deltaR, double p, double deltas){
	double resp;

	resp = pow((1.0 + deltas)/(1.0 + deltas + deltaR), p);

	return resp;
}

int main(int argc,char *argv[])
{
FILE *in, *out;
char infile[100], outfile[100];
int i, j, k, tmp, window;
int np, trashi, cont, nd, ng, Nk, *contk, ind, npart, *PDF;
float tmpf, trashf;
double L, X[3], den, *grid1, *grid2, *M1, *M2, *Pmm, *Pww, *PAA, *PMM, *K, *Kmean, n_bar1, n_bar2, kn, Normx, Normk, kmin, kmax, dk, kx, ky, kz, kmod, Lb, Mtot1, Mtot2, Rsmooth, R_times;
double Mdeltas, Mp;
fftw_complex *out1, *out2, *gridk1, *gridk2, *gridk3;
fftw_plan p1, p2;

if (argc != 7){
	printf("Wrong number of arguments.\n");
	printf("arg1: Preffix for the input and output files.\n");
	printf("arg2: Box size.\n");
	printf("arg3: Number of division by dimension in the grid.\n");
	printf("arg4: Window function to use to construct the density grid: NGP (0), CIC (1), Spherical (2) or Gaussian (3)\n");
	printf("arg5: Size of the smooth. (Only used for the spherical and Gaussian window functions)\n");
	printf("arg6: R_times. (Only used for the Gaussian window function)\n");
	exit(0);
}

/*Some parameters of simulation and grid*/
nd = atoi(argv[3]);		/*Number of division by dimensions in the grid*/
L = atof(argv[2]);		/*Needed in the auto power spectrum calculaitons*/
window = atoi(argv[4]);		/*Window function to be used to construct the density field*/
Rsmooth = atof(argv[5]);	/*Smoothing radius used in the spherical smooth*/
R_times = atof(argv[6]);	/*Number of smthooing radius to look in the gaussian window*/
ng = nd*nd*nd;			/*Total number of grid cells*/
Lb = L/nd;			/*Size of each cell*/
Normx = pow(1.0/L, 3.0/2.0);	/*Normalization for the Fourier transforms (from k to x)*/
Normk = pow(L/(nd*nd), 3.0/2.0);/*Normalization for the Fourier transforms (from x to k)*/
kn = 2.0*PI/L;			/*Fundamental mode for this box size*/
Nk = (int) nd/4;		/*Number of k bins for the final power spectrum*/

/*Define the parameters of the marked field*/
Mdeltas = 0.25;
Mp = 2.0;

/*Alloc the grids*/
grid1 = (double *)malloc(ng*sizeof(double));
grid2 = (double *)malloc(ng*sizeof(double)); 
M1 = (double *)malloc(ng*sizeof(double));
M2 = (double *)malloc(ng*sizeof(double)); 

/*Open the density grid file*/
sprintf(infile,"Density_Grid_%s.dat", argv[1]);
in = fopen(infile,"rb");
if (in == NULL) {
	printf("Unable to open %s\n", infile);
	exit(0);
}

/*Save the information of the particle file*/
printf("Reading the density grid\n");

for(i=0;i<ng;i++)
    fread(&grid1[i], sizeof(double), 1, in);
for(i=0;i<ng;i++)
    fread(&grid2[i], sizeof(double), 1, in);  
fclose(in);

/*Alloc the FFTW stuff*/
printf("Computing the density grid in Fourier space\n");
out1 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
out2 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));

p1 = fftw_plan_dft_r2c_3d(nd, nd, nd, grid1, out1, FFTW_ESTIMATE); 
p2 = fftw_plan_dft_r2c_3d(nd, nd, nd, grid2, out2, FFTW_ESTIMATE);

/*Compute the FFT of the grids*/
fftw_execute(p1);
fftw_execute(p2);

/*Take the mean over the two grids*/
gridk1 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd/2+1;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			tmp = nd*(nd/2 + 1)*i + (nd/2 + 1)*j + k;
			gridk1[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0*Normk;
			gridk1[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0*Normk;
		}
	}
}

/*Compute the mean and std and the density grid*/
kx = 0.0;
ky = 0.0;
for(i=0;i<ng;i++){
	kx += grid1[i];
	ky += grid1[i]*grid1[i];
}
kx = kx/ng;
ky = sqrt(ky/ng - kx*kx);
printf("(delta) Mean = %lf and std = %lf\n", kx, ky);

/*Construct the lognormal field*/
cont = 0;
for(i=0;i<ng;i++){
	if(grid1[i] <= -1.0 || grid2[i] <= -1.0){
		cont ++;
		grid1[i] = log(1e-10);
		grid2[i] = log(1e-10);
	}
	else{
		M1[i] = Marked(grid1[i], Mp, Mdeltas);
		grid1[i] = log(grid1[i] + 1.0);
		M2[i] = Marked(grid2[i], Mp, Mdeltas);
		grid2[i] = log(grid2[i] + 1.0);
	}
}
printf("There are %d cells with delta = -1\n", cont);

/*Compute the mean and std of the lognormal grid*/
kx = 0.0;
ky = 0.0;
kz = 0.0;
for(i=0;i<ng;i++){
	kx += grid1[i];
	kz += grid2[i];
	ky += grid1[i]*grid1[i];
}
kx = kx/ng;
kz = kz/ng;
ky = sqrt(ky/ng - kx*kx);
printf("(A) Mean = %lf (%lf) and std = %lf\n", kx, kz, ky);

/*Compute the overfield A
for(i=0;i<ng;i++){
	grid1[i] = grid1[i] - kx;
	grid2[i] = grid2[i] - kz;
}*/

/*Compute the FFT of the grids*/
fftw_execute(p1);
fftw_execute(p2);

/*Take the mean over the two grids*/
gridk2 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd/2+1;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			tmp = nd*(nd/2 + 1)*i + (nd/2 + 1)*j + k;
			gridk2[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0*Normk;
			gridk2[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0*Normk;
		}
	}
}

/*Compute the mean and std of the marked grid*/
kx = 0.0;
ky = 0.0;
kz = 0.0;
for(i=0;i<ng;i++){
	kx += M1[i];
	kz += M2[i];
	ky += M1[i]*M1[i];
}
kx = kx/ng;
kz = kz/ng;
ky = sqrt(ky/ng - kx*kx);
printf("(m) Mean = %lf (%lf) and std = %lf\n", kx, kz, ky);

/*Compute the overfield m*/
for(i=0;i<ng;i++){
	grid1[i] = M1[i];
	grid2[i] = M2[i];
}

/*Compute the FFT of the grids*/
fftw_execute(p1);
fftw_execute(p2);

/*Take the mean over the two grids*/
gridk3 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd/2+1;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			tmp = nd*(nd/2 + 1)*i + (nd/2 + 1)*j + k;
			gridk3[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0*Normk;
			gridk3[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0*Normk;
		}
	}
}

/*Free the density grids*/
free(grid1);
free(grid2);
free(M1);
free(M2);
fftw_free(out1);
fftw_free(out2);
fftw_free(p1);
fftw_free(p2);

/*Set the main values for the power spectra*/
printf("Computing the power spectrum\n");
K = (double *)malloc(Nk*sizeof(double));
Kmean = (double *)malloc(Nk*sizeof(double));
contk = (int *)malloc(Nk*sizeof(int));
Pmm = (double *)malloc(Nk*sizeof(double));
Pww = (double *)malloc(Nk*sizeof(double));
PAA = (double *)malloc(Nk*sizeof(double));
PMM = (double *)malloc(Nk*sizeof(double));

/*Set the final vectors*/
kmin = kn;
kmax = kn*nd/2.0;
dk = (kmax - kmin)/Nk;
for(i=0;i<Nk;i++){
	K[i] = kmin + (i + 0.5)*dk;
	Pmm[i] = 0.0;
	Pww[i] = 0.0;
	PAA[i] = 0.0;
	PMM[i] = 0.0;
	contk[i] = 0;
	Kmean[i] = 0.0;
}

/*Compute the power spectra*/
for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd/2+1;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			kmod = sqrt(kx*kx + ky*ky + kz*kz);
			ind = Indice(kmod, kmin, dk);

			if(ind < Nk && ind >= 0){
				tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;

				Pmm[ind] += (gridk1[tmp][0]*gridk1[tmp][0] + gridk1[tmp][1]*gridk1[tmp][1]);
				Pww[ind] += (gridk1[tmp][0]*gridk1[tmp][0] + gridk1[tmp][1]*gridk1[tmp][1])/pow(W(kx, ky, kz, Lb, Rsmooth, window), 2.0);
				PAA[ind] += (gridk2[tmp][0]*gridk2[tmp][0] + gridk2[tmp][1]*gridk2[tmp][1]);
				PMM[ind] += (gridk3[tmp][0]*gridk3[tmp][0] + gridk3[tmp][1]*gridk3[tmp][1]);
				Kmean[ind] += kmod;
				contk[ind] += 1;
			}
		}
	}
}

/*Take the mean*/
for(i=0;i<Nk;i++)
	if(contk[i]>0){
		Pmm[i] = Pmm[i]/contk[i];
		PAA[i] = PAA[i]/contk[i];
		PMM[i] = PMM[i]/contk[i];
		Pww[i] = Pww[i]/contk[i];
		Kmean[i] = Kmean[i]/contk[i];	
	}

/*Open the output file*/
sprintf(outfile,"Pmm_%s.dat", argv[1]);
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

/*Save the Pmm in the out file*/
for(i=0;i<Nk;i++)
	fprintf(out, "%lf %lf %lf %lf %lf %lf %d\n", K[i], Kmean[i], Pww[i], Pmm[i], PAA[i], PMM[i], contk[i]);
fclose(out);

/*Free the memory*/
free(K);
free(Kmean);
free(contk);
free(Pmm);
free(PAA);
free(PMM);
fftw_free(gridk1);
fftw_free(gridk2);
fftw_free(gridk3);

printf("Everything done\n");

return 0;
}
