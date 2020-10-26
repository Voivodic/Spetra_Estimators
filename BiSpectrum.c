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

/*Give the density to each grid using the NGP density assignment*/
void NGP(double *grid, double *pos, int nd, double L, double mass){
	double Ld =  L/nd;
	int xt[3];
	
	ind(pos, xt, Ld, nd);
	
	grid[xt[0]*nd*nd + xt[1]*nd + xt[2]] += mass;
}

/*Give the density to each grid using the CIC density assignment*/
void CIC(double *grid, double *pos, int nd, double L, double npart){
	double Ld =  L/nd;
	double dx[3], t[3];
	int xt[3], sign[3], i, j, k;
	
	ind(pos, xt, Ld, nd);
	
	for(i=0;i<3;i++){
		dx[i] = pos[i]/Ld - (xt[i] + 1.0/2.0);	

		if(dx[i]>=0.0){ 
			sign[i] = 1;
			t[i] = 1.0 - dx[i];
		}
		else{
			sign[i] = -1;
			dx[i] = -dx[i];
			t[i] = 1.0 - dx[i];
		}
	}

	grid[xt[0]*nd*nd + xt[1]*nd + xt[2]] += npart*t[0]*t[1]*t[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + xt[1]*nd + xt[2]] += npart*dx[0]*t[1]*t[2];
	grid[xt[0]*nd*nd + mod(xt[1],sign[1],nd)*nd + xt[2]] += npart*t[0]*dx[1]*t[2];
	grid[xt[0]*nd*nd + xt[1]*nd + mod(xt[2],sign[2],nd)] += npart*t[0]*t[1]*dx[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + mod(xt[1],sign[1],nd)*nd + xt[2]] += npart*dx[0]*dx[1]*t[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + xt[1]*nd + mod(xt[2],sign[2],nd)] += npart*dx[0]*t[1]*dx[2];
	grid[xt[0]*nd*nd + mod(xt[1],sign[1],nd)*nd + mod(xt[2],sign[2],nd)] += npart*t[0]*dx[1]*dx[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + mod(xt[1],sign[1],nd)*nd + mod(xt[2],sign[2],nd)] += npart*dx[0]*dx[1]*dx[2];
}

/*Give the density to each grid using a sphere*/
double Sphere(double *grid, double *pos, int nd, double L, double R, double mass){
	double Ld =  L/nd;
	double dx[3], M;
	int xt[3], i, j, k, l, times;
	
	ind(pos, xt, Ld, nd);
	times = floor(0.5 + R/Ld);
	
	M = 0.0;
	for(i=-times;i<=times;i++)
		for(j=-times;j<=times;j++)
			for(k=-times;k<=times;k++){
				dx[0] = pos[0] - (mod(xt[0], i, nd) + 1.0/2.0)*Ld;
				dx[1] = pos[1] - (mod(xt[1], j, nd) + 1.0/2.0)*Ld;
				dx[2] = pos[2] - (mod(xt[2], k, nd) + 1.0/2.0)*Ld;
	
				for(l=0;l<3;l++){	
					if(dx[l] > Ld*(times+1))	dx[l] = L - dx[l];
					if(dx[k] < -Ld*(times+1))	dx[l] = L + dx[l];
				}

				if(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] <= R*R){
					grid[mod(xt[0],i,nd)*nd*nd + mod(xt[1],j,nd)*nd + mod(xt[2],k,nd)] += mass;
					M += mass;
				}
			}
	return M;
}

/*Give the density to each grid using a sphere*/
double Exp(double *grid, double *pos, int nd, double L, double R, double R_times, double mass){
	double Ld =  L/nd;
	double dx[3], M, r2, w;
	int xt[3], i, j, k, l, times;
	
	ind(pos, xt, Ld, nd);
	times = floor(0.5 + R*R_times/Ld);
	
	M = 0.0;
	for(i=-times;i<=times;i++)
		for(j=-times;j<=times;j++)
			for(k=-times;k<=times;k++){
				dx[0] = pos[0] - (mod(xt[0], i, nd) + 1.0/2.0)*Ld;
				dx[1] = pos[1] - (mod(xt[1], j, nd) + 1.0/2.0)*Ld;
				dx[2] = pos[2] - (mod(xt[2], k, nd) + 1.0/2.0)*Ld;
	
				for(l=0;l<3;l++){	
					if(dx[l] > Ld*(times+1))	dx[l] = L - dx[l];
					if(dx[k] < -Ld*(times+1))	dx[l] = L + dx[l];
				}

				r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
				if(r2 <= R*R*R_times*R_times){
					w = exp(-0.5*r2/R/R);
					grid[mod(xt[0],i,nd)*nd*nd + mod(xt[1],j,nd)*nd + mod(xt[2],k,nd)] += mass*w;
					M += mass*w;
				}
			}
	return M;
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
int i, j, k, a, b, c, tmp, window;
int np, trashi, cont, nd, ng, Nk, ind, contk1, contk2, contk3, tmpb, tmpc, npart, *PDF;
float tmpf, trashf;
double L, X[3], den, *grid1, *grid2, *M1, *M2, *K, n_bar1, n_bar2, kn, Normx, Normk, kmin, kmax, dk, kx, ky, kz, kmod, Lb, Mtot1, Mtot2, Rsmooth, R_times, Bm, Bw, BA, BM, I, Pm, Pw, PA, PM, IP, Kmean1, Kmean2, Kmean3;
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
			gridk1[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0;
			gridk1[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0;
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
			gridk2[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0;
			gridk2[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0;
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
			gridk3[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0;
			gridk3[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0;
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

/*Alloc the final arrays*/
K = (double *)malloc(Nk*sizeof(double));

/*Set the final vectors*/
kmin = kn;
kmax = kn*nd/2.0;
dk = (kmax - kmin)/Nk;
for(i=0;i<Nk;i++)
	K[i] = kmin + dk*(i + 0.5);

/*Set the arrays for the FFT's*/
fftw_complex *in_Bm, *in_Bw, *in_BA, *in_BM, *in_I;
double *out_Bm1, *out_Bw1, *out_BA1, *out_BM1, *out_I1, *out_Bm2, *out_Bw2, *out_BA2, *out_BM2, *out_I2, *out_Bm3, *out_Bw3, *out_BA3, *out_BM3, *out_I3;

fftw_plan p_Bm1, p_Bw1, p_BA1, p_BM1, p_I1, p_Bm2, p_Bw2, p_BA2, p_BM2, p_I2, p_Bm3, p_Bw3, p_BA3, p_BM3, p_I3;

/*Alloc the arrays for the FFT's arrays*/
in_Bm = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
in_Bw = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
in_BA = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
in_BM = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
in_I = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));

out_Bm1 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_Bw1 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_BA1 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_BM1 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_I1 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_Bm2 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_Bw2 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_BA2 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_BM2 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_I2 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_Bm3 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_Bw3 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_BA3 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_BM3 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));
out_I3 = (double*) fftw_malloc(nd*nd*nd*sizeof(double));

p_Bm1 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_Bm, out_Bm1, FFTW_ESTIMATE); 
p_Bw1 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_Bw, out_Bw1, FFTW_ESTIMATE); 
p_BA1 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_BA, out_BA1, FFTW_ESTIMATE); 
p_BM1 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_BM, out_BM1, FFTW_ESTIMATE); 
p_I1 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_I, out_I1, FFTW_ESTIMATE); 
p_Bm2 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_Bm, out_Bm2, FFTW_ESTIMATE); 
p_Bw2 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_Bw, out_Bw2, FFTW_ESTIMATE); 
p_BA2 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_BA, out_BA2, FFTW_ESTIMATE); 
p_BM2 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_BM, out_BM2, FFTW_ESTIMATE);
p_I2 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_I, out_I2, FFTW_ESTIMATE); 
p_Bm3 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_Bm, out_Bm3, FFTW_ESTIMATE); 
p_Bw3 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_Bw, out_Bw3, FFTW_ESTIMATE); 
p_BA3 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_BA, out_BA3, FFTW_ESTIMATE); 
p_BM3 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_BM, out_BM3, FFTW_ESTIMATE);
p_I3 = fftw_plan_dft_c2r_3d(nd, nd, nd, in_I, out_I3, FFTW_ESTIMATE); 

/*Open the Bispectrum output file*/
sprintf(outfile,"Bmmm_%s.dat", argv[1]);
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

/**************************/
/* Compute the Bispectrum */
/**************************/
printf("Computing the Bispectrum.\n");
	
for(a=0;a<Nk;a++){
	printf("Doing the i = %d/%d configurations.\n", a+1, Nk);

	/*Compute the density grids for the first field*/
	Kmean1 = 0.0;
	contk1 = 0;
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
				tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;
				
				if(ind == a){
					Kmean1 += kmod;
					contk1 += 1;

					in_Bm[tmp][0] = gridk1[tmp][0];
					in_Bm[tmp][1] = gridk1[tmp][1];

                    in_Bw[tmp][0] = gridk1[tmp][0]/W(kx, ky, kz, Lb, Rsmooth, window);
					in_Bw[tmp][1] = gridk1[tmp][1]/W(kx, ky, kz, Lb, Rsmooth, window);

					in_BA[tmp][0] = gridk2[tmp][0];
					in_BA[tmp][1] = gridk2[tmp][1];

					in_BM[tmp][0] = gridk3[tmp][0];
					in_BM[tmp][1] = gridk3[tmp][1];

					in_I[tmp][0] = 1.0;
					in_I[tmp][1] = 0.0;
				}
				else{
					in_Bm[tmp][0] = 0.0;
					in_Bm[tmp][1] = 0.0;

					in_Bw[tmp][0] = 0.0;
					in_Bw[tmp][1] = 0.0;                    

					in_BA[tmp][0] = 0.0;
					in_BA[tmp][1] = 0.0;

					in_BM[tmp][0] = 0.0;
					in_BM[tmp][1] = 0.0;
	
					in_I[tmp][0] = 0.0;
					in_I[tmp][1] = 0.0;
				}				
			}
		}
	}

	/*Compute the FFTs*/
	fftw_execute(p_Bm1);
	fftw_execute(p_Bw1);    
	fftw_execute(p_BA1);
	fftw_execute(p_BM1);
	fftw_execute(p_I1);	

	tmpb = floor(a/2);
	for(b=tmpb;b<a;b++){
		/*Compute the density grids for the second field*/
		Kmean2 = 0.0;
		contk2 = 0;
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
					tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;
				
					if(ind == b){
						Kmean2 += kmod;
						contk2 += 1;

				    	in_Bm[tmp][0] = gridk1[tmp][0];
				    	in_Bm[tmp][1] = gridk1[tmp][1];

                        in_Bw[tmp][0] = gridk1[tmp][0]/W(kx, ky, kz, Lb, Rsmooth, window);
					    in_Bw[tmp][1] = gridk1[tmp][1]/W(kx, ky, kz, Lb, Rsmooth, window);

				    	in_BA[tmp][0] = gridk2[tmp][0];
				    	in_BA[tmp][1] = gridk2[tmp][1];

					    in_BM[tmp][0] = gridk3[tmp][0];
				    	in_BM[tmp][1] = gridk3[tmp][1];

						in_I[tmp][0] = 1.0;
						in_I[tmp][1] = 0.0;
					}
					else{
				    	in_Bm[tmp][0] = 0.0;
				    	in_Bm[tmp][1] = 0.0;

 				    	in_Bw[tmp][0] = 0.0;
				    	in_Bw[tmp][1] = 0.0;                       

				    	in_BA[tmp][0] = 0.0;
				    	in_BA[tmp][1] = 0.0;

					    in_BM[tmp][0] = 0.0;
					    in_BM[tmp][1] = 0.0;
	
						in_I[tmp][0] = 0.0;
						in_I[tmp][1] = 0.0;
					}				
				}
			}
		}

		/*Compute the FFTs*/
	    fftw_execute(p_Bm2);
        fftw_execute(p_Bw2);
    	fftw_execute(p_BA2);
	    fftw_execute(p_BM2);
	    fftw_execute(p_I2);	

		tmpc = a - b - 1;
		if(tmpc < 0) tmpc = 0;
		for(c=tmpc;c<b;c++){
			/*Compute the density grids for the second field*/
			Kmean3 = 0.0;
			contk3 = 0;
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
						tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;
				
						if(ind == c){
							Kmean3 += kmod;
							contk3 += 1;

				        	in_Bm[tmp][0] = gridk1[tmp][0];
				    	    in_Bm[tmp][1] = gridk1[tmp][1];

                            in_Bw[tmp][0] = gridk1[tmp][0]/W(kx, ky, kz, Lb, Rsmooth, window);
					        in_Bw[tmp][1] = gridk1[tmp][1]/W(kx, ky, kz, Lb, Rsmooth, window);

				    	    in_BA[tmp][0] = gridk2[tmp][0];
				    	    in_BA[tmp][1] = gridk2[tmp][1];

					        in_BM[tmp][0] = gridk3[tmp][0];
				    	    in_BM[tmp][1] = gridk3[tmp][1];

							in_I[tmp][0] = 1.0;
							in_I[tmp][1] = 0.0;
						}
						else{
				    	    in_Bm[tmp][0] = 0.0;
				        	in_Bm[tmp][1] = 0.0;

				    	    in_Bw[tmp][0] = 0.0;
				        	in_Bw[tmp][1] = 0.0;                            

				    	    in_BA[tmp][0] = 0.0;
				        	in_BA[tmp][1] = 0.0;

					        in_BM[tmp][0] = 0.0;
					        in_BM[tmp][1] = 0.0;
	
							in_I[tmp][0] = 0.0;
							in_I[tmp][1] = 0.0;
						}				
					}
				}
			}

			/*Compute the FFTs*/
	        fftw_execute(p_Bm3);
            fftw_execute(p_Bw3);
    	    fftw_execute(p_BA3);
	        fftw_execute(p_BM3);
	        fftw_execute(p_I3);	

			/*Compute the sum of the grids in real space and save teh Bispectrum*/
			Bm = 0.0;
            Bw = 0.0;
            BA = 0.0;
            BM = 0.0;
			I = 0.0;
			for(i=0;i<nd;i++)
       				for(j=0;j<nd;j++)
					    for(k=0;k<nd;k++){
					    	tmp = i*nd*nd + j*nd + k;

					    	Bm += out_Bm1[tmp]*out_Bm2[tmp]*out_Bm3[tmp];
                            Bw += out_Bw1[tmp]*out_Bw2[tmp]*out_Bw3[tmp];
					    	BA += out_BA1[tmp]*out_BA2[tmp]*out_BA3[tmp];
					    	BM += out_BM1[tmp]*out_BM2[tmp]*out_BM3[tmp];
					    	I += out_I1[tmp]*out_I2[tmp]*out_I3[tmp];
					    }
			Bm = Bm/I*pow(L, 6.0)/pow(ng, 3.0);
            Bw = Bw/I*pow(L, 6.0)/pow(ng, 3.0);
			BA = BA/I*pow(L, 6.0)/pow(ng, 3.0);
			BM = BM/I*pow(L, 6.0)/pow(ng, 3.0);

			fprintf(out, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", a+1, b+1, c+1, Kmean1/contk1, Kmean2/contk2, Kmean3/contk3, Bw, Bm, BA, BM, I/ng/ng);	
		}/*Close the loop in the third field*/
		
		/*Compute the sum of the grids in real space and save teh Bispectrum*/
		Bm = 0.0;
        Bw = 0.0;
        BA = 0.0;
        BM = 0.0;
		I = 0.0;
		for(i=0;i<nd;i++)
       			for(j=0;j<nd;j++)
				    for(k=0;k<nd;k++){
				    	tmp = i*nd*nd + j*nd + k;

				    	Bm += out_Bm1[tmp]*out_Bm2[tmp]*out_Bm2[tmp];
                        Bw += out_Bw1[tmp]*out_Bw2[tmp]*out_Bw2[tmp];
				    	BA += out_BA1[tmp]*out_BA2[tmp]*out_BA2[tmp];
				    	BM += out_BM1[tmp]*out_BM2[tmp]*out_BM2[tmp];                        
				    	I += out_I1[tmp]*out_I2[tmp]*out_I2[tmp];
				    }
		Bm = Bm/I*pow(L, 6.0)/pow(ng, 3.0);
        Bw = Bw/I*pow(L, 6.0)/pow(ng, 3.0);
		BA = BA/I*pow(L, 6.0)/pow(ng, 3.0);
		BM = BM/I*pow(L, 6.0)/pow(ng, 3.0);

		fprintf(out, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", a+1, b+1, b+1, Kmean1/contk1, Kmean2/contk2, Kmean2/contk2, Bw, Bm, BA, BM, I/ng/ng);
	}/*Close the loop in the second field*/

	/*Compute the sum of the grids in real space and save teh Bispectrum*/
	Bm = 0.0;
    Bw = 0.0;
    BA = 0.0;
    BM = 0.0;
	Pm = 0.0;
    Pw = 0.0;
    PA = 0.0;
    PM = 0.0;
	I = 0.0;
	IP = 0.0;
	for(i=0;i<nd;i++)
       		for(j=0;j<nd;j++)
			    for(k=0;k<nd;k++){
			    	tmp = i*nd*nd + j*nd + k;

			    	Bm += out_Bm1[tmp]*out_Bm1[tmp]*out_Bm1[tmp];
                    Bw += out_Bw1[tmp]*out_Bw1[tmp]*out_Bw1[tmp];
                    BA += out_BA1[tmp]*out_BA1[tmp]*out_BA1[tmp];
                    BM += out_BM1[tmp]*out_BM1[tmp]*out_BM1[tmp];
			    	I += out_I1[tmp]*out_I1[tmp]*out_I1[tmp];

			    	Pm += out_Bm1[tmp]*out_Bm1[tmp];
                    Pw += out_Bw1[tmp]*out_Bw1[tmp];
                    PA += out_BA1[tmp]*out_BA1[tmp];
                    PM += out_BM1[tmp]*out_BM1[tmp];
			    	IP += out_I1[tmp]*out_I1[tmp];					
		    	}
	Bm = Bm/I*pow(L, 6.0)/pow(ng, 3.0);
    Bw = Bw/I*pow(L, 6.0)/pow(ng, 3.0);
	BA = BA/I*pow(L, 6.0)/pow(ng, 3.0);
	BM = BM/I*pow(L, 6.0)/pow(ng, 3.0);

	Pm = Pm/IP*pow(L, 3.0)/pow(ng, 2.0);
    Pw = Pw/IP*pow(L, 3.0)/pow(ng, 2.0);
	PA = PA/IP*pow(L, 3.0)/pow(ng, 2.0);
	PM = PM/IP*pow(L, 3.0)/pow(ng, 2.0);

	fprintf(out, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", a+1, a+1, a+1, Kmean1/contk1, Kmean1/contk1, Kmean1/contk1, Bw, Bm, BA, BM, I/ng/ng);
	fprintf(out, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", a+1, -1, -1, Kmean1/contk1, Kmean1/contk1, Kmean1/contk1, Pw, Pm, PA, PM, I/ng);
}/*Close the loop in the first field*/

/*Free the grids*/
fftw_free(gridk1);
fftw_free(gridk2);
fftw_free(gridk3);

/*Free the FFTW3 arrays*/
fftw_free(in_Bm); fftw_free(in_Bw); fftw_free(in_BA); fftw_free(in_BM); fftw_free(in_I);
fftw_free(out_Bm1); fftw_free(out_Bm2); fftw_free(out_Bm3); fftw_free(out_Bw1); fftw_free(out_Bw2); fftw_free(out_Bw3); fftw_free(out_BA1); fftw_free(out_BA2); fftw_free(out_BA3); fftw_free(out_BM1); fftw_free(out_BM2); fftw_free(out_BM3);fftw_free(out_I1); fftw_free(out_I2); fftw_free(out_I3);
fftw_free(p_Bm1); fftw_free(p_Bm2); fftw_free(p_Bm3); fftw_free(p_Bw1); fftw_free(p_Bw2); fftw_free(p_Bw3); fftw_free(p_BA1); fftw_free(p_BA2); fftw_free(p_BA3); fftw_free(p_BM1); fftw_free(p_BM2); fftw_free(p_BM3); fftw_free(p_I1); fftw_free(p_I2); fftw_free(p_I3);

/*Free the memory*/
free(K);

/*Close the output files*/
fclose(out);

printf("Everything done\n");

return 0;
}