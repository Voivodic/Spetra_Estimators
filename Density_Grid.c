#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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
					if(dx[l] > L/2.0)	dx[l] = L - dx[l];
					if(dx[l] < -L/2.0)	dx[l] = L + dx[l];
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
					if(dx[l] > L/2.0)	dx[l] = L - dx[l];
					if(dx[l] < -L/2.0)	dx[l] = L + dx[l];
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

int main(int argc,char *argv[])
{
FILE *first, *out;
char firstfile[100], outfile[100];
int i, j, k, tmp, window;
int np, trashi, cont, nd, ng, Nk, *contk, ind, npart, *PDF;
float tmpf, trashf;
double L, X[3], den, *grid1, *grid2,  n_bar1, n_bar2, kn, Normx, Normk, kmin, kmax, dk, kx, ky, kz, kmod, Lb, Mtot1, Mtot2, Rsmooth, R_times;

if (argc != 8){
	printf("Wrong number of arguments.\n");
	printf("arg1: Name of the first file.\n");
	printf("arg2: Preffix for the output files.\n");
	printf("arg3: Box size.\n");
	printf("arg4: Number of division by dimension in the grid.\n");
	printf("arg5: Window function to use to construct the density grid: NGP (0), CIC (1), Spherical (2) or Gaussian (3)\n");
	printf("arg6: Size of the smooth. (Only used for the spherical and Gaussian window functions)\n");
	printf("arg7: R_times. (Only used for the Gaussian window function)\n");
	exit(0);
}

/*Get the name of all files*/
sprintf(firstfile,"%s", argv[1]);

/*Some parameters of simulation and grid*/
nd = atoi(argv[4]);		/*Number of division by dimensions in the grid*/
L = atof(argv[3]);		/*Needed in the auto power spectrum calculaitons*/
window = atoi(argv[5]);		/*Window function to be used to construct the density field*/
Rsmooth = atof(argv[6]);	/*Smoothing radius used in the spherical smooth*/
R_times = atof(argv[7]);	/*Number of smthooing radius to look in the gaussian window*/
ng = nd*nd*nd;			/*Total number of grid cells*/
Lb = L/nd;			/*Size of each cell*/
Normx = pow(1.0/L, 3.0/2.0);	/*Normalization for the Fourier transforms (from k to x)*/
Normk = pow(L/(nd*nd), 3.0/2.0);/*Normalization for the Fourier transforms (from x to k)*/
kn = 2.0*PI/L;			/*Fundamental mode for this box size*/
Nk = (int) nd/4;		/*Number of k bins for the final power spectrum*/

/*Alloc the grids*/
grid1 = (double *)malloc(ng*sizeof(double));
grid2 = (double *)malloc(ng*sizeof(double)); 

/*Initialize the grid quantities*/
for(i=0;i<ng;i++){
	grid1[i] = 0.0;
	grid2[i] = 0.0;
}

/*Open the first file*/
first = fopen(firstfile,"rb");
if (first == NULL) {
	printf("Unable to open %s\n",firstfile);
	exit(0);
}

/*Read the total number of particles and the box size in the matter calogue*/
fread(&np, sizeof(int), 1, first);
//fscanf(first, "%d", &np);

/*Save the information of the particle file*/
printf("Reading the catalogue and constructin the density grids\n");
Mtot1 = 0.0;
Mtot2 = 0.0;
for(i=0;i<np;i++){
	//if(i%1000000 == 0) 	printf("i = %d of %d with Mbar = %lf\n", (int)floor(i/1000000), (int)floor(np/1000000), Mtot1/((double) i));
	
	for(j=0;j<3;j++){
		fread(&tmpf, sizeof(float), 1, first);
		X[j] = (double) tmpf;
	}
	for(j=0;j<3;j++)
		fread(&tmpf, sizeof(float), 1, first);

	if(window == 0){
		NGP(grid1, X, nd, L, 1.0);
		Mtot1 += 1.0;
	}
	if(window == 1){
		CIC(grid1, X, nd, L, 1.0);
		Mtot1 += 1.0;
	}
	if(window == 2)
		Mtot1 += Sphere(grid1, X, nd, L, Rsmooth, 1.0);
	if(window == 3)
		Mtot1 += Exp(grid1, X, nd, L, Rsmooth, R_times, 1.0);


	for(j=0;j<3;j++)
		X[j] = cysumf(X[j], -Lb/2.0, L);

	if(window == 0){
		NGP(grid2, X, nd, L, 1.0);
		Mtot2 += 1.0;
	}
	if(window == 1){
		CIC(grid2, X, nd, L, 1.0);
		Mtot2 += 1.0;
	}
	if(window == 2)
		Mtot2 += Sphere(grid2, X, nd, L, Rsmooth, 1.0);
	if(window == 3)
		Mtot2 += Exp(grid2, X, nd, L, Rsmooth, R_times, 1.0);
}
fclose(first);

/*Compute the overdensity (delta)*/
n_bar1 = Mtot1/(nd*nd*nd);
n_bar2 = Mtot2/(nd*nd*nd);
for(i=0;i<ng;i++){
	grid1[i] =  grid1[i]/n_bar1 - 1.0;
	grid2[i] =  grid2[i]/n_bar2 - 1.0;
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
printf("Mean = %lf and std = %lf\n", kx, ky);

/*Save a slice of the density field*/
printf("Saving a density slice\n");
sprintf(outfile,"Density_Slice_%s.dat", argv[2]);
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

for(i=0;i<nd;i++){
	for(j=0;j<nd;j++)
		fprintf(out, "%lf ", grid1[(int)floor(nd/2)*nd*nd + i*nd + j]);
	fprintf(out, "\n");
}
fclose(out);

/*Compute the PDF of the lognormal density field*/
printf("Computing the PDF of the lognormal field\n");
PDF = (int *)malloc(2000*sizeof(int));
for(i=0;i<2000;i++)	PDF[i] = 0;

/*Construct the lognormal field*/
cont = 0;
for(i=0;i<ng;i++){
	if(grid1[i] <= -1.0 || grid2[i] <= -1.0){
		cont ++;
		kx = log(1e-10);
	}
	else{
		kx = log(grid1[i] + 1.0);
		tmp = (int) ((kx + 10.0)/0.01);
		if(tmp >= 0 && tmp < 2000)
			PDF[tmp] += 1;
	}
}
printf("There are %d cells with delta = -1\n", cont);

/*Save the PDF*/
sprintf(outfile,"PDF_%s.dat", argv[2]);
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}
for(i=0;i<2000;i++)
	fprintf(out, "%lf %d\n", -10.0 + (i+0.5)*0.01, PDF[i]);
fclose(out);

/*Save the density grids*/
sprintf(outfile,"Density_Grid_%s.dat", argv[2]);
out = fopen(outfile,"wb");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

for(i=0;i<ng;i++)
    fwrite(&grid1[i], sizeof(double), 1, out);
for(i=0;i<ng;i++)
    fwrite(&grid2[i], sizeof(double), 1, out);  
fclose(out);

/*Free the density grids*/
free(grid1);
free(grid2);

return 0;
}