/*
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * and reimplementation by Scott B. Baden, UCSD
 *
 * Modified and  restructured by Didem Unat, Koc University
 *
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>
using namespace std;

void checkCUDAError(const char *msg);

// Utilities
//

// Timer
// Make successive calls and take a difference to get the elapsed time.
static const double kMicro = 1.0e-6;
static const int BLOCKSIZE = 16;

double getTime()
{
    struct timeval TV;
    struct timezone TZ;

    const int RC = gettimeofday(&TV, &TZ);
    if(RC == -1) {
            cerr << "ERROR: Bad call to gettimeofday" << endl;
            return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()

// Allocate a 2D array
double **alloc2D(int m,int n){
   double **E;
   int nx=n, ny=m;
   E = (double**)malloc(sizeof(double*)*ny + sizeof(double)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++)
     E[j] = (double*)(E+ny) + j*nx;
   return(E);
}

double *flatten(double **array, int m, int n)
{
  double *a;
  a = (double*)malloc(sizeof(double)*(m+2)*(n+2));
  int i, j;
  for(j=0;j<=m + 1; j++){
    for (i = 0; i <= n + 1; i++) {
      a[(j * (n+2)) + i] = array[j][i];
    }
  }
  return a;
}

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
 double stats(double *E, int m, int n, double *_mx){
     double mx = -1;
     double l2norm = 0;
     int i, j;
     for (j=1; j<=m; j++)
       for (i=1; i<=n; i++) {
	   l2norm += E[(j * (n+2)) + i]*E[(j * (n+2)) + i];
	   if (E[(j * (n+2)) + i] > mx)
	       mx = E[(j * (n+2)) + i];
      }
     *_mx = mx;
     l2norm /= (double) ((m)*(n));
     l2norm = sqrt(l2norm);
     return l2norm;
 }
 void checkCUDAError(const char *msg)
 {
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err)
  {
  fprintf(stderr, "Cuda error: %s: %s.\n", msg,
  cudaGetErrorString( err) );
  exit(EXIT_FAILURE);
  }
 }

// External functions
extern "C" {
    void splot(double **E, double T, int niter, int m, int n);
}
void cmdLine(int argc, char *argv[], double& T, int& n, int& px, int& py, int& plot_freq, int& no_comm, int&num_threads);

__global__
void vecODEKernel(double* R, double* E, double epsilon, double M1, double M2, double dt, double kk, double a, double b, int n)
{
  int row = blockIdx.y*blockDim.y+threadIdx.y+1;
  int col = blockIdx.x*blockDim.x+threadIdx.x+1;

  if((row < n) && (col < n)) {
    row = row * (n+2);
    E[row + col] = E[row + col] -dt*(kk* E[row + col]*(E[row + col] - a)*(E[row + col]-1)+ E[row + col] *R[row + col]);

    R[row + col] = R[row + col] + dt*(epsilon+M1* R[row + col]/( E[row + col]+M2))*(-R[row + col]-kk* E[row + col]*(E[row + col]-b-1));
  }
}

__global__
void boundaryKernel(double *E_prev, int m, int n)
{
    int row = blockIdx.y*blockDim.y+threadIdx.y+1;
    int col = blockIdx.x*blockDim.x+threadIdx.x+1;
    row = row * (n+2);

    E_prev[row] = E_prev[row + 2];
    E_prev[row + n + 1] = E_prev[row + (n-1)];
    E_prev[col] = E_prev[col+2];
    E_prev[(m+1)*(n+2) + col] = E_prev[(m-1)*(n+2) + col];
}

__global__
void matAllKernel(double alpha, double* E, double* E_prev, double* R, int n, int m, double epsilon, double M1, double M2, double dt, double kk, double a, double b)
{
  int row = blockIdx.y*blockDim.y+threadIdx.y;
  int col = blockIdx.x*blockDim.x+threadIdx.x;
  int row_m = row * (n+2);

  // Mirror boundary setup
  if(col == 0 || col == (n+1)) {
    E_prev[row_m] = E_prev[row_m + 2];
    E_prev[row_m + n + 1] = E_prev[row_m + (n-1)];
  }
  if(row == 0 || row == (n+1)) {
    E_prev[col] = E_prev[col+2];
    E_prev[(m+1)*(n+2) + col] = E_prev[(m-1)*(n+2) + col];
  }

  __syncthreads();
  row = row + 1;
  col = col + 1;

  if((row < n) && (col < n)) {
    row = row * (n+2);
    //PDE
    E[row + col] = E_prev[row + col]+alpha*(E_prev[row + col + 1]+E_prev[row + col -1]-4*E_prev[row + col]+E_prev[row + col + (n+2)]+E_prev[row + col - (n+2)]);
    //ODE
    E[row + col] = E[row + col] -dt*(kk* E[row + col]*(E[row + col] - a)*(E[row + col]-1)+ E[row + col] *R[row + col]);
    R[row + col] = R[row + col] + dt*(epsilon+M1* R[row + col]/( E[row + col]+M2))*(-R[row + col]-kk* E[row + col]*(E[row + col]-b-1));
  }
}

void simulate (double* E,  double* E_prev,double* R,
	       const double alpha, const int n, const int m, const double kk,
	       const double dt, const double a, const double epsilon,
	       const double M1,const double  M2, const double b)
{
    dim3 DimBlock(BLOCKSIZE,BLOCKSIZE,1);
    dim3 DimGrid(ceil((double)n/DimBlock.x), ceil((double)n/DimBlock.y));

    matAllKernel<<<DimGrid, DimBlock>>>(alpha, E, E_prev, R, n, m, epsilon, M1, M2, dt, kk, a, b);
}


// Main program
int main (int argc, char** argv)
{
  /*
   *  Solution arrays
   *   E is the "Excitation" variable, a voltage
   *   R is the "Recovery" variable
   *   E_prev is the Excitation variable for the previous timestep,
   *      and is used in time integration
   */
  double **E, **R, **E_prev;

  // Various constants - these definitions shouldn't change
  const double a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;

  double T=1000.0;
  int m=200,n=200;
  int plot_freq = 0;
  int px = 1, py = 1;
  int no_comm = 0;
  int num_threads=1;

  cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);
  m = n;
  // Allocate contiguous memory for solution arrays
  // The computational box is defined on [1:m+1,1:n+1]
  // We pad the arrays in order to facilitate differencing on the
  // boundaries of the computation box
  E = alloc2D(m+2,n+2);
  E_prev = alloc2D(m+2,n+2);
  R = alloc2D(m+2,n+2);

  int i,j;
  // Initialization
  for (j=1; j<=m; j++)
    for (i=1; i<=n; i++)
      E_prev[j][i] = R[j][i] = 0;

  for (j=1; j<=m; j++)
    for (i=n/2+1; i<=n; i++)
      E_prev[j][i] = 1.0;

  for (j=m/2+1; j<=m; j++)
    for (i=1; i<=n; i++)
      R[j][i] = 1.0;

  double *Ef, *Rf, *E_prevf;

  Ef = flatten(E, m, n);
  Rf = flatten(R, m, n);
  E_prevf = flatten(E_prev, m, n);

  double dx = 1.0/n;

  // For time integration, these values shouldn't change
  double rp= kk*(b+1)*(b+1)/4;
  double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
  double dtr=1/(epsilon+((M1/M2)*rp));
  double dt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
  double alpha = d*dt/(dx*dx);

  cout << "Grid Size       : " << n << endl;
  cout << "Duration of Sim : " << T << endl;
  cout << "Time step dt    : " << dt << endl;
  cout << "Process geometry: " << px << " x " << py << endl;
  if (no_comm)
  cout << "Communication   : DISABLED" << endl;

  cout << endl;

  // Integer timestep number
  int niter=0;

  int size = ((n+2)*(m+2) * sizeof(double));
  double *d_E, *d_E_prev, *d_R;

  // allocate memory for the devices
  cudaMalloc((void **) &d_E, size);
  cudaMalloc((void **) &d_E_prev, size);
  cudaMalloc((void **) &d_R, size);
  checkCUDAError("Error allocating device memory arrays");

  // copy all arrays to device
  cudaMemcpy(d_R, Rf, size, cudaMemcpyHostToDevice);
  checkCUDAError("Unable to copy to device, R");

  cudaMemcpy(d_E_prev, E_prevf, size, cudaMemcpyHostToDevice);
  checkCUDAError("Unable to copy to device, E_prev");

  cudaMemcpy(d_E, Ef, size, cudaMemcpyHostToDevice);
  checkCUDAError("Unable to copy to device, E");

  // Simulated time is different from the integer timestep number
  // Simulated time
  double t = 0.0;

  // Start the timer
  double t0 = getTime();

  while (t<T) {

    t += dt;
    niter++;

    simulate(d_E, d_E_prev, d_R, alpha, n, m, kk, dt, a, epsilon, M1, M2, b);

    //swap current E with previous E
    double *tmp = d_E; d_E = d_E_prev; d_E_prev = tmp;

    if (plot_freq){
      int k = (int)(t/plot_freq);
      if ((t - k * plot_freq) < dt){
        splot(E,t,niter,m+2,n+2);
      }
    }
  }//end of while loop

  double time_elapsed = getTime() - t0;

  // copy back all arrays
  cudaMemcpy(E_prevf, d_E_prev, size, cudaMemcpyDeviceToHost);
  checkCUDAError("Unable to retrieve result from device, E_prev");
  cudaMemcpy(Rf, d_R, size, cudaMemcpyDeviceToHost);
  checkCUDAError("Unable to retrieve result from device, R");
  cudaMemcpy(Ef, d_E, size, cudaMemcpyDeviceToHost);
  checkCUDAError("Unable to retrieve result from device, E");

  // free memory
  cudaFree(d_R); cudaFree(d_E); cudaFree(d_E_prev);

  double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed ;
  double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;

  cout << "Number of Iterations        : " << niter << endl;
  cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
  cout << "Sustained Gflops Rate       : " << Gflops << endl;
  cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl;

  double mx;
  double l2norm = stats(E_prevf,m,n,&mx);
  cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;

  if (plot_freq){
    cout << "\n\nEnter any input to close the program and the plot..." << endl;
    getchar();
  }

  free (E);
  free (E_prev);
  free (R);

  return 0;
}
