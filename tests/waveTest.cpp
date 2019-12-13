#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <operators.h>
#include <cmath>
#include <cstdio>
#include <output.h>
#include <cubic.h>

const double EPSILON = 1.e-15;
const double sigma1st = 0.0;
const double sigma2nd = 0.0;

const double PI = 3.14159265358979;

// Our variables. We define phi as the amplitude of the wave, 
// pi = dphi/dt, and chi = dphi/dx.
const unsigned int U_PHI = 0;
const unsigned int U_PI  = 1;
const unsigned int U_CHI = 2;
const unsigned int VARS1 = 3;
const unsigned int VARS2 = 2;

void applyKODiss(Grid& grid, double **u, double **dudt, const double sigma, const unsigned int nvars);

// rhsSecondOrderWave {{{
/**
 * The righthand side for the wave equation split into a two-variable system
 * that is first order in time and second order in space.
 */
void rhsSecondOrderWave(Grid& grid, double **u, double **dudt){
  // We need at least five points on our grid. If not, throw an error and quit.
  if(grid.getSize() < 5){
    printf("Grid is too small. Need at least 5 points.\n");
    return;
  }
  // Define some variables we'll need.
  double stencil3[3] = {0.0, 0.0, 0.0};
  double stencil4[4] = {0.0, 0.0, 0.0, 0.0};
  double stencil5[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double dx = grid.getSpacing();
  int shp = grid.getSize();

  // Calculate the left boundary. The leftmost point needs to use a one-sided
  // derivative operator. The second leftmost point needs to use a centered
  // second-order operator. The third point can use the full five-point stencil.

  // Leftmost point.
  dudt[0][U_PHI] = u[0][U_PI];
  stencil4[0] = u[0][U_PHI];
  stencil4[1] = u[1][U_PHI];
  stencil4[2] = u[2][U_PHI];
  stencil4[3] = u[3][U_PHI];
  dudt[0][U_PI] = operators::dxx_2off(stencil4, dx);

  // Second leftmost point.
  dudt[1][U_PHI] = u[1][U_PI];
  stencil3[0] = u[0][U_PHI];
  stencil3[1] = u[1][U_PHI];
  stencil3[2] = u[2][U_PHI];
  dudt[1][U_PI] = operators::dxx_2(stencil3, dx);

  // Loop through all the interior points.
  for(int i = 2; i < shp - 2; i++){
    dudt[i][U_PHI] = u[i][U_PI];
    for(int j = 0; j < 5; j++){
      stencil5[j] = u[i - 2 + j][U_PHI];
    }
    dudt[i][U_PI] = operators::dxx_4(stencil5, dx);
  }
  /*for(int i = 2; i < shp - 2; i++){
    dudt[i][U_PHI] = u[i][U_PI];
    for(int j = 0; j < 3; j++){
      stencil3[j] = u[i - 1 + j][U_PHI];
    }
    dudt[i][U_PI] = operators::dxx_2(stencil3, dx);
  }*/

  // Second rightmost point.
  dudt[shp - 2][U_PHI] = u[shp - 2][U_PI];
  stencil3[0] = u[shp - 3][U_PHI];
  stencil3[1] = u[shp - 2][U_PHI];
  stencil3[2] = u[shp - 1][U_PHI];
  dudt[shp - 2][U_PI] = operators::dxx_2(stencil3, dx);

  // Rightmost point. Make sure to fill the stencil backward.
  dudt[shp - 1][U_PHI] = u[shp - 1][U_PI];
  stencil4[3] = u[shp - 1][U_PHI];
  stencil4[2] = u[shp - 2][U_PHI];
  stencil4[1] = u[shp - 3][U_PHI];
  stencil4[0] = u[shp - 4][U_PHI];
  dudt[shp - 1][U_PI] = operators::dxx_2off(stencil4, dx);

  applyKODiss(grid, u, dudt, sigma2nd, VARS2);
}
// }}}

// rhsFirstOrderWave {{{
/**
 * The righthand side for the wave equation split into a three-variable system
 * that is first order in space and time.
 */
void rhsFirstOrderWave(Grid& grid, double **u, double **dudt){
  // We need at least five points on our grid. If not, throw an error and quit.
  if(grid.getSize() < 5){
    printf("Grid is too small. Need at least 5 points.\n");
    return;
  }
  // Define some variables we'll need.
  double stencil3[3] = {0.0, 0.0, 0.0};
  double stencil5[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double dx = grid.getSpacing();
  int shp = grid.getSize();

  // Calculate the left boundary. The leftmost point needs to use a one-sided
  // derivative operator. The second leftmost needs to use the centered second-order
  // operator. The third point can use the full five-point stencil.

  // Leftmost point.
  dudt[0][U_PHI] = u[0][U_PI];
  stencil3[0] = u[0][U_CHI];
  stencil3[1] = u[1][U_CHI];
  stencil3[2] = u[2][U_CHI];
  dudt[0][U_PI] = operators::dx_2off(stencil3, dx);
  stencil3[0] = u[0][U_PI];
  stencil3[1] = u[1][U_PI];
  stencil3[2] = u[2][U_PI];
  dudt[0][U_CHI] = operators::dx_2off(stencil3, dx);

  // Second leftmost point.
  dudt[1][U_PHI] = u[1][U_PI];
  stencil3[0] = u[0][U_CHI];
  stencil3[1] = u[1][U_CHI];
  stencil3[2] = u[2][U_CHI];
  dudt[1][U_PI] = operators::dx_2(stencil3, dx);
  stencil3[0] = u[0][U_PI];
  stencil3[1] = u[1][U_PI];
  stencil3[2] = u[2][U_PI];
  dudt[1][U_CHI] = operators::dx_2(stencil3, dx);

  // Now loop through all the interior points.
  for(int i = 2; i < shp - 2; i++){
    dudt[i][U_PHI] = u[i][U_PI];

    for(int j = 0; j < 5; j++){
      stencil5[j] = u[i - 2 + j][U_CHI];
    }
    dudt[i][U_PI] = operators::dx_4(stencil5, dx);

    for(int j = 0; j < 5; j++){
      stencil5[j] = u[i - 2 + j][U_PI];
    }
    dudt[i][U_CHI] = operators::dx_4(stencil5, dx);
  }
  /*for(int i = 2; i < shp - 2; i++){
    dudt[i][U_PHI] = u[i][U_PI];

    for(int j = 0; j < 3; j++){
      stencil3[j] = u[i - 1 + j][U_CHI];
    }
    dudt[i][U_PI] = operators::dx_2(stencil3, dx);

    for(int j = 0; j < 3; j++){
      stencil3[j] = u[i - 1 + j][U_PI];
    }
    dudt[i][U_CHI] = operators::dx_2(stencil3, dx);
  }*/

  // Second rightmost point.
  dudt[shp - 2][U_PHI] = u[shp - 2][U_PI];
  stencil3[0] = u[shp - 3][U_CHI];
  stencil3[1] = u[shp - 2][U_CHI];
  stencil3[2] = u[shp - 1][U_CHI];
  dudt[shp - 2][U_PI] = operators::dx_2(stencil3, dx);
  stencil3[0] = u[shp - 3][U_PI];
  stencil3[1] = u[shp - 2][U_PI];
  stencil3[2] = u[shp - 1][U_PI];
  dudt[shp - 2][U_CHI] = operators::dx_2(stencil3, dx);

  // Rightmost point. Make sure to fill the stencils in reverse
  // order to flip the edge operator.
  dudt[shp - 1][U_PHI] = u[shp - 1][U_PI];
  stencil3[2] = u[shp - 3][U_CHI];
  stencil3[1] = u[shp - 2][U_CHI];
  stencil3[0] = u[shp - 1][U_CHI];
  dudt[shp - 1][U_PI] = operators::dx_2off(stencil3, dx);
  stencil3[2] = u[shp - 3][U_PI];
  stencil3[1] = u[shp - 2][U_PI];
  stencil3[0] = u[shp - 1][U_PI];
  dudt[shp - 1][U_CHI] = operators::dx_2off(stencil3, dx);

  // We'll need to do something for Kreiss-Oliger dissipation.
  applyKODiss(grid, u, dudt, sigma1st, VARS1);
}
// }}}

// applyKODiss {{{

/**
 * Apply Kreiss-Oliger dissipation to the righthand side of a dataset.
 */
void applyKODiss(Grid& grid, double **u, double **dudt, const double sigma, const unsigned int nvars){
  // The grid needs to have at least 7 points for this to work.
  if(grid.getSize() < 7){
    printf("Grid is too small to use Kreiss-Oliger dissipation.\n");
    return;
  }

  // Set up some quantities we'll be using.
  unsigned int shp = grid.getSize();
  double dx = grid.getSpacing();
  double stencil4[4] = {0.0, 0.0, 0.0, 0.0};
  double stencil5[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double stencil6[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double stencil7[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // Apply KO dissipation to the leftmost point.
  for(unsigned int m = 0; m < nvars; m++){
    for(unsigned int j = 0; j < 4; j++){
      stencil4[j] = u[j][m];
    }
    dudt[0][m] += sigma * operators::ko_dx_off1(stencil4, dx);
  }

  // Apply KO dissipation to the second leftmost point.
  for(unsigned int m = 0; m < nvars; m++){
    for(unsigned int j = 0; j < 5; j++){
      stencil5[j] = u[j][m];
    }
    dudt[1][m] += sigma * operators::ko_dx_off2(stencil5, dx);
  }

  // Apply KO dissipation to the third leftmost point.
  for(unsigned int m = 0; m < nvars; m++){
    for(unsigned int j = 0; j < 6; j++){
      stencil6[j] = u[j][m];
    }
    dudt[2][m] += sigma * operators::ko_dx_off3(stencil6, dx);
  }

  // Apply KO dissipation to the interior points.
  for(unsigned int i = 3; i < shp - 3; i++){
    for(unsigned int m = 0; m < nvars; m++){
      for(unsigned int j = 0; j < 7; j++){
        stencil7[j] = u[i - 3 + j][m];
      }
      dudt[i][m] += sigma * operators::ko_dx(stencil7, dx);
    }
  }

  // Apply KO dissipation to the third rightmost point.
  for(unsigned int m = 0; m < nvars; m++){
    for(unsigned int j = 0; j < 6; j++){
      stencil6[j] = u[shp - 1 - j][m];
    }
    dudt[shp - 3][m] += sigma * operators::ko_dx_off3(stencil6, dx);
  }

  // Apply KO dissipation to the second rightmost point.
  for(unsigned int m = 0; m < nvars; m++){
    for(unsigned int j = 0; j < 5; j++){
      stencil5[j] = u[shp - 1 - j][m];
    }
    dudt[shp - 2][m] += sigma * operators::ko_dx_off2(stencil5, dx);
  }

  // Apply KO dissipation to the rightmost point.
  for(unsigned int m = 0; m < nvars; m++){
    for(unsigned int j = 0; j < 4; j++){
      stencil4[j] = u[shp - 1 - j][m];
    }
    dudt[shp - 1][m] += sigma * operators::ko_dx_off1(stencil4, dx);
  }
}

// }}}

// applyBoundariesFirst {{{
void applyBoundariesFirst(Domain& domain, Grid& grid, double **data){
  // We want outflow conditions.
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();

  // Periodic conditions
  for(int i = 0; i < nb; i++){
    data[i][U_PHI] = data[shp - 1 - nb - i][U_PHI];
    data[i][U_PI] = data[shp - 1 - nb - i][U_PI];
    data[i][U_CHI] = data[shp - 1 - nb - i][U_CHI];

    data[shp - 1 - i][U_PHI] = data[nb + i][U_PHI];
    data[shp - 1 - i][U_PI] = data[nb + i][U_PI];
    data[shp - 1 - i][U_CHI] = data[nb + i][U_CHI];
  }
}
// }}}

// applyBoundariesSecond {{{
void applyBoundariesSecond(Domain& domain, Grid& grid, double **data){
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();

  // Periodic boundaries
  for(int i = 0; i < nb; i++){
    data[i][U_PHI] = data[shp - 1 - nb - i][U_PHI];
    data[i][U_PI] = data[shp - 1 - nb - i][U_PI];

    data[shp - 1 - i][U_PHI] = data[nb + i][U_PHI];
    data[shp - 1 - i][U_PI] = data[nb + i][U_PI];
  }
}
// }}}

// applyAdvectionBoundary {{{

/**
 * Apply an advection boundary condition to a variable.
 * @param domain - The Domain object this system is associated with.
 * @param grid   - The particular grid this system is being solved on.
 * @param data   - The individual point data for the system on this grid.
 * @param rhs    - This needs to be an array that contains the corresponding
 *                 at the boundary points. In other words, it should be of
 *                 dimension [2][nvars].
 */
void applyAdvectionBoundary(Domain& domain, Grid& grid, double **data, double **rhs, int nvars){
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();
  double dx = grid.getSpacing();

  // Apply the boundaries.
  for(int m = 0; m < nvars; m++){
    // Right boundary
    data[shp - nb][m] = -2.0*dx*rhs[1][m] + data[shp - nb - 2][m];
    data[shp - nb + 1][m] = -4.0*dx*rhs[1][m] + data[shp - nb - 3][m];

    // Left boundary
    data[nb - 1][m] = -2.0*dx*rhs[0][m] + data[nb + 1][m];
    data[nb - 2][m] = -4.0*dx*rhs[0][m] + data[nb + 2][m];
  }

  // Fix the rest of the ghost points.
  for(int i = 0; i < nb - 1; i++){
    for(int m = 0; m < nvars; m++){
      data[i][m] = data[nb - 1][m];
      data[shp - 1 - i][m] = data[shp - nb][m];
    }
  }
}

// }}}

// applyAdvectionBoundaryFirst {{{

void applyAdvectionBoundaryFirst(Domain& domain, Grid& grid, double **data, double **rhs){
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();
  double dx = grid.getSpacing();

  // Apply the boundaries
  // Right boundary
  //data[shp - nb][U_PHI] = -2.0*dx*rhs[1][U_PHI] + data[shp - nb - 2][U_PHI];
  //data[shp - nb][U_PI ] = -2.0*dx*rhs[1][U_PI ] + data[shp - nb - 2][U_PI ];
  //data[shp - nb][U_CHI] = -2.0*dx*rhs[1][U_CHI] + data[shp - nb - 2][U_CHI];
  data[shp - nb][U_PHI] = -dx*rhs[1][U_PHI] + data[shp - nb - 1][U_PHI];
  data[shp - nb][U_PI ] = -dx*rhs[1][U_PI ] + data[shp - nb - 1][U_PI ];
  data[shp - nb][U_CHI] = -dx*rhs[1][U_CHI] + data[shp - nb - 1][U_CHI];

  // Left boundary
  //data[nb - 1][U_PHI] = -2.0*dx*rhs[0][U_PHI] + data[nb + 1][U_PHI];
  //data[nb - 1][U_PI ] = -2.0*dx*rhs[0][U_PI ] + data[nb + 1][U_PI ];
  //data[nb - 1][U_CHI] = -2.0*dx*rhs[0][U_CHI] + data[nb + 1][U_CHI];
  data[nb - 1][U_PHI] = -dx*rhs[0][U_PHI] + data[nb][U_PHI];
  data[nb - 1][U_PI ] = -dx*rhs[0][U_PI ] + data[nb][U_PI ];
  data[nb - 1][U_CHI] = -dx*rhs[0][U_CHI] + data[nb][U_CHI];

  // Fix the rest of the ghost points.
  for(int i = 0; i < nb - 1; i++){
    for(int m = 0; m < VARS1; m++){
      data[i][m] = data[nb - 1][m];
      data[shp - 1 - i][m] = data[shp - nb][m];
    }
  }
}

// }}}

// calcErrorGrid {{{

/**
 * Calculate the error across the domain by comparing the expected solution to the actual. The
 * calculation excludes ghost points.
 */
double calcErrorGrid(Domain& domain, Grid& grid, double *data, double (*actual)(double, double), 
                     double t){
  // Get some characteristics of the grid so we can figure out where to integrate.
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();
  const double *bounds = grid.getBounds();
  double L = (bounds[1] - bounds[0]);
  double dx = grid.getSpacing();
  
  // Loop over all physical points on the grid, calculate the expected value at that point,
  // and calculate the error.
  const double *points = grid.getPoints();
  double sum = 0.0;
  for(int i = 0; i < shp - 2*nb; i++){
    double local = actual(points[i], t);
    sum += (local - data[i])*(local - data[i]);
  }
  sum = std::sqrt(sum*dx/L);

  return sum;
}

// }}}

// calcL2Norm {{{

double calcL2Norm(Domain& domain, Grid& grid, double *data){
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();
  const double *bounds = grid.getBounds();
  double L = (bounds[1] - bounds[0]);
  double dx = grid.getSpacing();

  // Loop over only the physical points on the grid, then calculate the L2 norm.
  double sum = 0.0;
  double err = 0.0;
  double old;
  double term;
  for(int i = 0; i < shp - 2*nb - 1; i++){
    //term = (data[i] + data[i+1])/2.0;
    //term = term * term;
    term = data[i]*data[i];
    old = sum;
    sum += term + err;
    err = term - ((sum - old) - err);
  }
  sum = std::sqrt(dx*sum/L);

  return sum;
}
// }}}

// exactGaussianWave {{{

/**
 * Calculate the exact wave assuming Gaussian initial conditions. Assume the center of the wave
 * is 0.5 on a domain [0,1].
 */
double exactGaussianWave(double x, double t){
  double xl = x + t;
  /*while(xl > 1.0){
    xl -= 1.0;
  }*/
  double xr = x - t;
  /*while(xr < 0.0){
    xr += 1.0;
  }*/
  // Calculate the leftgoing wave.
  double left = std::exp(-(xl - 0.5)*(xl - 0.5)*64.0);
  // Calculate the rightgoing wave.
  double right = std::exp(-(xr - 0.5)*(xr - 0.5)*64.0);

  // Get the combination of the two waves.
  return 0.5*(left + right);
}

// }}}

// exactSineWave {{{
/**
 * Calculate a sine wave with a wavelength of exactly 1.
 */
double exactSineWave(double x, double t){
  double xl = x + t;
  while(xl > 1.0){
    xl -= 1.0;
  }
  double xr = x - t;
  while(xr < 0.0){
    xr += 1.0;
  }
  // Calculate the leftgoing wave.
  double left = 0.5*std::sin(2.0*PI*xl);
  // Calculate the rightgoing wave.
  double right = 0.5*std::sin(2.0*PI*xr);

  // Get the combination of the two waves.
  return 0.5*(left + right);
}
// }}}

// exportToCSV {{{
void exportToCSV(char* file, double t, double* x, Grid& grid){
  FILE *fp;

  fp = fopen(file,"a");
  const double *points = grid.getPoints();
  unsigned int shp = grid.getSize();
  for(unsigned int i = 0; i < shp; i++){
    fprintf(fp,"%20e, %20e, %20e\n", t, points[i], x[i]);
  }
  fclose(fp);
}
// }}}

// writeToCSV {{{
void writeToCSV(char* file, double t, double x){
  FILE *fp;

  fp = fopen(file,"a");
  fprintf(fp,"%20e, %20e\n",t, x);
  fclose(fp);
}
// }}}

int main(int argc, char *argv[]){
  RK4 rk4 = RK4();

  // Construct a domain and a grid.
  Domain domain = Domain();
  domain.setGhostPoints(3);
  int N = 201;
  int frequency = 4;
  Grid grid = Grid(domain.getBounds(), N, domain.getGhostPoints());
  domain.setCFL(0.25);
  double dt = domain.getCFL() * grid.getSpacing();
  double ti = 0.0;
  double tf = 1.0;
  int nb = domain.getGhostPoints();
  int shp = grid.getSize();

  // Clear our error data files.
  fclose(fopen("error_1st.csv","w"));
  fclose(fopen("error_2nd.csv","w"));
  fclose(fopen("norm_1st.csv","w"));
  fclose(fopen("norm_2nd.csv","w"));

  // Memory allocation {{{
  // Allocate memory for all our new data sets.
  double **u = new double*[grid.getSize()];
  double **uint = new double*[grid.getSize()];
  for(int i = 0; i < grid.getSize(); i++){
    u[i] = new double[VARS1];
    uint[i] = new double[VARS1];
    u[i][U_PHI] = 0.0;
    u[i][U_PI ] = 0.0;
    u[i][U_CHI] = 0.0;
    uint[i][U_PHI] = 0.0;
    uint[i][U_PI ] = 0.0;
    uint[i][U_CHI] = 0.0;
  }
  double ***k = new double**[4];
  for(int i = 0; i < 4; i++){
    k[i] = new double*[grid.getSize()];
    for(int j = 0; j < grid.getSize(); j++){
      k[i][j] = new double[VARS1];
      k[i][j][U_PHI] = 0.0;
      k[i][j][U_PI ] = 0.0;
      k[i][j][U_CHI] = 0.0;
    }
  }
  double *phi = new double[N];
  double *pi = new double[N];
  double *chi = new double[N];
  double *r = new double[N];
  for(int i = 0; i < N; i++){
    r[i] = grid.getPoints()[i + nb];
  }
  // }}}

  // First-order system {{{
  // Initialize the data. To be consistent, we need to set both PHI and CHI, which
  // is defined as the spatial derivative of PHI.
  const double *points = grid.getPoints();
  for(int i = 0; i < grid.getSize(); i++){
    double x = points[i];
    u[i][U_PHI] = exactGaussianWave(x, 0.0);
    u[i][U_CHI] = -128.0*(x - 0.5)*u[i][U_PHI];
    //u[i][U_PI] = 1.0;
    //u[i][U_PHI] = exactSineWave(points[i], 0.0);
    //u[i][U_CHI] = 2*PI*exactSineWave(points[i]+0.25,0.0);
  }
  //applyBoundariesFirst(domain,grid,uint);

  // Perform the integration.
  printf("First-order integration.\n");
  double rhs_b1[VARS1] = {0.0, 0.0, 0.0};
  double rhs_b2[VARS1] = {0.0, 0.0, 0.0};
  double **rhs_b = new double*[2];
  rhs_b[0] = rhs_b1;
  rhs_b[1] = rhs_b2;
  int Nt = (int)(tf / dt);
  for(int n = 0; n < Nt; n++){
    double t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      rk4.calcStage(&rhsFirstOrderWave, u, uint, k[i], grid, dt, VARS1, i);
      //applyBoundariesFirst(domain,grid,uint);
      for(int m = 0; m < VARS1; m++){
        rhs_b[0][m] = k[i][nb][m];
        rhs_b[1][m] = k[i][shp - nb - 1][m];
      }
      applyAdvectionBoundary(domain, grid, uint, rhs_b, VARS1);
      //applyAdvectionBoundaryFirst(domain, grid, uint, rhs_b);
    }
    // For the advective boundary conditions, we need to save the old solution.
    for(int i = 0; i < grid.getSize(); i++){
      for(int m = 0; m < VARS1; m++){
        uint[i][m] = u[i][m];
      }
    }
    rk4.combineStages(k, u, grid, dt, VARS1);
    // For the advective boundary conditions, we need to get the approximate right-hand side.
    for(int m = 0; m < VARS1; m++){
      rhs_b[0][m] = (u[nb][m] - uint[nb][m]);
      rhs_b[1][m] = (u[shp - nb - 1][m] - uint[shp - nb - 1][m]);
    }
    //applyBoundariesFirst(domain,grid,u);
    applyAdvectionBoundary(domain, grid, u, rhs_b, VARS1);
    //applyAdvectionBoundaryFirst(domain, grid, u, rhs_b);
    // Print out the error in PHI.
    for(int i = 0; i < N; i++){
      phi[i] = u[i + nb][U_PHI];
      pi[i] = u[i + nb][U_PI];
      chi[i] = u[i + nb][U_CHI];
    }
    if((n + 1) % frequency == 0){
      output_data("Phi200", phi, r, N, t + dt);
      output_data("Pi200", pi, r, N, t + dt);
      output_data("Chi200", chi, r, N, t + dt);
    }
    double e = calcErrorGrid(domain, grid, phi, &exactGaussianWave, t + dt);
    //double e = calcErrorGrid(domain, grid, phi, &exactSineWave, t + dt);
    /*double e = calcErrorGrid(domain, grid, phi, 
      [](double x, double t){
        return t;
      },
      t + dt
    );*/
    double norm = calcL2Norm(domain, grid, phi);
    //printf("t = %g, Error = %g\n", t, e);
    writeToCSV("error_1st.csv", t + dt, e);
    writeToCSV("norm_1st.csv", t + dt, norm);
  }
  // }}}

  // Second-order system {{{
  // Initialize the data.
  for(int i = 0; i < grid.getSize(); i++){
    //u[i][U_PHI] = 0.0;
    u[i][U_PHI] = exactGaussianWave(points[i], 0.0);
    //u[i][U_PHI] = exactSineWave(points[i], 0.0);
    u[i][U_PI] = 0.0;
    //u[i][U_PI] = 1.0;
  }
  //applyBoundariesFirst(domain,grid,uint);

  // Perform the integration.
  printf("Second-order integration.\n");
  for(int n = 0; n < Nt; n++){
    double t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      rk4.calcStage(&rhsSecondOrderWave, u, uint, k[i], grid, dt, VARS2, i);
      for(int m = 0; m < VARS1; m++){
        rhs_b[0][m] = k[i][nb][m];
        rhs_b[1][m] = k[i][shp - nb - 1][m];
      }
      //applyBoundariesSecond(domain, grid, uint);
      applyAdvectionBoundary(domain, grid, uint, rhs_b, VARS2);
    }
    // For the advective boundary conditions, we need to save the old solution.
    for(int i = 0; i < grid.getSize(); i++){
      for(int m = 0; m < VARS1; m++){
        uint[i][m] = u[i][m];
      }
    }
    rk4.combineStages(k, u, grid, dt, VARS2);
    for(int m = 0; m < VARS1; m++){
      rhs_b[0][m] = (u[nb][m] - uint[nb][m]);
      rhs_b[1][m] = (u[shp - nb - 1][m] - uint[shp - nb - 1][m]);
    }
    applyAdvectionBoundary(domain, grid, u, rhs_b, VARS2);
    //applyBoundariesSecond(domain, grid, u);
    // Print out the error in PHI.
    for(int i = 0; i < N; i++){
      phi[i] = u[i + nb][U_PHI];
      pi[i] = u[i + nb][U_PI];
    }
    output_data("Phi2nd", phi, r, N, t);
    output_data("Pi2nd", pi, r, N, t);

    double e = calcErrorGrid(domain, grid, phi, &exactGaussianWave, t);
    //double e = calcErrorGrid(domain, grid, phi, &exactSineWave, t);
    double norm = calcL2Norm(domain, grid, phi);
    /*double e = calcErrorGrid(domain, grid, phi, 
      [](double x, double t){
        return t;
      },
      t
    );*/

    //printf("t = %g, Error = %g\n", t, e);
    writeToCSV("error_2nd.csv",t,e);
    writeToCSV("norm_2nd.csv", t, norm);
  }
  // }}}

  // Exact solution {{{
  for(int n = 0; n < Nt; n++){
    double t = n * dt;
    for(int i = 0; i < grid.getSize(); i++){
      phi[i] = exactGaussianWave(points[i], t);
      //phi[i] = exactSineWave(points[i], t);
      //phi[i] = t;
    }
    output_data("PhiExact",phi, r, grid.getSize(), t);
  }
  // }}}

  // Memory cleanup {{{
  // Release all the memory we allocated.
  delete[] rhs_b;
  delete[] r;
  delete[] chi;
  delete[] pi;
  delete[] phi;
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < grid.getSize(); j++){
      delete[] k[i][j];
    }
    delete[] k[i];
  }
  delete[] k;
  for(int i = 0; i < grid.getSize(); i++){
    delete[] uint[i];
    delete[] u[i];
  }
  delete[] uint;
  delete[] u;
  // }}}
  return 0;
}
