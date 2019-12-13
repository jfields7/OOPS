#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <operators.h>
#include <cmath>
#include <cstdio>
#include <output.h>
#include <firstorderwave.h>

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

int main(int argc, char* argv[]){
  // Construct a domain and a set of grids.
  Domain domain = Domain();
  domain.setGhostPoints(3);
  int N0 = 65;
  int N = N0 - 1;
  int nb = domain.getGhostPoints();
  
  int ngrids = 1;
  while (N > 2){
    N = N >> 1;
    ngrids++;
  }
  double dx_grid = (domain.getBounds()[1] - domain.getBounds()[0])/ngrids;
  for(int i = 0; i < ngrids; i++){
    double bounds[2] = {dx_grid * i, dx_grid * (i + 1)};
    domain.addGrid(bounds, ((N0 - 1) >> i) + 1);
  }

  // First, let's confirm that the grids are all placed correctly.
  for(auto it = domain.getGrids().begin(); it != domain.getGrids().end(); ++it){
    const double *bounds = it->getBounds();
    unsigned int shp = it->getSize();
    double dx = it->getSpacing();
    
    // Print out all the information we collected on this grid.
    printf("Grid:\n");
    printf("  Bounds  = [%g, %g]\n",bounds[0], bounds[1]);
    printf("  Size    = %d\n",shp);
    printf("  Spacing = %g\n", dx);
  }

  // Now we need to try to construct our ODE system.
  RK4 rk4 = RK4();
  FirstOrderWave ode = FirstOrderWave(domain, rk4);
  return 0;
}
