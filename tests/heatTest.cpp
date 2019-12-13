#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <operators.h>
#include <cmath>
#include <cstdio>
#include <output.h>

const double sigma = 1e-2;

// rhsHeat {{{

/**
 * A righthand side for the one-dimensional heat equation.
 */
void rhsHeat(const Grid& grid, double **u, double **dudt){
  // We need at least five points on our grid. If not, throw an error and
  // quit.
  if (grid.getSize() < 5){
    printf("Grid is too small. Need at least 5 points.\n");
    return;
  }

  // Define some variables we'll need.
  double stencil3[3] = {0.0, 0.0, 0.0};
  double stencil5[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double dx = grid.getSpacing();
  int shp = grid.getSize();

  // Calculate the left boundary.
  // Leftmost point.
  stencil3[0] = u[0][0];
  stencil3[1] = u[1][0];
  stencil3[2] = u[2][0];
  dudt[0][0] = sigma*operators::dxx_2off(stencil3, dx);

  // Second leftmost point.
  dudt[1][0] = sigma*operators::dxx_2(stencil3, dx);

  // Interior points.
  /*for (unsigned int i = 2; i < shp - 2; i++){
    for(int j = 0; j < 5; j++){
      stencil5[j] = u[i - 2 + j][0];
    }
    dudt[i][0] = sigma*operators::dxx_4(stencil5, dx);
  }*/
  for (unsigned int i = 2; i < shp - 2; i++){
    for(int j = 0; j < 3; j++){
      stencil3[j] = u[i - 1 + j][0];
    }
    dudt[i][0] = sigma*operators::dxx_2(stencil3, dx);
  }

  // Calculate the right boundary.
  // Rightmost point.
  stencil3[0] = u[shp - 1][0];
  stencil3[1] = u[shp - 2][0];
  stencil3[2] = u[shp - 3][0];
  dudt[shp - 1][0] = sigma*operators::dxx_2off(stencil3, dx);

  // Second rightmost point.
  dudt[shp - 2][0] = sigma*operators::dxx_2(stencil3, dx);
}
// }}}

// applyBoundaries {{{
void applyBoundaries(Domain& domain, Grid& grid, double **data){
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();

  // We need du/dx = 0 at the boundaries.
  //data[nb - 1][0] = data[nb + 1][0];
  //data[shp - nb][0] = data[shp - nb - 2][0];
  for (int i = 0; i < nb; i++){
    data[i][0] = data[nb + 1 + i][0];
    data[shp - 1 - i][0] = data[shp - nb - 2 - i][0];
  }
}
// }}}

// calcL2Norm {{{
double calcL2Norm(Domain& domain, Grid& grid, double *data){
  unsigned int nb = domain.getGhostPoints();
  unsigned int shp = grid.getSize();
  const double *bounds = grid.getBounds();
  double L = (bounds[1] - bounds[0]);
  double dx = grid.getSpacing();

  // Loop over only the physical grid points, then calculate the L2 norm.
  double sum = 0.0;
  double err = 0.0;
  double old;
  double term;
  for(unsigned int i = nb; i < shp - nb - 1; i++){
    term = data[i]*data[i];
    old = sum;
    sum += term + err;
    err = term - ((sum - old) - err);
  }
  sum = std::sqrt(dx*sum/L);
  
  return sum;
}
// }}}

// calcGaussian {{{
/**
 * Calculate the amplitude of a Gaussian at x.
 */
double calcGaussian(double x){
  return std::exp(-(x - 0.5)*(x - 0.5)*64.0);
}

// }}}

// writeToCSV {{{
void writeToCSV(char *file, double t, double x){
  FILE *fp;

  fp = fopen(file, "a");
  fprintf(fp, "%20e, %20e\n",t, x);
  fclose(fp);
}
// }}}

int main(int argc, char *argv[]){
  RK4 rk4 = RK4();

  // Construct a domain and a grid.
  Domain domain = Domain();
  domain.setGhostPoints(3);
  int N = 401;
  Grid grid = Grid(domain.getBounds(), N, domain.getGhostPoints());
  domain.setCFL(0.25);
  double dt = domain.getCFL() * grid.getSpacing();
  double ti = 0.0;
  double tf = 10.0;
  int nb = domain.getGhostPoints();
  int shp = grid.getSize();

  // Clear our norm data files.
  fclose(fopen("norm_heat.csv","w"));

  // Memory allocation {{{
  double **u = new double*[grid.getSize()];
  double **uint = new double*[grid.getSize()];
  for(int i = 0; i < grid.getSize(); i++){
    u[i] = new double[1];
    uint[i] = new double[1];
    u[i][0] = 0.0;
    uint[i][0] = 0.0;
  }
  double ***k = new double**[4];
  for(int i = 0; i < 4; i++){
    k[i] = new double*[grid.getSize()];
    for(int j = 0; j < grid.getSize(); j++){
      k[i][j] = new double[1];
      k[i][j][0] = 0.0;
    }
  }
  double *T = new double[grid.getSize()];
  double *r = new double[grid.getSize()];
  for(int i = 0; i < grid.getSize(); i++){
    r[i] = grid.getPoints()[i];
  }
  // }}}

  // Heat equation {{{
  for(int i = 0; i < grid.getSize(); i++){
    u[i][0] = calcGaussian(r[i]);
  }

  // Perform the integration
  int Nt = (int) (tf / dt);
  for(int n = 0; n < Nt; n++){
    double t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      rk4.calcStage(&rhsHeat, u, uint, k[i], grid, dt, 1, i);
      applyBoundaries(domain, grid, uint);
    }
    rk4.combineStages(k, u, grid, dt, 1);
    applyBoundaries(domain, grid, u);
    for(int i = 0; i < grid.getSize(); i++){
      T[i] = u[i][0];
    }
    output_data("T", T, r, grid.getSize(), t);
    double norm = calcL2Norm(domain, grid, T);
    writeToCSV("norm_heat.csv", t, norm);
  }
  // }}}

  // Memory cleanup {{{
  delete[] r;
  delete[] T;
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
