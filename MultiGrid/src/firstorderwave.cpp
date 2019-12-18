#include <firstorderwave.h>
#include <operators.h>
#include <iostream>
#include <waveparameters.h>
#include <cmath>

// FirstOrderWave {{{
FirstOrderWave::FirstOrderWave(Domain& d, Solver& s) : ODE(3, 1){
  if(d.getGhostPoints() < 2){
    std::cerr << "Warning: domain has fewer ghost points than expected. Solution may not behave correctly.\n";
  }
  domain = &d;
  solver = &s;

  // Set some default parameters.
  params = new WaveParameters();

  reallocateData();
}
// }}}

// ~FirstOrderWave {{{
FirstOrderWave::~FirstOrderWave(){
  delete params;
}
// }}}

// rhs {{{
void FirstOrderWave::rhs(const Grid& grid, double **u, double **dudt){
  // Check that the grid is actually big enough.
  if(grid.getSize() < 5){
    printf("Grid is too small. Need at least 5 points.\n");
    return;
  }
  // Define some variables we'l need.
  double stencil3[3] = {0.0, 0.0, 0.0};
  double stencil5[5] = {0.0, 0.0, 0.0};
  double dx = grid.getSpacing();
  int shp = grid.getSize();

  // Calculate the left boundary. The leftmost point needs to use a one-sided derivative operator.
  // The second leftmost needs to use the centered second-order operator. The third point can use the
  // full five-point stencil.

  // Leftmost point
  dudt[0][U_PHI] = u[0][U_PI];
  stencil3[0] = u[0][U_CHI];
  stencil3[1] = u[1][U_CHI];
  stencil3[2] = u[2][U_CHI];
  dudt[0][U_PI] = operators::dx_2off(stencil3, dx);
  stencil3[0] = u[0][U_PI];
  stencil3[1] = u[1][U_PI];
  stencil3[2] = u[2][U_PI];
  dudt[0][U_CHI] = operators::dx_2off(stencil3, dx);

  // Second leftmost point
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

  // Second rightmost point
  dudt[shp - 2][U_PHI] = u[shp - 2][U_PI];
  stencil3[0] = u[shp - 3][U_CHI];
  stencil3[1] = u[shp - 2][U_CHI];
  stencil3[2] = u[shp - 1][U_CHI];
  dudt[shp - 2][U_PI] = operators::dx_2(stencil3, dx);
  stencil3[0] = u[shp - 3][U_PI];
  stencil3[1] = u[shp - 2][U_PI];
  stencil3[2] = u[shp - 1][U_PI];
  dudt[shp - 2][U_CHI] = operators::dx_2(stencil3, dx);

  // Rightmost point
  dudt[shp - 1][U_PHI] = u[shp - 1][U_PI];
  stencil3[0] = u[shp - 3][U_CHI];
  stencil3[1] = u[shp - 2][U_CHI];
  stencil3[2] = u[shp - 1][U_CHI];
  dudt[shp - 1][U_PI] = operators::dx_2off(stencil3, dx);
  stencil3[0] = u[shp - 3][U_PI];
  stencil3[1] = u[shp - 2][U_PI];
  stencil3[2] = u[shp - 1][U_PI];
  dudt[shp - 1][U_CHI] = operators::dx_2off(stencil3, dx);

  // We'll need to do something for Kreiss-Oliger dissipation.
  applyKODiss(grid, u, dudt);
}
// }}}

// applyKODiss {{{
void FirstOrderWave::applyKODiss(const Grid& grid, double **u, double **dudt){
  WaveParameters *wp = (WaveParameters*) params;
  double koSigma = wp->getKOSigma();
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
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int j = 0; j < 4; j++){
      stencil4[j] = u[j][m];
    }
    dudt[0][m] += koSigma * operators::ko_dx_off1(stencil4, dx);
  }

  // Apply KO dissipation to the second leftmost point.
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int j = 0; j < 5; j++){
      stencil5[j] = u[j][m];
    }
    dudt[1][m] += koSigma * operators::ko_dx_off2(stencil5, dx);
  }

  // Apply KO dissipation to the third leftmost point.
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int j = 0; j < 6; j++){
      stencil6[j] = u[j][m];
    }
    dudt[2][m] += koSigma * operators::ko_dx_off3(stencil6, dx);
  }

  // Apply KO dissipation to the interior points.
  for(unsigned int i = 3; i < shp - 3; i++){
    for(unsigned int m = 0; m < nEqs; m++){
      for(unsigned int j = 0; j < 7; j++){
        stencil7[j] = u[i - 3 + j][m];
      }
      dudt[i][m] += koSigma * operators::ko_dx(stencil7, dx);
    }
  }

  // Apply KO dissipation to the third rightmost point.
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int j = 0; j < 6; j++){
      stencil6[j] = u[shp - 1 - j][m];
    }
    dudt[shp - 3][m] += koSigma * operators::ko_dx_off3(stencil6, dx);
  }

  // Apply KO dissipation to the second rightmost point.
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int j = 0; j < 5; j++){
      stencil5[j] = u[shp - 1 - j][m];
    }
    dudt[shp - 2][m] += koSigma * operators::ko_dx_off2(stencil5, dx);
  }

  // Apply KO dissipation to the rightmost point.
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int j = 0; j < 4; j++){
      stencil4[j] = u[shp - 1 - j][m];
    }
    dudt[shp - 1][m] += koSigma * operators::ko_dx_off1(stencil4, dx);
  }
}
// }}}

// applyBoundaries {{{
void FirstOrderWave::applyBoundaries(){
  unsigned int nb = domain->getGhostPoints();
  auto left_it = data.begin();
  auto right_it = --data.end();

  double **left = left_it->getData();
  double **right = right_it->getData();
  unsigned int nr = right_it->getGrid().getSize();

  // For now, we'll just apply fixed boundaries.
  for(int i = 0; i < nb; i++){
    // Set the left boundary.
    //left[i][U_PHI] = right[nr - 1 - nb - i][U_PHI];
    //left[i][U_PI] = right[nr - 1 - nb - i][U_PI];
    //left[i][U_CHI] = right[nr - 1 - nb - i][U_CHI];

    // Set the right boundary.
    //right[nr - 1 - i][U_PHI] = left[nb + i][U_PHI];
    //right[nr - 1 - i][U_PI] = left[nb + i][U_PI];
    //right[nr - 1 - i][U_CHI] = left[nb + i][U_CHI];

    left[i][U_PHI] = 0.0;
    left[i][U_PI] = 0.0;
    left[i][U_CHI] = 0.0;

    right[nr - 1 - i][U_PHI] = 0.0;
    right[nr - 1 - i][U_PI] = 0.0;
    right[nr - 1 - i][U_CHI] = 0.0;
  }
}
// }}}

// initData {{{
void FirstOrderWave::initData(){
  // Cast our parameters to wave parameters.
  WaveParameters *wp = (WaveParameters*) params;

  switch(wp->getInitialConditions()){
    case WaveParameters::GAUSSIAN:
      applyGaussian();
      break;
    default:
      std::cerr << "Warning: unrecognized initial conditions selected. Defaulting to Gaussian.\n";
      applyGaussian();
      break;
  }
}
// }}}

// applyGaussian {{{
void FirstOrderWave::applyGaussian(){
  // First, let's identify the center.
  double x0 = 0.5*(domain->getBounds()[1] + domain->getBounds()[0]);

  // Next, let's loop through every grid and start assigning points.
  for(auto it = data.begin(); it != data.end(); ++it){
    const double *x = it->getGrid().getPoints();
    unsigned int nx = it->getGrid().getSize();
    double **u = it->getData();
    for(unsigned int i = 0; i < nx; i++){
      double val = std::exp(-(x[i] - x0)*(x[i] - x0)*64.0);
      u[i][U_PHI] = val;
      u[i][U_PI ] = 0.0;
      u[i][U_CHI] = -128.0*(x[i] - x0)*val;
    }
  }
}
// }}}
