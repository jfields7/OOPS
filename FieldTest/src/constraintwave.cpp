#include <constraintwave.h>
#include <operators.h>
#include <iostream>
#include <fieldtestparameters.h>
#include <cmath>
#include <memory>
#include <solverdata.h>

// ConstraintWave {{{
ConstraintWave::ConstraintWave(Domain &d, Solver &s) : ODE(3, 1){
  domain = &d;
  solver = &s;

  addField("Evolution", nEqs, true, true);
  addField("Constraint", 1, false, true);

  reallocateData();
}
// }}}

// ~ConstraintWave {{{
ConstraintWave::~ConstraintWave(){
}
// }}}

// rhs {{{
void ConstraintWave::rhs(std::shared_ptr<FieldMap>& fieldMap){
  unsigned int nb = domain->getGhostPoints();
  // Check that the grid is actually big enough.
  const Grid& grid = fieldMap->getGrid();
  if(grid.getSize() < 1){
    printf("Grid is too small. Need at least 5 points.\n");
    return;
  }
  double **dudt = fieldMap->getSolverField("Evolution")->getCurrentRHS();
  double **u = fieldMap->getSolverField("Evolution")->getIntermediateData();

  // Define some variables we'll need.
  double stencil3[3] = {0.0, 0.0, 0.0};
  double dx = grid.getSpacing();
  int shp = grid.getSize();

  // Interior points
  for(int i = nb + 1; i < shp - nb - 1; i++){
    dudt[U_PHI][i] = u[U_PI][i];

    for(int j = 0; j < 3; j++){
      stencil3[j] = u[U_CHI][i - 1 + j];
    }
    dudt[U_PI][i] = operators::dx_2(stencil3, dx);

    for(int j = 0; j < 3; j++){
      stencil3[j] = u[U_PI][i - 1 + j];
    }
    dudt[U_CHI][i] = operators::dx_2(stencil3, dx);
  }

  // Boundary points.
  for(unsigned int m = 0; m < nEqs; m++){
    for(unsigned int i = 0; i < 3; i++){
      stencil3[i] = u[m][nb + i];
    }
    dudt[m][nb] = operators::dx_2off(stencil3, dx);
    for(unsigned int i = 0; i < 3; i++){
      stencil3[i] = u[m][shp - nb - 1 - i];
    }
    dudt[m][shp - nb - 1] = operators::dx_2off(stencil3, dx);
  }
}
// }}}

// initData {{{
void ConstraintWave::initData(){
  double x0 = 0.5*(domain->getBounds()[1] + domain->getBounds()[0]);
  double amp = 1.0;
  for(auto it = fieldData.begin(); it != fieldData.end(); ++it){
    auto evol = (**it)["Evolution"];
    const double *x = evol->getGrid().getPoints();
    unsigned int nx = evol->getGrid().getSize();
    double **u = evol->getData();
    for(unsigned int i = 0; i < nx; i++){
      double val = amp*std::exp(-(x[i] - x0)*(x[i] - x0)*64.0);
      u[U_PHI][i] = val;
      u[U_PI ][i] = 0.0;
      u[U_CHI][i] = -128.0*(x[i] - x0)*val;
    }
  }
}
// }}}

// setParameters {{{
void ConstraintWave::setParameters(FieldTestParameters *p){
  params = p;
}
// }}}

// getParameters {{{
FieldTestParameters* ConstraintWave::getParameters(){
  return params;
}
// }}}
