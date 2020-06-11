#include <domain.h>
#include <grid.h>
#include <rkck.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <wave.h>
#include <polynomialinterpolator.h>
#include <waveparser.h>
#include <waveparameters.h>

int main(int argc, char* argv[]){
  if(argc < 2){
    std::cout << "Usage: ./Wave <parameter file>\n";
    return 0;
  }

  WaveParameters params;
  WaveParser parser;
  parser.updateParameters(argv[1], &params);

  // Construct our domain and a grid to fit on it.
  Domain domain = Domain();
  int N = params.getGridPoints();
  double bounds[2] = {0.0};
  bounds[0] = domain.getBounds()[0];
  bounds[1] = domain.getBounds()[1];
  domain.addGrid(bounds, N);

  // Set up our ODE system.
  RKCK rkck = RKCK();
  PolynomialInterpolator interpolator = PolynomialInterpolator(4);
  Wave ode = Wave(domain, rkck);
  ode.setInterpolator(&interpolator);
  ode.setParameters(&params);
  ode.initData();

  double ti = 0.0;
  double tf = 5.0;
  double dx = domain.getGrids().begin()->getSpacing();
  double dt = domain.getCFL()*dx;
  double dtmin = params.getMinCFL()*dx;
  double dtmax = params.getMaxCFL()*dx;
  rkck.setErrorTolerance(params.getErrorTolerance());
  //ode.dumpCSV("Evolution","phi00000.csv", 0, 0);
  ode.outputSDFField("Evolution","Phi",0, 0);
  for(double t = ti; t < tf; t+=dt){ 
    if(t != ti){
      dt = fmin(dtmax, fmax(dtmin, rkck.getRecommendedStepSize()));
    }
    ode.evolveStep(dt);

    //char buffer[12];
    //sprintf(buffer, "phi%05d.csv",i+1);
    ode.outputSDFField("Evolution","Phi",t, 0);
    //printf("Current dt = %g\n", dt);
  }

  return 0;
}
