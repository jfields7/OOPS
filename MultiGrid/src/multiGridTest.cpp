#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <operators.h>
#include <cmath>
#include <cstdio>
#include <output.h>
#include <firstorderwave.h>
#include <waveparameters.h>
#include <cubicinterpolator.h>
#include <polynomialinterpolator.h>
#include <waveparser.h>

const double EPSILON = 1.e-15;
const double sigma1st = 0.0;
const double sigma2nd = 0.0;

const double PI = 3.14159265358979;

int main(int argc, char* argv[]){
  // Construct a domain and a set of grids.
  Domain domain = Domain();
  domain.setGhostPoints(3);
  domain.setCFL(0.25);
  int N0 = 65;
  int N = N0 - 1;
  int nb = domain.getGhostPoints();
  
  int ngrids = 1;
  while (N > 4*nb){
    N = N >> 1;
    if(N >= 4*nb){
      ngrids++;
    }
  }
  // Debugging only
  //ngrids = 1;

  double dx_grid = (domain.getBounds()[1] - domain.getBounds()[0])/ngrids;
  for(int i = 0; i < ngrids; i++){
    double bounds[2] = {dx_grid * i, dx_grid * (i + 1)};
    //double bounds[2] = {1.0 - dx_grid * (i + 1), 1.0 - dx_grid * i};
    domain.addGrid(bounds, ((N0 - 1) >> i) + 1);
    //domain.addGrid(bounds, N0);
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
  WaveParser parser = WaveParser();
  WaveParameters params;
  //CubicInterpolator interpolator = CubicInterpolator();
  PolynomialInterpolator interpolator = PolynomialInterpolator(4);
  FirstOrderWave ode = FirstOrderWave(domain, rk4);
  ode.setInterpolator(&interpolator);
  //WaveParameters *wp = (WaveParameters*) ode.getParameters();
  //wp->setInitialConditions(WaveParameters::GAUSSIAN);
  parser.updateParameters("wave.par",&params);
  ode.setParameters(&params);
  ode.initData();

  double ti = 0.0;
  double tf = 5.0;
  double dt = domain.getCFL()*(--domain.getGrids().end())->getSpacing();
  //double dt = domain.getCFL()*domain.getGrids().begin()->getSpacing();
  unsigned int M = (tf - ti)/dt;
  //ode.dump_csv("phi00000.csv", 0, 0);
  ode.output_frame("Phi", 0, 0);
  //ode.dump_csv("chi00000.csv", 0, 2);
  for(unsigned int i = 0; i < M; i++){
    double t = (i + 1)*dt;
    ode.evolveStep(dt);

    char buffer[12];
    sprintf(buffer, "phi%05d.csv",i + 1);
    //ode.dump_csv(buffer, t, 0);
    ode.output_frame("Phi", t, 0);
    //sprintf(buffer, "chi%05d.csv",i + 1);
    //ode.dump_csv(buffer, t, 2);
  }

  return 0;
}
