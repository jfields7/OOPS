#include <constraintwave.h>
#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <fieldtestparser.h>
#include <fieldtestparameters.h>
#include <polynomialinterpolator.h>

int main(int argc, char* argv[]){
  Domain domain = Domain();
  domain.setGhostPoints(0);
  domain.setCFL(0.25);
  int N0 = 101;
  double bounds[2] = {0.0, 1.0};
  domain.addGrid(bounds, N0);

  RK4 rk4 = RK4();
  //FieldTestParser parser = FieldTestParser();
  FieldTestParameters params;
  PolynomialInterpolator interpolator = PolynomialInterpolator(4);
  ConstraintWave ode = ConstraintWave(domain, rk4);
  ode.setInterpolator(&interpolator);
  ode.setParameters(&params);
  ode.initData();

  double ti = 0.0;
  double tf = 5.0;
  double dt = domain.getCFL()*(domain.getGrids().begin())->getSpacing();
  unsigned int M = (tf - ti)/dt;
  //ode.output_frame("Phi", 0, 0);
  //ode.outputSDFField("Evolution","Phi",0,0);
  ode.dumpCSV(std::string("Evolution"), std::string("Phi00000.csv"), 0, 0);
  for(unsigned int i = 0; i < M; i++){
    double t = (i + 1)*dt;
    ode.evolveStep(dt);

    //ode.output_frame("Phi", t, 0);
    //ode.outputSDFField("Evolution","Phi",t,0);
    char name[32];
    sprintf(name,"Phi%05d.csv",i+1);
    ode.dumpCSV(std::string("Evolution"), std::string(name), t, 0);
  }

  return 0;
}
