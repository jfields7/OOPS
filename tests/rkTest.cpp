#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <cmath>
#include <cstdio>

const double EPSILON = 1e-15;

/**
 * The righthand side for a simple harmonic oscillator.
 */
void rhsSHO(Grid& grid, double **u, double **dudt){
  for(int i = 0; i < grid.getSize(); i++){
    dudt[i][0] = u[i][1];
    dudt[i][1] = -u[i][0];
  }
}

/**
 * The righthand side for a very stiff ODE.
 */
void rhsStiff(Grid& grid, double **u, double **dudt){
  for(int i = 0; i < grid.getSize(); i++){
    dudt[i][0] = 2.0*u[i][0] - std::exp(u[i][1]);
    dudt[i][1] = 1.0;
  }
}

/**
 * Calculate the relative error.
 */
double calcError(double exp, double act){
  return fabs((act - exp)/(exp + EPSILON));
}

/**
 * Setup the initial data.
 */
void initializeData(double **u, double **uint, double ***k){
  /*u[0][0] = 1.0;
  u[0][1] = 0.0;
  u[1][0] = 0.0;
  u[1][1] = 1.0;*/
  u[0][0] = 1.0/3.0;
  u[0][1] = 0.0;
  u[1][0] = 1.0/3.0;
  u[1][1] = 0.0;

  uint[0][0] = 0.0;
  uint[0][1] = 0.0;
  uint[1][0] = 0.0;
  uint[1][1] = 0.0;

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 2; j++){
      k[i][j][0] = 0.0;
      k[i][j][1] = 0.0;
    }
  }
}

int main(int argc, char *argv[]){
  // We don't need to construct a domain; this is just a time integrator.
  RK4 rk4 = RK4();
  double bounds[2] = {0.0, 1.0};
  Grid grid = Grid(bounds, 2, 0);


  // Set up information about the simulation.
  //double dt = 0.1;
  double dt = 1.e-3;
  double t = 0.0;
  double tf = 5.0;

  // Set up initial conditions and initialize all the data we'll need.
  double **u = new double*[2];
  u[0] = new double[2];
  u[1] = new double[2];
  u[0][0] = 1.0/3.0;
  u[0][1] = 0.0;
  u[1][0] = 1.0/3.0;
  u[1][1] = 0.0;
  double **uint = new double*[2];
  uint[0] = new double[2];
  uint[1] = new double[2];
  uint[0][0] = 0.0;
  uint[0][1] = 0.0;
  uint[1][0] = 0.0;
  uint[1][1] = 0.0;
  double ***k = new double**[4];
  for(int i = 0; i < 4; i++){
    k[i]  = new double*[2];
    for(int j = 0; j < 2; j++){
      k[i][j] = new double[2];
      k[i][j][0] = 0.0;
      k[i][j][1] = 0.0;
    }
  }

  unsigned int N;
  // Perform the integration.
  N = tf / dt;
  for(unsigned int n = 0; n < N; n++){
    t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      //rk4.calcStage(&rhsSHO,u,uint,k[i], grid, dt, 2, i);
      rk4.calcStage(&rhsStiff,u,uint,k[i], grid, dt, 2, i);
    }
    rk4.combineStages(k, u, grid, dt, 2);
  }

  //double exp = cos(tf);
  double exp = std::exp(-tf)/3.0;
  printf("Initial amplitude test.\n");
  printf("  Expected: %g\n",exp);
  printf("  Actual: %g\n",u[0][0]);
  printf("  Error: %g\n",calcError(exp,u[0][0]));
  double c1 = u[0][0];

  /*exp = sin(tf);
  printf("Initial velocity test.\n");
  printf("  Expected: %g\n",exp);
  printf("  Actual: %g\n",u[1][0]);
  printf("  Error: %g\n",calcError(exp,u[1][0]));*/
  exp = tf;
  printf("Time test.\n");
  printf("  Expected: %g\n",tf);
  printf("  Actual: %g\n",u[0][1]);
  printf("  Error: %g\n",calcError(exp,u[0][1]));

  // Do a convergence test. Reset the initial data.
  // The dt = 0.05 test.
  initializeData(u,uint,k);
  dt = 5.e-4;
  N = tf / dt;
  for(unsigned int n = 0; n < N; n++){
    t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      //rk4.calcStage(&rhsSHO, u, uint, k[i], grid, dt, 2, i);
      rk4.calcStage(&rhsStiff, u, uint, k[i], grid, dt, 2, i);
    }
    rk4.combineStages(k, u, grid, dt, 2);
  }
  double c2 = u[0][0];

  // The dt = 0.025 test.
  initializeData(u, uint, k);
  dt = 2.5e-4;
  N = tf / dt;
  for(unsigned int n = 0; n < N; n++){
    t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      //rk4.calcStage(&rhsSHO, u, uint, k[i], grid, dt, 2, i);
      rk4.calcStage(&rhsStiff, u, uint, k[i], grid, dt, 2, i);
    }
    rk4.combineStages(k, u, grid, dt, 2);
  }
  double c3 = u[0][0];

  // The dt = 0.0125 test.
  initializeData(u, uint, k);
  dt = 1.e-4;
  N = tf / dt;
  for(unsigned int n = 0; n < N; n++){
    t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      //rk4.calcStage(&rhsSHO, u, uint, k[i], grid, dt, 2, i);
      rk4.calcStage(&rhsStiff, u, uint, k[i], grid, dt, 2, i);
    }
    rk4.combineStages(k, u, grid, dt, 2);
  }
  double c4 = u[0][0];

  initializeData(u, uint, k);
  dt = 1.e-5;
  N = tf / dt;
  for(unsigned int n = 0; n < N; n++){
    t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      //rk4.calcStage(&rhsSHO, u, uint, k[i], grid, dt, 2, i);
      rk4.calcStage(&rhsStiff, u, uint, k[i], grid, dt, 2, i);
    }
    rk4.combineStages(k, u, grid, dt, 2);
  }
  double c5 = u[0][0];

  initializeData(u, uint, k);
  dt = 1.e-6;
  N = tf / dt;
  for(unsigned int n = 0; n < N; n++){
    t = n * dt;
    for(int i = 0; i < rk4.getNStages(); i++){
      //rk4.calcStage(&rhsSHO, u, uint, k[i], grid, dt, 2, i);
      rk4.calcStage(&rhsStiff, u, uint, k[i], grid, dt, 2, i);
    }
    rk4.combineStages(k, u, grid, dt, 2);
  }
  double c6 = u[0][0];

  // Print the results of the convergence test.
  /*exp = cos(tf);
  double e1, e2, e3, e4, e5, e6;
  e1 = calcError(exp,c1);
  e2 = calcError(exp,c2);
  e3 = calcError(exp,c3);
  e4 = calcError(exp,c4);
  e5 = calcError(exp,c5);
  e6 = calcError(exp,c6);
  double p12 = log(e2/e1)/log(0.1);
  double p23 = log(e3/e2)/log(0.1);
  double p34 = log(e4/e3)/log(0.1);
  double p45 = log(e5/e4)/log(0.1);
  double p56 = log(e6/e5)/log(0.1);
  double avg = (p12 + p23 + p34 + p45 + p56)/5.0;
  printf("Convergence test results:\n");
  printf("  Expected Amplitude: %g\n",cos(tf));
  printf("  Actual Amplitudes:\n");
  printf("    dt = 1e-1 -> u = %g\n",c1);
  printf("    dt = 1e-2 -> u = %g\n",c2);
  printf("    dt = 1e-3 -> u = %g\n",c3);
  printf("    dt = 1e-4 -> u = %g\n",c4);
  printf("    dt = 1e-5 -> u = %g\n",c5);
  printf("    dt = 1e-6 -> u = %g\n",c6);
  printf("  1e-1 -> 1e-2: p = %g\n",p12);
  printf("  1e-2 -> 1e-3: p = %g\n",p23);
  printf("  1e-3 -> 1e-4: p = %g\n",p34);
  printf("  1e-4 -> 1e-5: p = %g\n",p45);
  printf("  1e-5 -> 1e-6: p = %g\n",p56);
  printf("  Average: p = %g\n",avg);*/
  //double p123 = log(fabs(c1 - c2)/fabs(c2 - c3))/log(2.0);
  double p123 = fabs(c1 - c2)/fabs(c2 - c3);
  printf("Estimated order of convergence: %g\n",p123);

  // Release all the memory we've used.
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 2; j++){
      delete[] k[i][j];
    }
    delete[] k[i];
  }
  delete[] k;
  delete[] uint[1];
  delete[] uint[0];
  delete[] uint;
  delete[] u[1];
  delete[] u[0];
  delete[] u;
}
