#include <rk4.h>
#include <ode.h>

Result RK4::calcStage(void (*rhs)(const Grid&, double**,double**), double *data0[], double *dataint[], 
                      double *dest[], const Grid& grid, double dt, const unsigned int vars, 
                      unsigned int stage){
  int shp = grid.getSize();
  switch(stage){
    case 0:
      rhs(grid, data0, dest);
      for(int m = 0; m < vars; m++){
        for(int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + 0.5*dest[m][i]*dt;
        }
      }
      return SUCCESS;
      break;
    case 1:
      rhs(grid, dataint, dest);
      for(int m = 0; m < vars; m++){
        for(int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + 0.5*dest[m][i]*dt;
        }
      }
      return SUCCESS;
      break;
    case 2:
      rhs(grid, dataint, dest);
      for(int m = 0; m < vars; m++){
        for(int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + dest[m][i]*dt;
        }
      }
      return SUCCESS;
      break;
    case 3:
      rhs(grid, dataint, dest);
      return SUCCESS;
      break;
    default:
      return INVALID_STAGE;
      break;
  }
}

Result RK4::calcStage(ODE *ode, double *data0[], double *dataint[], double *dest[],
                      const Grid& grid, double dt, unsigned int stage){
  unsigned int shp = grid.getSize();
  unsigned int vars = ode->getNEqs();
  switch(stage){
    case 0:
      ode->rhs(grid, data0, dest);
      for(int m = 0; m < vars; m++){
        for(int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + 0.5*dest[m][i]*dt;
        }
      }
      return SUCCESS;
      break;
    case 1:
      ode->rhs(grid, dataint, dest);
      for(int m = 0; m < vars; m++){
        for(int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + 0.5*dest[m][i]*dt;
        }
      }
      return SUCCESS;
      break;
    case 2:
      ode->rhs(grid, dataint, dest);
      for(int m = 0; m < vars; m++){
        for(int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + dest[m][i]*dt;
        }
      }
      return SUCCESS;
      break;
    case 3:
      ode->rhs(grid, dataint, dest);
      return SUCCESS;
      break;
    default:
      return INVALID_STAGE;
      break;
  }
}

Result RK4::combineStages(double **data[], double *dest[], const Grid& grid, double dt, int vars){
  int shp = grid.getSize();
  // Combine all the stages according to the RK4 method.
  double **k1 = data[0];
  double **k2 = data[1];
  double **k3 = data[2];
  double **k4 = data[3];
  double error;
  double term;
  double old;
  for(int m = 0; m < vars; m++){
    for(int i = 0; i < shp; i++){
      dest[m][i] = dest[m][i] + ((k1[m][i] + 2.0*k2[m][i]) + (2.0*k3[m][i] + k4[m][i]))*dt/6.0;
    }
  }
}
