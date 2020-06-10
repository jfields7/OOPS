#include <rkck.h>
#include <ode.h>
#include <fieldmap.h>
#include <cmath>

Result RKCK::setStageTime(double srcTime, double &destTime, double dt, unsigned int stage){
  if(stage >= nStages){
    return INVALID_STAGE;
  }
  destTime = srcTime + tc[stage]*dt;
  return SUCCESS;
}

Result RKCK::calcStage(ODE *ode, std::shared_ptr<FieldMap>& fieldMap, double dt, unsigned int stage){
  unsigned int shp = fieldMap->getGrid().getSize();
  unsigned int vars;
  if(stage >= nStages){
    return INVALID_STAGE;
  }
  ode->rhs(fieldMap);
  for(auto field : fieldMap->getSolverFields()){
    auto solverData = field.second;
    double **data0 = solverData->getData();
    double **dataint = solverData->getIntermediateData();
    //double **dest = solverData->getCurrentRHS();
    double ***work = solverData->getWorkData();
    unsigned int vars = solverData->getEqCount();
    if(kc[stage][0] != 0.0){
      for(unsigned int m = 0; m < vars; m++){
        for(unsigned int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i];
          for(unsigned int k = 0; k <= stage; k++){
            dataint[m][i] += kc[stage][k]*work[k][m][i]*dt;
          }
        }
      }
    }
  }

  return SUCCESS;
}

Result RKCK::combineStages(std::shared_ptr<FieldMap>& fieldMap, double dt){
  unsigned int shp = fieldMap->getGrid().getSize();
  double err = 0.0;
  double temp = 0.0;
  for(auto field : fieldMap->getSolverFields()){
    auto solverData = field.second;
    double **data0 = solverData->getData();
    double **dataint = solverData->getIntermediateData();
    double ***work = solverData->getWorkData();
    double **k1 = work[0];
    double **k2 = work[1];
    double **k3 = work[2];
    double **k4 = work[3];
    double **k5 = work[4];
    double **k6 = work[5];
    unsigned int vars = solverData->getEqCount();
    for(unsigned int m = 0; m < vars; m++){
      for(unsigned int i = 0; i < shp; i++){
        dataint[m][i] = data0[m][i] + ((c5[0]*k1[m][i] + c5[2]*k3[m][i]) + (c5[3]*k4[m][i] + c5[5]*k6[m][i]))*dt;
        temp = (dataint[m][i] - data0[m][i]) - ((c4[0]*k1[m][i] + c4[2]*k3[m][i]) + (c4[3]*k4[m][i]) +
                (c4[4]*k5[m][i] + c4[5]*k6[m][i]))*dt;
        err = fmax(temp, err);
      }
    }
  }
  dtrec = 0.9*dt*fmin(fmax(sqrt(tol/(2.0*fabs(err))),0.3),2.0);

  return SUCCESS;
}
