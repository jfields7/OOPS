#include <rk4.h>
#include <ode.h>
#include <fieldmap.h>

Result RK4::setStageTime(double srcTime, double &destTime, double dt, unsigned int stage){
  if(stage >= nStages){
    return INVALID_STAGE;
  }
  destTime = srcTime + tc[stage]*dt;
  return SUCCESS;
}

Result RK4::calcStage(ODE *ode, std::shared_ptr<FieldMap>& fieldMap, double dt, unsigned int stage){
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
    double **dest = solverData->getCurrentRHS();
    unsigned int vars = solverData->getEqCount();
    if(kc[stage] != 0.0){
      for(unsigned int m = 0; m < vars; m++){
        for(unsigned int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i] + kc[stage]*dest[m][i]*dt;
        }
      }
    }
  }

  return SUCCESS;
}

Result RK4::combineStages(std::shared_ptr<FieldMap>& fieldMap, double dt){
  unsigned int shp = fieldMap->getGrid().getSize();
  for(auto field : fieldMap->getSolverFields()){
    auto solverData = field.second;
    double **data0 = solverData->getData();
    double **dataint = solverData->getIntermediateData();
    double ***work = solverData->getWorkData();
    double **k1 = work[0];
    double **k2 = work[1];
    double **k3 = work[2];
    double **k4 = work[3];
    unsigned int vars = solverData->getEqCount();
    for(unsigned int m = 0; m < vars; m++){
      for(unsigned int i = 0; i < shp; i++){
        dataint[m][i] = data0[m][i] + ((k1[m][i] + 2.0*k2[m][i]) + (2.0*k3[m][i] + k4[m][i]))*dt/6.0;
      }
    }
  }

  return SUCCESS;
}
