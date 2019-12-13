#include <ode.h>

ODE::ODE(const unsigned int n, const unsigned int id) : nEqs(n), pId(id){
  domain = nullptr;
  params = nullptr;
  solver = nullptr;
}

ODE::~ODE(){
  data.clear();
}

Result ODE::reallocateData(){
  data.clear();

  std::set<Grid> grids = domain->getGrids();

  for(auto it = grids.begin(); it != grids.end(); ++it){
    data.emplace(nEqs, solver->getNStages(), *it);
  }

  return SUCCESS;
}

Result ODE::setDomain(Domain *d){
  domain = d;

  return reallocateData();
}

Result ODE::setSolver(Solver *s){
  solver = s;

  return reallocateData();
}

Result ODE::evolveStep(double dt){
  // Loop over every stage for the solver.
  for(unsigned int i = 0; i < solver->getNStages(); i++){
    // Loop over every data set in the domain.
    for(auto it = data.begin(); it != data.end(); ++it){
      solver->calcStage(rhs, it->getData(), it->getIntermediateData(), (it->getWorkData())[i],
                       it->getGrid(), dt, nEqs, i);
    }

    // Perform the grid exchange and apply boundary conditions.
    performGridExchange();
    applyBoundaries();
  }

  for(auto it = data.begin(); it != data.end(); ++it){
    solver->combineStages(it->getWorkData(), it->getData(), it->getGrid(), dt, nEqs);
  }
  performGridExchange();
  applyBoundaries();

  return SUCCESS;
}

void ODE::performGridExchange(){
  // TODO
}

Result ODE::setParameters(Parameters *p){
  if(p->getId() == pId){
    params = p;
    return SUCCESS;
  }
  else{
    return UNRECOGNIZED_PARAMS;
  }
}
