#include <ode.h>
#include <cubic.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <output.h>
#include <memory>

ODE::ODE(const unsigned int n, const unsigned int id) : nEqs(n), pId(id){
  domain = nullptr;
  solver = nullptr;
  max_dx = 0.0;
  interpolator = nullptr;
  time = 0.0;
  for(unsigned int i = 0; i < nEqs; i++){
    evolutionIndices.push_back(i);
  }
}

ODE::~ODE(){
  data.clear();
}

// reallocateData {{{
Result ODE::reallocateData(){
  data.clear();

  //std::set<Grid> grids = domain->getGrids();

  /*for(auto it = domain->getGrids().begin(); it != domain->getGrids().end(); ++it){
    max_dx = (it->getSpacing() > max_dx) ? it->getSpacing() : max_dx;
    data.emplace(nEqs, solver->getNStages(), *it);
  }*/

  for(auto it = domain->getGrids().begin(); it != domain->getGrids().end(); ++it){
    max_dx = (it->getSpacing() > max_dx) ? it->getSpacing() : max_dx;
    //fieldData.emplace(std::unique_ptr<FieldMap>(new FieldMap(*it, fieldList)));
    fieldData.emplace(std::make_shared<FieldMap>(*it, fieldList));
  }

  return SUCCESS;
}
// }}}

// addField {{{
Result ODE::addField(std::string name, unsigned int eqs, bool isEvolved, bool isComm) {
  if(fieldList.find(name) != fieldList.end()){
    return FIELD_EXISTS;
  }
  else{
    unsigned int stages = (isEvolved ? solver->getNStages() : 0);
    fieldList.insert({name, FieldInfo(name, eqs, stages, isComm)});
  }
  return SUCCESS;
}
// }}}

// removeField {{{
Result ODE::removeField(std::string name){
  if(fieldList.find(name) == fieldList.end()){
    return UNRECOGNIZED_FIELD;
  }
  else{
    fieldList.erase(name);
  }
  return SUCCESS;
}
// }}}

// setDomain {{{
Result ODE::setDomain(Domain *d){
  domain = d;

  return reallocateData();
}
// }}}

// setSolver {{{
Result ODE::setSolver(Solver *s){
  solver = s;

  return reallocateData();
}
// }}}

// setInterpolator {{{

Result ODE::setInterpolator(Interpolator* interp){
  interpolator = interp;
  return SUCCESS;
}

// }}}

// evolveStep {{{
Result ODE::evolveStep(double dt){
  // Keep track of the current time as well as the original time.
  double old_time = time;

  // Loop over every stage for the solver.
  solver->setStageTime(old_time, time, dt, 0);
  for(auto fields : fieldData){
    // Copy the current solution into the intermediate variables.
    for(auto fieldpair : fields->getSolverFields()){
      auto field = fieldpair.second;
      field->setCurrentStage(0);
      double **data0 = field->getData();
      double **dataint = field->getIntermediateData();
      unsigned int vars = field->getEqCount();
      unsigned int shp = field->getGrid().getSize();
      for(unsigned int m = 0; m < vars; m++){
        for(unsigned int i = 0; i < shp; i++){
          dataint[m][i] = data0[m][i];
        }
      }
    }
    // Calculate the first stage.
    solver->calcStage(this, fields, dt, 0);
  }
  // Perform the grid exchange and apply boundary conditions.
  doAfterStage();
  performGridExchange();
  doAfterExchange();
  applyBoundaries();
  doAfterBoundaries();

  for(unsigned int i = 1; i < solver->getNStages(); i++){
    solver->setStageTime(old_time, time, dt, i);
    // Loop over every data set in the domain.
    for(auto fields : fieldData){
      for(auto fieldpair : fields->getSolverFields()){
        fieldpair.second->setCurrentStage(i);
      }
      solver->calcStage(this, fields, dt, i);
    }

    // Perform the grid exchange and apply boundary conditions.
    doAfterStage();
    performGridExchange();
    doAfterExchange();
    applyBoundaries();
    doAfterBoundaries();
  }

  time = old_time + dt;

  for(auto fields : fieldData){
    solver->combineStages(fields, dt);
  }
  doAfterStage();
  performGridExchange();
  doAfterExchange();
  applyBoundaries();
  doAfterBoundaries();

  // Swap all the intermediate data into the solution array.
  for(auto fields : fieldData){
    // Copy the intermediate solution into the real solution.
    for(auto fieldpair : fields->getSolverFields()){
      auto field = fieldpair.second;
      double **data0 = field->getData();
      double **dataint = field->getIntermediateData();
      unsigned int vars = field->getEqCount();
      unsigned int shp = field->getGrid().getSize();
      for(unsigned int m = 0; m < vars; m++){
        for(unsigned int i = 0; i < shp; i++){
          data0[m][i] = dataint[m][i];
        }
      }
    }
  }

  return SUCCESS;
}
// }}}

// performGridExchange {{{
void ODE::performGridExchange(){
  // If the domain doesn't have multiple grids, no transfer is necessary.
  // Similarly, if no interpolator is available, there's no reason to interpolate.
  if(fieldData.size() < 2 || interpolator == nullptr){
    return;
  }
  // Otherwise, we need to exchange the ghost points between all the grids.
  /*for(auto it = data.begin(); it != (--data.end()); ++it){
    auto it2 = it;
    ++it2;
    exchangeGhostPoints(*it, *it2);
  }*/
  for(auto it = fieldData.begin(); it != (--fieldData.end()); ++it){
    auto it2 = it;
    ++it2;
    exchangeGhostPoints(*it, *it2);
  }
}
// }}}

// exchangeGhostPoints {{{
//void ODE::exchangeGhostPoints(const SolverData& data1, const SolverData& data2){
void ODE::exchangeGhostPoints(const std::shared_ptr<FieldMap>& data1, const std::shared_ptr<FieldMap>& data2){
  double **u1;
  double **u2;
  double dx1 = data1->getGrid().getSpacing();
  double dx2 = data2->getGrid().getSpacing();
  unsigned int nb = domain->getGhostPoints();
  unsigned int shp1 = data1->getGrid().getSize();
  for(auto field : fieldList){
    // Save some information that we'll be referencing a lot to save typing.
    FieldInfo &info = field.second;
    if(!info.isComm){
      continue;
    }
    if(info.nStages == 0){
      u1 = (*data1)[field.first]->getData();
      u2 = (*data2)[field.first]->getData();
    }
    else{
      u1 = data1->getSolverField(field.first)->getIntermediateData();
      u2 = data2->getSolverField(field.first)->getIntermediateData();
    }
    // Check if the spacing is essentially the same. If so, we just need to swap the points.
    if(fabs((dx1 - dx2)/dx1) < 1e-14){
      for(unsigned int m = 0; m < nEqs; m++){
        for(unsigned int i = 0; i < nb; i++){
          // Copy the first physical points on u2 into the ghost points of u1.
          u1[m][shp1 - nb + i] = u2[m][nb + i];
          // Copy the last physical points of u1 into the ghost points of u2.
          u2[m][i] = u1[m][shp1 - 2*nb + i];
        }
      }
    }
    // Otherwise, let's perform some interpolation. We first need to enforce a 2x difference
    // between the grids and complain if that's not the case.
    else if(fabs((dx1*2.0 - dx2)/dx2) < 1e-14){
      // The left grid is finer than the right grid.
      // Note the nb + 1. We also need to copy the one physical point so that they stay
      // consistent. We use the one from the more refined grid.
      for(unsigned int m = 0; m < nEqs; m++){
        for(unsigned int i = 0; i < nb + 1; i++){
          // Copy the exact point from u1 into the ghost region of u2.
          u2[m][i] = u1[m][shp1 - 1 - 3*nb + 2*i];
        }
      }
      interpolateLeft(data1, data2);
    }
    else if(fabs((dx1 - 2.0*dx2)/dx2) < 1e-14){
      // The right grid is finer than the left grid.
      //for(unsigned int i = 0; i < nb; i++){
      // Note the nb + 1. We also need to copy the one physical point so that they stay
      // consistent. We use the one from the more refined grid.
      for(unsigned int m = 0; m < nEqs; m++){
        for(unsigned int i = 0; i < nb + 1; i++){
          // Copy the exact point from u2 into the ghost region of u1.
          u1[m][shp1 - nb - 1 + i] = u2[m][nb + 2*i];
        }
      }
      interpolateRight(data1, data2);
    }
    else{
      std::cerr << "Grids are not spaced correctly; spacing between neighboring grids" << "\n"
                << "must be equal or a factor of two different.\n";
    }
  }
}
// }}}

// interpolateLeft {{{
void ODE::interpolateLeft(const std::shared_ptr<FieldMap>& datal, const std::shared_ptr<FieldMap>& datar){
  double **ul;
  double **ur;
  unsigned int nb = domain->getGhostPoints();
  unsigned int shpl = datal->getGrid().getSize();
  double *stencil = interpolator->getStencil();
  unsigned int nStart = interpolator->getStencilSize()/2;
  for(auto field : fieldList){
    FieldInfo &info = field.second;
    if(!info.isComm){
      continue;
    }
    if(info.nStages == 0){
      ul = (*datal)[field.first]->getData();
      ur = (*datar)[field.first]->getData();
    }
    else{
      ul = datal->getSolverField(field.first)->getIntermediateData();
      ur = datar->getSolverField(field.first)->getIntermediateData();
    }
    for(unsigned int n = 0; n < info.nEqs; n++){
      for(unsigned int i = 0; i < nb; i++){
        // If we're on an even ghost point, we don't need to interpolate. This corresponds to an
        // odd index because counting starts at 0, hence the screwy math.
        if( i & 1 == 1){
          ul[n][shpl - nb + i] = ur[n][nb + (i + 1)/2];
        }
        else{
          // Fill in the stencil. Use the left points from the finer grid and the right
          // points from the coarser grid. This makes sure that we don't run out of points.
          for(int m = -nStart; m < 0; m++){
            stencil[m + nStart] = ul[n][shpl - nb + i + 1 + 2*m];
          }
          for(int m = 0; m < nStart; m++){
            stencil[m + nStart] = ur[n][nb + i/2 + 1 + m];
          }
          ul[n][shpl - nb + i] = interpolator->interpolate();
        }
      }
    }
  }
}
// }}}

// interpolateRight {{{
void ODE::interpolateRight(const std::shared_ptr<FieldMap>& datal, const std::shared_ptr<FieldMap>& datar){
  // FIXME: Right now, we assume cubic interpolation. This needs to be changed to be general.
  double **ul;
  double **ur;
  unsigned int nb = domain->getGhostPoints();
  unsigned int shpl = datal->getGrid().getSize();
  double *stencil = interpolator->getStencil();
  unsigned int nStart = interpolator->getStencilSize()/2;
  for(auto field : fieldList){
    const FieldInfo& info = field.second;
    if(!info.isComm){
      continue;
    }
    if(info.nStages == 0){
      ul = (*datal)[field.first]->getData();
      ur = (*datar)[field.first]->getData();
    }
    else{
      ul = datal->getSolverField(field.first)->getIntermediateData();
      ur = datar->getSolverField(field.first)->getIntermediateData();
    }
    for(unsigned int n = 0; n < info.nEqs; n++){
      for(unsigned int i = 0; i < nb; i++){
        // If we're on an even ghost point, we don't need to interpolate. This corresponds to an
        // odd index because counting starts at 0, hence the screwy math.
        if(i & 1 == 1){
          ur[n][nb - 1 - i] = ul[n][shpl - nb - 1 - (i + 1)/2];
        }
        else{
          // Fill in the stencil. Use the right points from the finer grid and the left
          // points from the coarser grid. This makes sure that we don't run out of points.
          for(int m = -nStart; m < 0; m++){
            stencil[m + nStart] = ul[n][shpl - nb - 1 + m - i/2];
          }
          for(int m = 0; m < nStart; m++){
            stencil[m + nStart] = ur[n][nb + 2*m - i];
          }
          ur[n][nb - 1 - i] = interpolator->interpolate();
        }
      }
    }
  }
}
// }}}

// output_field {{{
void ODE::outputSDFField(std::string field, char* name, double t, unsigned int var){
  unsigned int nb = domain->getGhostPoints();
  unsigned int size = 0;
  double *x = nullptr;
  for(auto it = fieldData.begin(); it != fieldData.end(); ++it){
    auto evol = (**it)[field];
    unsigned int shp = evol->getGrid().getSize();
    if(shp > size){
      if(x != nullptr){
        delete[] x;
      }
      size = shp*2;
      x = new double[size];
    }
    const double *points = evol->getGrid().getPoints();
    for(unsigned int i = 0; i < shp; i++){
      x[i] = points[i];
    }

    double **u = evol->getData();
    output::output_data(name, u[var] + nb, x + nb, shp - 2*nb, t);
  }
  delete[] x;
}
// }}}

// dump_frame {{{
void ODE::dumpCSV(std::string field, char* name, double t, unsigned int var){
  static bool first_call = 1;
  FILE *f;
  if(first_call){
    f = fopen(name, "w");
    //first_call = 0;
  }
  else{
    f = fopen(name, "a");
  }

  int nb = domain->getGhostPoints();
  for(auto it = fieldData.begin(); it != fieldData.end(); ++it){
    auto evol = (**it)[field];
    const Grid& grid = (*it)->getGrid();
    unsigned int shp = grid.getSize();
    for(int i = nb; i < shp - nb; i++){
      //fprintf(f,"%08g,   %08g,   %08g\n", t, r[i], data[i][var]);
      fprintf(f,"%08g,   %08g,   %0g\n", t,
              grid.getPoints()[i],
              evol->getData()[var][i]);
    }
  }

  fclose(f);
}
// }}}

// getTime {{{
double ODE::getTime(){
  return time;
}
// }}}

// setTime {{{
void ODE::setTime(double t){
  time = t;
}
