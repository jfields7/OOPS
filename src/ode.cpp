#include <ode.h>
#include <cubic.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <output.h>

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
  auxiliaryFields.clear();
  auxiliaryFieldNames.clear();
  data.clear();
}

// reallocateData {{{
Result ODE::reallocateData(){
  data.clear();

  //std::set<Grid> grids = domain->getGrids();

  for(auto it = domain->getGrids().begin(); it != domain->getGrids().end(); ++it){
    max_dx = (it->getSpacing() > max_dx) ? it->getSpacing() : max_dx;
    data.emplace(nEqs, solver->getNStages(), *it);
  }
  for(const auto &field : auxiliaryFieldNames){
    auxiliaryFields[field].clear();
    for(auto it = domain->getGrids().begin(); it != domain->getGrids().end(); ++it){
      auxiliaryFields[field].emplace(nEqs, *it);
    }
  }

  return SUCCESS;
}
// }}}

// addAuxiliaryField {{{
Result ODE::addAuxiliaryField(std::string name) {
  if(auxiliaryFields.find(name) != auxiliaryFields.end()){
    return FIELD_EXISTS;
  }
  else{
    auxiliaryFields[name] = std::set<ODEData>();
    auxiliaryFieldNames.insert(name);
    return SUCCESS;
  }
}
// }}}

// removeAuxiliaryField {{{
Result ODE::removeAuxiliaryField(std::string name){
  if(auxiliaryFields.find(name) == auxiliaryFields.end()){
    return UNRECOGNIZED_FIELD;
  }
  else{
    auxiliaryFields.erase(name);
    auxiliaryFieldNames.erase(name);
  }
}
// }}}

// getAuxiliaryField {{{
std::set<ODEData>* ODE::getAuxiliaryField(std::string name){
  if(auxiliaryFields.find(name) != auxiliaryFields.end()){
    return &auxiliaryFields[name];
  }
  else{
    return nullptr;
  }
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
  for(unsigned int i = 0; i < solver->getNStages(); i++){
    solver->setStageTime(old_time, time, dt, i);
    // Loop over every data set in the domain.
    for(auto it = data.begin(); it != data.end(); ++it){
      //solver->calcStage(rhs, it->getData(), it->getIntermediateData(), (it->getWorkData())[i],
      //                 it->getGrid(), dt, nEqs, i);
      solver->calcStage(this, it->getData(), it->getIntermediateData(), (it->getWorkData())[i],
                        it->getGrid(), dt, i);
    }

    // Perform the grid exchange and apply boundary conditions.
    doAfterStage(true);
    performGridExchange();
    doAfterExchange(true);
    applyBoundaries(true);
    doAfterBoundaries(true);
  }

  time = old_time + dt;

  for(auto it = data.begin(); it != data.end(); ++it){
    solver->combineStages(it->getWorkData(), it->getData(), it->getGrid(), dt, evolutionIndices);
  }
  doAfterStage(false);
  performGridExchange();
  doAfterExchange(false);
  applyBoundaries(false);
  doAfterBoundaries(false);

  return SUCCESS;
}
// }}}

// performGridExchange {{{
void ODE::performGridExchange(){
  // If the domain doesn't have multiple grids, no transfer is necessary.
  if(data.size() < 2){
    return;
  }
  // Otherwise, we need to exchange the ghost points between all the grids.
  for(auto it = data.begin(); it != (--data.end()); ++it){
    auto it2 = it;
    ++it2;
    exchangeGhostPoints(*it, *it2);
  }
}
// }}}

// exchangeGhostPoints {{{
void ODE::exchangeGhostPoints(const SolverData& data1, const SolverData& data2){
  // Save some information that we'll be referencing a lot to save typing.
  double **u1 = data1.getData();
  double **u2 = data2.getData();
  double dx1 = data1.getGrid().getSpacing();
  double dx2 = data2.getGrid().getSpacing();
  unsigned int nb = domain->getGhostPoints();
  unsigned int shp1 = data1.getGrid().getSize();
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
// }}}

// interpolateLeft {{{
void ODE::interpolateLeft(const SolverData& datal, const SolverData& datar){
  // FIXME: Right now, we assume cubic interpolation. This needs to be changed to be general.
  double **ul = datal.getData();
  double **ur = datar.getData();
  unsigned int nb = domain->getGhostPoints();
  unsigned int shpl = datal.getGrid().getSize();
  double *stencil = interpolator->getStencil();
  unsigned int nStart = interpolator->getStencilSize()/2;
  for(unsigned int n = 0; n < nEqs; n++){
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
// }}}

// interpolateRight {{{
void ODE::interpolateRight(const SolverData& datal, const SolverData& datar){
  // FIXME: Right now, we assume cubic interpolation. This needs to be changed to be general.
  double **ul = datal.getData();
  double **ur = datar.getData();
  unsigned int nb = domain->getGhostPoints();
  unsigned int shpl = datal.getGrid().getSize();
  double *stencil = interpolator->getStencil();
  unsigned int nStart = interpolator->getStencilSize()/2;
  for(unsigned int n = 0; n < nEqs; n++){
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
// }}}

// output_frame {{{
void ODE::output_frame(char* name, double t, unsigned int var){
  
  unsigned int nb = domain->getGhostPoints();
  unsigned int size = 0;
  double *x = nullptr;
  for(auto it = data.begin(); it != data.end(); ++it){
    unsigned int shp = it->getGrid().getSize();
    // Unfortunately, SDF requires all arrays passed in
    // to be mutable, but getPoints() only returns a const
    // array. We could simply change getPoints() to return
    // mutable arrays, but that's not a sound design principle.
    // Therefore, we'll just bite the bullet and allocate memory
    // to manually copy the array. This may change if too many
    // people complain. We allocate more memory than we need
    // to limit how often this has to be rebuilt during the
    // output procedure.
    if(shp > size){
      if(x != nullptr){
        delete[] x;
      }
      size = shp*2;
      x = new double[size];
    }
    const double *points = it->getGrid().getPoints();
    for(unsigned int i = 0; i < shp; i++){
      x[i] = points[i];
    }

    double **u = it->getData();
    output::output_data(name, u[var] + nb, x + nb, shp - 2*nb, t);
  }
  delete[] x;
}
// }}}

// dump_frame {{{
void ODE::dump_csv(char* name, double t, unsigned int var){
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
  for(auto it = data.begin(); it != data.end(); ++it){
    //const double *r = it->getGrid().getPoints();
    //double **data = it->getData();
    unsigned int shp = it->getGrid().getSize();
    for(int i = nb; i < shp - nb; i++){
      //fprintf(f,"%08g,   %08g,   %08g\n", t, r[i], data[i][var]);
      fprintf(f,"%08g,   %08g,   %0g\n", t,
              it->getGrid().getPoints()[i],
              it->getData()[var][i]);
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
