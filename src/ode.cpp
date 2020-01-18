#include <ode.h>
#include <cubic.h>
#include <cmath>
#include <iostream>
#include <cstdio>

ODE::ODE(const unsigned int n, const unsigned int id) : nEqs(n), pId(id){
  domain = nullptr;
  params = nullptr;
  solver = nullptr;
  max_dx = 0.0;
  interpolator = nullptr;
}

ODE::~ODE(){
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
  // Loop over every stage for the solver.
  for(unsigned int i = 0; i < solver->getNStages(); i++){
    // Loop over every data set in the domain.
    for(auto it = data.begin(); it != data.end(); ++it){
      //solver->calcStage(rhs, it->getData(), it->getIntermediateData(), (it->getWorkData())[i],
      //                 it->getGrid(), dt, nEqs, i);
      solver->calcStage(this, it->getData(), it->getIntermediateData(), (it->getWorkData())[i],
                        it->getGrid(), dt, i);
    }

    // Perform the grid exchange and apply boundary conditions.
    performGridExchange();
    applyBoundaries(true);
  }

  for(auto it = data.begin(); it != data.end(); ++it){
    solver->combineStages(it->getWorkData(), it->getData(), it->getGrid(), dt, nEqs);
  }
  performGridExchange();
  applyBoundaries(false);

  // Loop over all the data sets.
  /*for(auto it = data.begin(); it != data.end(); ++it){
    // Figure out how many time steps we need to take for this grid.
    unsigned int N = (unsigned int) round(max_dx/it->getGrid().getSpacing());
    double dt_i = dt/N;
    for(int t = 0; t < N; t++){
      for(unsigned int i = 0; i < solver->getNStages(); i++){
        solver->calcStage(this, it->getData(), it->getIntermediateData(), (it->getWorkData())[i],
                          it->getGrid(), dt_i, i);
        applyBoundaries();
      }

      solver->combineStages(it->getWorkData(), it->getData(), it->getGrid(), dt_i, nEqs);
      applyBoundaries();
    }
  }

  performGridExchange();*/

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
    for(unsigned int i = 0; i < nb; i++){
      for(unsigned int j = 0; j < nEqs; j++){
        // Copy the first physical points on u2 into the ghost points of u1.
        u1[shp1 - nb + i][j] = u2[nb + i][j];
        // Copy the last physical points of u1 into the ghost points of u2.
        u2[i][j] = u1[shp1 - 2*nb + i][j];
      }
    }
  }
  // Otherwise, let's perform some interpolation. We first need to enforce a 2x difference
  // between the grids and complain if that's not the case.
  if(fabs((dx1*2.0 - dx2)/dx2) < 1e-14){
    // The left grid is finer than the right grid.
    // Note the nb + 1. We also need to copy the one physical point so that they stay
    // consistent. We use the one from the more refined grid.
    for(unsigned int i = 0; i < nb + 1; i++){
      for(unsigned int j = 0; j < nEqs; j++){
        // Copy the exact point from u1 into the ghost region of u2.
        u2[i][j] = u1[shp1 - 1 - 3*nb + 2*i][j];
      }
    }
    interpolateLeft(data1, data2);
  }
  else if(fabs((dx1 - 2.0*dx2)/dx2) < 1e-14){
    // The right grid is finer than the left grid.
    //for(unsigned int i = 0; i < nb; i++){
    // Note the nb + 1. We also need to copy the one physical point so that they stay
    // consistent. We use the one from the more refined grid.
    for(unsigned int i = 0; i < nb + 1; i++){
      for(unsigned int j = 0; j < nEqs; j++){
        // Copy the exact point from u2 into the ghost region of u1.
        u1[shp1 - nb - 1 + i][j] = u2[nb + 2*i][j];
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
  for(unsigned int i = 0; i < nb; i++){
    for(unsigned int j = 0; j < nEqs; j++){
      // If we're on an even ghost point, we don't need to interpolate. This corresponds to an
      // odd index because counting starts at 0, hence the screwy math.
      if( i & 1 == 1){
        ul[shpl - nb + i][j] = ur[nb + (i + 1)/2][j];
      }
      else{
        /*ul[shpl - nb + i][j] = interp::cubicInterpCenter(ul[shpl - nb + i - 3][j], 
                                                         ul[shpl - nb + i - 1][j],
                                                         ur[nb + i/2 + 1][j],
                                                         ur[nb + i/2 + 2][j]);*/
        // Fill in the stencil. Use the left points from the finer grid and the right
        // points from the coarser grid. This makes sure that we don't run out of points.
        for(int m = -nStart; m < 0; m++){
          stencil[m + nStart] = ul[shpl - nb + i + 1 + 2*m][j];
        }
        for(int m = 0; m < nStart; m++){
          stencil[m + nStart] = ur[nb + i/2 + 1 + m][j];
        }
        ul[shpl - nb + i][j] = interpolator->interpolate();
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
  for(unsigned int i = 0; i < nb; i++){
    for(unsigned int j = 0; j < nEqs; j++){
      // If we're on an even ghost point, we don't need to interpolate. This corresponds to an
      // odd index because counting starts at 0, hence the screwy math.
      if(i & 1 == 1){
        ur[nb - 1 - i][j] = ul[shpl - nb - 1 - (i + 1)/2][j];
      }
      else{
        /*ur[nb - 1 - i][j] = interp::cubicInterpCenter(ul[shpl - nb - 3 - i/2][j],
                                                      ul[shpl - nb - 2 - i/2][j],
                                                      ur[nb - i][j],
                                                      ur[nb + 2 - i][j]);*/
        // Fill in the stencil. Use the right points from the finer grid and the left
        // points from the coarser grid. This makes sure that we don't run out of points.
        for(int m = -nStart; m < 0; m++){
          stencil[m + nStart] = ul[shpl - nb - 1 + m - i/2][j];
        }
        for(int m = 0; m < nStart; m++){
          stencil[m + nStart] = ur[nb + 2*m - i][j];
        }
        ur[nb - 1 - i][j] = interpolator->interpolate();
      }
    }
  }
}
// }}}

// setParameters {{{
Result ODE::setParameters(Parameters *p){
  if(p->getId() == pId){
    params = p;
    return SUCCESS;
  }
  else{
    return UNRECOGNIZED_PARAMS;
  }
}
// }}}

// output_frame {{{
void ODE::output_frame(char* name, double t, unsigned int var){
  
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
              it->getData()[i][var]);
    }
  }

  fclose(f);
}
// }}}
