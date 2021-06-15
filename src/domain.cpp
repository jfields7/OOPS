#include <domain.h>
#include <iostream>

Domain::Domain(){
  x_bounds[0] = 0.0;
  x_bounds[1] = 1.0;
  d_cfl = 0.5;
  nghosts = 3;
}

Domain::~Domain(){
  grids.clear();
}

void Domain::setCFL(double cfl){
  d_cfl = cfl;
}

double Domain::getCFL() const{
  return d_cfl;
}

void Domain::setBounds(double bounds[2]){
  grids.clear();
  // Make sure that the bounds are ordered correctly.
  // Otherwise, flip them.
  if(bounds[0] < bounds[1]){
    x_bounds[0] = bounds[0];
    x_bounds[1] = bounds[1];
  }
  else{
    x_bounds[0] = bounds[1];
    x_bounds[1] = bounds[0];
  }
}

const double* Domain::getBounds() const{
  return x_bounds;
}

void Domain::setGhostPoints(unsigned int n){
  if(nghosts == n){
    return;
  }
  nghosts = n;
  grids.clear();
}

/*unsigned int Domain::getGhostPoints() const{
  return nghosts;
}*/

std::set<Grid>& Domain::getGrids(){
  return grids;
}

Result Domain::addGrid(double bounds[2], unsigned int n){
  // Check that the bounds are valid.
  if(bounds[0] < bounds[1]){
    if(bounds[0] < x_bounds[0] || bounds[1] > x_bounds[1]){
      return OUT_OF_BOUNDS;
    }
  }
  //grids.insert(Grid(bounds, n, nghosts));
  grids.emplace(bounds, n, nghosts);

  return SUCCESS;
}

void Domain::clearGrids(){
  grids.clear();
}
