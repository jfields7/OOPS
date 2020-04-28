#include <grid.h>
#include <iostream>
#include <new>
#include <cmath>

// Constructor
Grid::Grid(const double bounds[2], unsigned int n, unsigned int nghosts){
  // Make sure the bounds are actually ordered correctly.
  // If they aren't, swap them so they are.
  if(bounds[0] < bounds[1]){
    x_bounds[0] = bounds[0];
    x_bounds[1] = bounds[1];
  }
  else{
    x_bounds[0] = bounds[1];
    x_bounds[1] = bounds[0];
  }

  // Now we need to figure out our grid spacing. Assume it's a cell-edge grid.
  dx = (x_bounds[1] - x_bounds[0])/(double)(n-1);

  // Next, we need to construct the grid itself. Add the ghost points to the grid.
  nx = n + 2*nghosts;
  // Allocate memory for the new grid and calculate all its positions.
  try{
    points = new double[nx];
    for(int i = 0; i < nx; i++){
      points[i] = x_bounds[0] + (i - (int)nghosts)*dx;
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for grid.\n";
    nx = 0;
    dx = 0;
    points = NULL;
  }
}

// Copy constructor
Grid::Grid(const Grid& other){
  x_bounds[0] = other.getBounds()[0];
  x_bounds[1] = other.getBounds()[1];
  dx = other.getSpacing();
  nx = other.getSize();
  try{
    points = new double[nx];
    for(int i = 0; i < nx; i++){
      points[i] = other.getPoints()[i];
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for grid.\n";
    nx = 0;
    dx = 0;
    points = NULL;
  }
}

// Destructor
Grid::~Grid(){
  // Clear all the memory from the points.
  delete[] points;
}

// Build a new grid with new parameters.
Result Grid::rebuildGrid(const double bounds[2], unsigned int n, unsigned int nghosts){
  // Delete the old grid. If the construction failed, it should be NULL anyway, so
  // this shouldn't fail.
  delete[] points;

  // Make sure the bounds are actually ordered correctly.
  // Swap them if they aren't.
  if(bounds[0] < bounds[1]){
    x_bounds[0] = bounds[0];
    x_bounds[1] = bounds[1];
  }
  else{
    x_bounds[0] = bounds[1];
    x_bounds[1] = bounds[0];
  }

  // Now set up the grid spacing. Assume it's a cell-edge grid.
  dx = (x_bounds[1] - x_bounds[0])/(double)(n-1);

  // Construct the new grid.
  nx = n + 2*nghosts;
  // Try to allocate memory for a new grid.
  try{
    points = new double[nx];
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for grid.\n";
    n = 0;
    dx = 0;
    points = NULL;
    return BAD_ALLOC;
  }

  // Rebuild the points.
  for(int i = 0; i < nx; i++){
    points[i] = x_bounds[0] + (i-(int)nghosts)*dx;
  }
  
  return SUCCESS;
}

// Get all the points.
const double* Grid::getPoints() const{
  return points;
}

// Get the size of the grid.
const unsigned int Grid::getSize() const{
  return nx;
}

// Get the grid spacing.
const double Grid::getSpacing() const{
  return dx;
}

// Get the bounds of the grid.
const double* Grid::getBounds() const{
  return x_bounds;
}

// Compare this grid to another by its bounds.
bool Grid::operator < (const Grid& g) const{
  const double *bounds = g.getBounds();
  if(x_bounds[0] < bounds[0] && x_bounds[1] < bounds[1]){
    return true;
  }
  return false;
}

// Find out if two grids share a boundary.
Boundary Grid::whichNeighbor(const Grid& g) const{
  const double *bounds = g.getBounds();
  // Because numerical error sometimes creeps in, we can't compare the bounds directly and
  // expect them to be equal. Instead, we take the smallest dx from the grids and compare
  // the boundaries using that.
  double buffer = 0.5*fmin(dx, g.getSpacing());
  if(fabs(x_bounds[1] - bounds[0]) < buffer){
    return RIGHT;
  }
  if(fabs(x_bounds[0] - bounds[1]) < buffer){
    return LEFT;
  }
  return NONE;
}
