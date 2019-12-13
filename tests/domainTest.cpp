#include <domain.h>
#include <grid.h>
#include <cstdio>

// Unit tests for the domain class.

bool testAddGrid(Domain& domain, double bounds[2], unsigned int n){
  Result r = domain.addGrid(bounds,n);
  const double *x_bounds = domain.getBounds();
  // If the bounds are clearly not inside the domain, the grid add should fail.
  switch(r){
    case BAD_ALLOC:
      printf("There was an error while trying to allocate memory for the grid.\n");
      return false;
      break;
    case OUT_OF_BOUNDS:
      if (x_bounds[0] <= bounds[0] && x_bounds[1] >= bounds[1]){
        printf("Valid grid was not added to the domain.\n");
        return false;
      }
      else{
        // There's no need to test any further if the grid wasn't added.
        return true;
      }
      break;
    case SUCCESS:
      if( !(x_bounds[0] <= bounds[0] && x_bounds[1] >= bounds[1])){
        printf("Invalid grid was added to the domain.\n");
        return false;
      }
      break;
  }

  // Now verify that the grid is actually in the set by constructing a grid directly and
  // comparing the two.
  Grid g = Grid(bounds, n, domain.getGhostPoints());

  std::set<Grid>& grids = domain.getGrids();
  if(grids.size() != 1){
    printf("Grid set has wrong size.\n");
    return false;
  }

  std::set<Grid>::iterator it = grids.begin();
  const double* ret_bounds = it->getBounds();
  if(ret_bounds[0] != bounds[0] || ret_bounds[1] != bounds[1]){
    printf("Grid has incorrect bounds.\n");
    return false;
  }

  if(it->getSpacing() != g.getSpacing()){
    printf("Grid has incorrect spacing.\n");
    return false;
  }

  return true;
}

void printTest(char *name, bool success){
  if(success){
    printf("\033[1;32m%s test passed.\033[0m\n\n",name);
  }
  else{
    printf("\033[1;31m%s test failed.\033[0m\n\n",name);
  }
}

int main(int argc, char *argv[]){
  Domain domain = Domain();
  double bounds[2] = {0.0, 0.5};
  int n = 100;

  printTest("Add Valid Grid", testAddGrid(domain,bounds,n));

  return 0;
}
