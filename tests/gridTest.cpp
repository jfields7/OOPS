// Test whether or not the Grid class is working properly.
#include <grid.h>
#include <cstdio>

bool testGrid(Grid& grid, const double bounds[], int n, int nghosts){
  double dx = (bounds[1] - bounds[0])/(n-1);
  // Check that the grid spacing is correct.
  if(dx != grid.getSpacing()){
    printf("\033[1;31mIncorrect grid spacing.\033[0m\n");
    printf("  Expected: %g\n",dx);
    printf("  Actual: %g\n",grid.getSpacing());
    return false;
  }
  // Check that the grid count is correct.
  if(grid.getSize() != n + 2*nghosts){
    printf("\033[1;31mIncorrect point count.\033[0m\n");
    printf("  Expected: %d\n",n + 2*nghosts);
    printf("  Actual: %d\n",grid.getSize());
    return false;
  }
  // Make sure the bounds match.
  if((grid.getBounds())[0] != bounds[0]){
    printf("\033[1;31mLeft bound is incorrect.\033[0m\n");
    return false;
  }
  if((grid.getBounds())[1] != bounds[1]){
    printf("\033[1;31mRight bound is incorrect.\033[0m\n");
    return false;
  }
  // Check that the grid itself is correct.
  const double *points = grid.getPoints();
  for(int i = 0; i < n + 2*nghosts; i++){
    if(points[i] != (i-nghosts)*dx){
      printf("\033[1;31mIncorrect grid.\033[0m\n");
      printf("  Index %d\n",i);
      printf("  Expected: %g\n",(i-nghosts)*dx);
      printf("  Actual: %g\n",points[i]);
      return false;
    }
  }

  return true;
}

bool lessThan(Grid& a, Grid& b){
  const double *boundsA = a.getBounds();
  const double *boundsB = b.getBounds();

  bool result = (boundsA[0] < boundsB[0] && boundsA[1] < boundsB[1]);
  return result == (a < b);
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
  double bounds[] = {0.0, 1.0};
  Grid grid = Grid(bounds, 100, 3);
  // Test that the constructor works.
  printTest("Constructor", testGrid(grid,bounds,100,3));
  
  // Test that we can rebuild the grid with a different set of parameters.
  bounds[0] = 0.0;
  bounds[1] = 10.0;
  grid.rebuildGrid(bounds, 200, 2);
  printTest("Rebuild", testGrid(grid,bounds,200,2));

  // Now we need to do a comparison test.
  double bounds2[] = {10.0, 12.0};
  Grid grid2 = Grid(bounds2, 50, 2);
  printTest("Less Than Comparison", lessThan(grid,grid2));
  bounds2[0] = -1.0;
  bounds2[1] = 0.0;
  grid2.rebuildGrid(bounds2, 50, 2);
  printTest("Greater Than Comparison", lessThan(grid,grid2));
  bounds[0] = 0.25;
  bounds[1] = 0.50;
  printTest("Nested Comparison", lessThan(grid,grid2));

  // Do a neighbor test.
  bounds2[0] = 10.0;
  bounds2[1] = 12.0;
  grid2.rebuildGrid(bounds2, 50, 2);
  printTest("Right Neighbor Test", grid.whichNeighbor(grid2) == RIGHT);
  bounds2[0] = -2.0;
  bounds2[1] = 0.0;
  grid2.rebuildGrid(bounds2, 50, 2);
  printTest("Left Neighbor Test", grid.whichNeighbor(grid2) == LEFT);
  bounds2[0] = 11.0;
  bounds2[1] = 12.0;
  grid2.rebuildGrid(bounds2, 50, 2);
  printTest("Not Neighboring Test - Right", grid.whichNeighbor(grid2) == NONE);
  bounds2[0] = -2.0;
  bounds2[1] = -1.0;
  grid2.rebuildGrid(bounds2, 50, 2);
  printTest("Not Neighboring Test - Left", grid.whichNeighbor(grid2) == NONE);

  return 0;
}
