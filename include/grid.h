#ifndef GRID_H
#define GRID_H

#include "types.h"

/*************************************************************************
 *
 * Class name: Grid
 * Author: Jacob Fields
 * Date Modified: 10-17-2019
 * 
 * Description: This defines a 1d grid with fixed spacing.
 *
 * Usage: The grid has no default constructor. Construct a new grid with
 *        Grid(double bounds[2], int n) 
 *
 *
 *************************************************************************/

class Grid{
  private:
    // Grid spacing.
    double dx;
    // Number of grid points.
    unsigned int nx;
    // The start and end of the grid.
    double x_bounds[2];
    // The grid points themselves.
    double *points;

  public:
    /**
     * A grid constructor. There is no default constructor.
     * @param bounds An ordered pair describing the left- and rightmost
     *   bounds of the grid. This should not include ghost points.
     * @param n The number of points (excluding ghost points) in the grid.
     * @param nghosts The number of ghost points on the grid.
     */
    Grid(const double bounds[2], unsigned int n, unsigned int nghosts);

    /**
     * A copy constructor. Because of the memory allocated for our points,
     * we need to define a copy constructor that performs a deep copy.
     */
    Grid(const Grid& other);

    /**
     * A grid destructor. It frees any memory used.
     */
    ~Grid();

    /**
     * Free the memory on the current grid, then rebuild it with the new
     *   parameters.
     * @param bounds An ordered pair describing the left- and rightmost
     *   bounds of the grid. This should not include ghost points.
     * @param n The number of points (excluding ghost points) in the grid.
     * @nghosts The number of ghost points on the grid.
     * @returns A Result enumerator indicating success or what the
     *   specific error is.
     */
    Result rebuildGrid(double bounds[2], unsigned int n, unsigned int nghosts);

    /**
     * Get the data points on the grid, including ghost points.
     * @returns A const double* containing all the points.
     */
    const double* getPoints() const;

    /**
     * Get the number of grid points. This does include ghost points.
     * @returns An int containing the number of grid points.
     */
    unsigned int getSize() const;

    /**
     * Get the spacing dx of the grid.
     * @returns A double for the grid spacing.
     */
    double getSpacing() const;

    /**
     * Get the bounds (without ghost points) of this grid.
     * @returns An ordered pair containing the left and right bounds.
     */
     const double* getBounds() const;

    /**
     * Compare this grid to another grid by the bounds (i.e., is the right
     *   bound less than or equal to the other grid's left bound).
     * @returns true if this grid is to the left of the other grid.
     */
    bool operator < (const Grid& g) const;

    /**
     * Compare this grid to another grid and determine if they are neighbors.
     * @param g The other grid to compare this one to.
     * @return Which boundary, if any, is shared.
     */
    Boundary whichNeighbor(const Grid& g) const;
};

#endif
