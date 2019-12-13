#ifndef DOMAIN_H
#define DOMAIN_H

#include "types.h"
#include "grid.h"
#include <set>

/*************************************************************************
 *
 * Class name: Domain
 * Author: Jacob Fields
 * Date Modified: 10-17-2019
 * 
 * Description: This defines the 1d domain used for the solver. It should
 *   maintain a list of all the grids currently in use by this processor.
 *
 *************************************************************************/

class Domain{
  private:
    // The bounds of the problem.
    double x_bounds[2];

    // The CFL to use on the grid.
    double d_cfl;

    // The number of ghost points permitted on the grids
    unsigned int nghosts;

    // An ordered set of all the grids in use.
    std::set<Grid> grids;
    
  public:

    /**
     * The default constructor. It sets up the grid with a set of default
     *   parameters, i.e., a region from 0 to 1 with 3 ghost points and a 
     *   CFL of 0.5.
     */
    Domain();

    /**
     * The destructor. It clears all the grids from the domain. If your 
     *   code is designed properly, they should automatically be deleted.
     *   I am not responsible for your memory leaks.
     */
    ~Domain();

    /**
     * Set the CFL for this domain.
     * @param cfl The CFL to use.
     */
    void setCFL(double cfl);

    /**
     * Get the CFL for this domain.
     * @returns A double indicating the CFL.
     */
    double getCFL() const;

    /**
     * Set the bounds for this domain.
     * WARNING: This will clear any grids still on the domain.
     * @param bounds An ordered pair
     */
    void setBounds(double bounds[2]);

    /**
     * Get the bounds for this domain.
     * @returns An ordered pair specifying the bounds.
     */
    const double* getBounds() const;

    /**
     * Set the number of ghost points to use on this domain.
     * @param n The number of ghost points.
     */
    void setGhostPoints(unsigned int n);
    
    /**
     * Get the number of ghost points used by the domain.
     * @returns An unsigned integer indicating the number of ghost points.
     */
    unsigned int getGhostPoints() const;

    /**
     * Get all the grids used by this domain.
     * @returns A reference to an ordered set of Grid objects.
     */
    std::set<Grid>& getGrids();

    /**
     * Add a new grid to the domain in the specified region with the 
     *   specified number of points.
     * @param bounds The boundaries of the new grid.
     * @param n The number of points to use on the grid.
     * @returns A Result indicating success or failure.
     */
    Result addGrid(double bounds[2], unsigned int n);

    /**
     * Clear all the grids on the domain. See the comment on the destructor.
     */
    void clearGrids();

};

#endif
