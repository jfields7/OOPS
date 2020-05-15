#ifndef SOLVER_H
#define SOLVER_H
#include "grid.h"
#include <vector>

class ODE;

/**************************************************************************************
 *
 * Class name: Solver
 * Author: Jacob Fields
 * Date Modified: 10-23-2019
 *
 * Description: An abstract class that can be used to create a new ODE solver for
 *              initial value problems. Consequently, it should be independent of the
 *              number of dimensions of the problem.
 *
 * Usage: 
 *
 **************************************************************************************/

class Solver{
  protected:
    /**
     * The number of stages used by this solver.
     */
    const int nStages;
    
  public:
    /**
     * A very basic constructor that just sets nStages.
     */
    Solver(const int n) : nStages(n){};

    /**
     *
     */
    virtual Result setStageTime(double srcTime, double &destTime, double dt, unsigned int stage) = 0;

    /**
     * Calculate a single stage of the solver with vital information being provided by an ODE object. 
     * Although the grid and work data are all included in the ODE object, we need to specify them anyway
     * so that we only operate on a single grid.
     */
    virtual Result calcStage(ODE *ode, double *data0[], double *dataint[], double *dest[], const Grid& grid,
                             double dt, unsigned int stage) = 0;

    /**
     * Combine all the calculations from every stage.
     * @param data The data coming into the solver. This should be an array of all
     *             the stages, so it should contain nStages * shp points.
     * @param dest The location for the final result to be stored. This is assumed
     *             to contain the original data and is so utilized in the calculation.
     * @param grid The grid this dataset belongs to.
     * @param dt The size of time step used for this solver.
     * @param evolutionIndices The indices for the equations that should be evolved.
     * @return A Result enum designating what the outcome of the calculation was.
     */
    virtual Result combineStages(double **data[], double *dest[], const Grid& grid, double dt, 
      const std::vector<unsigned int>& evolutionIndices) = 0;

    /**
     * Find out how many stages the solver has.
     * @return An integer containing the number of stages used in the solver.
     */
    inline int getNStages() const{
      return nStages;
    }
};
#endif
