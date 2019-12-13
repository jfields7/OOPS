#ifndef SOLVER_H
#define SOLVER_H

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
     * Calculate a single stage of the solver. Some methods, such as conservation laws,
     * require specialized behavior after each stage of the calculation, so we force
     * this to be split.
     * @param rhs A function representing the righthand side of the ODE.
     * @param data0 The original data coming into the solver.
     * @param dataint The intermediate data from the previous stage. This is updated
     *                during each RK stage to contain the intermediate value of the
     *                solution that will be used in the next stage.
     * @param dest The location for each stage's increment/righthand side.
     * @param grid The grid this dataset belongs to.
     * @param dt The size of time step to use for this solver.
     * @param vars The number of variables in the system.
     * @param stage Which stage of the solver should be calculated.
     * @return If the operation was successful or why it failed.
     */
    virtual Result calcStage(void (*rhs)(const Grid&, double**,double**), double *data0[], double *dataint[],
                             double *dest[], const Grid& grid, double dt, const unsigned int vars, 
                             unsigned int stage) = 0;

    /**
     * Combine all the calculations from every stage.
     * @param data The data coming into the solver. This should be an array of all
     *             the stages, so it should contain nStages * shp points.
     * @param dest The location for the final result to be stored. This is assumed
     *             to contain the original data and is so utilized in the calculation.
     * @param grid The grid this dataset belongs to.
     * @param dt The size of time step used for this solver.
     * @param vars The number of variables in the system.
     * @return A Result enum designating what the outcome of the calculation was.
     */
    virtual Result combineStages(double **data[], double *dest[], const Grid& grid, double dt, const int vars) = 0;

    /**
     * Find out how many stages the solver has.
     * @return An integer containing the number of stages used in the solver.
     */
    inline int getNStages() const{
      return nStages;
    }
};
#endif
