#ifndef ODE_H
#define ODE_H

#include "types.h"
#include "domain.h"
#include "solver.h"
#include "rk4.h"
#include "solverdata.h"
#include "parameters.h"
#include <set>

/**
 * An abstract class representing a system of ODEs. Despite the name, which
 * is *technically* accurate, it can and should be used for implementing
 * PDEs, too.
 */
class ODE{
  protected:
    /**
     * The number of independent equations in this system.
     */
    const unsigned int nEqs;

    /**
     * The type of parameter expected by this system.
     */
    const unsigned int pId;

    /**
     * The parameters for this system of equations.
     */
    Parameters *params;

    /**
     * The domain to solve this ODE on.
     */
    Domain *domain;

    /**
     * The data for this domain.
     */
    std::set<SolverData> data;

    /**
     * The solver to use for this system of ODEs.
     */
    Solver *solver;

    /**
     * The maximum grid spacing on the domain in question.
     */
    double max_dx;

    /**
     * Because everything is dynamically allocated, copying an ODE object
     * is a baaaaaaad idea. We make the copy constructor private to prevent 
     * unintended disaster.
     */
    ODE(const ODE& other) :nEqs(0), pId(0){};

    /**
     * Apply boundary conditions to the data. The default version assumes that
     * it is a true set of ODEs, so no boundaries are necessary.
     */
    virtual void applyBoundaries() {};

    /**
     * Clear and reallocate all the data.
     */
    Result reallocateData();

    /**
     * Exchange the ghost points between the grids. If the grids are different
     * sizes, this involves interpolation.
     */
    void performGridExchange();

    /**
     * Swap the ghost points between two data sets. Again, if the grids are different
     * sizes, we perform interpolation.
     */
    void exchangeGhostPoints(const SolverData &data1, const SolverData &data2);

    /**
     * Interpolate from data2 to get the ghost points for data1.
     * @param datal - The left data, assumed to be finer.
     * @param datar - The right data, assumed to be coarser.
     */
    void interpolateLeft(const SolverData &datal, const SolverData &datar);

    /**
     * Interpolate from data1 to get the ghost points for data2.
     * @param datal - The left data, assumed to be coarser.
     * @param datar - The right data, assumed to be finer.
     */
    void interpolateRight(const SolverData &datal, const SolverData &datar);
      
  public:
    /**
     * Since this is an abstract class, the only thing the constructor needs
     * to do is set the number of equations in the ODE and the id for the
     * expected parameter id.
     */
    /*ODE(const unsigned int n, const unsigned int id) : nEqs(n), pId(id) {
      params = Parameters();
      domain = Domain();
      solver = RK4();
    };*/
    ODE(const unsigned int n, const unsigned int id);

    /**
     * A simple destructor. It just clears the set of data.
     */
    ~ODE();

    /**
     * Set the domain for this problem. This will clear all existing data and
     * reinitialize it for the new domain.
     * @param domain - the new domain to use.
     * @returns SUCCESS or an error message, usually BAD_ALLOC.
     */
    Result setDomain(Domain *domain);

    /**
     * Set the solver for this problem. This will clear all existing data and
     * reinitialize it for the new solver.
     * @param solver - the new solver to use.
     * @returns SUCCESS or an error message, usually BAD_ALLOC.
     */
    Result setSolver(Solver *solver);

    /**
     * The evolution function for a single step. We make this virtual because
     * different ODEs might require different things to happen in between each
     * stage of the solver. A default version is included that simply loops over
     * each stage of the solver and applies boundary conditions.
     * @param dt - The time step to evolve.
     * @returns failure or success during the evolution.
     */
    virtual Result evolveStep(double dt);

    /**
     * The initial condition function for setting up the data based on the
     * Parameter object assigned to the class.
     */
    virtual void initData() = 0;

    /**
     * The righthand side routine for the ODE solver. This, of course, should be
     * overwritten in the descendent ODE.
     * @param grid - The specific grid to perform the calculation on.
     * @param data - A 2d array of data containing the current data for the 
     *               system on the Grid.
     * @param dudt - A 2d array of data to store the righthand side data for
     *               this Grid object.
     */
    virtual void rhs(const Grid& grid, double** data, double** dudt) = 0;

    /**
     * Set the Parameters object for this object. The Parameters id must
     * match what the ODE object has been set to recognize or it returns
     * an error.
     * @param p - The parameter object to associate with this ODE.
     * @returns SUCCESS or UNRECOGNIZED_PARAMS
     */
    Result setParameters(Parameters *p);

    /**
     * Get the number of equations in this system.
     */
    inline unsigned int getNEqs() const {
      return nEqs;
    }

    /**
     * Get the domain this problem is solved on.
     */
    inline Domain *getDomain(){
      return domain;
    }

    /**
     * Get the parameters for this system of ODEs.
     */
    inline Parameters *getParameters(){
      return params;
    }

    /**
     * Output a frame of one variable in the ODE to the specified .sdf file.
     */
    void output_frame(char *name, double t, unsigned int var);

    /**
     * Dump all of the current data to a .csv file.
     */
    void dump_csv(char *name, double t, unsigned int var);
};

#endif
