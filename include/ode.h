#ifndef ODE_H
#define ODE_H

#include "types.h"
#include "domain.h"
#include "solver.h"
#include "rk4.h"
#include "odedata.h"
#include "solverdata.h"
#include "parameters.h"
#include "interpolator.h"
#include <set>
#include <map>
#include <string>
#include <vector>
#include <fieldmap.h>

// FieldInfo {{{
/**
 * A struct to help us keep track of all the fields when we reallocate
 * data.
 */
struct FieldInfo{
  std::string name;
  unsigned int nEqs;
  unsigned int nStages;
  FieldInfo(std::string n, unsigned int eqs, unsigned int stages){
    name = n;
    nEqs = eqs;
    nStages = stages;
  }
  FieldInfo(const FieldInfo& other){
    name = std::string(other.name);
    nEqs = other.nEqs;
    nStages = other.nStages;
  }
};
// }}}

/**
 * An abstract class representing a system of ODEs. Despite the name, which
 * is *technically* accurate, it can and should be used for implementing
 * PDEs, too.
 */
class ODE{
  public:

  private:
    /**
     * A list of all the fields and their info.
     */
    std::map<std::string, FieldInfo> fieldList;
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
     * The domain to solve this ODE on.
     */
    Domain *domain;

    /**
     * The data for this domain.
     */
    std::set<SolverData> data;

    /**
     * The data for all fields on this domain.
     */
    std::set<std::unique_ptr<FieldMap>> fieldData;

    /**
     * The indices of the variables that should be evolved.
     */
    std::vector<unsigned int> evolutionIndices;

    /**
     * The solver to use for this system of ODEs.
     */
    Solver *solver;

    /**
     * The interpolation method to use for multiresolution grids.
     */
    Interpolator *interpolator;
    

    /**
     * The maximum grid spacing on the domain in question.
     */
    double max_dx;

    /**
     * The evolution time.
     */
    double time;

    /**
     * Because everything is dynamically allocated, copying an ODE object
     * is a baaaaaaad idea. We make the copy constructor private to prevent 
     * unintended disaster.
     */
    ODE(const ODE& other) :nEqs(0), pId(0){};

    /**
     * Apply boundary conditions to the data. The default version assumes that
     * it is a true set of ODEs, so no boundaries are necessary.
     * @param intermediate - Whether we're modifying the intermediate dataset
     *                       or the original dataset.
     */
    virtual void applyBoundaries(bool intermediate) {};

    /**
     * A function that does nothing by default but can be overwritten in a base
     * class to do stuff after the solver is called but before the grid data is
     * exchanged.
     * @param intermediate - Whether we're working with an intermediate solver
     *                       stage or the final combining step.
     */
    virtual void doAfterStage(bool intermediate){};

    /**
     * A function that does nothing by default but can be overwritten in a base
     * class to do stuff after the grid data is exchanged but before the
     * boundary conditions are applied.
     * @param intermediate - Whether we're working with an intermediate solver
     *                       stage or the final combining step.
     */
    virtual void doAfterExchange(bool intermediate){};

    /**
     * A function that does nothing by default but can be overwritten in a base
     * class to do stuff after the boundaries are applied.
     * @param intermediate - Whether we're working with an intermediate solver
     *                       stage or the final combining step.
     */
    virtual void doAfterBoundaries(bool intermediate){};

    /**
     * Add a new field to the ODE object. The reallocateData() function
     * needs to be called for any memory to be updated.
     */
    Result addField(std::string name, unsigned int eqs, bool isEvolved);

    /**
     * Remove an field from the ODE object.
     */
    Result removeField(std::string name);

    /**
     * Check if a particular field exists.
     */
    bool hasField(std::string name);

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
     * Set the interpolation method for this problem.
     * @param interp - the new interpolator to use.
     * @returns SUCCESS or an error message, usually BAD_ALLOC.
     */
    Result setInterpolator(Interpolator *interp);

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
     * Output a frame of one variable in the ODE to the specified .sdf file.
     */
    void output_frame(char *name, double t, unsigned int var);

    /**
     * Output a frame of one variable in an ODE field to the specified .sdf file.
     */
    void output_field(std::string field, char* name, double t, unsigned int var);

    /**
     * Dump all of the current data to a .csv file.
     */
    void dump_csv(char *name, double t, unsigned int var);

    /**
     * Get the evolution time.
     * @returns the current time in the evolution.
     */
    double getTime();

    /**
     * Set the evolution time.
     * @param t - The evolution time that the ODE should be at.
     */
    void setTime(double t);

};

#endif
