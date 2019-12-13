#ifndef SOLVERDATA_H
#define SOLVERDATA_H

#include "types.h"
#include "grid.h"

class SolverData{
  private:
  // Private members {{{
  /**
   * The number of independent equations in the system.
   */
  unsigned int nEq;

  /**
   * The number of stages used by the solver.
   */
  unsigned int nStages;

  /**
   * The raw data for the solution.
   */
  double **data;

  /**
   * The intermediate data for the solution.
   */
  double **data_int;

  /**
   * The storage for the results from the ODE solver.
   */
  double ***work;

  /**
   * The grid this data belongs to.
   */
  Grid& grid;
  // }}}

  /**
   * A copy constructor for SolverData. We make this private so that
   * the object cannot be copied.
   */
  //SolverData(const SolverData& other);

  public:
  /**
   * The constructor for SolverData.
   * @param eqCount - the number of equations in the system.
   * @param stages  - the number of stages used by the solver.
   * @param grid    - the Grid this data belongs to.
   */
  SolverData(unsigned int eqCount, unsigned int stages, Grid& grid);

  /**
   * A destructor for SolverData that deletes all the memory that was allocated.
   */
  ~SolverData();

  // Inline getter functions {{{
  /**
   * Get the array containing the data. This is stored as [points][vars].
   */
  inline double** getData() const{
    return data;
  }

  /**
   * Get the array containing the intermediate data. This is stored as [points][vars].
   */
  inline double** getIntermediateData() const{
    return data_int;
  }

  /**
   * Get the array containing the work data. This is stored as [stage][points][vars].
   */
  inline double*** getWorkData() const{
    return work;
  }

  /**
   * Get the number of equations for the ODE this represents.
   */
  inline unsigned int getEqCount() const{
    return nEq;
  }

  /**
   * Get the grid this data corresponds to.
   */
  inline const Grid& getGrid() const{
    return grid;
  }
  // }}}

  /**
   * Compare this SolverData object to another SolverData object to find out how
   * they are spatially related. In other words, this compares the associated grids
   * for each data object.
   * @returns true if this data is to the left of the other data.
   */
  bool operator < (const SolverData& data) const;
};

#endif
