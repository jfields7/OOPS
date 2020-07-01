#ifndef ODEDATA_H
#define ODEDATA_H

#include "types.h"
#include "grid.h"

class ODEData {
  protected:
    /**
     * A copy constructor for ODEData. We make this protected so that
     * the object cannot be copied.
     */
    ODEData(const ODEData& other);

    // Protected members {{{
    /**
     * The number of independent equations in this system.
     */
    unsigned int nEq;

    /**
     * The number of extra lines to allocate.
     */
    unsigned int lines;

    /**
     * Raw data for an ODE.
     */
    double **data;

    /**
     * Empty arrays to use to get a line of data.
     */
    double ***line;

    /**
     * The Grid this data belongs to.
     */
    const Grid& mGrid;
    // }}}
  public:
  /**
   * The constructor for ODEData.
   * @param eqCount - the number of equations in this system.
   * @param grid    - the Grid this data belongs to.
   * @param lines   - the number of extra lines to store.
   */
  ODEData(unsigned int eqCount, const Grid& grid, unsigned int lines=0);

  /**
   * A destructor for ODEData that deletes all the memory that was allocated.
   */
  virtual ~ODEData();

  // Inline getter functions {{{
  /**
   * Get the array containing the data. This is stored as [vars][points].
   */
  inline double** getData() const{
    return data;
  }

  /**
   * Get an empty array for a line of data. This is stored as [vars][points].
   */
  inline double** getLine(unsigned int i) const{
    return line[i];
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
    return mGrid;
  }
  // }}}

  /**
   * Compare this ODEData object to another ODEData object to find out how
   * they are spatially related. In other words, this compares the associated grids
   * for each data object.
   * @returns true if this data is to the left of the other data.
   */
  bool operator < (const ODEData& data) const;
};

#endif
