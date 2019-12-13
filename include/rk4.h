#ifndef RK4_H
#define RK4_H

#include "types.h"
#include "grid.h"
#include "domain.h"
#include "solver.h"

/****************************************************************************************
 *
 * Class name: RK4
 * Author: Jacob Fields
 * Date Modified: 10-23-2019
 * 
 * Description: This defines an ODE solver that implements the classic
 *              4th-order Runge-Kutta algorithm.
 *
 * Usage:
 *
 ****************************************************************************************/

class RK4 : public Solver{
  private:
  public:
    /**
     * The default constructor for an RK4 object.
     */
    RK4() : Solver(4) {};

    /**
     * An RK4 destructor. It frees any memory used.
     */
    ~RK4() {};

    // Inherited methods
    virtual Result calcStage(void (*rhs)(const Grid&,double**,double**), double *data0[], double *dataint[], 
                             double *dest[], const Grid& grid, double dt, const unsigned int vars, int stage);
    virtual Result combineStages(double **data[], double *dest[], Grid& grid, double dt, const int vars);
};

#endif
