#ifndef RK4_H
#define RK4_H

#include "types.h"
#include "grid.h"
#include "domain.h"
#include "solver.h"
#include <memory>

class ODE;
class FieldMap;

/****************************************************************************************
 *
 * Class name: RK4
 * Author: Jacob Fields
 * Date Modified: 12-16-2019
 * 
 * Description: This defines an ODE solver that implements the classic
 *              4th-order Runge-Kutta algorithm.
 *
 * Usage:
 *
 ****************************************************************************************/

class RK4 : public Solver{
  private:
    /**
     * The coefficients for each stage.
     */
    const double kc[4];
    const double tc[4];

  public:
    /**
     * The default constructor for an RK4 object.
     */
    RK4() : Solver(4),
            kc{0.5, 0.5, 1.0, 0.0},
            tc{0.0, 0.5, 0.5, 1.0} {
    };

    /**
     * An RK4 destructor. It frees any memory used.
     */
    virtual ~RK4() {};

    // Inherited methods
    virtual Result setStageTime(double srcTime, double &destTime, double dt, unsigned int stage);
    virtual Result calcStage(ODE *ode, std::shared_ptr<FieldMap>& fieldMap, double dt, unsigned int stage);
    virtual Result combineStages(std::shared_ptr<FieldMap>& fieldMap, double dt);
};

#endif
