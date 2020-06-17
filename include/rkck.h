#ifndef RKCK_H
#define RKCK_H

#include "types.h"
#include "grid.h"
#include "domain.h"
#include "solver.h"
#include <memory>

class ODE;
class FieldMap;

/****************************************************************************************
 *
 * Class name: RKCK
 * Author: Jacob Fields
 * Date Modified: 6-10-2020
 * 
 * Description: This defines an ODE solver that implements the the Cash-Karp
 *              5th-order adaptive integrator.
 *
 * Usage:
 *
 ****************************************************************************************/

class RKCK : public Solver{
  private:
    /**
     * The coefficients for each stage.
     */
    const double kc[6][6];
    const double tc[6];
    const double c4[6];
    const double c5[6];

    /**
     * The recommended step size.
     */
    double dtrec;

    /**
     * The error tolerance.
     */
    double tol;

  public:
    /**
     * The default constructor for an RKCK object.
     */
    RKCK() : Solver(6),
             kc{{ 1.0/5.0     ,  0          ,  0            ,  0               ,  0           ,  0},
                { 3.0/40.0    ,  9.0/40.0   ,  0            ,  0               ,  0           ,  0},
                { 3.0/10.0    , -9.0/10.0   ,  6.0/5.0      ,  0               ,  0           ,  0},
                {-11.0/54.0   ,  5.0/2.0    , -70.0/27.0    ,  35.0/27.0       ,  0           ,  0},
                { 1631.0/55296,  175.0/512.0,  575.0/13824.0,  44275.0/110592.0,  253.0/4096.0,  0},
                { 0           ,  0          ,  0            ,  0               ,  0           ,  0}},
             tc{0.0, 0.2, 0.3, 0.6, 1.0, 0.875},
             c5{37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0},
             c4{2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0} {
      dtrec = 1.0;
      tol = 1e-6;
    };

    /**
     * An RKCK destructor. It frees any memory used.
     */
    virtual ~RKCK(){};

    /**
     * Inherited methods
     */
    virtual Result setStageTime(double srcTime, double &destTime, double dt, unsigned int stage);
    virtual Result calcStage(ODE *ode, std::shared_ptr<FieldMap>& fieldMap, double dt, unsigned int stage);
    virtual Result combineStages(std::shared_ptr<FieldMap>& fieldMap, double dt);

    /**
     * Get the recommended step size.
     */
    inline double getRecommendedStepSize() const{
      return dtrec;
    }

    /**
     * Get the local error tolerance.
     */
    inline double getErrorTolerance() const{
      return tol;
    }

    /**
     * Set the local error tolerance.
     */
    inline void setErrorTolerance(double t){
      tol = t;
    }

    /**
     * Reset the recommended time step.
     */
    inline void resetRecommendedStepSize(double dt){
      dtrec = dt;
    }
};

#endif
