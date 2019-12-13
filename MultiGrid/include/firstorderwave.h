#ifndef FIRST_ORDER_WAVE_H
#define FIRST_ORDER_WAVE_H

#include <ode.h>

class FirstOrderWave : public ODE {
  private:
    // Variable labels.
    const unsigned int U_PHI = 0;
    const unsigned int U_PI = 1;
    const unsigned int U_CHI = 2;

    void applyKODiss(Grid& grid, double **u, double **dudt);

    void applyGaussian();
  protected:
    virtual void rhs(Grid& grid, double **u, double **dudt);
    virtual void applyBoundaries();

  public:
    FirstOrderWave(Domain& d, Solver& s);
    virtual ~FirstOrderWave();

    virtual void initData();
};

#endif
