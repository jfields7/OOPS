#ifndef FIRST_ORDER_WAVE_H
#define FIRST_ORDER_WAVE_H

#include <ode.h>
#include <waveparameters.h>

class FirstOrderWave : public ODE {
  private:
    // Variable labels.
    static const unsigned int U_PHI = 0;
    static const unsigned int U_PI = 1;
    static const unsigned int U_CHI = 2;

    void applyKODiss(const Grid& grid, double **u, double **dudt);

    void applyGaussian();

    WaveParameters *params;
  protected:
    virtual void applyBoundaries(bool intermediate);

    virtual void rhs(const Grid& grid, double **u, double **dudt);

  public:
    FirstOrderWave(Domain& d, Solver& s);
    virtual ~FirstOrderWave();

    virtual void initData();

    void setParameters(WaveParameters* p);
    WaveParameters* getParameters();
};

#endif
