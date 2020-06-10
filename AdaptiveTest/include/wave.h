#ifndef WAVE_H
#define WAVE_H

#include <ode.h>
#include <waveparameters.h>

class Wave : public ODE {
  private:
    // Variable labels
    static const unsigned int U_PHI = 0;
    static const unsigned int U_PI = 1;
    static const unsigned int U_CHI = 2;

    WaveParameters *params;

  protected:
    virtual void applyBoundaries();

    virtual void rhs(std::shared_ptr<FieldMap>& fieldMap);

  public:
    Wave(Domain& d, Solver& s);
    virtual ~Wave();

    virtual void initData();

    void setParameters(WaveParameters *p);
    WaveParameters* getParameters();
};

#endif
