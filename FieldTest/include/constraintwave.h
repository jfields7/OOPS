#ifndef CONSTRAINT_WAVE_H
#define CONSTRAINT_WAVE_H

#include <ode.h>
#include <fieldtestparameters.h>

class ConstraintWave : public ODE{
  private:
    // Variable labels
    static const unsigned int U_PHI = 0;
    static const unsigned int U_PI = 1;
    static const unsigned int U_CHI = 2;

    // Constraint labels
    static const unsigned int C_PHI = 0;

    FieldTestParameters *params;
  protected:
    virtual void rhs(std::shared_ptr<FieldMap>& fieldMap);

  public:
    ConstraintWave(Domain& d, Solver& s);
    virtual ~ConstraintWave();

    virtual void initData();

    void setParameters(FieldTestParameters* p);
    FieldTestParameters* getParameters();
};

#endif
