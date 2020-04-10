#ifndef NORECON_H 
#define NORECON_H 

#include "reconstruction.h"
#include <cmath>

class NoRecon : Reconstruction {
  private:
  public:
    NoRecon();
    virtual ~NoRecon();

    virtual Result reconstruct(const int n, const double* const RESTRICT u,
                               double* const RESTRICT ul, double* const RESTRICT ur);
};

#endif
