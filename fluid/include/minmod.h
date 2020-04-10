#ifndef MINMOD_H
#define MINMOD_H

#include "reconstruction.h"
#include <cmath>

class MinmodRecon : Reconstruction {
  private:
    /**
     * Perform a minmod operation: if x and y have the same sign, return
     * the one with the smallest absolute value. If they have different
     * signs, return zero. The operation is symmetric, so
     * minmod(x,y) = minmod(y,x).
     * @param x The first value for comparison.
     * @param y The second value for comparison.
     * @return The x or y with the smallest absolute value or zero.
     */
    inline double minmod(double x, double y){
      return 0.5*(copysign(1.0,x) + copysign(1.0,y)) * fmin(fabs(x), fabs(y));
    }
  public:
    MinmodRecon();
    virtual ~MinmodRecon();

    virtual Result reconstruct(const int n, const double* const RESTRICT u,
                               double* const RESTRICT ul, double* const RESTRICT ur);
};

#endif
