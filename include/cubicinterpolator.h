#ifndef CUBIC_INTERPOLATOR_H
#define CUBIC_INTERPOLATOR_H

#include <interpolator.h>

/**
 * An object that performs a cubic (four-point) interpolation.
 */
class CubicInterpolator : public Interpolator{
  public:
    /**
     * Construct a cubic interpolator object.
     */
    CubicInterpolator();

    /**
     * Destroy a cubic interpolator object.
     */
    virtual ~CubicInterpolator();

    /**
     * Perform a centered cubic interpolation on the stencil.
     */
    virtual double interpolate();
};

#endif
