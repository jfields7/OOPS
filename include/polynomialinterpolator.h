#ifndef POLYNOMIAL_INTERPOLATOR_H
#define POLYNOMIAL_INTERPOLATOR
#include <interpolator.h>

/**
 * An object that performs an n-point centered polynomial interpolation.
 * Note: the number of points for interpolation can only be selected in the
 * constructor. Changing the number of points requires building a completely
 * new interpolator.
 */
class PolynomialInterpolator : public Interpolator{
  private:
    /**
     * The interpolation weights calculated using the method presented in Holmstrom 1996.
     */
    double *weights;

    /**
     * Calculate a new interpolation weight at a given index.
     */
    double calculateWeight(int l);

  public:
    /**
     * Construct a polynomial interpolator object.
     */
    PolynomialInterpolator(const unsigned int p);

    /**
     * Destroy a polynomial interpolator object.
     */
    virtual ~PolynomialInterpolator();

    /**
     * Perform a centered n-point interpolation.
     */
    virtual double interpolate();
};

#endif
