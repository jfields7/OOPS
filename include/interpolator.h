#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

/**
 * An abstract class representing a generic interpolation scheme.
 */
class Interpolator{
  protected:
    /**
     * The stencil size for this interpolation scheme.
     */
    const unsigned int nStencil;

    /**
     * The stencil for this interpolator. This is allocated
     * dynamically when the scheme is created and destroyed
     * when the destructor is called.
     */
    double *stencil;
  public:
    /**
     * The constructor for this class sets the stencil size
     * and allocates memory.
     */
    Interpolator(const unsigned int n);

    /**
     * The destructor. It just clears the stencil memory.
     */
    ~Interpolator();

    /**
     * Get the stencil so it can be modified. 
     * @return the stencil of size nStencil.
     */
    inline double* getStencil(){
      return stencil;
    }

    /**
     * Get the stencil size.
     * @return the stencil size
     */
    inline unsigned int getStencilSize(){
      return nStencil;
    }

    /**
     * An overridable interpolation scheme. Ideally, the 
     * child object overriding it should be a centered scheme,
     * interpolating between indices nStencil/2 - 1 and
     * nStencil/2. This also means that interpolation schemes
     * should ideally have even stencils (linear, cubic, etc.).
     * @return the interpolated point.
     */
    virtual double interpolate() = 0;
};

#endif
