#ifndef PARAMS_H
#define PARAMS_H

/**
 * A class representing parameters. ODEs and PDEs with special parameters for
 * their initial conditions should take a parameter object descended from this
 * object. ODEs and PDEs not requiring special parameters should just use this
 * Parameters object.
 * Regarding descendent Parameter objects, they should specify a unique id that
 * can be used to recognize that particular set of parameters. This can be used
 * as a safety mechanism by an ODE object to guarantee that it doesn't try to
 * cast a Parameters object that doesn't have all the necessary data.
 */
class Parameters{
  public:
    enum InterpolationScheme{
      CUBIC,
    };

  protected:
    /**
     * The id specifying what kind of parameter set this is.
     */
    const unsigned int pId;

    /**
     * The kind of interpolation we can use for grids.
     */
    InterpolationScheme interp;

    /**
     * Create a parameters object with the specified id. This is protected so
     * that only a descendent of the Parameters object can change the id.
     */
    Parameters(unsigned int id) : pId(id), interp(InterpolationScheme::CUBIC){};
  public:
    /**
     * Create a parameters object with a null id.
     */
    Parameters() : pId(0), interp(InterpolationScheme::CUBIC){};
    
    /**
     * Get the id for this particular paramters object.
     */
    inline unsigned int getId() const{
      return pId;
    }

    /**
     * Get the interpolation scheme we should be using.
     */
    inline InterpolationScheme getInterpolation() const{
      return interp;
    }

    /**
     * Set the interpolation scheme we should use.
     */
    inline void setInterpolation(InterpolationScheme scheme){
      interp = scheme;
    }
};

#endif
