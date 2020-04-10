#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <types.h>

/**********************************************************************************
 *
 * Class: Reconstruction
 * Author: Jacob Fields
 * Date Created: 5 Mar 2020
 * 
 * Description: This class describes the behavior for a general reconstruction
 *              routine for use in a Riemann solver.
 *
 *********************************************************************************/

class Reconstruction{
  public:
    /**
     * The constructor for the reconstruction routine.
     */
    Reconstruction();
    /**
     * The destructor for the reconstruction routine.
     */
    virtual ~Reconstruction();
    
    /**
     * Calculate the reconstructed right and left states ur and ul for a Riemann
     * problem based on an input set u. As a note on conventions, left and right
     * is determined based on the position relative to the interface at u[i-1/2].
     * Therefore, ur[i] is on the right of this interface (placing it in u[i], but
     * on the left side), and ul[i] is to the left of this interface (on the right
     * of u[i-1]).
     * @param n The size of the input state.
     * @param u The input state for the reconstruction.
     * @param ul The reconstructed left state.
     * @param ur The reconstructed right state.
     * @returns A Result enum specifying success or failure.
     */
    virtual Result reconstruct(const int n, const double* const RESTRICT u,
                               double* const RESTRICT ul, double* const RESTRICT ur)=0;
};

#endif
