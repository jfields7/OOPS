#ifndef OPERATORS_H
#define OPERATORS_H

/**************************************************************************************
 * 
 * File name: operators.h
 * Author: Jacob Fields
 * Date Modified: 11-8-2019
 *
 * Description: This defines some simple mathematical operators we'll need, mostly
 *  derivatives.
 *
 *************************************************************************************/

namespace operators{
  // Derivative Operators {{{
  /* A fourth-order centered difference first derivative operator.
   * @param u:  The five-point stencil to use to calculate the derivative. The value
   *            u[2] is assumed to be the middle.
   * @param dx: The spacing between points. This is assumed to be uniform.
   * @return du/dx[2] as a double.
   */
  inline double dx_4(const double u[5], const double dx){
    return (u[0] - 8.0*u[1] + 8.0*u[3] - u[4])/(12.0*dx);
  }

  /* A fourth-order centered difference second derivative operator.
   * @param u:  The five-point stencil to use to calculate the derivative. The value
   *            u[2] is assumed to be the middle.
   * @param dx: The spacing between points. This is assumed to be uniform.
   * @return d^2u/dx^2[2] as a double.
   */
  inline double dxx_4(const double u[5], const double dx){
    return (- u[0] + 16.0*u[1] - 30.0*u[2] + 16.0*u[3] - u[4])/(12.0*dx*dx);
  }

  /* A second-order centered difference first derivative operator.
   * @param u:  The three-point stencil to use to calculate the derivative. The value
   *            u[1] is assumed to be the middle.
   * @param dx: The spacing between points. This is assumed to be uniform.
   * @return du/dx[1] as a double.
   */
  inline double dx_2(const double u[3], const double dx){
    return (-u[0] + u[2])/(2.0*dx);
  }

  /* A second-order forward difference first derivative operator.
   * @param u:  The three-point stencil use to calculate the derivative. The value
   *            u[0] is assumed to be the leftmost point.
   * @param dx: The spacing between points. This is assumed to be uniform.
   * @return du/dx[0] as a double.
   */
  inline double dx_2off(const double u[3], const double dx){
    return (-3.0*u[0] + 4.0*u[1] - u[2])/(2.0*dx);
  }

  /* A second-order centered difference second derivative operator.
   * @param u:  The three-point stencil to use to calculate the derivative. The value
   *            u[1] is assumed to be the middle.
   * @param dx: The spacing between points. This is assumed to be uniform.
   * @return d^2u/dx^2[1] as a double.
   */
  inline double dxx_2(const double u[3], const double dx){
    return (u[0] - 2.0*u[1] + u[2])/(dx*dx);
  }

  /* A second-order forward difference second derivative operator.
   * @param u:  The four-point stencil to use to calculate the derivative. The value
   *            u[0] is assumed to be the leftmost point.
   * @param dx: The spacing between points. This is assumed to be uniform.
   * @return d^2u/dx^2[0] as a double.
   */
  inline double dxx_2off(const double u[4], const double dx){
    return (2.0*u[0] - 5.0*u[1] + 4.0*u[2] - u[3])/(dx*dx);
  }
  // }}}

  // Kreiss-Oliger dissipation {{{
  /* A fourth-order centered difference fourth derivative operator for Kreiss-Oliger
   * dissipation.
   * @param u:  The seven-point stencil containing the points to use to calculate the
   *            derivative operator.
   * @param dx: The grid spacing to use, which is assumed to be uniform in the stencil.
   * @return d^4u[3]/dx^4 as a double.
   */
  inline double ko_dx(const double u[7], const double dx){
    return - (-u[0] + 6.0*u[1] - 15.0*u[2] + 20.0*u[3] - 15.0*u[4] + 6.0*u[5] - u[6])/(64.0*dx);
  }

  inline double ko_dx_off1(const double u[4], const double dx){
    return (-u[0] + 3.0*u[1] - 3.0*u[2] + u[3])*48.0/(59.0*64.0*dx);
  }

  inline double ko_dx_off2(const double u[5], const double dx){
    return (3.0*u[0] - 10.0*u[1] + 12.0*u[2] - 6.0*u[3] + u[4])*48.0/(43.0*64.0*dx);
  }

  inline double ko_dx_off3(const double u[6], const double dx){
    return (-3.0*u[0] + 12.0*u[1] - 19.0*u[2] + 15.0*u[3] - 6.0*u[4] + u[5])*48.0/(49.0*64.0*dx);
  }

  /**
   * A lower-order (2nd?) centered different operator for Kreiss-Oliger dissipation.
   * @param u: The five-point stencil containing the points to use to calculate the
   *           derivative operator.
   * @param dx: The grid spacing to use, which is assumed to be uniform in the stencil.
   * @return: The KO dissipation at u[2] as a double.
   */
  inline double ko_dx_2(const double u[5], const double dx){
    return (u[0] - 4.0*u[1] + 6.0*u[2] - 4.0*u[3] + u[4])*dx;
  }

    
  // }}}
};

#endif
