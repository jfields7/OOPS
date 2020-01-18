#include <cubicinterpolator.h>

// CubicInterpolator {{{
CubicInterpolator::CubicInterpolator() : Interpolator(4) {
}
// }}}

// ~CubicInterpolator {{{
CubicInterpolator::~CubicInterpolator() {
}
// }}}

// interpolate {{{
double CubicInterpolator::interpolate(){
  return (9.0*(stencil[1] + stencil[2]) - (stencil[0] + stencil[3]))/16.0;
}
// }}}
