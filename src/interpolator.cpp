#include <interpolator.h>

Interpolator::Interpolator(const unsigned int n) : nStencil(n) {
  stencil = new double[n];
}

Interpolator::~Interpolator(){
  delete[] stencil;
}
