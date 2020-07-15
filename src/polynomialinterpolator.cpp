#include <polynomialinterpolator.h>

// PolynomialInterpolator {{{
PolynomialInterpolator::PolynomialInterpolator(const unsigned int p) : Interpolator(p) {
  // Allocate memory for the weights.
  weights = new double[p];

  // Build the weights.
  for(unsigned int n = 0; n < p; n++){
    weights[n] = calculateWeight(n - p/2 + 1);
  }
}
// }}}

// ~PolynomialInterpolator {{{
PolynomialInterpolator::~PolynomialInterpolator(){
  delete[] weights;
}
// }}}

// calculateWeight {{{
double PolynomialInterpolator::calculateWeight(int l){
  int p = (unsigned int)nStencil;
  double result = 1.0;
  // If we've got an odd number, the weight should be negative
  if((p/2 + l - 1) & 1){
    result = -1.0;
  }

  for(int i = 0; i < p; i++){
    result *= ((i - p/2) + 0.5);
  }

  // Calculate the denominator.
  double den = l - 0.5;
  for(int i = 1; i <= p/2 + l - 1; i++){
    den *= i;
  }
  for(int i = 1; i <= p/2 - l; i++){
    den *= i;
  }

  return result/den;
}
// }}}

// interpolate {{{
double PolynomialInterpolator::interpolate(){
  double result = 0.0;
  for(int i = 0; i < nStencil; i++){
    result += weights[i]*stencil[i];
  }
  return result;
}
// }}}
