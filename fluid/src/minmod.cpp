#include <minmod.h>

// Minmod {{{
MinmodRecon::MinmodRecon(){

}
// }}}

// ~Minmod {{{
MinmodRecon::~MinmodRecon(){

}
// }}}

// reconstruct {{{
Result MinmodRecon::reconstruct(const int n, const double* const RESTRICT u,
                           double* const RESTRICT ul, double* const RESTRICT ur){
  ur[0] = u[0];
  ul[0] = u[0];

  for (int i = 1; i < n-1; i++){
    double df1 = (u[i] - u[i-1]);
    double df2 = (u[i+1] - u[i]);

    double slp = minmod(df1, df2);
    ul[i+1] = u[i] + 0.5 * slp;
    ur[i]   = u[i] - 0.5 * slp;
  }

  ul[n-1] = u[n-1];
  ur[n-1] = u[n-1];
}
// }}}
