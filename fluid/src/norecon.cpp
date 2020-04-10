#include <norecon.h>

// NoRecon {{{
NoRecon::NoRecon(){

}
// }}}

// ~NoRecon {{{
NoRecon::~NoRecon(){

}
// }}}

// reconstruct {{{
Result NoRecon::reconstruct(const int n, const double* const RESTRICT u,
                            double* const RESTRICT ul, double* const RESTRICT ur){
  ul[0] = u[0];

  for(int i = 0; i < n-1; i++){
    ul[i+1] = u[i];
    ur[i] = u[i];
  }

  ur[n-1] = u[n-1];
}
// }}}
