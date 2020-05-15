#include <odedata.h>
#include <new>
#include <iostream>

// ODEData {{{
ODEData::ODEData(unsigned int eqCount, const Grid& grid) : mGrid(grid){
  nEq = eqCount;

  unsigned int nx = mGrid.getSize();

  // Try to allocate memory for the array.
  try{
    data = new double*[eqCount];
    for(unsigned int i = 0; i < eqCount; i++){
      data[i] = new double[nx];
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for solver data.\n";
    nEq = 0;
    data = NULL;
  }
}
// }}}

// Copy constructor {{{
ODEData::ODEData(const ODEData& other): mGrid(other.getGrid()){
  std::cout << "ODEData copy constructor: This shouldn't be getting called, but it is.\n";
}
// }}}

// ~ODEData {{{
ODEData::~ODEData(){
  for(unsigned int i = 0; i < nEq; i++){
    delete[] data[i];
  }
  delete[] data;
}
// }}}

// {{{
bool ODEData::operator < (const ODEData& data) const{
  return mGrid < data.getGrid();
}
// }}}
