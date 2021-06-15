#include <odedata.h>
#include <new>
#include <iostream>

// ODEData {{{
ODEData::ODEData(unsigned int eqCount, const Grid& grid, unsigned int lines) : mGrid(grid){
  nEq = eqCount;
  this->lines = lines;

  unsigned int nx = mGrid.getSize();

  // Try to allocate memory for the array.
  if(lines > 0){
    line = new double**[lines];
  }
  try{
    data = new double*[eqCount];
    for(unsigned int i = 0; i < lines; i++){
      line[i] = new double*[eqCount];
      for(unsigned int m = 0; m < eqCount; m++){
        line[i][m] = new double[nx];
      }
    }
    for(unsigned int i = 0; i < eqCount; i++){
      data[i] = new double[nx];
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for solver data.\n";
    nEq = 0;
    data = NULL;
    line = NULL;
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
  for(unsigned int i = 0; i < lines; i++){
    for(unsigned int m = 0; m < nEq; m++){
      delete[] line[i][m];
    }
    delete[] line[i];
  }
  delete[] line;
  delete[] data;
}
// }}}

// {{{
bool ODEData::operator < (const ODEData& data) const{
  return mGrid < data.getGrid();
}
// }}}
