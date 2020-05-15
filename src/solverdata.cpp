#include <solverdata.h>
#include <new>
#include <iostream>

// {{{
SolverData::SolverData(unsigned int eqCount, unsigned int nStages, const Grid& grid): ODEData(eqCount, grid){
  this->nStages = nStages;
  unsigned int nx = mGrid.getSize();
  // Try to allocate memory for the arrays.
  try{
    data_int = new double*[eqCount];
    for(int i = 0; i < eqCount; i++){
      data_int[i] = new double[nx];
    }
    work = new double**[nStages];
    for(int i = 0; i < nStages; i++){
      work[i] = new double*[eqCount];
      for(int j = 0; j < eqCount; j++){
        work[i][j] = new double[nx];
      }
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for solver data.\n";
    nEq = 0;
    data_int = NULL;
    work = NULL;
  }
}
// }}}

// {{{
SolverData::SolverData(const SolverData& other): ODEData(other){
  std::cout << "SolverData copy constructor: This shouldn't be getting called, but it is.\n";
}
// }}}

// {{{
SolverData::~SolverData(){
  for(int i = 0; i < nStages; i++){
    for(int j = 0; j < nEq; j++){
      delete[] work[i][j];
    }
    delete[] work[i];
  }
  for(int i = 0; i < nEq; i++){
    delete[] data_int[i];
  }
  delete[] work;
  delete[] data_int;
}
// }}}

// {{{
bool SolverData::operator < (const SolverData& data) const{
  return mGrid < data.getGrid();
}
// }}}
