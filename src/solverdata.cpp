#include <solverdata.h>
#include <new>
#include <iostream>

SolverData::SolverData(unsigned int eqCount, unsigned int nStages, const Grid& grid): mGrid(grid){
  this->nStages = nStages;
  nEq = eqCount;

  unsigned int nx = mGrid.getSize();
  
  // Try to allocate memory for the arrays.
  try{
    data = new double*[nx];
    data_int = new double*[nx];
    for(int i = 0; i < nx; i++){
      data[i] = new double[eqCount];
      data_int[i] = new double[eqCount];
    }
    work = new double**[nStages];
    for(int i = 0; i < nStages; i++){
      work[i] = new double*[nx];
      for(int j = 0; j < nx; j++){
        work[i][j] = new double[eqCount];
      }
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for solver data.\n";
    nEq = 0;
    data = NULL;
    data_int = NULL;
    work = NULL;
  }
}

SolverData::SolverData(const SolverData& other): mGrid(other.getGrid()){

}

SolverData::~SolverData(){
  for(int i = 0; i < nStages; i++){
    for(int j = 0; j < mGrid.getSize(); j++){
      delete[] work[i][j];
    }
    delete[] work[i];
  }
  for(int i = 0; i < mGrid.getSize(); i++){
    delete[] data_int[i];
    delete[] data[i];
  }
  delete[] work;
  delete[] data_int;
  delete[] data;
}

bool SolverData::operator < (const SolverData& data) const{
  return mGrid < data.getGrid();
}
