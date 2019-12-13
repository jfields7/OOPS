#include <solverdata.h>
#include <new>

SolverData::SolverData(unsigned int eqCount, unsigned int nStages, Grid& grid){
  this->grid = grid;
  this->nStages = nStages;
  nEq = eqCount;

  unsigned int nx = grid.getSize();
  
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

SolverData::SolverData(const SolverData& other){

}

SolverData::~SolverData(){
  for(int i = 0; i < nStages; i++){
    for(int j = 0; j < grid.getSize(); i++){
      delete[] work[i][j];
    }
    delete[] work[i];
  }
  for(int i = 0; i < grid.getSize(); i++){
    delete[] data_int[i];
    delete[] data[i];
  }
  delete[] work;
  delete[] data_int;
  delete[] data;
}

bool SolverData::operator < (const SolverData& data) const{
  return grid < data.getGrid();
}
