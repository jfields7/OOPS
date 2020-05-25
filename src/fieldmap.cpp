#include <fieldmap.h>
#include <odedata.h>
#include <solverdata.h>
#include <grid.h>
#include <stdexcept>
#include <typeinfo>
#include <typeindex>
#include <memory>
#include <ode.h>

// Constructor {{{
FieldMap::FieldMap(const Grid& grid) : mGrid(grid){
  
}

FieldMap::FieldMap(const Grid& grid, std::map<std::string, FieldInfo>& fields) : mGrid(grid){
  for(auto i : fields){
    FieldInfo info = i.second;
    addField(info.name, info.nEqs, info.nStages);
  }
}
// }}}

// Destructor {{{
FieldMap::~FieldMap(){
  fields.clear();
}
// }}}

// addField {{{
void FieldMap::addField(std::string name, unsigned int nEqs, unsigned int stages=0){
  // FIXME: This should really throw an exception if the field already exists.
  if(hasField(name)){
    return;
  }

  if(stages > 0){
    //fields[name] = SolverData(nEqs, stages, mGrid);
    //fields.emplace(name, std::shared_ptr<SolverData>(new SolverData(nEqs, stages, mGrid)));
    solverFields.emplace(name, std::make_shared<SolverData>(nEqs, stages, mGrid));
    fields.insert({name, solverFields.at(name)});
  }
  else{
    //fields.emplace(name, nEqs, mGrid);
    //fields.emplace(name, std::shared_ptr<ODEData>(new ODEData(nEqs, mGrid)));
    fields.emplace(name, std::make_shared<ODEData>(nEqs, mGrid));
  }
}
// }}}

// removeField {{{
void FieldMap::removeField(std::string name){
  if(!hasField(name)){
    throw std::out_of_range("Cannot find field.");
  }
  else{
    fields.erase(name);
  }
}
// }}}

// hasField {{{
bool FieldMap::hasField(std::string name){
  return (fields.find(name) != fields.end());
}
// }}}

// isSolverField {{{
bool FieldMap::isSolverField(std::string name){
  // Check that the field actually exists.
  if(!hasField(name)){
    throw std::out_of_range("Cannot find field.");
    return false;
  }

  return (std::type_index(typeid(fields.at(name))) == std::type_index(typeid(SolverData)));
}
// }}}

// operator < {{{
bool FieldMap::operator< (const FieldMap& other){
  return mGrid < other.getGrid();
}
// }}}
