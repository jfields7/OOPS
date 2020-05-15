#ifndef PARAM_PARSER_H
#define PARAM_PARSER_H
#include "paramreader.h"
#include "parameters.h"
#include <string.h>

/**
 * A class that will parse a parameter file. It can function in two ways: it can either generate
 * a new Parameters object from a file, or it can take an existing Parameters object and update its values
 * from a file. It has a unique id that much match the id for the Parameters objects it modifies or creates.
 */
class ParamParser{
  protected:
    /**
     * The id specifying what kind of parameters this parses.
     */
    const unsigned int pId;
    /**
     * The ParamReader this object uses.
     */
    ParamReader reader;

  public:
    /**
     * The constructor for a ParamParser.
     */
    ParamParser(unsigned int id) : pId(id){
      reader = ParamReader();
    }

    /**
     * Update an existing Parameters object from a file.
     */
    virtual void updateParameters(std::string fname, Parameters *params)=0;

    /**
     * Check if the Parameters object matches the id.
     */
    inline bool checkId(Parameters* params) const{
      return params->getId() == pId;
    }

    /**
     * Get this parser's id.
     */
    inline unsigned int getId() const{
      return pId;
    }

};

#endif
