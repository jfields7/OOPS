#ifndef WAVE_PARSER_H
#define WAVE_PARSER_H

#include <paramparser.h>
#include <waveparameters.h>

class WaveParser : ParamParser{
  public:
    WaveParser() : ParamParser(1){}

    virtual ~WaveParser(){}

    virtual Parameters& getParameters(std::string fname){
      WaveParameters params = WaveParameters();
      updateParameters(fname, &params);
      return params;
    }

    virtual void updateParameters(std::string fname, Parameters* params);
};

#endif
