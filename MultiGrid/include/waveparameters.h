#ifndef WAVE_PARAMETERS_H
#define WAVE_PARAMETERS_H

#include <parameters.h>

class WaveParameters : public Parameters {
  public:
    enum InitialConditions{
      GAUSSIAN,
    };

    WaveParameters() : Parameters(1){
      ics = GAUSSIAN;
      koSigma = 0.1;
    }

    inline void setInitialConditions(InitialConditions i){
      ics = i;
    }

    inline InitialConditions getInitialConditions(){
      return ics;
    }

    inline void setKOSigma(double sigma){
      koSigma = sigma;
    }

    inline double getKOSigma(){
      return koSigma;
    }
  private:
    InitialConditions ics;

    double koSigma;
};

#endif
