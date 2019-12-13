#ifndef CUBIC_H
#define CUBIC_H

namespace interp{
  double cubicInterpCenter(double ym2, double ym1, double yp1, double yp2);

  double cubicInterpLeft(double ym1, double yp1, double yp2, double yp3);

  double cubicInterpRight(double ym3, double ym2, double ym1, double yp1);
}

#endif
