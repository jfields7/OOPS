#include <cubic.h>

double interp::cubicInterpCenter(double ym2, double ym1, double yp1, double yp2){
  return (9.0*(ym1 + yp1) - ym2 - yp2)/16.0;
}

double interp::cubicInterpLeft(double ym1, double yp1, double yp2, double yp3){
  return (5.0*(ym1 - yp2) + 15.0*yp1 + yp3)/16.0;
}

double interp::cubicInterpRight(double ym3, double ym2, double ym1, double yp1){
  return (5.0*(yp1 - ym2) + 15.0*ym1 + ym3)/16.0;
}
