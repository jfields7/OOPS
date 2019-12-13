#ifndef OUTPUT_H
#define OUTPUT_H
#include "sdf.h"

void output_data(char* name, double v[], double r[], int size, double time){
  gft_out_full(name,time,&size,"r",1,r,v);
}

#endif
