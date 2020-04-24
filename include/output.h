#ifndef OUTPUT_H
#define OUTPUT_H
#include <oopsconfig.h>
#include <cstdio>

#ifdef USE_SDF
#include "sdf.h"
#endif

void output_data(char* name, double v[], double r[], int size, double time){
  #ifdef USE_SDF
  gft_out_full(name,time,&size,"r",1,r,v);
  #else
  static bool first = true;
  if(first){
    printf("[0;34mSDF output currently disabled. Further warnings will be suppressed.\n[0m");
  }
  #endif
}

#endif
