#ifndef  MATERIAL_H
#define  MATERIAL_H

// INCLUDES
#include <iostream>
#include <cmath>

class Material
{ public:
  // double c[9] = {1044000.0,0.0,-22730,0.0,0.0,336.0,124.0,-2.47,0.0};
  double c[9] = {1000000.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double relaxation_time;
  double stiff_ratio;

  void set_material(double _stiff_ratio, double _relaxation_time);

};
#endif
