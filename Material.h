#ifndef  MATERIAL_H
#define  MATERIAL_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class
#include <cblas>

#include "UELMAT_ShapeFunctions.h"

using namespace Eigen;

class Material
{
  double c[9] = {1044000.0,0.0,-22730,0.0,0.0,336.0,124.0,-2.47,0.0};
  double relaxation_time;
  double stiff_ratio;

 public:
  void set_material(double _stiff_ratio, double _relaxation_time);

};
#endif
