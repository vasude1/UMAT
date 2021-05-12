// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
// #include <cblas.h>

// OWN INCLUDES
#include "Material_ps.h"
//#include "ArithmeticSequence.h"

// Namespaces
using namespace Eigen;
using namespace std;


void Material::set_material(double _stiff_ratio, double _relaxation_time){
    // c = {1044000.0,0.0,-22730,0.0,0.0,336.0,0.0,0.0,0.0};
    stiff_ratio = _stiff_ratio;
    relaxation_time = _relaxation_time;
    // cblas_dscal ( 10,  stiff_ratio, c,c );
    for(int i=0;i<9;++i){
      c[i] *= stiff_ratio;
      // std::cout << c[i] << '\n';
    }

};
