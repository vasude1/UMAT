// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
// #include <cblas.h>

// OWN INCLUDES
#include "Material.h"

// Namespaces
using namespace Eigen;
using namespace std;


void Material::set_material(double _stiff_ratio, double _relaxation_time){
    stiff_ratio = _stiff_ratio;
    relaxation_time = _relaxation_time;
    // cblas_dscal ( 10,  stiff_ratio, c,c );
    for(int i=0;i<9;++i){
      c[i] *= stiff_ratio;
      // std::cout << c[i] << '\n';
    }

};
