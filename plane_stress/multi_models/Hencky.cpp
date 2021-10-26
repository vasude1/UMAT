#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Hencky.h"
using namespace Eigen;
using namespace std;



void Hencky::set_material(double _stiff_ratio, double _relaxation_time){
    stiff_ratio = _stiff_ratio;
    relaxation_time = _relaxation_time;
    c *= stiff_ratio;
};

double Hencky::coeff(){
return c;
}

double Hencky::rel_time(){
return relaxation_time;
}

Vector2d Hencky::compute_derivative(Vector2d& lambda){
    epsilon = 0.5*log(lambda.array());
    dW_d(0)=2.0*c/lambda(0)*(2.0*epsilon(0)+epsilon(1));
    dW_d(1)=2.0*c/lambda(1)*(epsilon(0)+2.0*epsilon(1));
    return dW_d;
};

double Hencky::compute_strain_energy(){
    psi=0.0;
    psi += 4*c*(sq(epsilon(0))+sq(epsilon(1))+epsilon(0)*epsilon(1));
    return psi;
};

Matrix2d& Hencky::compute_stiff_principal(Vector2d& lambda){
    Stiff_principal(0,0) = 8*c;
    Stiff_principal(0,1) = 4*c;
    Stiff_principal(1,1) = 8*c;
    Stiff_principal(1,0) = 4*c;
    return Stiff_principal;
};
