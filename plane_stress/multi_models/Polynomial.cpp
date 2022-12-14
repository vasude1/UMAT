#define EIGEN_NO_DEBUG

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Polynomial.h"
using namespace Eigen;
using namespace std;



void Polynomial::set_material(double _stiff_ratio, double _relaxation_time){
  // std::cout << "Here 1" << '\n';
    stiff_ratio = _stiff_ratio;
    // std::cout << "Here 2" << '\n';
    relaxation_time = _relaxation_time;
    // std::cout << "Here 222" << '\n';
    for(i=0;i<9;++i){
      // std::cout << i << '\n';
      c[i] *= stiff_ratio;
    }
    // std::cout << "Here 3" << '\n';
  };

double Polynomial::coeff(){
return c[0];
}

double Polynomial::rel_time(){
return relaxation_time;
}

void Polynomial::compute_invar(Vector2d C){
    c33 = 1.0/C(0)/C(1);
    invar[0] = C(0)+ C(1) + c33;
    invar[1] = C(0)*C(1)+C(1)*c33+c33*C(0);
    dI1_d[0] = 1.0-1.0/C(0)/C(0)/C(1);
    dI1_d[1] = 1.0-1.0/C(0)/C(1)/C(1);
    dI2_d[0] = C(1)-1.0/C(0)/C(0);
    dI2_d[1] = C(0)-1.0/C(1)/C(1);
};

Vector2d Polynomial::compute_derivative(Vector2d lambda){
    compute_invar(lambda);
    dW_dI[0]= c[0]+2.0*c[2]*(invar[0]-3)+3.0*c[5]*(invar[0]-3)*(invar[0]-3)+2.0*c[6]*(invar[0]-3)*(invar[1]-3);
    dW_dI[1]= c[6]*(invar[0]-3)*(invar[0]-3);
    dW_d(0)=dW_dI[0]*dI1_d[0]+dW_dI[1]*dI2_d[0];
    dW_d(1)=dW_dI[0]*dI1_d[1]+dW_dI[1]*dI2_d[1];
    return dW_d;
};

double Polynomial::compute_strain_energy(){
    psi=0.0;
    psi += c[0]*(invar[0]-3)+c[1]*(invar[1]-3)+c[2]*sq(invar[0]-3)+c[3]*(invar[0]-3)*(invar[1]-3);
    psi += c[4]*sq(invar[1]-3)+c[5]*cu(invar[0]-3)+c[6]*sq(invar[0]-3)*(invar[1]-3);
    psi += c[7]*(invar[0]-3)*sq(invar[1]-3)+c[8]*cu(invar[1]-3);
    return psi;
};

Matrix2d Polynomial::compute_stiff_principal(Vector2d lambda){
    d2W_dI1I1= 2*c[2] + 6*c[5]*(invar[0]-3)+2.0*c[6]*(invar[1]-3);
    d2W_dI1I2= 2.0*c[6]*(invar[0]-3);
    d2W_dI2I1= 2.0*c[6]*(invar[0]-3);
    d2W_dI2I2= 0.0;

    d2AI1_dAdA = 1.0+1.0/lambda(0)/lambda(0)/lambda(1);
    d2AI1_dAdB = 1.0/lambda(1)/lambda(1)/lambda(0);

    d2BI1_dBdB = 1.0+1.0/lambda(1)/lambda(1)/lambda(0);
    d2BI1_dBdA = 1.0/lambda(0)/lambda(0)/lambda(1);

    d2AI2_dAdA = lambda(1)+1.0/lambda(0)/lambda(0);
    d2AI2_dAdB = lambda(0);

    d2BI2_dBdB = lambda(0)+1.0/lambda(1)/lambda(1);
    d2BI2_dBdA = lambda(1);

    Stiff_principal(0,0) = 2.0*dW_dI[0]*d2AI1_dAdA+2.0*lambda(0)*dI1_d[0]*(d2W_dI1I1*dI1_d[0]+d2W_dI1I2*dI2_d[0]);
    Stiff_principal(0,0) +=2.0*dW_dI[1]*d2AI2_dAdA+2.0*lambda(0)*dI2_d[0]*(d2W_dI2I1*dI1_d[0]+d2W_dI2I2*dI2_d[0]);
    Stiff_principal(0,0) *= 2*lambda(0);
    Stiff_principal(0,1) = 2.0*dW_dI[0]*d2AI1_dAdB+2.0*lambda(0)*dI1_d[0]*(d2W_dI1I1*dI1_d[1]+d2W_dI1I2*dI2_d[1]);
    Stiff_principal(0,1) +=2.0*dW_dI[1]*d2AI2_dAdB+2.0*lambda(0)*dI2_d[0]*(d2W_dI2I1*dI1_d[1]+d2W_dI2I2*dI2_d[1]);
    Stiff_principal(0,1) *= 2*lambda(1);

    Stiff_principal(1,1) = 2.0*dW_dI[0]*d2BI1_dBdB+2.0*lambda(1)*dI1_d[1]*(d2W_dI1I1*dI1_d[1]+d2W_dI1I2*dI2_d[1]);
    Stiff_principal(1,1) +=2.0*dW_dI[1]*d2BI2_dBdB+2.0*lambda(1)*dI2_d[1]*(d2W_dI2I1*dI1_d[1]+d2W_dI2I2*dI2_d[1]);
    Stiff_principal(1,1) *= 2*lambda(1);
    Stiff_principal(1,0) = 2.0*dW_dI[0]*d2BI1_dBdA+2.0*lambda(1)*dI1_d[1]*(d2W_dI1I1*dI1_d[0]+d2W_dI1I2*dI2_d[0]);
    Stiff_principal(1,0) +=2.0*dW_dI[1]*d2BI2_dBdA+2.0*lambda(1)*dI2_d[1]*(d2W_dI1I2*dI1_d[0]+d2W_dI2I2*dI2_d[0]);
    Stiff_principal(1,0) *= 2*lambda(0);
    return Stiff_principal;
};
