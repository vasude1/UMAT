#ifndef  POLYNOMIAL_H
#define  POLYNOMIAL_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class


using namespace Eigen;
using namespace std;

class Polynomial
{
  double c[9] = {1044000.0,0.0,-22730.0,0.0,0.0,336.0,124.0,0.0,0.0};
  double relaxation_time,stiff_ratio,c33;
  Vector2d dW_d;
  double dW_dI[2],invar[2],dI1_d[2],dI2_d[2];
  double d2W_dI1I1,d2W_dI1I2,d2W_dI2I1,d2W_dI2I2;
  double psi;
  double d2AI1_dAdA,d2AI1_dAdB,d2BI1_dBdB,d2BI1_dBdA;
  double d2AI2_dAdA,d2AI2_dAdB,d2BI2_dBdB,d2BI2_dBdA;
  Matrix2d Stiff_principal;

public:
  void set_material(double, double);
  double coeff();
  double rel_time();
  void compute_invar(Vector2d&);
  Vector2d compute_derivative(Vector2d&);
  double compute_strain_energy();
  Matrix2d& compute_stiff_principal(Vector2d&);
  template <typename number>
  number sq( number n ) {
    return (n*n);
   }

  template <typename number>
  number cu( number n ) {
    return (n*n*n);
   }
};

#endif
