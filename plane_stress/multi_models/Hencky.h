#ifndef  HENCKY_H
#define  HENCKY_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class


using namespace Eigen;
using namespace std;
// Strain energy form
// PSI = 4c(e_1^2+e_2^2+e_1e_2)
// e_i s are the logarithmic strains
class Hencky
{
  double c = 1E7;
  double relaxation_time,stiff_ratio;
  Vector2d dW_d,epsilon;
  double psi;
  Matrix2d Stiff_principal;

public:
  void set_material(double, double);
  double coeff();
  double rel_time();
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
