#ifndef  VISCOUS_BRANCH_H
#define  VISCOUS_BRANCH_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class


using namespace Eigen;
using namespace std;

#include "Material_ps.h"

typedef Matrix<double, 6, 6> Matrix6d;

class viscous_branch: public Material
{
  public:
  Matrix2d C_i,be_tr,be,be_inverse;
  Vector2d res_principal,eigs,eigs_tr;
  Vector2d depsilon;

  Matrix2d Stiff_principal,Tangent_principal,C_alg;

  Matrix2d tau;
  Vector2d tau_principal,epsilon, epsilon_tr;
  Matrix3d mat_tan;
  Matrix3d mat_tan_vol;
  Matrix3d Rotation_mat;
  Matrix3d Rotation_mat_transp;

  double lambda_A,lambda_B,lambda_C,eta;
  Vector2d v0,v1,v2,v0_,v1_,v2_;
  double invar[2];
  double derivative[2];
  double second_derivative[2][2];
  double I1, I2, dW_dI1, dW_dI2, d2W_dI1I1, d2W_dI1I2, d2W_dI2I2;
  double delta_t;

  viscous_branch();
  void set_inelastic_strain(MatrixXd);
  void compute_be_tr(MatrixXd);
  void compute_invar(MatrixXd, bool);
  void compute_derivative();
  void compute_second_derivative();
  void compute_stress_tau_principal();
  void compute_stress_tau_principal_hencky();
  MatrixXd compute_stress_tau();
  double compute_strain_energy(MatrixXd);
  double compute_strain_energy_hencky(MatrixXd);
  double compute_dissipation();
  void compute_tangent_principal();
  void compute_tangent_principal_hencky();
  void compute_tangent_nr_principal();
  void compute_residual_principal();
  void compute_residual_principal_hencky();
  void mat_tan_principal(bool);
  MatrixXd rotate_mat_tan(bool);
  MatrixXd mat_tan_volu();
  MatrixXd update_intervar_newton_principal(MatrixXd);
  MatrixXd update_intervar_newton_principal_hencky(MatrixXd);


  friend void convert_matrix_to_vector(MatrixXd&, MatrixXd&);
  friend void convert_vector_to_matrix(MatrixXd&, MatrixXd&);
};

template <typename number>
number sq( number n ) {
  return (n*n);
 }

template <typename number>
number cu( number n ) {
  return (n*n*n);
 }
 
#endif
