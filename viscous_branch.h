#ifndef  VISCOUS_BRANCH_H
#define  VISCOUS_BRANCH_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class


using namespace Eigen;
using namespace std;

#include "Material.h"

class viscous_branch: public Material
{
  public:
  Matrix3d C_i,be_tr,be;
  Vector3d res_principal,eigs,eigs_tr;

  Matrix3d Stiff_principal,Tangent_principal;

  Matrix3d tau;
  Vector3d tau_principal,epsilon, epsilon_tr;
  MatrixXd mat_tan = MatrixXd::Zero(6,6);
  Matrix3d mat_tan_rotated;
  MatrixXd Rotation_mat = mat_tan;

  double lambda_A,lambda_B,lambda_C,eta;
  Vector3d v0,v1,v2;
  double invar[2];
  double derivative[2];
  double second_derivative[2][2];
  double I1, I2, dW_dI1, dW_dI2, d2W_dI1I1, d2W_dI1I2, d2W_dI2I2;
  double delta_t;


  void set_inelastic_strain(MatrixXd);
  void compute_be_tr(MatrixXd&);
  void compute_invar(MatrixXd, bool);
  void compute_derivative();
  void compute_second_derivative();
  void compute_stress_tau_principal();
  MatrixXd compute_stress_tau(MatrixXd);
  double compute_strain_energy(MatrixXd);
  double compute_dissipation();
  void compute_tangent_principal();
  void compute_tangent_nr_principal();
  void compute_residual_principal();
  void mat_tan_principal();
  MatrixXd rotate_mat_tan();
  MatrixXd update_intervar_newton_principal(MatrixXd&);

  friend void convert_matrix_to_vector(MatrixXd&, MatrixXd&);
  friend void convert_vector_to_matrix(MatrixXd&, MatrixXd&);
};
#endif
