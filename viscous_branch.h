#ifndef  VISCOUS_BRANCH_H
#define  VISCOUS_BRANCH_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

#include "Material.h"


using namespace Eigen;
using namespace std;

class viscous_branch: public Material
{
  public:
  MatrixXd C_i(3,3);
  MatrixXd be_tr(3,3);
  MatrixXd be(3,3);
  MatrixXd res(3,3);
  MatrixXd res_principal(3,1);
  MatrixXd eigs(3,1);
  MatrixXd be_col(6,1);
  MatrixXd dbe_col(6,1);
  MatrixXd res_col(6,1);
  MatrixXd Stiff(6,6);
  MatrixXd Stiff_principal(3,3);
  MatrixXd Tangent(6,6);
  MatrixXd Tangent_principal(3,3);
  MatrixXd tau(3,3);
  MatrixXd tau_principal(3,1);
  MatrixXd epsilon(3,1), epsilon_tr(3,1);
  MatrixXd mat_tan = MatrixXd::Zero(6,6);
  MatrixXd mat_tan_rotated(3,3);
  MatrixXd Rotation_mat(6,6);

  double lambda_A,lambda_B,lambda_C;
  Vector3d v0,v1,v2;
  float invar[2];
  float derivative[2];
  float second_derivative[2,2];
  float I1, I2, dW_dI1, dW_dI2, d2W_dI1I1, d2W_dI1I2, d2W_dI2I2;


  void set_inelastic_strain(MatrixXd _C_i);
  void compute_be_tr(MatrixXd F);
  void update_intervar_newton(MatrixXd be_tr);
  void compute_residual(MatrixXd be_trial, MatrixXd be, MatrixXd res);
  void compute_tangent(MatrixXd Stiff, MatrixXd be);
  void compute_tangent_nr(MatrixXd be);
  void compute_tangent_nr_principal(MatrixXd be);
  void compute_invar(MatrixXd C);
  void compute_derivative();
  void compute_stress_pk2(MatrixXd C,MatrixXd S);
  MatrixXd compute_stress_tau(MatrixXd be,MatrixXd tau);

  friend void convert_matrix_to_vector(MatrixXd Matrix, MatrixXd Vector);
  friend void convert_vector_to_matrix(MatrixXd Vector, MatrixXd Matrix);
};
#endif
