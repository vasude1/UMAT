#ifndef  VISCOUS_BRANCH_H
#define  VISCOUS_BRANCH_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

using namespace Eigen;
using namespace std;

class viscous_branch
{
  Matrix2d C_i,be_tr,be,be_inverse;
  Vector2d res_principal,eigs,eigs_tr;
  Vector2d depsilon;
  Matrix2d Stiff_principal,Tangent_principal,C_alg,pressure_tangent;
  Matrix2d tau;
  Vector2d tau_principal,epsilon, epsilon_tr,tau_dev;
  Matrix3d mat_tan,Rotation_mat,Rotation_mat_transp;
  double eta,p;
  Vector2d v0,v1,v0_,v1_,dW_d;
  double gamma_0;
  double delta_t;
public:
  viscous_branch();
  template <typename model>
  void set_viscous_para(double,model*);
  void set_inelastic_strain(MatrixXd);
  void compute_be_tr(MatrixXd&);
  void compute_stress_tau_principal(Vector2d);
  MatrixXd compute_stress_tau();
  double compute_dissipation();
  void compute_tangent_pressure();
  template <typename model>
  void compute_tangent_nr_principal(model*);
  template <typename model>
  void compute_residual_principal(model*);
  void mat_tan_principal();
  MatrixXd rotate_mat_tan();
  template <typename model>
  MatrixXd update_intervar_newton_principal(MatrixXd&,model*);
};

#endif
