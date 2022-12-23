#ifndef  DASHPOT_H
#define  DASHPOT_H
#define EIGEN_NO_DEBUG
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

using namespace Eigen;
using namespace std;

class dashpot
{
  Matrix2d C_i,be_tr,be,be_inverse;
  Vector2d res_principal,eigs,eigs_tr;
  Vector2d depsilon;
  Matrix2d Stiff_principal,Tangent_principal,C_alg,pressure_tangent,tot_tangent;
  Matrix2d tau;
  Vector2d tau_principal,epsilon, epsilon_tr,tau_dev;
  Matrix3d mat_tan,Rotation_mat,Rotation_mat_transp;
  double eta,p;
  Vector2d v0,v1,v0_,v1_,dW_d;
  double gamma_0,m,tau_v,tau_hat,rt_to=sqrt(2.0);
  double eps_max,delta_t,to_rt_to=2.0*sqrt(2.0);
  int count;
  EigenSolver<Matrix2d> es;
public:
  dashpot();
  template <typename model>
  void set_viscous_para(double,model*);
  void set_inelastic_strain(Matrix2d);
  void compute_be_tr(Matrix2d);
  void compute_stress_tau_principal(Vector2d&);
  Matrix2d compute_stress_tau();
  double compute_dissipation();
  void compute_tangent_pressure();
  template <typename model>
  void compute_tangent_nr_principal(model*);
  template <typename model>
  void compute_residual_principal(model*);
  void mat_tan_principal();
  Matrix3d rotate_mat_tan();
  template <typename model>
  Matrix2d update_intervar_newton_principal(Matrix2d,model*);
};

#endif
