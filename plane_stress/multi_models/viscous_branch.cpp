// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

#include "viscous_branch.h"
#include "viscous_branch.t.hpp"

viscous_branch::viscous_branch(){
mat_tan.setZero();
}

void viscous_branch::set_inelastic_strain(MatrixXd _C_i){
  C_i = _C_i;
};

void viscous_branch::compute_be_tr(MatrixXd& F){
  be_tr = F*(C_i.inverse()*F.transpose());
  EigenSolver<Matrix2d> es(be_tr);
  v0_ = es.eigenvectors().col(0).real()(seq(0,1),0);
  v1_ = es.eigenvectors().col(1).real()(seq(0,1),0);

  v0 = es.eigenvectors().row(0).real()(0,seq(0,1));
  v1 = es.eigenvectors().row(1).real()(0,seq(0,1));

  eigs_tr(0) = es.eigenvalues()(0,0).real();
  eigs_tr(1) = es.eigenvalues()(1,0).real();

  Rotation_mat(0,0) = v0(0)*v0(0);
  Rotation_mat(1,0) = v1(0)*v1(0);
  Rotation_mat(2,0) = v0(0)*v1(0);

  Rotation_mat(0,1) = v0(1)*v0(1);
  Rotation_mat(1,1) = v1(1)*v1(1);
  Rotation_mat(2,1) = v0(1)*v1(1);

  Rotation_mat(0,2) = 2.0*v0(0)*v0(1);
  Rotation_mat(1,2) = 2.0*v1(0)*v1(1);
  Rotation_mat(2,2) = v0(0)*v1(1)+v0(1)*v1(0);

  Rotation_mat_transp = Rotation_mat.transpose();
};

void viscous_branch::compute_stress_tau_principal(Vector2d dW_d){
  tau_principal(0) = 2.0*dW_d(0)*eigs(0);
  tau_principal(1) = 2.0*dW_d(1)*eigs(1);
  p = -(tau_principal(0)+tau_principal(1))/3.0;
};

MatrixXd viscous_branch::compute_stress_tau(){
  tau = tau_principal(0)*v0_*v0_.transpose()+tau_principal(1)*v1_*v1_.transpose();
  return tau;
};

double viscous_branch::compute_dissipation(){
  double psi = 0.0;
  psi = tau_dev.transpose()*tau_dev+p*p;
  psi /= 2.0*eta;
  return psi;
};

void viscous_branch::compute_tangent_pressure(){
  pressure_tangent(0,0) = -1.0/3.0*(Stiff_principal(0,0)+Stiff_principal(1,0));
  pressure_tangent(0,1) = -1.0/3.0*(Stiff_principal(0,1)+Stiff_principal(1,1));
  pressure_tangent(1,0) = pressure_tangent(0,0);
  pressure_tangent(1,1) = pressure_tangent(0,1);
};

void viscous_branch::mat_tan_principal(){
    C_alg =  Stiff_principal*(Tangent_principal.inverse());
    mat_tan(seq(0,1),seq(0,1)).array() = C_alg - 2.0*Matrix2d(tau_principal.asDiagonal());
    if (eigs_tr(0)==eigs_tr(1)) {
      mat_tan(2,2) = 0.5*(C_alg(0,0)-C_alg(1,0))-tau_principal(0);
    }
    else {
      mat_tan(2,2) = (tau_principal(0)*eigs_tr(1) - tau_principal(1)*eigs_tr(0))/(eigs_tr(0) - eigs_tr(1));
    }
};

MatrixXd viscous_branch::rotate_mat_tan(){
  mat_tan_principal();
  mat_tan = Rotation_mat*(mat_tan*Rotation_mat_transp);
  return mat_tan;
};
