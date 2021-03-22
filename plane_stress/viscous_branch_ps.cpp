// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "Material_ps.h"
#include "viscous_branch_ps.h"


// Namespaces
using namespace Eigen;
using namespace std;


void viscous_branch::set_inelastic_strain(MatrixXd _C_i){
  C_i = _C_i;
  eta = relaxation_time*c[0];
};

void viscous_branch::compute_be_tr(MatrixXd F){
  be_tr = F*(C_i.inverse()*F.transpose());
  EigenSolver<Matrix2d> es(be_tr);
  v0_ = es.eigenvectors().col(0).real()(seq(0,1),0);
  v1_ = es.eigenvectors().col(1).real()(seq(0,1),0);

  v0 = es.eigenvectors().row(0).real()(0,seq(0,1));
  v1 = es.eigenvectors().row(1).real()(0,seq(0,1));

  eigs_tr(0) = es.eigenvalues()(0,0).real();
  eigs_tr(1) = es.eigenvalues()(1,0).real();

  //
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

void viscous_branch::compute_invar(MatrixXd C, bool a){
  if(a){
    double c33 = 1.0/(C(1,1)*C(0,0) - C(1,0)*C(0,1));
    invar[0] = C(0,0) + C(1,1) + c33;
    invar[1] = -0.5*(C(0,0)*C(0,0)+C(1,1)*C(1,1)+c33*c33+2.0*C(1,0)*C(1,0) - invar[0]*invar[0]);
  }
  if(!a){
    double c33 = 1.0/C(0)/C(0);
    invar[0] = C(0)+ C(1) + c33;
    invar[1] = C(0)*C(1) + C(1)*c33 + C(0)*c33;
  }

};

void viscous_branch::compute_derivative(){
  I1 = invar[0];
  I2 = invar[1];
  derivative[0]= c[0];
  derivative[1]= c[1];
};

void viscous_branch::compute_second_derivative(){
  I1 = invar[0];
  I2 = invar[1];
  second_derivative[0][0]= 2*c[2] + 6*c[5]*(I1-3) + 2*c[6]*(I2-3);
  second_derivative[0][1]= c[3] + 2*c[6]*(I1-3)+ 2*c[7]*(I2-3);
  second_derivative[1][0]= second_derivative[0][1];
  second_derivative[1][1]= 2*c[4] + 6*c[8]*(I2-3)+ 2*c[7]*(I1-3);

  // // std::cout << second_derivative[0][0] << second_derivative[0][1] <<second_derivative[1][0]<<second_derivative[1][1] << '\n';

};

void viscous_branch::compute_stress_tau_principal(){
  lambda_A = eigs(0);
  lambda_B = eigs(1);

  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  tau_principal(0) = (dW_dI1*(1.0-1.0/lambda_A/lambda_A/lambda_B) + dW_dI2*(lambda_B-1.0/lambda_A/lambda_A))*lambda_A;
  tau_principal(1) = (dW_dI1*(1.0-1.0/lambda_A/lambda_B/lambda_B) + dW_dI2*(lambda_A-1.0/lambda_B/lambda_B))*lambda_B;

  tau_principal = 2*tau_principal;
};

void viscous_branch::compute_stress_tau_principal_hencky(){

  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  tau_principal(0) = 2*dW_dI1*(epsilon(0)+epsilon(1));
  tau_principal(1) = 2*dW_dI1*(epsilon(0)+epsilon(1));

};

MatrixXd viscous_branch::compute_stress_tau(){
  tau = tau_principal(0)*v0_*v0_.transpose()+tau_principal(1)*v1_*v1_.transpose();
  return tau;
};

double viscous_branch::compute_strain_energy(MatrixXd _be){
  double psi = 0.0;
  compute_invar(_be,1);
  psi += c[0]*(I1-3);
  return psi;
};

double viscous_branch::compute_strain_energy_hencky(MatrixXd _be){
  double psi = 0.0;
  psi += epsilon(0)*epsilon(0)+epsilon(1)*epsilon(1)+(epsilon(1)-epsilon(0))*(epsilon(1)-epsilon(0));
  return 0.5*dW_dI1*psi;
};

double viscous_branch::compute_dissipation(){
  double psi = 0.0;
  psi += tau_principal.transpose()*tau_principal;
  psi /= eta;
  return psi;
};

void viscous_branch::compute_tangent_principal(){
  compute_second_derivative();
  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  Stiff_principal(0,0) = dW_dI1*(1.0+1.0/lambda_A/lambda_A/lambda_B);
  Stiff_principal(0,0) *= 2*lambda_A;
  Stiff_principal(0,1) = 2*dW_dI1/lambda_A/lambda_B;

  Stiff_principal(1,0) = 2*dW_dI1/lambda_A/lambda_B;
  Stiff_principal(1,1) = dW_dI1*(1.0+1.0/lambda_A/lambda_B/lambda_B);;
  Stiff_principal(1,1) *= 2*lambda_B;

  Stiff_principal *= 2.0;

};


void viscous_branch::compute_tangent_principal_hencky(){
  compute_second_derivative();
  dW_dI1 = derivative[0];
  Stiff_principal(seq(0,1),seq(0,1)) = 2*dW_dI1*MatrixXd::Ones(2,2);

};

void viscous_branch::compute_tangent_nr_principal(){
  Tangent_principal = delta_t/2.0/eta*Stiff_principal+MatrixXd::Identity(2,2);
};



void viscous_branch::compute_residual_principal(){

  compute_invar(eigs,0);
  compute_derivative();
  compute_stress_tau_principal();

  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  res_principal = epsilon + delta_t/2.0/eta*tau_principal - epsilon_tr;
};


void viscous_branch::compute_residual_principal_hencky(){

  compute_invar(eigs,0);
  compute_derivative();
  compute_stress_tau_principal_hencky();

  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  res_principal = epsilon + delta_t/2.0/eta*tau_principal - epsilon_tr;
};

void viscous_branch::mat_tan_principal(bool elastic){
    //
    C_alg =  Stiff_principal*(Tangent_principal.inverse());
    mat_tan(seq(0,1),seq(0,1)).array() = C_alg - 2.0*Matrix2d(tau_principal.asDiagonal());

    if (eigs_tr(0)==eigs_tr(1)) {
      mat_tan(2,2) = 0.5*(C_alg(0,0)-C_alg(1,0))-tau_principal(0);
    }
    else {
      mat_tan(2,2) = (tau_principal(0)*eigs_tr(1) - tau_principal(1)*eigs_tr(0))/(eigs_tr(0) - eigs_tr(1));
    }

};


MatrixXd viscous_branch::rotate_mat_tan(bool elastic){
  mat_tan_principal(elastic);
  mat_tan = Rotation_mat*(mat_tan*Rotation_mat_transp);
  return mat_tan;

};

MatrixXd viscous_branch::update_intervar_newton_principal(MatrixXd F){
  eigs = eigs_tr;
  epsilon_tr = 0.5*log(eigs_tr.array());
  epsilon = epsilon_tr;
    do {

      compute_residual_principal();
      compute_tangent_principal();
      compute_tangent_nr_principal();
      depsilon = -1.0*Tangent_principal.householderQr().solve(res_principal);
      epsilon += depsilon;
      eigs = exp(2.0*epsilon.array());
    }
    while (res_principal.norm() > 1E-4);
  be = eigs(0)*v0_*v0_.transpose()+eigs(1)*v1_*v1_.transpose();
  be_inverse = 1.0/eigs(0)*v0_*v0_.transpose()+1.0/eigs(1)*v1_*v1_.transpose();
  C_i = F.transpose()*(be_inverse*F);
  return C_i;
};


MatrixXd viscous_branch::update_intervar_newton_principal_hencky(MatrixXd F){
  eigs = eigs_tr;
  epsilon_tr = 0.5*log(eigs_tr.array());
  epsilon = epsilon_tr;


  compute_residual_principal_hencky();
  compute_tangent_principal_hencky();
  compute_tangent_nr_principal();
  depsilon = -1.0*Tangent_principal.householderQr().solve(res_principal);
  epsilon += depsilon;
  eigs = exp(2.0*epsilon.array());
  compute_residual_principal_hencky();
  compute_tangent_principal_hencky();
  compute_tangent_nr_principal();

  be = eigs(0)*v0_*v0_.transpose()+eigs(1)*v1_*v1_.transpose();
  be_inverse = 1.0/eigs(0)*v0_*v0_.transpose()+1.0/eigs(1)*v1_*v1_.transpose();
  C_i = F.transpose()*(be_inverse*F);

  return C_i;
};
