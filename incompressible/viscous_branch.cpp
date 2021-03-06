// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "Material.h"
#include "viscous_branch.h"


// Namespaces
using namespace Eigen;
using namespace std;


void viscous_branch::set_inelastic_strain(MatrixXd _C_i){
  C_i = _C_i;
  eta = relaxation_time*c[0];
};

void viscous_branch::compute_be_tr(MatrixXd F){
  be_tr = F*(C_i.inverse()*F.transpose());
  EigenSolver<Matrix3d> es(be_tr);
  v0_ = es.eigenvectors().col(0).real()(seq(0,2),0);
  v1_ = es.eigenvectors().col(1).real()(seq(0,2),0);
  v2_ = es.eigenvectors().col(2).real()(seq(0,2),0);

  v0 = es.eigenvectors().row(0).real()(0,seq(0,2));
  v1 = es.eigenvectors().row(1).real()(0,seq(0,2));
  v2 = es.eigenvectors().row(2).real()(0,seq(0,2));

  eigs_tr(0) = es.eigenvalues()(0,0).real();
  eigs_tr(1) = es.eigenvalues()(1,0).real();
  eigs_tr(2) = es.eigenvalues()(2,0).real();

  //
  Rotation_mat(0,0) = v0(0)*v0(0);
  Rotation_mat(1,0) = v1(0)*v1(0);
  Rotation_mat(2,0) = v2(0)*v2(0);
  Rotation_mat(3,0) = v0(0)*v1(0);
  Rotation_mat(4,0) = v1(0)*v2(0);
  Rotation_mat(5,0) = v2(0)*v0(0);

  Rotation_mat(0,1) = v0(1)*v0(1);
  Rotation_mat(1,1) = v1(1)*v1(1);
  Rotation_mat(2,1) = v2(1)*v2(1);
  Rotation_mat(3,1) = v0(1)*v1(1);
  Rotation_mat(4,1) = v1(1)*v2(1);
  Rotation_mat(5,1) = v2(1)*v0(1);

  Rotation_mat(0,2) = v0(2)*v0(2);
  Rotation_mat(1,2) = v1(2)*v1(2);
  Rotation_mat(2,2) = v2(2)*v2(2);
  Rotation_mat(3,2) = v0(2)*v1(2);
  Rotation_mat(4,2) = v1(2)*v2(2);
  Rotation_mat(5,2) = v2(2)*v0(2);

  Rotation_mat(0,3) = 2.0*v0(0)*v0(1);
  Rotation_mat(1,3) = 2.0*v1(0)*v1(1);
  Rotation_mat(2,3) = 2.0*v2(0)*v2(1);
  Rotation_mat(3,3) = v0(0)*v1(1)+v0(1)*v1(0);
  Rotation_mat(4,3) = v1(0)*v2(1)+v1(1)*v2(0);
  Rotation_mat(5,3) = v2(0)*v0(1)+v2(1)*v0(0);

  Rotation_mat(0,4) = 2.0*v0(1)*v0(2);
  Rotation_mat(1,4) = 2.0*v1(1)*v1(2);
  Rotation_mat(2,4) = 2.0*v2(1)*v2(2);
  Rotation_mat(3,4) = v0(1)*v1(2)+v0(2)*v1(1);
  Rotation_mat(4,4) = v1(1)*v2(2)+v1(2)*v2(1);
  Rotation_mat(5,4) = v2(1)*v0(2)+v2(2)*v0(1);

  Rotation_mat(0,5) = 2.0*v0(2)*v0(0);
  Rotation_mat(1,5) = 2.0*v1(2)*v1(0);
  Rotation_mat(2,5) = 2.0*v2(2)*v2(0);
  Rotation_mat(3,5) = v0(2)*v1(0)+v0(0)*v1(2);
  Rotation_mat(4,5) = v1(2)*v2(1)+v1(0)*v2(2);
  Rotation_mat(5,5) = v2(2)*v0(0)+v2(0)*v0(2);
  Rotation_mat_transp = Rotation_mat.transpose();

};

void viscous_branch::compute_invar(MatrixXd C, bool a){
  if(a){
    invar[0] = C(0,0) + C(1,1) + C(2,2);
    invar[1] = -0.5*(C(0,0)*C(0,0)+C(1,1)*C(1,1)+C(2,2)*C(2,2)+2.0*(C(1,0)*C(1,0)+C(1,2)*C(1,2)+C(0,2)*C(0,2)) - invar[0]*invar[0]);
  }
  if(!a){
    invar[0] = C(0)+ C(1) + C(2);
    invar[1] = C(0)*C(1) + C(1)*C(2) + C(0)*C(2);
  }

};

void viscous_branch::compute_derivative(){
  I1 = invar[0];
  I2 = invar[1];
  derivative[0]= c[0] + 2*c[2]*(I1-3) +c[3]*(I2-3) + 3*c[5]*(I1-3)*(I1-3) + 2*c[6]*(I1-3)*(I2-3) + c[7]*(I2-3)*(I2-3);
  derivative[1]= c[1] + 2*c[4]*(I2-3) +c[3]*(I1-3) + 3*c[8]*(I2-3)*(I2-3) + 2*c[7]*(I2-3)*(I1-3) + c[6]*(I1-3)*(I1-3);
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
  lambda_C = eigs(2);

  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  tau_principal(0) = (dW_dI1 + dW_dI2*(I1-lambda_A))*lambda_A;
  tau_principal(1) = (dW_dI1 + dW_dI2*(I1-lambda_B))*lambda_B;
  tau_principal(2) = (dW_dI1 + dW_dI2*(I1-lambda_C))*lambda_C;
  tau_principal = 2*tau_principal;
};

void viscous_branch::compute_stress_tau_principal_hencky(){

  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  tau_principal(0) = dW_dI1*epsilon(0);
  tau_principal(1) = dW_dI1*epsilon(1);
  tau_principal(2) = dW_dI1*epsilon(2);

};

MatrixXd viscous_branch::compute_stress_tau(){
  tau = tau_principal(0)*v0_*v0_.transpose()+tau_principal(1)*v1_*v1_.transpose()+tau_principal(2)*v2_*v2_.transpose();
  return tau;
};

double viscous_branch::compute_strain_energy(MatrixXd _be){
  double psi = 0.0;
  compute_invar(_be,1);
  psi += c[0]*(I1-3)+ c[1]*(I2-3);
  psi += c[2]*pow(I1-3,2.0)+c[3]*(I1-3)*(I2-3)+c[4]*pow(I2-3,2.0);
  psi += c[5]*pow(I1-3,3.0)+c[6]*pow(I1-3,2.0)*(I2-3);
  psi += c[7]*(I1-3)*pow(I2-3,2.0)+c[8]*pow(I2-3,3.0);
  return psi;
};

double viscous_branch::compute_strain_energy_hencky(MatrixXd _be){
  double psi = 0.0;
  psi += epsilon.squaredNorm();
  return 0.5*dW_dI1*psi;
};

double viscous_branch::compute_dissipation(){
  double psi = 0.0;
  double temp = 1.0/3.0*tau_principal.sum();
  psi += tau_principal.transpose()*(tau_principal-temp*VectorXd::Ones(3));
  psi /= eta;
  return psi;
};

void viscous_branch::compute_tangent_principal(){
  compute_second_derivative();
  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  Stiff_principal(0,0) = dW_dI1;
  Stiff_principal(0,0) *= 2*lambda_A;
  Stiff_principal(0,1) = 0.0;
  Stiff_principal(0,2) = 0.0;

  Stiff_principal(1,0) = 0.0;
  Stiff_principal(1,1) = dW_dI1;
  Stiff_principal(1,1) *= 2*lambda_B;
  Stiff_principal(1,2) = 0.0;

  Stiff_principal(2,0) = 0.0;
  Stiff_principal(2,1) = 0.0;
  Stiff_principal(2,2) = dW_dI1;
  Stiff_principal(2,2) *= 2*lambda_C;

  Stiff_principal *= 2.0;

};


void viscous_branch::compute_tangent_principal_hencky(){
  compute_second_derivative();
  dW_dI1 = derivative[0];
  Stiff_principal = dW_dI1*MatrixXd::Identity(3,3);
};

void viscous_branch::compute_tangent_nr_principal(){
  Tangent_principal = delta_t/2.0/eta*Stiff_principal+MatrixXd::Identity(3,3);
  for(int i=0;i<3;++i){
      Tangent_principal(seq(0,2),i).array() -= delta_t/2.0/eta/3.0*Stiff_principal(seq(0,2),i).sum();
  }
};



void viscous_branch::compute_residual_principal(){

  compute_invar(eigs,0);
  compute_derivative();
  compute_stress_tau_principal();

  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  res_principal = epsilon + delta_t/2.0/eta*(tau_principal - 1.0/3.0*tau_principal.sum()*VectorXd::Ones(3)) - epsilon_tr;
};


void viscous_branch::compute_residual_principal_hencky(){

  compute_invar(eigs,0);
  compute_derivative();
  compute_stress_tau_principal_hencky();

  I1 = invar[0];
  I2 = invar[1];
  dW_dI1 = derivative[0];
  dW_dI2 = derivative[1];

  res_principal = epsilon + delta_t/2.0/eta*(tau_principal - 1.0/3.0*tau_principal.sum()*VectorXd::Ones(3)) - epsilon_tr;
};

void viscous_branch::mat_tan_principal(bool elastic){
    //
    C_alg =  Stiff_principal*(Tangent_principal.inverse());
    mat_tan(seq(0,2),seq(0,2)).array() = C_alg - 2.0*Matrix3d(tau_principal.asDiagonal());

    if (eigs_tr(0)==eigs_tr(1)) {
      mat_tan(3,3) = 0.5*(C_alg(0,0)-C_alg(1,0))-tau_principal(0);
    }
    else {
      mat_tan(3,3) = (tau_principal(0)*eigs_tr(1) - tau_principal(1)*eigs_tr(0))/(eigs_tr(0) - eigs_tr(1));
    }

    if (eigs_tr(1)==eigs_tr(2)) {
      mat_tan(4,4) = 0.5*(C_alg(1,1)-C_alg(2,1))-tau_principal(1);
    }
    else {
      mat_tan(4,4) = (tau_principal(1)*eigs_tr(2) - tau_principal(2)*eigs_tr(1))/(eigs_tr(1) - eigs_tr(2));
    }

    if (eigs_tr(0)==eigs_tr(2)) {
      mat_tan(5,5) = 0.5*(C_alg(2,2)-C_alg(2,0))-tau_principal(2);
    }
    else {
      mat_tan(5,5) = (tau_principal(2)*eigs_tr(0) - tau_principal(0)*eigs_tr(2))/(eigs_tr(2) - eigs_tr(0));
    }

};

MatrixXd viscous_branch::mat_tan_volu(){
    double temp[3];
    for(int i=0;i<3;++i){
        temp[i] = -1.0*(C_alg(2,i)-1.0/3.0*(C_alg(seq(0,2),i).sum()));
    }
    double temp_sum = temp[0]+ temp[1]+temp[2];
    for(int i=0;i<3;++i){
        mat_tan_vol(seq(0,2),i).array() = temp[i]-1.0/3.0*temp_sum;;
    }
    mat_tan_vol = Rotation_mat*mat_tan_vol*(Rotation_mat.transpose());
    // std::cout << mat_tan_vol << '\n';
    mat_tan_vol = Rotation_mat*mat_tan_vol*(Rotation_mat_transp);
    MatrixXd mat_tan_vol_transp = mat_tan_vol.transpose();
    // std::cout << "mat tan volu" << '\n';
    // std::cout << (mat_tan_vol) << '\n';
    // return 0.5*(mat_tan_vol+mat_tan_vol_transp);
    return mat_tan_vol;
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
  be = eigs(0)*v0_*v0_.transpose()+eigs(1)*v1_*v1_.transpose()+eigs(2)*v2_*v2_.transpose();
  be_inverse = 1.0/eigs(0)*v0_*v0_.transpose()+1.0/eigs(1)*v1_*v1_.transpose()+1.0/eigs(2)*v2_*v2_.transpose();
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

  be = eigs(0)*v0_*v0_.transpose()+eigs(1)*v1_*v1_.transpose()+eigs(2)*v2_*v2_.transpose();
  be_inverse = 1.0/eigs(0)*v0_*v0_.transpose()+1.0/eigs(1)*v1_*v1_.transpose()+1.0/eigs(2)*v2_*v2_.transpose();
  C_i = F.transpose()*(be_inverse*F);

  return C_i;
};

void convert_matrix_to_vector(MatrixXd& Matrix, MatrixXd& Vector){
  Vector(0) = Matrix(0,0);
  Vector(1) = Matrix(1,1);
  Vector(2) = Matrix(2,2);
  Vector(3) = Matrix(0,1);
  Vector(4) = Matrix(1,2);
  Vector(5) = Matrix(0,2);
};

void convert_vector_to_matrix(MatrixXd& Vector, MatrixXd& Matrix){
   Matrix(0,0) = Vector(0);
   Matrix(1,1) = Vector(1);
   Matrix(2,2) = Vector(2);
   Matrix(0,1) = Vector(3);
   Matrix(1,2) = Vector(4);
   Matrix(0,2) = Vector(5);
   Matrix(1,0) = Vector(3);
   Matrix(2,1) = Vector(4);
   Matrix(2,0) = Vector(5);
};

//
// MatrixXd viscous_branch::mat_tan_volu(){
//     double temp[3];
//     for(int i=0;i<3;++i){
//         mat_tan_vol(seq(0,2),i).array() = -1.0*(C_alg(2,i)-1.0/3.0*(C_alg(seq(0,2),i).sum()));
//     }
//     mat_tan_vol = Rotation_mat*mat_tan_vol*(Rotation_mat_transp);
//     MatrixXd mat_tan_vol_transp = mat_tan_vol.transpose();
//     return 0.5*(mat_tan_vol+mat_tan_vol_transp);
// };
