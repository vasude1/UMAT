//--------------------------------*/

#include <stdio.h>
#include <iostream>
#include <aba_for_c.h>
#include <Eigen/Dense> // Eigen class
#include <list>
#include <string>

using namespace Eigen;
using namespace std;

#include "Parser.h"
#include "/home/vasudevan/PhD/Code/UMAT/Material.h"
#include "/home/vasudevan/PhD/Code/UMAT/viscous_branch.h"

extern "C"
//
void FOR_NAME(umat)(double* STRESS,double* STATEV,double* DDSDDE,double* SSE,double* SPD,double* SCD,
double* RPL,double* DDSDDT,double* DRPLDE,double* DRPLDT,double* STRAN,double* DSTRAN,
double* TIME,double* DTIME,double* TEMP,double* DTEMP,double* PREDEF,double* DPRED,char* CMNAME,int& NDI,int& NSHR,int& NTENS,
int& NSTATV,double* PROPS,int& NPROPS,double* COORDS,double* DROT,double* PNEWDT,double* CELENT,
double* DFGRD0,double* DFGRD1,int& NOEL,int& NPT,int& LAYER,int& KSPT,int& JSTEP,int& KINC)
{
  // std::cout << "Begin" << '\n';
// PLANE STRESS ONLY
  int number_branches = 10;
  MatrixXd F =  MatrixXd::Zero(3,3);
	MatrixXd F_old =  MatrixXd::Zero(3,3);
  parse_deformation_variables(DFGRD0,DFGRD1,F,F_old);
  MatrixXd C = F.transpose()*F;
  MatrixXd b = F*F.transpose();

  // std::cout << "C: " << '\n';
  // std::cout << C << '\n';
  double stiff_ratio[number_branches];
  double relaxation_time[number_branches];
  parse_material_properties(PROPS,stiff_ratio,relaxation_time,number_branches);


  // Matrix3d *ivar = Matrix3d[10];
  Matrix3d ivar[number_branches];

  parse_internal_variables(STATEV,ivar,number_branches);


	MatrixXd C_old = F_old.transpose()*F_old;


	// viscous_branch* branch = viscous_branch[10];
  viscous_branch branch[number_branches];

  viscous_branch hyper_branch;
  MatrixXd tau = MatrixXd::Zero(3,3);
  Matrix3d temp_tau;
  MatrixXd tangent = MatrixXd::Zero(6,6);


  *SSE = 0.0;
  *SCD = 0.0;
  //
  // for(int i=0;i<number_branches;++i)
  // {
  //     // std::cout << i << '\n';
  //     (branch+i)->delta_t = *DTIME;
  //     (branch+i)->Material::set_material(stiff_ratio[i],relaxation_time[i]);
  //     // std::cout << *(ivar+i) << '\n';
  //     (branch+i)->set_inelastic_strain(*(ivar+i));
  //     (branch+i)->compute_be_tr(F);
  //     *(ivar+i)=(branch+i)->update_intervar_newton_principal(F);
  //     temp_tau = (branch+i)->compute_stress_tau(branch[i].be);
  //     tau += temp_tau;
  //     tangent += (branch+i)->rotate_mat_tan(0);
  //     *SSE += (branch+i)->compute_strain_energy(branch[i].be);
  //     *SCD += (branch+i)->compute_dissipation();
  //     // std::cout << *(ivar+i) << '\n';
  //
  // }
  // std::cout << "Elastic Branch" << '\n';
  hyper_branch.Material::set_material(1.0,1E7);
  hyper_branch.set_inelastic_strain(MatrixXd::Identity(3,3));
  hyper_branch.compute_be_tr(F);
  hyper_branch.update_intervar_newton_principal(F);
  tau += hyper_branch.compute_stress_tau(b);
  tangent += hyper_branch.rotate_mat_tan(1);
  *SSE += hyper_branch.compute_strain_energy(b);

  MatrixXd tau_temp = tau;
  deviatoric_projector(tau);
  double p = tau(2,2);
  tau += -p*MatrixXd::Identity(3,3);
  deviatoric_projector_tangent(p,tau_temp,tangent);
  // tangent = (tangent + tangent.transpose());

  // throw;
  return_stress(STRESS,tau);
  return_tangent(DDSDDE,tangent);
  return_internalvar(STATEV,ivar,number_branches);

return;
};
