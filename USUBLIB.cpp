//--------------------------------*/

#include <stdio.h>
#include <iostream>
#include <aba_for_c.h>
#include <Eigen/Dense> // Eigen class
#include <list>
#include <string>
#include <chrono>

using namespace Eigen;
using namespace std;
using namespace std::chrono;

typedef Matrix<double, 6, 6> Matrix6d;

#include "Parser.h"
#include "/home/vasudevan/PhD/Code/UMAT/perfect_incomp/Material.h"
#include "/home/vasudevan/PhD/Code/UMAT/perfect_incomp/viscous_branch.h"

extern "C"
//
void FOR_NAME(umat)(double* STRESS,double* STATEV,double* DDSDDE,double* SSE,double* SPD,double* SCD,
double* RPL,double* DDSDDT,double* DRPLDE,double* DRPLDT,double* STRAN,double* DSTRAN,
double* TIME,double* DTIME,double* TEMP,double* DTEMP,double* PREDEF,double* DPRED,char* CMNAME,int& NDI,int& NSHR,int& NTENS,
int& NSTATV,double* PROPS,int& NPROPS,double* COORDS,double* DROT,double* PNEWDT,double* CELENT,
double* DFGRD0,double* DFGRD1,int& NOEL,int& NPT,int& LAYER,int& KSPT,int& JSTEP,int& KINC)
{
  // std::cout << NOEL << '\n';
  // std::cout << DDSDDE << '\n';
  // std::cout << "Begin" << '\n';
// PLANE STRESS ONLY
  int number_branches = 0;
  int left=0;
  MatrixXd F =  MatrixXd::Zero(3,3);
	MatrixXd F_old =  MatrixXd::Zero(3,3);
  parse_deformation_variables(DFGRD0,DFGRD1,F,F_old);
  MatrixXd C = F.transpose()*F;
  MatrixXd b = F*F.transpose();
  MatrixXd C_inv = C.inverse();
  // MatrixXd C_old = F_old.transpose()*F_old;

  double stiff_ratio[number_branches];
  double relaxation_time[number_branches];
  parse_material_properties(PROPS,stiff_ratio,relaxation_time,number_branches);

  Matrix3d ivar[number_branches];
  parse_internal_variables(STATEV,ivar,number_branches);

  viscous_branch branch[number_branches+1];

  Matrix3d tau[number_branches+1];
  Matrix6d tangent[number_branches+1];
  Matrix6d tangent_volu[number_branches+1];
  double _SSE[number_branches+1];
  double _SCD[number_branches+1];


  Matrix3d tau_final = MatrixXd::Zero(3,3);
  Matrix6d tangent_final = MatrixXd::Zero(6,6);
  Matrix6d tangent_volu_final = MatrixXd::Zero(6,6);
  // auto start = high_resolution_clock::now();

  //
  // #pragma omp parallel num_threads(number_branches+1)
  // {
  // #pragma omp for
  for(int i=left;i<number_branches+1;++i)
  {
      if(i==number_branches)
      {
        (branch+i)->delta_t = *DTIME;
        (branch+i)->Material::set_material(1.0,1E7);
        (branch+i)->set_inelastic_strain(MatrixXd::Identity(3,3));
        (branch+i)->compute_be_tr(F);
        (branch+i)->update_intervar_newton_principal(F);
        *(tau+i) = (branch+i)->compute_stress_tau();
        *(tangent+i) = (branch+i)->rotate_mat_tan(1);
        *(tangent_volu+i) = (branch+i)->mat_tan_volu();
        *(_SSE+i) = (branch+i)->compute_strain_energy(branch[i].be);
        *(_SCD+i) = (branch+i)->compute_dissipation();
      }
      else
      {
        (branch+i)->delta_t = *DTIME;
        (branch+i)->Material::set_material(stiff_ratio[i],relaxation_time[i]);
        (branch+i)->set_inelastic_strain(*(ivar+i));
        (branch+i)->compute_be_tr(F);
        // *(ivar+i)=(branch+i)->update_intervar_newton_principal_hencky(F);
        (branch+i)->update_intervar_newton_principal(F);
        *(tau+i) = (branch+i)->compute_stress_tau();
        *(tangent+i) = (branch+i)->rotate_mat_tan(0);
        *(tangent_volu+i) = (branch+i)->mat_tan_volu();
        *(_SSE+i) = (branch+i)->compute_strain_energy_hencky(branch[i].be);
        *(STATEV+10*i+9) = *(_SCD+i) = (branch+i)->compute_dissipation();
      }
  }
  // }
  // auto stop = high_resolution_clock::now();
  // auto duration = duration_cast<microseconds>(stop - start);
  // cout << "duration = " << duration.count() << endl;

  add_all(&tau_final,&tangent_final,&tangent_volu_final,SSE,SCD,tau,tangent,tangent_volu,_SSE,_SCD,number_branches+1,left);

  MatrixXd tau_temp = tau_final;
  deviatoric_projector(&tau_final);
  double p = -1.0*tau_final(2,2);
  tau_final += p*MatrixXd::Identity(3,3);
  deviatoric_projector_tangent(p,tau_temp,&tangent_final);
  tangent_final = tangent_final + tangent_volu_final;
  return_stress(STRESS,tau_final);
  return_tangent(DDSDDE,tangent_final);
  return_internalvar(STATEV,ivar,number_branches);

return;
};
