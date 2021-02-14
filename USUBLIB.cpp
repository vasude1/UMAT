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
#include "Material.h"
#include "viscous_branch.h"

extern "C"
//
void FOR_NAME(umat)(double* STRESS,double* STATEV,double* DDSDDE,double* SSE,double* SPD,double* SCD,
double* RPL,double* DDSDDT,double* DRPLDE,double* DRPLDT,double* STRAN,double* DSTRAN,
double* TIME,double* DTIME,double* TEMP,double* DTEMP,double* PREDEF,double* DPRED,char* CMNAME,int& NDI,int& NSHR,int& NTENS,
int& NSTATV,double* PROPS,int& NPROPS,double* COORDS,double* DROT,double* PNEWDT,double* CELENT,
double* DFGRD0,double* DFGRD1,int& NOEL,int& NPT,int& LAYER,int& KSPT,int& JSTEP,int& KINC)
{

// PLANE STRESS ONLY

  MatrixXd F =  MatrixXd::Zero(3,3);
	MatrixXd F_old =  MatrixXd::Zero(3,3);
  parse_deformation_variables(DFGRD0,DFGRD1,F,F_old);
  MatrixXd C = F.transpose()*F;
  MatrixXd b = F*F.transpose();

  double stiff_ratio[10];
  double relaxation_time[10];
  parse_material_properties(PROPS,stiff_ratio,relaxation_time);

  Matrix3d *ivar =new Matrix3d [10];
  parse_internal_variables(STATEV,ivar);

	MatrixXd C_old = F_old.transpose()*F_old;

	viscous_branch* branch = new viscous_branch[10];
  viscous_branch hyper_branch;
  MatrixXd tau = MatrixXd::Zero(3,3);
  MatrixXd tangent = MatrixXd::Zero(6,6);

  *SSE = 0.0;
  *SCD = 0.0;
  for(int i=0;i<10;++i)
  {
      (branch+i)->delta_t = *DTIME;
      (branch+i)->set_material(stiff_ratio[i],relaxation_time[i]);
      (branch+i)->set_inelastic_strain(*(ivar+i));
      (branch+i)->compute_be_tr(F);
      *(ivar+i)=(branch+i)->update_intervar_newton_principal(F);
      tau += (branch+i)->compute_stress_tau(branch[i].be);
      tangent += (branch+i)->rotate_mat_tan();
      *SSE += (branch+i)->compute_strain_energy(branch[i].be);
      *SCD += (branch+i)->compute_dissipation();
  }

  hyper_branch.set_material(1.0,1E7);
  hyper_branch.set_inelastic_strain(MatrixXd::Identity(3,3));
  hyper_branch.compute_be_tr(F);
  hyper_branch.update_intervar_newton_principal(F);
  tau += hyper_branch.compute_stress_tau(b);
  tangent += hyper_branch.rotate_mat_tan();
  *SSE += hyper_branch.compute_strain_energy(b);

  double p = tau(2,2);

  deviatoric_projector_tangent(p,tau,tangent);
  deviatoric_projector(tau);
  tau += -p*MatrixXd::Identity(3,3);

  return_stress(STRESS,tau);
  return_tangent(DDSDDE,tangent);
  return_internalvar(STATEV,ivar);



return;
};
