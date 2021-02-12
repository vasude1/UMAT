//--------------------------------*/

#include <stdio.h>
#include <iostream>
#include <aba_for_c.h>
#include <Eigen/Dense> // Eigen class
#include <list>
#include <string>

#include "Parser.h"


using namespace Eigen;
using namespace std;

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

  MatrixXd *ivar =new Eigen::MatrixXd(3,3)[10];
  parse_internal_variables(STATEV,ivar);

	MatrixXd C_old = F_old.transpose()*F_old;

	branch = new viscous_branch[10];
  viscous_branch hyper_branch;
  MatrixXd tau = MatrixXd::Zero(3,3);
  MatrixXd tangent = MatrixXd::Zero(6,6);
  *SSE = 0.0;
  *SCD = 0.0;
  for(int i=0;i<10;++i)
  {
      branch[i]->set_material(stiff_ratio[i],relaxation_time[i]);
      branch[i]->set_inelastic_strain(*(ivar+i));
      branch[i]->compute_be_tr(F);
      *(ivar+i)=branch[i]->update_intervar_newton_principal(F);
      tau += branch[i]->compute_stress_tau(branch[i].be);
      tangent += branch[i]->rotate_mat_tan();
      *SSE += branch[i]->compute_strain_energy(branch[i].be);
      *SCD += branch[i]->compute_dissipation();
  }
  hyper_branch->set_inelastic_strain(MatrixXd::Identity(3,3));
  tau += hyper_branch->compute_stress_tau(b);
  tangent += hyper_branch->rotate_mat_tan();
  *SSE += hyper_branch->compute_strain_energy(b);

  tangent = deviatoric_projector_tangent(tau,tangent);
  tau = deviatoric_projector(tau);

  float p = tau(2,2);
  tau += -p*MatrixXd::Identity(3,3);





return;
};
