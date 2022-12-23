//--------------------------------*/
#define EIGEN_NO_DEBUG
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

#include "Parser_ps.h"
#include "path_to/elastic_branch.h"
#include "path_to/elastic_branch.t.hpp"
#include "path_to/dashpot.h"
#include "path_to/dashpot.t.hpp"

#include "path_to/Polynomial.h"
#include "path_to/Hencky.h"

template <class model>
struct branch{dashpot dash_pot; model spring;};

template <class model>
struct elas_branch{elastic_branch dash_pot; model spring;};

//
extern "C"
//
void FOR_NAME(umat)(double* STRESS,double* STATEV,double* DDSDDE,double* SSE,double* SPD,double* SCD,
double* RPL,double* DDSDDT,double* DRPLDE,double* DRPLDT,double* STRAN,double* DSTRAN,
double* TIME,double* DTIME,double* TEMP,double* DTEMP,double* PREDEF,double* DPRED,char* CMNAME,int& NDI,int& NSHR,int& NTENS,
int& NSTATV,double* PROPS,int& NPROPS,double* COORDS,double* DROT,double* PNEWDT,double* CELENT,
double* DFGRD0,double* DFGRD1,int& NOEL,int& NPT,int& LAYER,int& KSPT,int& JSTEP,int& KINC)
{
  int number_branches =7;
  int left=0;
  Matrix2d F =  MatrixXd::Zero(2,2);
	Matrix2d F_old =  MatrixXd::Zero(2,2);
  parse_deformation_variables(DFGRD0,DFGRD1,F,F_old);
  double stiff_ratio[number_branches];
  double relaxation_time[number_branches];
  parse_material_properties(PROPS,stiff_ratio,relaxation_time,number_branches);

  Matrix2d ivar[number_branches];
  parse_internal_variables(STATEV,ivar,number_branches);

  elas_branch<Polynomial> hyper;
  branch<Polynomial> viscous[number_branches];

  Matrix2d tau[number_branches+1];
  Matrix3d tangent[number_branches+1];
  double _SSE[number_branches+1];
  double _SCD[number_branches];

  Matrix2d tau_final = MatrixXd::Zero(2,2);
  Matrix3d tangent_final = MatrixXd::Zero(3,3);


  // #pragma omp parallel num_threads(1)
  // {
  // #pragma omp for
  for(int i=left;i<number_branches+1;++i)
  {
      if(i==number_branches)
      {
        (hyper.dash_pot).compute_be_tr(F);
        (hyper.dash_pot).update_intervar_newton_principal(F,&(hyper.spring));
        *(tau+i) = (hyper.dash_pot).compute_stress_tau();
        *(tangent+i)= (hyper.dash_pot).rotate_mat_tan();
        *(STATEV+5*i) = *(_SSE+i) = (hyper.spring).compute_strain_energy();
      }
      else
      {
        ((viscous+i)->spring).set_material(stiff_ratio[i],relaxation_time[i]);
        ((viscous+i)->dash_pot).set_viscous_para(*DTIME,&((viscous+i)->spring));
        ((viscous+i)->dash_pot).set_inelastic_strain(*(ivar+i));
        ((viscous+i)->dash_pot).compute_be_tr(F);
        *(ivar+i)=((viscous+i)->dash_pot).update_intervar_newton_principal(F,&((viscous+i)->spring));
         *(tau+i)= ((viscous+i)->dash_pot).compute_stress_tau();
         *(tangent+i)= ((viscous+i)->dash_pot).rotate_mat_tan();
        *(STATEV+5*i+4) = *(_SSE+i) = ((viscous+i)->spring).compute_strain_energy();
        *(_SCD+i) = ((viscous+i)->dash_pot).compute_dissipation();

      }
    }
  // }
	double SCD_rate;
  add_all(&tau_final,&tangent_final,SSE,&SCD_rate,tau,tangent,_SSE,_SCD,number_branches+1,left);
	*SCD += SCD_rate*(*DTIME);
  deviatoric_projector_tangent(tau_final,&tangent_final);
  return_stress(STRESS,tau_final);
  return_tangent(DDSDDE,tangent_final);
  return_internalvar(STATEV,ivar,number_branches);
return;
};
