
void parse_deformation_variables(double* DFGRD0,double* DFGRD1,MatrixXd& F,MatrixXd& F_old){
  F_old(0,0) = DFGRD0[0];
  F_old(1,0) = DFGRD0[1];
  F_old(0,1) = DFGRD0[3];
  F_old(1,1) = DFGRD0[4];
  F_old(2,2) = 1.0/(F_old(0,0)*F_old(1,1) - F_old(0,1)*F_old(1,0));


  F(0,0) = DFGRD1[0];
  F(1,0) = DFGRD1[1];
  F(0,1) = DFGRD1[3];
  F(1,1) = DFGRD1[4];
  F(2,2) = 1.0/(F(0,0)*F(1,1) - F(0,1)*F(1,0));
};

void parse_material_properties(double* props, double* stiff_ratio, double* relax, int number_branches){
  for(int i=0;i<number_branches;++i){
    stiff_ratio[i]=props[2*i];
    relax[i]=props[2*i+1];
  }
};

void parse_internal_variables(double* STATEV,Matrix3d* ivar, int number_branches){
  for(int i=0;i<number_branches;++i){
    *(ivar+i) << *(STATEV+9*i),*(STATEV+9*i+1),*(STATEV+9*i+2),
                *(STATEV+9*i+3),*(STATEV+9*i+4),*(STATEV+9*i+5),
                *(STATEV+9*i+6),*(STATEV+9*i+7),*(STATEV+9*i+8);
   *(ivar+i) += MatrixXd::Identity(3,3);
   // cout<<*(ivar+i)<<endl;
  }
};

void deviatoric_projector(MatrixXd& tau){
  tau -= 1.0/3.0*(tau(0,0)+tau(1,1)+tau(2,2))*MatrixXd::Identity(3,3);
};


void deviatoric_projector_tangent(double p,MatrixXd& tau,MatrixXd& tangent){
  MatrixXd tangent_temp = tangent;
  tangent.col(0).array() -= (tangent_temp.col(0).array() + tangent_temp.col(1).array() + tangent_temp.col(2).array())/3.0;
  tangent.col(2).array() -= (tangent_temp.col(0).array() + tangent_temp.col(1).array() + tangent_temp.col(2).array())/3.0;
  tangent.col(3).array() -= (tangent_temp.col(0).array() + tangent_temp.col(1).array() + tangent_temp.col(2).array())/3.0;

  tangent.row(0).array() -= (tangent_temp.row(0).array() + tangent_temp.row(1).array() + tangent_temp.row(2).array())/3.0;
  tangent.row(2).array() -= (tangent_temp.row(0).array() + tangent_temp.row(1).array() + tangent_temp.row(2).array())/3.0;
  tangent.row(3).array() -= (tangent_temp.row(0).array() + tangent_temp.row(1).array() + tangent_temp.row(2).array())/3.0;

  tangent(seq(0,2),seq(0,2)).array() += tangent_temp(seq(0,2),seq(0,2)).sum()/9.0;

  VectorXd tau_row(6);
  double tr_tau = tau(0,0)+tau(1,1)+tau(2,2);
//
  tau_row(0) = tau(0,0)-tr_tau/3.0;
  tau_row(1) = tau(1,1)-tr_tau/3.0;
  tau_row(2) = tau(2,2)-tr_tau/3.0;
  tau_row(3) = tau(0,1);
  tau_row(4) = tau(1,2);
  tau_row(5) = tau(0,2);

  tangent(seq(0,5),0) -= 2.0/3.0*tau_row;
  tangent(seq(0,5),1) -= 2.0/3.0*tau_row;
  tangent(seq(0,5),2) -= 2.0/3.0*tau_row;

  tangent(0,seq(0,5)) -= 2.0/3.0*tau_row;
  tangent(1,seq(0,5)) -= 2.0/3.0*tau_row;
  tangent(2,seq(0,5)) -= 2.0/3.0*tau_row;
  //

  tr_tau *= 2.0/3.0;

  tangent(0,0) += 0.66*tr_tau;
  tangent(0,1) += -0.33*tr_tau;
  tangent(0,2) += -0.33*tr_tau;

  tangent(1,0) += -0.33*tr_tau;
  tangent(1,1) += 0.66*tr_tau;
  tangent(1,2) += -0.33*tr_tau;

  tangent(2,0) += -0.33*tr_tau;
  tangent(2,1) += -0.33*tr_tau;
  tangent(2,2) += 0.66*tr_tau;

  tangent(3,3) += 0.5*tr_tau;
  tangent(4,4) += 0.5*tr_tau;
  tangent(5,5) += 0.5*tr_tau;
  //
  // for(int i =0; i<6;++i){
  //   if(i<3){
  //     tangent(i,i) += 2*p;
  //   }
  //   else{
  //     tangent(i,i) += p;
  //   }
  // }
  // tangent(seq(0,2),seq(0,2)).array() -= p;

  tau_row(0) = tau(0,0);
  tau_row(1) = tau(1,1);
  tau_row(2) = tau(2,2);
  tau_row(3) = tau(0,1);
  tau_row(4) = tau(1,2);
  tau_row(5) = tau(0,2);

  tangent(seq(0,5),0) += 0.5*tau_row;
  tangent(seq(0,5),1) += 0.5*tau_row;
  tangent(seq(0,5),2) += 0.5*tau_row;

  tangent(0,seq(0,5)) += 0.5*tau_row;
  tangent(1,seq(0,5)) += 0.5*tau_row;
  tangent(2,seq(0,5)) += 0.5*tau_row;

  tangent(0,0) += 0.5*tau_row(0);
  tangent(1,1) += 0.5*tau_row(1);
  tangent(3,1) += 0.5*tau_row(3);
  tangent(4,2) += 0.5*tau_row(4);
  tangent(5,2) += 0.5*tau_row(5);
  tangent(2,2) += 0.5*tau_row(2);
  tangent(0,3) += 0.5*tau_row(3);
  tangent(1,4) += 0.5*tau_row(4);
  tangent(3,4) += 0.5*tau_row(5);
  tangent(0,5) += 0.5*tau_row(5);


  tangent(0,0) += 0.5*tau_row(0);
  tangent(3,0) += 0.5*tau_row(3);
  tangent(5,0) += 0.5*tau_row(5);
  tangent(1,1) += 0.5*tau_row(1);
  tangent(4,1) += 0.5*tau_row(4);
  tangent(3,3) += 0.5*tau_row(2);
  tangent(1,3) += 0.5*tau_row(3);
  tangent(4,3) += 0.5*tau_row(5);
  tangent(2,4) += 0.5*tau_row(4);
  tangent(2,5) += 0.5*tau_row(5);




};

void return_stress(double* stress, MatrixXd tau){
  *(stress+0) = tau(0,0);
  *(stress+1) = tau(1,1);
  *(stress+2) = tau(0,1);
};

void return_tangent(double* DDSDDE, MatrixXd tangent){
  *(DDSDDE+0) = tangent(0,0);
  *(DDSDDE+1) = tangent(1,0);
  *(DDSDDE+2) = tangent(4,0);
  *(DDSDDE+3) = tangent(0,1);
  *(DDSDDE+4) = tangent(1,1);
  *(DDSDDE+5) = tangent(4,1);
  *(DDSDDE+6) = tangent(0,4);
  *(DDSDDE+7) = tangent(1,4);
  *(DDSDDE+8) = tangent(4,4);
};

void return_internalvar(double* STATEV,Matrix3d* ivar, int number_branches){
  Matrix3d temp;
  for(int i=0;i<number_branches;++i){
    temp = *(ivar+i)-MatrixXd::Identity(3,3);

    *(STATEV+9*i) = temp(0,0);
    *(STATEV+9*i+1) = temp(0,1);
    *(STATEV+9*i+2) = temp(0,2);

    *(STATEV+9*i+3) = temp(1,0);
    *(STATEV+9*i+4) = temp(1,1);
    *(STATEV+9*i+5) = temp(1,2);

    *(STATEV+9*i+6) = temp(2,0);
    *(STATEV+9*i+7) = temp(2,1);
    *(STATEV+9*i+8) = temp(2,2);

  }
};
