
void parse_deformation_variables(double* DFGRD0,double* DFGRD1,MatrixXd F,MatrixXd F_old){
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

void parse_material_properties(double* props, double* stiff_ratio, double* relax){
  for(int i=0;i<10;++i){
    stiff_ratio[i]=props[2*i];
    relax[i]=props[2*i+1];
  }
};

void parse_internal_variables(double* STATEV,MatrixXd* ivar){
  for(int i=0;i<10;++i){
    *(ivar+i) << *(STATEV+5*i),*(STATEV+5*i+1),0.0,
                *(STATEV+5*i+2),*(STATEV+5*i+3),0.0,
                0.0,0.0,*(STATEV+5*i+4);
  }
};

MatrixXd deviatoric_projector(MatrixXd tau){
  tau -= 1.0/3.0*(tau(0,0)+tau(1,1)+tau(2,2))*MatrixXd::Identity(3,3);
  return tau;
};


MatrixXd deviatoric_projector_tangent(tau,tangent){
  MatrixXd tangent_temp = tangent;
  tangent.col(0).array() -= (tangent_temp.col(0).array() + tangent_temp.col(1).array() + tangent_temp.col(2).array())/3.0;
  tangent.col(2).array() -= (tangent_temp.col(0).array() + tangent_temp.col(1).array() + tangent_temp.col(2).array())/3.0;
  tangent.col(3).array() -= (tangent_temp.col(0).array() + tangent_temp.col(1).array() + tangent_temp.col(2).array())/3.0;

  tangent.row(0).array() -= (tangent_temp.row(0).array() + tangent_temp.row(1).array() + tangent_temp.row(2).array())/3.0;
  tangent.row(2).array() -= (tangent_temp.row(0).array() + tangent_temp.row(1).array() + tangent_temp.row(2).array())/3.0;
  tangent.row(3).array() -= (tangent_temp.row(0).array() + tangent_temp.row(1).array() + tangent_temp.row(2).array())/3.0;

  tangent(seq(0,2),seq(0,2)).array() += tangent_temp(seq(0,2),seq(0,2)).sum()/9.0;

  VectorXd tau_row(6);

  tau_row(0) = tau(0,0);
  tau_row(1) = tau(1,1);
  tau_row(2) = tau(2,2);
  tau_row(3) = tau(0,1);
  tau_row(4) = tau(1,2);
  tau_row(5) = tau(0,2);

  tangent(seq(0,5),0) -= 2.0/3.0*tau_row.array();
  tangent(seq(0,5),1) -= 2.0/3.0*tau_row.array();
  tangent(seq(0,5),2) -= 2.0/3.0*tau_row.array();

  tangent(0,seq(0,5)) -= 2.0/3.0*tau_row.array();
  tangent(1,seq(0,5)) -= 2.0/3.0*tau_row.array();
  tangent(2,seq(0,5)) -= 2.0/3.0*tau_row.array();





};
