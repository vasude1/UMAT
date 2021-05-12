
void parse_deformation_variables(double* DFGRD0,double* DFGRD1,MatrixXd& F,MatrixXd& F_old){

  F(0,0) = *DFGRD1;
  F(1,0) = *(DFGRD1+1);
  F(0,1) = *(DFGRD1+3);
  F(1,1) = *(DFGRD1+4);
};

void parse_material_properties(double* props, double* stiff_ratio, double* relax, int number_branches){
  for(int i=0;i<number_branches;++i){
    stiff_ratio[i]=props[2*i];
    relax[i]=props[2*i+1];
  }
};

void parse_internal_variables(double* STATEV,Matrix2d* ivar, int number_branches){
  for(int i=0;i<number_branches;++i){
    *(ivar+i) << *(STATEV+5*i),*(STATEV+5*i+1),*(STATEV+5*i+2),
                *(STATEV+5*i+3);
   *(ivar+i) += MatrixXd::Identity(2,2);
   // cout<<*(ivar+i)<<endl;
  }
};


void deviatoric_projector_tangent(Matrix2d tau,Matrix3d* tangent){
// //
  VectorXd tau_row(3);
  //
  tau_row<<tau(0,0),tau(1,1),tau(0,1);
  // //
  (*tangent)(0,0) += 2*tau_row(0);
  (*tangent)(1,1) += 2*tau_row(1);
  (*tangent)(2,2) += 0.5*(tau_row(0)+tau_row(1));

  (*tangent)(2,0) += tau_row(2);
  (*tangent)(2,1) += tau_row(2);

  (*tangent)(0,2) += tau_row(2);
  (*tangent)(1,2) += tau_row(2);
//
};

void return_stress(double* stress, MatrixXd tau){
  *(stress) = tau(0,0);
  *(stress+1) = tau(1,1);
  *(stress+2) = tau(0,1);
};

void return_tangent(double* DDSDDE, MatrixXd tangent){
  *(DDSDDE+0) = tangent(0,0);
  *(DDSDDE+1) = tangent(1,0);
  *(DDSDDE+2) = tangent(2,0);
  *(DDSDDE+3) = tangent(0,1);
  *(DDSDDE+4) = tangent(1,1);
  *(DDSDDE+5) = tangent(2,1);
  *(DDSDDE+6) = tangent(0,2);
  *(DDSDDE+7) = tangent(1,2);
  *(DDSDDE+8) = tangent(2,2);
};

void return_internalvar(double* STATEV,Matrix2d* ivar, int number_branches){
  Matrix2d temp;
  for(int i=0;i<number_branches;++i){
    temp = *(ivar+i)-MatrixXd::Identity(2,2);

    *(STATEV+5*i) = temp(0,0);
    *(STATEV+5*i+1) = temp(0,1);

    *(STATEV+5*i+2) = temp(1,0);
    *(STATEV+5*i+3) = temp(1,1);

  }
  double total_SE=0.0;
  for(int i=0;i<number_branches;++i){
    total_SE += *(STATEV+5*i+4);
  }
 
 total_SE += *(STATEV+5*number_branches);
 *(STATEV+5*number_branches+1) = total_SE/(*(STATEV+5*number_branches));

   for(int i=0;i<number_branches;++i){
    *(STATEV+5*i+4) /= *(STATEV+5*number_branches) ;
  }
	
};


void add_all(Matrix2d* tau_final,Matrix3d* tangent_final,double* SSE,double* SCD,Matrix2d* tau,
          Matrix3d* tangent,double* _SSE,double* _SCD ,int branches,int left){
  *SSE = 0.0;
  *SCD = 0.0;

  for(int i=left;i< branches; ++i)
  {
    *tau_final = *tau_final + *(tau+i);
    *tangent_final = *tangent_final + *(tangent+i);
    *SSE += *(_SSE+i);
    *SCD += *(_SCD+i);
  }
};
