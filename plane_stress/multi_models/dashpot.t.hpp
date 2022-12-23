template <typename model>
void dashpot::set_viscous_para(double dt,model* material){
  gamma_0 = 1.0; //1.0/60.0
  eta = 2.0*(material->rel_time())*(material->coeff());
  tau_hat = rt_to*eta;
  delta_t = dt;
  m=1.0;
};

template <typename model>
void dashpot::compute_tangent_nr_principal(model* material){
  Stiff_principal=material->compute_stiff_principal(eigs);
  compute_tangent_pressure();
  tot_tangent = Stiff_principal+pressure_tangent;
  Tangent_principal = gamma_0*delta_t/rt_to/pow(tau_hat,1)*pow(tau_v,m-1)*tot_tangent;
  Tangent_principal += MatrixXd::Identity(2,2);
  // tau_v = (tau_v<1.0)?1.0:tau_v;
  // if(tau_v<1){
  //   tau_v=1.0; }
  // Tangent_principal(0,0) += gamma_0*delta_t/pow(tau_hat,1)*(m-1)*pow(tau_v,m-3)/to_rt_to*tau_dev(0)*(tau_dev(0)*tot_tangent(0,0)+tau_dev(1)*tot_tangent(1,0)+p*pressure_tangent(0,0));
  // Tangent_principal(0,1) += gamma_0*delta_t/pow(tau_hat,1)*(m-1)*pow(tau_v,m-3)/to_rt_to*tau_dev(0)*(tau_dev(0)*tot_tangent(0,1)+tau_dev(1)*tot_tangent(1,1)+p*pressure_tangent(0,1));
  // Tangent_principal(1,0) += gamma_0*delta_t/pow(tau_hat,1)*(m-1)*pow(tau_v,m-3)/to_rt_to*tau_dev(1)*(tau_dev(0)*tot_tangent(0,0)+tau_dev(1)*tot_tangent(1,0)+p*pressure_tangent(0,0));
  // Tangent_principal(1,1) += gamma_0*delta_t/pow(tau_hat,1)*(m-1)*pow(tau_v,m-3)/to_rt_to*tau_dev(1)*(tau_dev(0)*tot_tangent(0,1)+tau_dev(1)*tot_tangent(1,1)+p*pressure_tangent(0,1));
};

template <typename model>
void dashpot::compute_residual_principal(model* material){
  dW_d = material->compute_derivative(eigs);
  compute_stress_tau_principal(dW_d);
  tau_dev = tau_principal.array()+p;
  tau_v = 1.0;//sqrt(tau_dev.transpose()*tau_dev+p*p)/rt_to;
  res_principal = epsilon + gamma_0*delta_t/rt_to/pow(tau_hat,1)*pow(tau_v,m-1)*tau_dev - epsilon_tr;
};

template <typename model>
Matrix2d dashpot::update_intervar_newton_principal(Matrix2d F,model* material){
  eigs = eigs_tr;
  epsilon_tr = 0.5*log(eigs_tr.array());
  epsilon = epsilon_tr;
  count=0;
  compute_residual_principal(material);
  compute_tangent_nr_principal(material);
  // eps_max=epsilon_tr.maxCoeff();
    while (res_principal.norm() > 1E-6 && count<10){ //do
      // compute_residual_principal(material);
      // compute_tangent_nr_principal(material);
      depsilon = Tangent_principal.colPivHouseholderQr().solve(res_principal);
      epsilon = epsilon - depsilon;
      eigs = exp(2.0*epsilon.array());
      compute_residual_principal(material);
      compute_tangent_nr_principal(material);
      ++count;
    }
    //while (res_principal.norm() > 1E-6 && count<10);
  // compute_residual_principal(material);
  // compute_tangent_nr_principal(material);
  be_inverse = 1.0/eigs(0)*v0_*v0_.transpose()+1.0/eigs(1)*v1_*v1_.transpose();
  C_i = F.transpose()*(be_inverse*F);
  return C_i;
};
