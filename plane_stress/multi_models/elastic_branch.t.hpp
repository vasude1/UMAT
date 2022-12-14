
template <typename model>
void elastic_branch::compute_tangent_nr_principal(model* material){
  Stiff_principal=material->compute_stiff_principal(eigs_tr);
  compute_tangent_pressure();
  tot_tangent = Stiff_principal+pressure_tangent;
};

template <typename model>
void elastic_branch::compute_residual_principal(model* material){
  dW_d = material->compute_derivative(eigs_tr);
  compute_stress_tau_principal(dW_d);
  tau_dev = tau_principal.array()+p;
};

template <typename model>
void elastic_branch::update_intervar_newton_principal(Matrix2d F,model* material){
  compute_residual_principal(material);
  compute_tangent_nr_principal(material);
};
