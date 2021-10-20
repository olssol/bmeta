data {
  int<lower = 1> ns;

  int<lower = 0>  n_pla[ns, 2];
  int<lower = 0>  n_trt[ns, 2];
  real<lower = 0> tau_max;
  real<lower = 0> tau2_eta;
}

parameters {
  real<lower = 0, upper = 1>       pi_c[ns];
  vector[ns]                       eta_i;
  real                             eta;
  real<lower = 0, upper = tau_max> tau;
}

transformed parameters {
  real<lower = 0, upper = 1> pi_t[ns];
  real<lower = 0>            tau2;

  {
    real lpi_c;
    for (i in 1:ns) {
      lpi_c   = log(pi_c[i] + 0.0001);
      pi_t[i] = exp(eta_i[i] + lpi_c) > 1 ? 0.9999999 : exp(eta_i[i] + lpi_c);
    }
  }

  tau2 = tau^2;
}

model {
  pi_c  ~ uniform(0,  1);
  eta   ~ normal(0,   tau2_eta);
  eta_i ~ normal(eta, tau2);
  tau   ~ uniform(0,  tau_max);


  // likelihood
  for (i in 1:ns) {
    n_pla[i, 1] ~ binomial(n_pla[i, 2], pi_c[i]);
    n_trt[i, 1] ~ binomial(n_trt[i, 2], pi_t[i]);
  }
}
