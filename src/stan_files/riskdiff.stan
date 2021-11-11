data {
  int<lower = 1>  ns;
  int<lower = 0>  n_pla[ns, 2];
  int<lower = 0>  n_trt[ns, 2];
  real<lower = 0> tau_max;
}

parameters {
  real<lower = 0,  upper = 1>       pi_c[ns];
  vector[ns]                        eta_i;
  real<lower = -1, upper = 1>       eta;
  real<lower = 0,  upper = tau_max> tau2;
}

transformed parameters {
  real<lower = 0, upper = 1> pi_t[ns];

  {
    real pt;
    for (i in 1:ns) {
      pt      = pi_c[i] + eta_i[i];
      pt      = pt > 1 ? 0.99999999 : pt;
      pi_t[i] = pt < 0 ? 0.00000001 : pt;
    }
  }
}

model {
  pi_c  ~ uniform(0,  1);
  eta   ~ uniform(-1, 1);
  eta_i ~ normal(eta, tau2);
  tau2  ~ uniform(0,  tau_max);

  // likelihood
  for (i in 1:ns) {
    n_pla[i, 1] ~ binomial(n_pla[i, 2], pi_c[i]);
    n_trt[i, 1] ~ binomial(n_trt[i, 2], pi_t[i]);
  }
}
