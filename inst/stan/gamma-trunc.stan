functions {
    real[] fstar_ode(real t, real[] z, real[] theta, data real[] x_r, int[] x_i) {
        int N = x_i[1]; // number of records

        real e[N] = theta[1:N];
        real param1 = theta[N+1];
        real param2 = theta[N+2];

        real upper_bound = x_r[1];
        real r = x_r[2];

        real dzdt[N];
        real tstar[N];

        for (k in 1:N) {
            tstar[k] = upper_bound - e[k];
            dzdt[k] = exp(gamma_lcdf(tstar[k]*(1.0-t) | param1, param2)) * r * tstar[k] * exp(-r*tstar[k]*t) / (1.0 - exp(-r*tstar[k]*t));
        }

        return dzdt;
    }
}

data {
    int<lower = 0> N; // number of records
    vector<lower = 0>[N] E_L;
    vector<lower = 0>[N] E_R;
    vector<lower = 0>[N] S_L;
    vector<lower = 0>[N] S_R;
    real<lower = 0> r;
    real<lower = 0> upper_bound;
}

transformed data {
    int X_i[1] = {N};

    real X_r[2] = {upper_bound, r};
}

parameters {
    real<lower = 0> mean_SI;
    real<lower = 0> sd_SI;

    vector<lower = 0, upper = 1>[N] s_raw;
    vector<lower = 0, upper = 1>[N] e_raw;
}

transformed parameters {
    real<lower = 0> param1 = square(mean_SI/sd_SI);
    real<lower = 0> param2 = mean_SI/square(sd_SI);

    vector<lower = min(S_L), upper = max(S_R)>[N] s;
    vector<lower = min(E_L), upper = max(E_R)>[N] e;

    vector[N] Z;

    {
        real theta[N+2];
        real Z0[N];

        s = S_L + (S_R - S_L) .* s_raw;
        for (k in 1:N)
            if (E_R[k] > s[k])
                e[k] = E_L[k] + (s[k] - E_L[k]) * e_raw[k];
            else
                e[k] = E_L[k] + (E_R[k] - E_L[k]) * e_raw[k];

        for (k in 1:N) {
            Z0[k] = 0.0;
            theta[k] = e[k];
        }
        theta[N+1] = param1;
        theta[N+2] = param2;

        Z = to_vector(to_array_1d(integrate_ode_rk45(fstar_ode, Z0, 0.001, {1.0}, theta, X_r, X_i)));
    }
}

model {
    mean_SI ~ normal(5.0, 10.0);
    sd_SI ~ cauchy(0, 5.0);

    e_raw ~ normal(0.5, 1.0);
    s_raw ~ normal(0.5, 1.0);

    target += gamma_lpdf(s - e | param1, param2) - log(Z);
}

generated quantities {
    vector[N] log_likelihood;
    for (k in 1:N)
        log_likelihood[k] = gamma_lpdf(s[k] - e[k] | param1, param2) - log(Z[k]);
}
