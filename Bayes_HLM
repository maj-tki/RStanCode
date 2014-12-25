//varying intercept and varying slope

data {
	int<lower=1> N;                     //# of observations
	int<lower=1> NS;                    //# of subjects
	int<lower=1> NI;                    //# of items
	int<lower=1, upper=NS> subject[N];  //vector of subject IDs
	int<lower=1, upper=NI> item[N];     //vector of item IDs
	vector[N] cond1;                    //dummy coded predictor1 (level 2)
	vector[N] cond2;                    //dummy coded predictor2 (level 3)
	vector[N] rt;                       //predicted
	vector[3] zero;                     //vector of 0s
}

parameters {
	vector[3] beta;                     //vector of betas: intercept + slopes for dummy variables 1 & 2
	vector[3] u[NS];                    //subject random effects
	vector[3] w[NI];                    //item random effects
	real<lower=0> sigma_e;              //SD of residuals
	vector<lower=0>[3] sigma_u;         //SD of subject random effects
	vector<lower=0>[3] sigma_w;         //SD of item random effects
	corr_matrix[3] omega_u;             //correlation matrix for subject intercept and slopes
	corr_matrix[3] omega_w;             //correlation matrix for item intercept and slopes
}

transformed parameters{
	cov_matrix[3] Sigma_u;              //VCV matrix for subject raneff
	cov_matrix[3] Sigma_w;              //VCV matrix for item raneff
	for(r in 1:3){ 
		for(c in 1:3){
			Sigma_u[r, c] <- sigma_u[r] * sigma_u[c] * omega_u[r, c];   
			Sigma_w[r, c] <- sigma_w[r] * sigma_w[c] * omega_w[r, c];
		}
	}
}

model{
	real mu[N];
	for(j in 1:NS){
		u[j] ~ multi_normal(zero, Sigma_u);
	}
	for(k in 1:NI){
		w[k] ~ multi_normal(zero, Sigma_w);
	}
	
	sigma_u ~ gamma(1.5, 1.0E-4);
	sigma_w ~ gamma(1.5, 1.0E-4);
	sigma_e ~ gamma(1.5, 1.0E-4);

	omega_u ~lkj_corr(2.0);
	omega_w ~lkj_corr(2.0);

	for(i in 1:N){
		mu[i] <- (beta[1] + u[subject[i], 1] + w[item[i], 1]) +
				 (beta[2] + u[subject[i], 2] + w[item[i], 2])*cond1[i] + 
				 (beta[3] + u[subject[i], 3] + w[item[i], 3])*cond2[i];
	}
	rt ~ normal(mu, sigma_e);
}
