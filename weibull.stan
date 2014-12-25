data{
	int<lower=1> N;
	int<lower=1> NS;
	int<lower=1> NI;
	int<lower=1, upper=NS> subject[N];
	int<lower=1, upper=NI> item[N];
	vector[N] cond1;
	vector[N] cond2;
	vector<lower=0>[N] rt;
}

parameters{
	vector[3] beta;
	vector[3] u[NS];
	vector[3] w[NI]; 
	vector[N] r;
	//real<lower=0> sigma_e; 
	vector<lower=0>[3] sigma_u;
	vector<lower=0>[3] sigma_w;
	vector<lower=0>[N] shape;
}

model{
	sigma_u ~ gamma(1.5, 1.0E-4);
	sigma_w ~ gamma(1.5, 1.0E-4);
	r ~ gamma(1, 0.001);

	for(i in 1:N){
		shape[i] <- exp(-((beta[1] + u[subject[i], 1] + w[item[i], 1]) +
				 	(beta[2] + u[subject[i], 2] + w[item[i], 2])*cond1[i] + 
				 	(beta[3] + u[subject[i], 3] + w[item[i], 3])*cond2[i])/r[i]);
	}
	for(j in 1:NS){
		u[j, 1] ~ normal(0, sigma_u[1]);
		u[j, 2] ~ normal(0, sigma_u[2]);
		u[j, 3] ~ normal(0, sigma_u[3]);
	}

	for(k in 1:NI){
		w[k, 1] ~ normal(0, sigma_w[1]);
		w[k, 2] ~ normal(0, sigma_w[2]);
		w[k, 3] ~ normal(0, sigma_w[3]);
	}

	rt ~ weibull(r, shape);



}