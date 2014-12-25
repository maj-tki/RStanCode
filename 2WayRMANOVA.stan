data {
	int<lower=1> N;
	int<lower=1> NS;
	int<lower=1> NI;
	int<lower=1, upper=NS> subject[N];
	int<lower=1, upper=NI> item[N];
	int<lower=1, upper=2> x1[N];
	int<lower=1, upper=2> x2[N];
	vector[N] rt;
}

parameters {
	matrix[2, 2] a1a2;
	matrix[2, 2] aS[NS];
	matrix[2, 2] aI[NI];
	real<lower=0> sigma_e; 
	real<lower=0> sigma_aS;
	real<lower=0> sigma_aI;
	real<lower=0> sigma_a1a2;	
}


model{
	real mu[N];
	sigma_aS ~ gamma(1.5, 1.0E-4);
	sigma_aI ~ gamma(1.5, 1.0E-4);
	sigma_a1a2 ~ gamma(1.5, 1.0E-4);
	sigma_e ~ gamma(1.5, 1.0E-4);

	for(i in 1:2){
		for(j in 1:2){
			a1a2[i, j] ~ normal(0, sigma_a1a2);
		}
	}

	for(j in 1:NS){
		for(i in 1:2){
			for(k in 1:2){
				aS[j, i, k] ~ normal(0, sigma_aS);
			}
		}
	}
	for(j in 1:NI){
		for(i in 1:2){
			for(k in 1:2){
				aI[j, i, k] ~ normal(0, sigma_aI);
			}
		}
	}

	for(i in 1:N){
	  mu[i] <- a1a2[x1[i], x2[i]] + aS[subject[i], x1[i], x2[i]] + aI[item[i], x1[i], x2[i]];
	}
	rt ~ normal(mu, sigma_e);
}
