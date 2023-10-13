data{
	int N; /* # of populations */
	vector[N] x; /* locations */
	int<lower=0> y[N]; /* number of stripe alleles */
	int<lower=0> n[N]; /* number of alleles */
}

parameters{
	real<lower=-10, upper=10> c; /* int */
	real<lower=-10, upper=10> beta; /* slope */
}

transformed parameters{
	vector<lower=0, upper=1>[N] p; /* vector of stripe freq. */
	real w; /* cline wdith */
	
	w = 4/beta;
	for(i in 1:N){
		p[i] = inv_logit(c + beta * x[i]);
	}
}

model{

	for(i in 1:N){
		/* increment likelihood */
		target += binomial_lpmf(y[i] | n[i], p[i]);
	}
	/* priors */
	target += normal_lpdf(c | 0, 10);
	target += normal_lpdf(beta | 0, 10);
}

