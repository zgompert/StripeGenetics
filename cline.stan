data{
	int N; /* # of populations */
	vector[N] x; /* locations */
	int<lower=0> y[N]; /* number of stripe alleles */
	int<lower=0> n[N]; /* number of alleles */
}

parameters{
	real<lower=-5, upper=5> c; /* cline center */
	real w; /* cline wdith */
}

transformed parameters{
	vector<lower=0, upper=1>[N] p; /* vector of stripe freq. */
	for(i in 1:N){
		p[i] = (1 + tanh(2 * (x[i] - c)/w))/2;
	}
}

model{

	for(i in 1:N){
		/* increment likelihood */
		target += binomial_lpmf(y[i] | n[i], p[i]);
	}
	/* priors */
	target += normal_lpdf(c | 0, 2);
	target += normal_lpdf(w | 0, 10);
}

