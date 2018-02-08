DS.mode.map.nnu <-
function(y.i, se.i, g.par, LP.par, B){
#####################################################
# INPUTS
#  yi, ni			success/trials of single sample
#  g.par		alpha and beta for parametric prior, generate
#  B 				length of theta range (from 0 to 1)
#  LP.par			vector of LP means
# OUTPUTS
#  ds.mode			Modes of posterior
	post.mu.i <- lambda.i(se.i, g.par[2]) * g.par[1] +
					(1-lambda.i(se.i, g.par[2]))* y.i
	post.tau2.i <- (1-lambda.i(se.i, g.par[2]))*se.i^2 #output is VARIANCE
	theta.rng <- seq(g.par[1] - 3*sqrt(g.par[2]),
						g.par[1] + 3*sqrt(g.par[2]), length.out = 1000)
	base.den <- dnorm(theta.rng, post.mu.i, sd = sqrt(post.tau2.i))
	##u <- pbeta(theta.rng ,g.par[1], g.par[2]) 
	Leg.G <- LP.basis.norm(theta.rng, g.par, length(LP.par))
	adj <-1+ (Leg.G %*% (LP.par))
	ds.den <- base.den * adj
	ds.mode <- theta.rng[which.max(ds.den)]
	return(ds.mode)
}