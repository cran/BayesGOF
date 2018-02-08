DS.mode.map.pgu <- 
function(y.i, g.par, LP.par, B){
#####################################################
# INPUTS
#  yi, ni			success/trials of single sample
#  g.par		alpha and beta for parametric prior, generate
#  B 				length of theta range (from 0 to 1)
#  LP.par			vector of LP means
# OUTPUTS
#  ds.mode			Modes of posterior
	post.par1.i <- y.i + g.par[1]
	post.par2.i <- g.par[2]/(1+g.par[2])
	theta.rng <- seq(1, 100, length.out = B)
	base.den <- dgamma(theta.rng, post.par1.i, scale = post.par2.i)
	u <- seq(1/B,1-1/B, length.out = B) 
	Leg.G <- LP.basis.gamma(theta.rng, g.par, length(LP.par))
	adj <-1+ (Leg.G %*% (LP.par))
	ds.den <- base.den * adj
	ds.mode <- theta.rng[which.max(ds.den)]
	return(ds.mode)
}