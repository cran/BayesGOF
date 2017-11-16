DS.mode.map <-
function(y.i, n.i, g.par, LP.par, B){
#####################################################
# INPUTS
#  yi, ni			success/trials of single sample
#  g.par		alpha and beta for parametric prior, generate
#  B 				length of theta range (from 0 to 1)
#  c.vec			vector of c-values
# OUTPUTS
#  ds.mode			Modes of posterior
	post.par1.i <- y.i + g.par[1]
	post.par2.i <- n.i - y.i + g.par[2]
	theta.rng <- seq(1/B,1-(1/B), length.out = B)
	base.den <- dbeta(theta.rng, post.par1.i, post.par2.i)
	u <- pbeta(theta.rng ,g.par[1], g.par[2]) 
	Leg.G <- LP.basis.beta(u,c(1,1),length(LP.par))
	adj <-1+ (Leg.G %*% (LP.par))
	ds.den <- base.den * adj
	ds.mode <- theta.rng[which.max(ds.den)]
	return(ds.mode)
}