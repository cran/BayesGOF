DS.PostMean.pge <- 
function(x.i, e.i, g.par, u, LP.par){
	fam <- "Poisson"
	B <- 250
	post.par1.i <- x.i + g.par[1]
	post.par2.i <- g.par[2]/(1+e.i*g.par[2])
	#u <- seq(1/B,1-(1/B), length.out = B) 
	weight <- weight.fun.univ(u, g.par[1],g.par[2],post.par1.i,post.par2.i, family = fam)
	leg.mat <- LP.basis.beta(u, c(1,1), length(LP.par))
	post.mu <- post.par1.i * post.par2.i
	ConMean.ind <- (post.mu + ConMean.prt.2.pg(LP.par, weight, leg.mat, u, g.par[1], g.par[2]) ) / EXP.denom(LP.par,weight,leg.mat)
	return(ConMean.ind)	
}
