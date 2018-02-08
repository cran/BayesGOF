DS.PostMean.pgu <-
function(x.i, g.par, u, LP.par){
		fam <- "Poisson"
		post.par1.i <- x.i + g.par[1]
		post.par2.i <- g.par[2]/(1+g.par[2])
		weight <- weight.fun.univ(u, g.par[1],g.par[2],post.par1.i,post.par2.i, family = fam)
		leg.mat <- LP.basis.beta(u, c(1,1), length(LP.par))
		if(sum(LP.par^2) == 0){
				ConMean.ind <- post.par1.i*post.par2.i
				} else {
				post.mu <- post.par1.i*post.par2.i
				ConMean.ind <- (post.mu + ConMean.prt.2.pg(LP.par, weight, leg.mat, u, g.par[1], g.par[2]) ) / EXP.denom(LP.par,weight,leg.mat)
				}
			return(ConMean.ind)
		}
