DS.PostMean.nnu <-
function(y.i, se.i, g.par, u, LP.par){
		fam <- "Normal"
		post.mu.i <- lambda.i(se.i, g.par[2]) * g.par[1] +
					(1-lambda.i(se.i, g.par[2]))* y.i
		post.tau2.i <- (1-lambda.i(se.i, g.par[2]))*se.i^2 #output is VARIANCE
		weight <- weight.fun.univ(u, g.par[1], g.par[2], post.mu.i, post.tau2.i, family = fam)
		leg.mat <- LP.basis.beta(u,c(1,1), length(LP.par))
		if(sum(LP.par^2) == 0){
				ConMean.ind <- post.mu.i
				} else {
				ConMean.ind <- ( post.mu.i + ConMean.prt.2.nn(LP.par, weight, leg.mat, u, g.par[1], g.par[2]) ) / EXP.denom(LP.par,weight,leg.mat)
				}
			return(ConMean.ind)
	}
