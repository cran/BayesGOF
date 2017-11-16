DS.PostMean <-
function(y.i, n.i, g.par, u, LP.par, method = c("LP-Score", "Theta")){
		method = match.arg(method)
		post.par1.i <- y.i + g.par[1]
		post.par2.i <- n.i - y.i + g.par[2]
		weight <- weight.fun.beta(u, g.par[1],g.par[2],post.par1.i,post.par2.i)
		leg.mat <- LP.basis.beta(u,c(1,1), length(LP.par))
		#leg.mat2 <- Leg.poly(u, length(LP.par))
		switch(method,
		"LP-Score" = {
			vec.map <- (EXP.num.1(weight, leg.mat) + 
						EXP.num.2(LP.par, weight, leg.mat) ) / 
						EXP.denom(LP.par, weight, leg.mat)
			return(vec.map)
			},
		"Theta" = {
			if(sum(LP.par^2) == 0){
				ConMean.ind <- (post.par1.i)/(post.par1.i + post.par2.i)
				} else {
				post.mu <- post.par1.i /(post.par1.i + post.par2.i)
				ConMean.ind <- ( post.mu + PP.prt.2(LP.par, weight, leg.mat, u, g.par[1], g.par[2]) ) / EXP.denom(LP.par,weight,leg.mat)
				}
			return(ConMean.ind)
		}
		)
	}
