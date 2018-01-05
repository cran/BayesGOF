EXP.score <- 
function(u, LP.par, weight.vec, leg.mat){
		#leg.mat <- LP.basis.beta(u,c(1,1), length(LP.par))
		vec.map <- (EXP.num.1(weight.vec, leg.mat) + 
						EXP.num.2(LP.par, weight.vec, leg.mat) ) / 
						EXP.denom(LP.par, weight.vec, leg.mat)
			return(vec.map)
}