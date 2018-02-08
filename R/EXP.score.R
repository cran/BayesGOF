EXP.score <- 
function(u, LP.par, weight.vec, leg.mat){
		vec.map <- (EXP.num.1(weight.vec, leg.mat) + 
						EXP.num.2(LP.par, weight.vec, leg.mat) ) / 
						EXP.denom(LP.par, weight.vec, leg.mat)
			return(vec.map)
}