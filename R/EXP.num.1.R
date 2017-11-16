EXP.num.1 <-
function(weight, leg.mat){
		Ep.d <- ( t(leg.mat) %*% weight )/ nrow(leg.mat) ##Vector length m
		return(as.vector(Ep.d))
		}
