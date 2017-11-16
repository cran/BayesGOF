EXP.denom <-
function(c.vec, weight, leg.mat){
		Ep.d <- NULL
		Ep.d <- ( t(leg.mat) %*% weight )/ length(weight)
		denom <- as.vector(1 + as.vector(Ep.d) %*% c.vec) #scalar
		return(denom)
		}
