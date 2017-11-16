EXP.num.2 <-
function(c.vec, weight, leg.mat){
		Ep.d <- NULL
		for(j in 1:length(c.vec)){
			inner.part <- apply(leg.mat, 2, function(x) x*leg.mat[,j])
			v <-(t(inner.part)%*% weight)/nrow(leg.mat)
			Ep.d[j] <- as.vector(v) %*% c.vec
			}
		return(Ep.d)
}
