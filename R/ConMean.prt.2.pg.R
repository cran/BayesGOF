ConMean.prt.2.pg <- 
function(c.vec, weight, leg.mat, u, par1, par2){
	#### theta is G inv u = qgamma(u,muhat,tauhat)
	inner.part <- qgamma(u, par1, scale = par2)*leg.mat
	v <-(t(inner.part)%*% weight)/nrow(leg.mat)
	num.2 <- as.vector(v) %*% c.vec
	return(num.2)
}
