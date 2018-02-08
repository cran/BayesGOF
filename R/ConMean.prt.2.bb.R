ConMean.prt.2.bb <-
function(c.vec, weight, leg.mat, u, par1, par2){
	#### theta is G inv u = qbeta(u,par1, par2)
	inner.part <- qbeta(u, par1, par2)*leg.mat
	v <-(t(inner.part)%*% weight)/nrow(leg.mat)
	num.2 <- as.vector(v) %*% c.vec
	return(num.2)
}