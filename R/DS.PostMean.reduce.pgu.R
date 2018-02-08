DS.PostMean.reduce.pgu <-
function(DS.GF.obj){
#####################################################
# INPUTS
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
#  ConMean.vec		vector of conditional mean of entire sample
	tbl <- table(DS.GF.obj$obs.data)
	cnt.vec <- as.integer(unlist(dimnames(tbl)))
	if(sum(DS.GF.obj$LP.par^2) == 0){
		ConMean.vec <- (DS.GF.obj$g.par[1]+cnt.vec)*(DS.GF.obj$g.par[2])/(1+DS.GF.obj$g.par[2])
		} else {
		B <- dim(DS.GF.obj$prior.fit)[1]
		u <- seq(1/B,1-(1/B), length.out = B) ##unit interval, 0 to 1
		ConMean.vec <- NULL
		ConMean.vec <- sapply(cnt.vec, function(x) DS.PostMean.pgu(x, DS.GF.obj$g.par, u, DS.GF.obj$LP.par))
	}
	return(ConMean.vec)
}
