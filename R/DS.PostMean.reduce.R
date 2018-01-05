DS.PostMean.reduce <-
function(DS.GF.obj){
#####################################################
# INPUTS
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
#  ConMean.vec		vector of conditional mean of entire sample
	y.i <- DS.GF.obj$obs.data[,1]
	n.i <- DS.GF.obj$obs.data[,2]
	if(sum(DS.GF.obj$LP.par^2) == 0){
		ConMean.vec <- (DS.GF.obj$g.par[1]+y.i)/(DS.GF.obj$g.par[1]+DS.GF.obj$g.par[2]+n.i)
		} else {
		B <- dim(DS.GF.obj$prior.data)[1]
		u <- seq(1/B,1-(1/B), length.out = B) ##unit interval, 0 to 1
		ConMean.vec <- NULL
		ConMean.vec <- apply(DS.GF.obj$obs.data, 1, function(x) DS.PostMean(x[1], x[2], DS.GF.obj$g.par, u, DS.GF.obj$LP.par))
	}
	return(ConMean.vec)
}
