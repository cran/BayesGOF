DS.PostMean.reduce.nnu <-
function(DS.GF.obj){
#####################################################
# INPUTS
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
#  ConMean.vec		vector of conditional mean of entire sample
	y.vec <- DS.GF.obj$obs.data[,1]
	se.vec <- DS.GF.obj$obs.data[,2]
	if(sum(DS.GF.obj$LP.par^2) == 0){
		ConMean.vec <- lambda.i(se.vec, DS.GF.obj$g.par[2]) * DS.GF.obj$g.par[1] +
					(1-lambda.i(se.vec, DS.GF.obj$g.par[2]))* y.vec
		} else {
		B <- dim(DS.GF.obj$prior.fit)[1]
		u <- seq(1/B,1-(1/B), length.out = B) ##unit interval, 0 to 1
		ConMean.vec <- NULL
		ConMean.vec <- apply(DS.GF.obj$obs.data, 1, function(x) DS.PostMean.nnu(x[1], x[2], DS.GF.obj$g.par, u, DS.GF.obj$LP.par))
	}
	return(ConMean.vec)
}
