DS.mode.reduce.nnu <-
function(DS.GF.obj){
#####################################################
# INPUTS
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
#  dc.mode.vec		vector of modes for entire sample
	y.vec <- DS.GF.obj$obs.data[,1]
	se.vec <- DS.GF.obj$obs.data[,2]
	if(sum(DS.GF.obj$LP.par^2) == 0){
		ds.mode.vec <- lambda.i(se.vec, DS.GF.obj$g.par[2]) * DS.GF.obj$g.par[1] +
					(1-lambda.i(se.vec, DS.GF.obj$g.par[2]))* y.vec
		} else {
		B <- dim(DS.GF.obj$prior.fit)[1]
		ds.mode.vec <- NULL
		for(i in 1:length(y.vec)){
			ds.mode.vec[i] <- DS.mode.map.nnu(y.vec[i], se.vec[i], DS.GF.obj$g.par, DS.GF.obj$LP.par, B)
			}
		}
	return(ds.mode.vec)
}