DS.mode.reduce <-
function(DS.GF.obj){
#####################################################
# INPUTS
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
#  dc.mode.vec		vector of modes for entire sample
	y.vec <- DS.GF.obj$obs.data[,1]
	n.vec <- DS.GF.obj$obs.data[,2]
	if(sum(DS.GF.obj$LP.par^2) == 0){
		dc.mode.vec <- ( (DS.GF.obj$g.par[1]+y.vec)-1)/
					   ( (DS.GF.obj$g.par[1]+y.vec) + (DS.GF.obj$g.par[2]+n.vec-y.vec)-2 )
		} else {
		B <- dim(DS.GF.obj$prior.data)[1]
		ds.mode.vec <- NULL
		for(i in 1:length(y.vec)){
			ds.mode.vec[i] <- DS.mode.map(y.vec[i], n.vec[i], DS.GF.obj$g.par, DS.GF.obj$LP.par, B)
			}
		}
	return(ds.mode.vec)
}
