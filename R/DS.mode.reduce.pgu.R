DS.mode.reduce.pgu<- 
function(DS.GF.obj){
#####################################################
# INPUTS
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
#  ds.mode.vec		vector of modes for entire sample
	tbl <- table(DS.GF.obj$obs.data)
	cnt.vec <- as.integer(unlist(dimnames(tbl)))
	ds.mode.vec <- NULL
	for(i in 1:length(cnt.vec)){
		ds.mode.vec[i] <- DS.mode.map.pgu(cnt.vec[i], DS.GF.obj$g.par, DS.GF.obj$LP.par, B=100)
		}
	return(ds.mode.vec)
}