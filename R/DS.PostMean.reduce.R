DS.PostMean.reduce <-
function(DS.GF.obj){
	fam = DS.GF.obj$fam
	switch(fam,
		"Normal" = {
			DS.PostMean.reduce.nnu(DS.GF.obj)
		 },
		 "Binomial" = {
			DS.PostMean.reduce.bbu(DS.GF.obj)
		 },
		 "Poisson" = {
			DS.PostMean.reduce.pgu(DS.GF.obj)				
			 }
		)
	}
