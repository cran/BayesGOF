DS.PostMode.reduce <-
function(DS.GF.obj){
	fam = DS.GF.obj$fam
	switch(fam,
		"Normal" = {
			DS.mode.reduce.nnu(DS.GF.obj)
		 },
		 "Binomial" = {
			DS.mode.reduce.bbu(DS.GF.obj)
		 },
		 "Poisson" = {
			DS.mode.reduce.pgu(DS.GF.obj)				
			 }
		)
	}