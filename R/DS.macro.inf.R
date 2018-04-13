DS.macro.inf <-
function(DS.GF.obj, num.modes =1, method = c("mean","mode"),
						iters = 25, e.vec = NULL){
	fam = DS.GF.obj$fam
	switch(fam,
		"Normal" = {
			DS.macro.inf.nnu(DS.GF.obj, num.modes , iters , method)
		 },
		 "Binomial" = {
			DS.macro.inf.bbu(DS.GF.obj, num.modes , iters , method)
		 },
		 "Poisson" = {
			if(is.null(e.vec) == TRUE){
				DS.macro.inf.pgu(DS.GF.obj, num.modes , iters , method)
				} else {
				DS.macro.inf.pge(DS.GF.obj, num.modes, iters, e.vec = e.vec, method)
				}
			 }
		)
	}