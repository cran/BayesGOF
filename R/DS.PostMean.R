DS.PostMean <-
function(y.0, n.0 = NULL, g.par, LP.par, B = 250, 
		 family = c("Normal","Binomial", "Poisson"), e.0 = NULL){
	fam = match.arg(family)
	u <- seq(1/B, 1-1/B, length.out = B)
	switch(fam,
		"Normal" = {
			DS.PostMean.nnu(y.i = y.0, se.i=n.0, g.par = g.par, u=u,LP.par = LP.par)
		 },
		 "Binomial" = {
			DS.PostMean.bbu(y.i = y.0, n.i=n.0,g.par = g.par, u=u,LP.par = LP.par)
		 },
		 "Poisson" = {
			if(is.null(e.0) == TRUE){
				DS.PostMean.pgu(x.i = y.0, g.par = g.par , u = u, LP.par = LP.par)
				} else {
				DS.PostMean.pge(x.i = y.0, e.i = e.0, g.par = g.par, u = u, LP.par = LP.par)
				}
				}
		)
	}