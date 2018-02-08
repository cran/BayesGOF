DS.PostMode <-
function(y.0, n.0 = NULL, g.par, LP.par, B =250, 
		 family = c("Normal","Binomial", "Poisson")){
	fam = match.arg(family)
	switch(fam,
		"Normal" = {
			DS.mode.map.nnu(y.i = y.0, se.i = n.0, g.par = g.par, LP.par = LP.par, B = B)
		 },
		 "Binomial" = {
			DS.mode.map.bbu(y.i = y.0, n.i=n.0, g.par = g.par, LP.par = LP.par, B = B)
		 },
		 "Poisson" = {
			DS.mode.map.pgu(y.i = y.0, g.par = g.par , LP.par = LP.par, B = B)
			}
			 
		)
	}
