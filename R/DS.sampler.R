DS.sampler <-
function(k, g.par, LP.par, con.prior = c("Normal", "Beta", "Gamma")){
	fam = match.arg(con.prior)
	switch(fam,
		"Normal" = {
			rDS.nnu(k, g.par, LP.par)
		 },
		 "Beta" = {
			rDS.bbu(k, g.par, LP.par)
		 },
		 "Gamma" = {
			rDS.pgu(k, g.par, LP.par)
			 }
		)
	}