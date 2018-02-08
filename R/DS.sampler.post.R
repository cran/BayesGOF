DS.sampler.post <-
function(k, g.par, LP.par, y.0, n.0 = NULL, con.prior = c("Normal", "Beta", "Gamma")){
	fam = match.arg(con.prior)
	switch(fam,
		"Normal" = {
			rDS.post.nnu(k, g.par, LP.par, y.i = y.0, se.i = n.0 )
		 },
		 "Beta" = {
			rDS.post.bbu(k, g.par, LP.par, y.i = y.0, n.i = n.0)
		 },
		 "Gamma" = {
			rDS.post.pgu(k, g.par, LP.par, x.i = y.0)
			 }
		)
	}