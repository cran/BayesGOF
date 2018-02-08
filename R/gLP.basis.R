gLP.basis <-
function(x, g.par, m, ind = NULL, con.prior = c("Normal", "Beta", "Gamma")){
#######################################
## INPUTS
##	x			set of values OR single value between (0,1)
##  g.par		Parameters for the beta parametric prior G
##  m			selected m value
##  ind			Extracts the jth column of matrix; ow entire matrix 
##  con.prior	Normal, Beta, or Gamma basis
##	
## OUTPUTS
##  TY      functional values for the first m legendre polynomials
##			evaluated over x
#######################################
	fam = match.arg(con.prior)
	switch(fam,
		"Normal" = {
			LP.basis.norm(x, g.par, m, ind)
		 },
		 "Beta" = {
			LP.basis.beta(x, g.par, m, ind)
		 },
		 "Gamma" = {
			LP.basis.gamma(x, g.par, m, ind)
			 }
		)
	}