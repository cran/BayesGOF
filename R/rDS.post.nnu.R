rDS.post.nnu <-
function(k, g.par, LP.par, y.i, se.i){
# Generates samples from the posterior using U-function
# accept/reject sampling
# INPUTS
#  k				number of desired samples
#  g.par			vector consisting of parameters 
#					for g: c(alpha, beta)
#  LP.par			vector of LP coefficients: c(c_1,..,c_m)
#  y.i, se.i 		new posterior value
# OUTPUTS
#  x.samp			vector of samples with length k
	B = 250
	m <- length(LP.par)
	post.mu.i <- lambda.i(se.i, g.par[2]) * g.par[1] +
					(1-lambda.i(se.i, g.par[2]))* y.i
	post.tau2.i <- (1-lambda.i(se.i, g.par[2]))*se.i^2 #output is VARIANCE
	x.samp <- NULL
	for(i in 1:k){
		DG.val <- 0
		DU.val <- 1
		unit.vec <- seq(1/B, 1-1/B, length = B)
		d.u <- 1+LP.basis.beta(unit.vec, c(1,1), m)%*%LP.par 
		d.max <- max(d.u)
		while(DG.val < DU.val){
			y <- rnorm(1, post.mu.i, sd = sqrt(post.tau2.i)) 
			Leg.G <- LP.basis.norm(y, g.par, m)
			DG.val <- 1+Leg.G%*%LP.par 
			DU.val <- runif(1) * d.max 
			}
		x.samp[i] <- y
		}
	return(x.samp)
	}
