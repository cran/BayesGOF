rDS.posterior <-
function(k, g.par, LP.par, y.i, n.i){
# Generates samples from the prior using U-function
# accept/reject sampling
# INPUTS
#  k				number of desired samples
#  g.par			vector consisting of parameters 
#					for g: c(alpha, beta)
#  LP.par			vector of LP coefficients: c(c_1,..,c_m)
# OUTPUTS
#  x.samp			vector of samples with length k
	B = 250
	m <- length(LP.par)
	post.alpha <- g.par[1] + y.i
	post.beta <- n.i - y.i + g.par[2]
	x.samp <- NULL
	cd.samp <- NULL
	for(i in 1:k){
		DG.val <- 0
		DU.val <- 1
		unit.vec <- seq(1/B, 1-1/B, length = B)
		d.u <- 1+LP.basis.beta(unit.vec, c(1,1), m)%*%LP.par 
		d.max <- max(d.u)
		while(DG.val < DU.val){
			y <- rbeta(1,post.alpha,post.beta) 
			Leg.G <- LP.basis.beta(y, g.par, m)
			DG.val <- 1+Leg.G%*%LP.par 
			DU.val <- runif(1) * d.max 
			}
		x.samp[i] <- y
		}
	return(x.samp)
	}
