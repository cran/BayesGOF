DS.micro.inf <-
function(DS.GF.obj, y.i, n.i){
#####################################################
# INPUTS
#  y.0			number of successes in sample
#  n.0			number of total trials in sample
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
		out <- list()
		alpha_0 <- DS.GF.obj$g.par[1]
		beta_0 <- DS.GF.obj$g.par[2]
		c.vec <- DS.GF.obj$LP.par
		post.alph.i <- y.i + alpha_0
		post.beta.i <- n.i - y.i + beta_0
		##calculate denominator/standardization constant
		B=250
		u.int <- seq(1/B, 1-(1/B), length.out = B)
		leg.mat.den <- LP.basis.beta(u.int, c(1,1), length(c.vec))
		wght.den <- weight.fun.beta(u.int,alpha_0,beta_0,
					post.alph.i,post.beta.i)
		##Get non-standard LP adjusted density
		eb.pos.den <- dbeta(u.int, post.alph.i, post.beta.i)
		#u.theta <- pbeta(u.int, alpha_0, beta_0)
		Leg.G <- LP.basis.beta(u.int, c(alpha_0, beta_0), length(c.vec))
		nstd.dc.pos <- eb.pos.den*(1+Leg.G%*%c.vec)
		denom <- EXP.denom(c.vec, wght.den, leg.mat.den)
		ds.pos <- nstd.dc.pos/denom
		out$post.data <- data.frame(theta.vals = u.int,
									parm.pos = eb.pos.den,
									ds.pos = ds.pos)				
		out$post.data[out$post.data<0]<-0
		out$post.mean <- DS.PostMean(y.i, n.i, DS.GF.obj$g.par, u.int, DS.GF.obj$LP.par, method = "Theta")
		out$post.mode <- DS.mode.map(y.i, n.i, DS.GF.obj$g.par, DS.GF.obj$LP.par, B)
		out$study <- c(y.i, n.i)
		class(out) <- "DS_GF_micro"
		return(out)
		}
