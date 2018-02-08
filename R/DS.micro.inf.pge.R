DS.micro.inf.pge <- 
function(DS.GF.obj, y.i, e.i){
#####################################################
# INPUTS
#  y.i				number of expected count
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
		fam <- "Poisson"
		out <- list()
		post.alph.i <- y.i + DS.GF.obj$g.par[1]
		post.beta.i <- DS.GF.obj$g.par[2]/(1+e.i*DS.GF.obj$g.par[2])
		##calculate denominator/standardization constant
		B=250
		u.int <- seq(1/B, 1-(1/B), length.out = B)
		theta.vals <- seq(0,max(DS.GF.obj$obs.data), length.out = B)
		eb.pos.den <- dgamma(theta.vals, post.alph.i, scale = post.beta.i)
		#determine LP density
		leg.mat.den <- LP.basis.beta(u.int, c(1,1), length(DS.GF.obj$LP.par))
		wght.den <- weight.fun.univ(u.int, DS.GF.obj$g.par[1], DS.GF.obj$g.par[2],
					post.alph.i, post.beta.i, family = fam)
		##Get non-standard LP adjusted density
		Leg.G <- LP.basis.gamma(theta.vals, c(DS.GF.obj$g.par[1], DS.GF.obj$g.par[2]), length(DS.GF.obj$LP.par))
		ds.pos <- eb.pos.den*(1+Leg.G%*%DS.GF.obj$LP.par) / EXP.denom(DS.GF.obj$LP.par, wght.den, leg.mat.den)
		out$post.fit <- data.frame(theta.vals = theta.vals,
									parm.pos = eb.pos.den,
									ds.pos = ds.pos)				
		out$post.fit[out$post.fit<0]<-0
		out$post.mean <- DS.PostMean.pge(y.i, e.i, DS.GF.obj$g.par, u.int, DS.GF.obj$LP.par)
		out$post.mode <- theta.vals[which.max(ds.pos)]
		out$cnt.val <- y.i
		class(out) <- "DS_GF_micro"
		return(out)
		}