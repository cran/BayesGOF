DS.micro.inf.bbu <-
function(DS.GF.obj, y.i, n.i){
#####################################################
# INPUTS
#  y.i			number of successes in new study
#  n.i			number of total trials in new study
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
		fam <- "Binomial"
		out <- list()
		#find posterior parameters
		post.alph.i <- y.i + DS.GF.obj$g.par[1]
		post.beta.i <- n.i - y.i + DS.GF.obj$g.par[2]
		##determine parametric density
		B=250
		u.int <- seq(1/B, 1-(1/B), length.out = B)
		eb.pos.den <- dbeta(u.int, post.alph.i, post.beta.i)
		##determine LP density
		#Leg(u) for weight function
		leg.mat.den <- LP.basis.beta(u.int, c(1,1), length(DS.GF.obj$LP.par))
		wght.den <- weight.fun.univ(u.int, DS.GF.obj$g.par[1], DS.GF.obj$g.par[2],
					post.alph.i,post.beta.i, family = fam)
		Leg.G <- LP.basis.beta(u.int, c(DS.GF.obj$g.par[1], DS.GF.obj$g.par[2]), length(DS.GF.obj$LP.par))
		ds.pos <- eb.pos.den*(1+Leg.G%*%DS.GF.obj$LP.par) / EXP.denom(DS.GF.obj$LP.par, wght.den, leg.mat.den)
		out$post.fit <- data.frame(theta.vals = u.int,
									parm.pos = eb.pos.den,
									ds.pos = ds.pos)				
		out$post.fit[out$post.fit<0]<-0
		out$post.mean <- DS.PostMean.bbu(y.i, n.i, DS.GF.obj$g.par, u.int, DS.GF.obj$LP.par)
		out$post.mode <- DS.mode.map.bbu(y.i, n.i, DS.GF.obj$g.par, DS.GF.obj$LP.par, B)
		out$study <- c(y.i, n.i)
		class(out) <- "DS_GF_micro"
		return(out)
		}
