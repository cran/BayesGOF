DS.micro.inf.nnu <-
function(DS.GF.obj, y.i, se.i){
#####################################################
# INPUTS
#  y.i			mean of new study
#  se.i			standard error of new study
#  DS.GF.obj		AutoBayes-DataCorrect object
# OUTPUTS
		fam <- "Normal"
		out <- list()
		#find posterior parameters
		post.mu.i <- lambda.i(se.i, DS.GF.obj$g.par[2]) * DS.GF.obj$g.par[1] +
					(1-lambda.i(se.i, DS.GF.obj$g.par[2]))* y.i
		post.tau2.i <- (1-lambda.i(se.i, DS.GF.obj$g.par[2]))*se.i^2 #output is VARIANCE
		##determine parametric density
		B=250
		theta.vals <- seq(DS.GF.obj$g.par[1] - 3*sqrt(DS.GF.obj$g.par[2]),
						DS.GF.obj$g.par[1] + 3*sqrt(DS.GF.obj$g.par[2]), length.out = 1000)
		u.int <- seq(1/B, 1-(1/B), length.out = B)
		eb.pos.den <- dnorm(theta.vals, post.mu.i, sd = sqrt(post.tau2.i))
		##determine LP density
		#Leg(u) for weight function
		leg.mat.den <- LP.basis.beta(u.int, c(1,1), length(DS.GF.obj$LP.par))
		wght.den <- weight.fun.univ(u.int, DS.GF.obj$g.par[1], DS.GF.obj$g.par[2],
					post.mu.i, post.tau2.i, family = fam)
		Leg.G <- LP.basis.norm(theta.vals, DS.GF.obj$g.par, length(DS.GF.obj$LP.par))
		ds.pos <- eb.pos.den*(1+Leg.G%*%DS.GF.obj$LP.par) / EXP.denom(DS.GF.obj$LP.par, wght.den, leg.mat.den)
		out$post.fit <- data.frame(theta.vals = theta.vals,
									parm.pos = eb.pos.den,
									ds.pos = ds.pos)				
		out$post.fit$ds.pos[out$post.fit$ds.pos<0]<-0
		out$post.mean <- DS.PostMean.nnu(y.i, se.i, DS.GF.obj$g.par, u.int, DS.GF.obj$LP.par)
		out$post.mode <- DS.mode.map.nnu(y.i, se.i, DS.GF.obj$g.par, DS.GF.obj$LP.par, B)
		out$study <- c(y.i, se.i)
		class(out) <- "DS_GF_micro"
		return(out)
		}