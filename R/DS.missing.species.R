DS.missing.species <- 
function(DS.GF.obj, t.ratio, find.se = FALSE, se.reps = 50){
	#Expectation functions
	MS.expect.new <- function(theta,t.val){
				exp(-theta)*(1-exp(-theta*t.val))/(1-exp(-theta))
				}
	species.val <- function(DS.object, t.ratio){
		if(DS.object$fam != "Poisson"){
			print("Not a DS prior from Poisson family.")
		} else {
			u.int <- DS.object$UF.data$UF.x
			new.theta <- qgamma(u.int, shape = DS.object$g.par[1], scale = DS.object$g.par[2])
			exp.new.theta <- MS.expect.new(new.theta, t.ratio)
			ms.new <- length(DS.object$obs.data)*sintegral(u.int, exp.new.theta*DS.object$UF.data$UF.y)$int
			return(ms.new)
			}
		}
	out <- list()
	#Find value
	out$ms.val <- species.val(DS.GF.obj, t.ratio)
	out$t.ratio <- t.ratio
	#Find standard error
	if(find.se == FALSE){
		return(out)
		} else {
		val.vec <- NULL
		for(i in 1:se.reps){
			par.g <- c(NA,NA)
			while(!is.finite(par.g[1])==TRUE | !is.finite(par.g[2])==TRUE){
				y.new <- rPPD.ds(DS.GF.obj,1, pred.type = "prior")$first.set
				possibleError <- tryCatch(par.g <- gMLE.pg(y.new, start.par = DS.GF.obj$g.par),
									  error = function(e) e )
				er.check <- !inherits(possibleError, "error")
				if(er.check == TRUE){
					par.g <- gMLE.pg(y.new, start.par = DS.GF.obj$g.par)
					} else {
					par.g <- c(NA, NA)
					}
				}
			new.LPc <- DS.prior.pgu(y.new, max.m = DS.GF.obj$m.val, 
								start.par = par.g)
			val.vec[i] <- species.val(new.LPc, t.ratio)
			}
		out$boot.ms.val <- val.vec
		out$boot.se <- sd(val.vec)
		return(out)
		}
}