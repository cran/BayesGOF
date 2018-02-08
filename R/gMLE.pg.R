gMLE.pg <-
function(cnt.vec, e.vec = NULL, start.par = c(1,1)){
	if(is.null(e.vec) == TRUE){ 
	e.vec <- rep(1,length(cnt.vec))
	} else {
	 e.vec <- e.vec}
	logL <- function(par, x, e) {
		f <- choose(x + par[1] - 1, x) * (par[2] / (e + par[2]))^par[1] *(e / (e+par[2]))^x 
		- sum(log(f))
		}
	suppressWarnings(par.g <- optim(start.par, logL, x = cnt.vec, e = e.vec, hessian = TRUE)$par)
	par.g[2] <- 1/par.g[2]
	return(par.g)
}