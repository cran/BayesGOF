Reduce.LP.coef <-
function(wght.mat, LP.par, u, leg.mat){
	L.mat <- matrix(0, nrow = length(LP.par), ncol = dim(wght.mat)[2])
	L.mat <- apply(wght.mat, 2, function(x) EXP.score(u, LP.par, x, leg.mat))
	if(!is.null(ncol(L.mat))){
	new.c.vec <- rowMeans(L.mat)
		} else {
		new.c.vec <- mean(L.mat)
		}
	return(new.c.vec)
}