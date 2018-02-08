plot.DS_GF_micro <-
function(x, ...){
	plot(x$post.fit$theta.vals, x$post.fit$ds.pos,  
		type = "l", lwd = 2, col = "red3",
		xlab = expression(theta), ylab = "", font.main = 1,
		cex.lab=1.45, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ...)
	title(ylab = expression(paste(hat(pi)(theta[i]~"|"~y[i]))), line = 2.3, cex.lab=1.45)
	lines(x$post.fit$theta.vals, x$post.fit$parm.pos,
		col = "blue", lty = "dashed", lwd = 2)
	}
