plot.DS_GF_micro <-
function(x, x.lim = c(0,1), ...){
	plot(x$post.data$theta.vals, x$post.data$ds.pos, xlim = x.lim, 
		main = "Posterior Density Estimate", type = "l", lwd = 2, col = "red3",
		xlab = expression(theta), ylab = "", font.main = 1,
		cex.lab=1.45, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
	title(ylab = expression(paste(hat(pi)(theta[i]~"|"~y[i]))), line = 2.3, cex.lab=1.45)
	lines(x$post.data$theta.vals, x$post.data$parm.pos,
		col = "blue", lty = "dashed", lwd = 2)
	}

