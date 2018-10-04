plot.DS_FB_obj <-
function(x, plot.type = "posterior",...){
	if(plot.type == "posterior"){
	par(mar=c(5,5,4,2)+0.3) #changes left margin to make large labels fit
	plot(x$post.fit$theta.vals, x$post.fit$parm.post, 
		type = "l", yaxs = "i", lwd = 2, col = "blue", lty = "dashed", 
		xlab = expression(theta), ylab = "", font.main = 1,
		cex.lab=1.45, cex.axis=1.5, cex.main=1.75, cex.sub=1.5, ...)
	title(ylab = expression(paste(hat(pi)(theta[i]~"|"~y[i]))), line = 2.3, cex.lab=1.45)
	lines(x$post.fit$theta.vals, x$post.fit$ds.post, col = "red3", lwd = 2)
	lines(x$post.fit$theta.vals, x$post.fit$finite.post, col = "darkgreen", lwd = 2)
	} else {
		par(mar=c(5,5,4,2)+0.3) #changes left margin to make large labels fit
		plot(x$prior.fit$theta.vals, x$prior.fit$parm.prior, 
			type = "l", yaxs = "i", lwd = 2, col = "blue", lty = "dashed", 
			xlab = expression(theta), ylab = "", font.main = 1,
			cex.lab=1.45, cex.axis=1.5, cex.main=1.75, cex.sub=1.5, ...)
		title(ylab = expression(paste(hat(pi)(theta[i]~"|"~y[i]))), line = 2.3, cex.lab=1.45)
		lines(x$prior.fit$theta.vals, x$prior.fit$ds.prior, col = "red3", lwd = 2)
		lines(x$prior.fit$theta.vals, x$prior.fit$finite.prior, col = "darkgreen", lwd = 2)
	}
}

