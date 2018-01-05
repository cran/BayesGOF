plot.DS_GF <-
function(x, plot.type = c("Ufunc","mDev","DSg"), ...) {
	method = match.arg(plot.type)
	switch(method,
		"Ufunc" = {
			plot(x$UF.dat$UF.x, x$UF.dat$UF.y, 
				xlim = c(0,1),
				#ylim = c(0,2.5), 
				main = "",
				type = "l",
				lwd = 2,
				xlab = "", ylab = "", font.main = 1,
				cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25,...)
			abline(h = 1, col = "red", lwd = 2, cex = 3, lty = "dashed")
			title(ylab = expression(paste(hat(d))), line = 2.3, cex.lab=1.25)
			},
		"mDev" = {
			if(x$sm.crit == "AIC"){
			plot(x$dev.df$m, x$dev.df$dev, 
				main = "", 
				type = "l", lwd = 2, col = "darkorchid",
				xlab = "m", ylab = "Deviance: AIC", font.main = 1,
				cex.lab=1.5, cex.axis=2, cex.main=2, cex.sub=1.5,...)
				} else {
			plot(x$dev.df$m, x$dev.df$dev, 
				main = "",
				type = "l", lwd = 2, col = "darkorchid",
				xlab = "m", ylab = "Deviance: BIC", font.main = 1,
				cex.lab=1.5, cex.axis=2, cex.main=2, cex.sub=1.5,...)
				}
			},
		"DSg" = {
			if(sum(x$LP.par^2)==0){
				plot(x$prior.data$theta.vals, x$prior.data$parm.prior, 
					main = "",
					type = "l", lwd = 2, col = "blue", lty = "dashed",
					xlab = "", #expression(theta) 
					ylab = "", font.main = 1,
					cex.lab=1.25, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ...)
				title(ylab = expression(paste(hat(pi)(theta))), line = 2.3, cex.lab=1.25)
				legend("bottom", legend = "g",
					xpd = TRUE, horiz = TRUE, inset = c(0, -.18),
					col = "blue", lty = "dashed", 
					lwd = 2, bty = "n") 
			} else {
				x.name <- deparse(substitute(x))
				plot(x$prior.data$theta.vals, x$prior.data$ds.prior, 
					main = "",
					type = "l", lwd = 2, col = "red3",
					xlab = "", ylab = "", font.main = 1,
					cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25,...)
				title(ylab = expression(paste(hat(pi)(theta))), line = 2.3, cex.lab=1.25)
				lines(x$prior.data$theta.vals, x$prior.data$parm.prior, 
					  col = "blue", lwd = 2, cex = 3, lty = "dashed")
				legend("bottom", legend = c(expression(hat(pi)),"g"),
						xpd = TRUE, horiz = TRUE, inset = c(0, -.18),
						col = c("red","blue"), lty = c("solid","dashed"), 
						lwd = 1, bty = "n")
				}
			}
	)
}