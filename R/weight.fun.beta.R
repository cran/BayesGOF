weight.fun.beta <-
function(u, alpha_0, beta_0, alpha.post.i, beta.post.i){
	if(alpha.post.i <= 0 | beta.post.i <=0){
		ratio = rep(0,length(u))
		return(ratio)
		} else {
		v1 <- dbeta(qbeta(u,shape1 = alpha_0, shape2 = beta_0), 
				shape1 = alpha.post.i, shape2 = beta.post.i)
		v2 <- dbeta(qbeta(u,shape1 = alpha_0, shape2 = beta_0), 
				shape1 = alpha_0, shape2 = beta_0)
		ratio <- v1/v2
		return(ratio)
		}
	}
