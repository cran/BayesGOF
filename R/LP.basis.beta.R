LP.basis.beta <- 
function(x, g.par, m, ind = NULL){
#######################################
## INPUTS
##	x		set of values OR single value between (0,1)
##  g.par	Parameters for the beta parametric prior G
##  m		selected m value
##  j		Extracts the jth column of matrix; ow entire matrix 
##	
## OUTPUTS
##  TY      functional values for the first m legendre polynomials
##			evaluated over x
#######################################
	u <- pbeta(x, g.par[1], g.par[2])
	poly <-  slegendre.polynomials(m,normalized=TRUE)  
	TY <- matrix(NA,length(u),m)
	for(j in 1:m) TY[,j] <- predict(poly[[j+1]],u)
	if(is.numeric(ind) == FALSE){
	return(TY)
	}else{
	return(TY[,ind])
	}
} 