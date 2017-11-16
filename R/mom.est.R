mom.est <-
function(success, trials){
						start.mu <- mean(success/trials)
						start.var <- var(success/trials)
						start.a <- ((1 - start.mu) / start.var - 1 / start.mu) * start.mu ^ 2
						start.b <- start.a * (1 / start.mu - 1)
						start.vals <- c(start.a, start.b)
					return(start.vals)
	}
