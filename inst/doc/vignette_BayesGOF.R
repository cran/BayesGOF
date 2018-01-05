## ---- message  = FALSE---------------------------------------------------
library(BayesGOF)
set.seed(8697)
data(rat)
###Use MLE to determine starting values
rat.start <- BetaBinoMLE(rat$y, rat$n)$estimate

## ------------------------------------------------------------------------
rat.ds <- DS.prior(rat, max.m = 6, rat.start)

## ---- fig.height = 5.5, fig.width = 5.5, fig.align = 'center'------------
plot(rat.ds, plot.type = "Ufunc")

## ------------------------------------------------------------------------
rat.ds

## ----fig.height = 5.5, fig.width = 5.5,fig.align = 'center'--------------
plot(rat.ds, plot.type = "DSg")

## ---- warning = FALSE----------------------------------------------------
rat.macro.md <- DS.macro.inf(rat.ds, num.modes = 2 , iters = 25, method = "mode") 

## ------------------------------------------------------------------------
rat.macro.md

## ---- fig.height = 5.5, fig.width = 5.5,fig.align = 'center'-------------

plot(rat.macro.md)

## ---- fig.align = 'center'-----------------------------------------------
rat.y71.micro <- DS.micro.inf(rat.ds, y.i = 4, n.i = 14)
rat.y71.micro

## ----fig.height = 5.5, fig.width = 5.5,fig.align = 'center'--------------
plot(rat.y71.micro)

