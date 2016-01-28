library(spatstat)
poissonProcess <- rpoispp(6, win = owin(c(0, 10), c(0, 10)))
thomasProcess <- rThomas(kappa = 2, scale = 0.2, mu = 3, win = owin(c(0, 10), c(0, 10)))
superProcess <- superimpose(poissonProcess, thomasProcess)
plot(superProcess)

genRetProbs <- function(X, rho){
  rho/density(X)
}

genRetProbs(superProcess, 6)
letsTry <- rthin(X = superProcess, genRetProbs(X = superProcess, 6))
plot(density(letsTry))

letsTryL <- Lest(letsTry, correction = "border")
myEnvelope <- envelope(letsTry, Lest, nsim = 1000)

plot(x = 0, y = 0, xlim = c(0, 2.5), ylim = c(min(myEnvelope$lo - myEnvelope$r),
                                              max(myEnvelope$hi - myEnvelope$r)), type = "n")
lines(letsTryL$bord - letsTryL$r ~ letsTryL$r, type = 'l', col = "tomato")
lines(myEnvelope$lo - myEnvelope$r ~ myEnvelope$r, lty =2, col = "steelblue")
lines(myEnvelope$hi - myEnvelope$r ~ myEnvelope$r, lty = 2, col = "steelblue")

##

fit <- kppm(superProcess, ~1, "Thomas")

plot(envelope(fit, Lest, nsim = 39))

#If we thin again, do things get better? Answer: Nope, we just get into
#a different issue. This is because then simple retention probabilities generator
#assumes a constant rho. However, the thinning probabilities get smaller than 6
#I might be able to fix this by doing a max between the two things. I'll test that
#out some vectorisation stuff in the console beforehand.
letsTry <- rthin(X = letsTry, genRetProbs(X = letsTry, 6))
plot(density(letsTry))

letsTryL <- Lest(letsTry, correction = "border")
myEnvelope <- envelope(letsTry, Lest, nsim = 1000)

plot(letsTryL$bord - letsTryL$r ~ letsTryL$r, type = 'l')
lines(myEnvelope$lo - myEnvelope$r ~ myEnvelope$r, lty =2)
lines(myEnvelope$hi - myEnvelope$r ~ myEnvelope$r, lty = 2)
