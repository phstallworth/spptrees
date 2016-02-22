library(spatstat)
library(spptrees)

lambdaSN_1 <- function(x, y, shot){
  n <- shot$n
  shotx <- shot$x
  shoty <- shot$y
  mySum <- rep(0, length(x))
  for(j in 1:length(mySum)){
    for(i in 1:n){
      mySum[j] <- mySum[j] + twod_dnorm(c(x[j], y[j]), c(shotx[i], shoty[i]), sigma = matrix(c(0.05, 0, 0, 0.05), nrow = 2, byrow = T) * 40)
    }
  }
  mySum
}
init <- rpoispp(5)
sn <- rpoispp(lambdaSN_1, shot = init)
plot(sn)
points(init, col = "red")

