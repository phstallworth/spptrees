library(spatstat)
library(spptrees)

lambdaSN_1 <- function(x, y, shot){
  n <- shot$n
  shotx <- shot$x
  shoty <- shot$y
  mySum <- rep(0, length(x))
  for(j in 1:length(mySum)){
    for(i in 1:n){
      mySum[j] <- mySum[j] + 5 * twod_dnorm(c(x[j], y[j]), c(shotx[i], shoty[i]), sigma = matrix(c(0.0005, 0, 0, 0.0005), nrow = 2, byrow = T))
    }
  }
  mySum
}

rshotnoise <- function(shot, noise, window = owin()){
  init <- rpoispp(shot)
  full <- rpoispp(noise, win = window,  shot = init,)
  full
}

x <- rep(NA, length(y))
for(i in 1:101){
  x[(101*(i-1) + 1):(101*i)] <- y[i]
}
y <- rep(seq(0, 1, by = 0.01), 101)

init <- rpoispp(10)
sn <- rpoispp(lambdaSN_1, shot = init)
intens <- lambdaSN_1(x, y, init)
contour(matrix(intens, nrow = 101, byrow = T), col = "green")
points(sn, pch = 16, col = "hotpink")
points(init, col = "forestgreen", pch = 12)

sn <- rshotnoise(10, lambdaSN_1)
plot(sn)
