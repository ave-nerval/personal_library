## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(profvis)
library(microbenchmark)
library(compiler)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# generiranje podatkov
set.seed(1234)
m <- 1000
n <- 100
X <- matrix(rnorm(m * n, mean = 1, sd = 1), nrow = n)
skupina <- rep(1:2, each = n/2)


## ----eval=F----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## # preverimo zapis funkcije
## stats:::t.test
## methods(t.test)
## stats:::t.test.formula
## stats:::t.test.default
## 
## # profiling
## profvis({
##   for(i in 1:m) t.test(X[ , i] ~ skupina)$stat
## })


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fun.start <- function(X, m, skupina){
  for(i in 1:m) t.test(X[ , i] ~ skupina)$stat
}

fun.fast <- function(X, skupina){
  n <- nrow(X)
  pod1 <- X[skupina==1,]
  pod2 <- X[skupina==2,]
  n12 <- n/2
  mean.1 <- colMeans(pod1)
  mean.2 <- colMeans(pod2)
  s2<-(colSums((pod1-rep(mean.1, each=n12))^2) + colSums((pod2-rep(mean.2, each=n12))^2))/(n-2)
  t.values<-(mean.1-mean.2)/(sqrt(s2*(1/n12+1/n12)))
  t.values
} # 23.3


fun.fast.compiled <- compiler::cmpfun(fun.fast)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# primerjava rezultatov
res.fun.old <- numeric(m) # rezerviramo prostor za rezultat
for(i in 1:m) res.fun.old[i]=t.test(X[ , i] ~ skupina)$stat
res.fun.fast <- fun.fast(X,skupina)

all.equal(res.fun.old, res.fun.fast)




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
microbenchmark(
  fun.start(X=X,m=m, skupina=skupina),
  fun.fast(X=X, skupina=skupina)
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# primerjava s system.time
system.time(for(i in 1:100){
  for(i in 1:m) t.test(X[ , i] ~ skupina)$stat
  })

system.time(for(i in 1:100) fun.fast(X,skupina))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
t_krit <- qt(0.975, df=n-2)
sum(abs(res.fun.fast) > t_krit)/m


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1234)
m <- 1000
n <- 100
X <- matrix(rbinom(m * n, size=1, prob=0.5), nrow = n)
skupina <- rep(1:2, each = n/2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# zacetna funkcija: chisq.test
fun.start.2 <- function(X,m, skupina) {
  for(i in 1:m) chisq.test(X[ , i], skupina)$stat
}


## ----eval=F----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## chisq.test
## profvis(for(i in 1:m) chisq.test(X[ , i], skupina)$stat)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# osnovna funkcija
res_slow_chi2 <- numeric(m)
for(i in 1:m) res_slow_chi2[i] = chisq.test(X[ , i], skupina, correct=F)$stat


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# hitrejsa funkcija 1
fun.fast.1 <- function(X, skupina, n){

  x <- table(X, skupina)
  sum.row <- rowSums(x)
  sum.col <- colSums(x)
  E <- outer(sum.row, sum.col, "*")/n
  stat <- sum((x - E)^2/E)
}

# hitrejsa funkcija 2
fun.fast.2 <- function(X, skupina, n){

  # priprava frekvenčne tabele
  t1 <- sum(X[skupina==1])
  t2 <- n/2-t1
  t3 <- sum(X)-t1
  t4 <- n/2-t3
  x <- matrix(data=c(t2,t1,t4,t3), nrow=2, ncol=2)
  
  # izračun pričakovanih frekvenc
  sum.row <- rowSums(x)
  sum.col <- colSums(x)
  E <- outer(sum.row, sum.col, "*")/n
  
  # izračun testne statistike
  stat <- sum((x - E)^2/E)
}

# hitrejsa funkcija 3
fun.fast.3 <- function(X, skupina, n){
  
  nn <- n/2
  # priprava frekvenčne tabele
  t1 <- sum(X[skupina==1])
  t2 <- nn-t1
  t3 <- sum(X)-t1
  t4 <- nn-t3
  x <- matrix(data=c(t2,t1,t4,t3), nrow=2, ncol=2)
  
  # izračun pričakovanih frekvenc
  sum.row <- rowSums(x)
  sum.col <- colSums(x)
  ee <- rbind((t2+t4)/2, (t1+t3)/2)
  E <- cbind(ee,ee)
  # izračun testne statistike
  stat <- sum((x - E)^2/E)
}

fun.fast.4 <- compiler::cmpfun(fun.fast.2)


X2 <- as.data.frame(X)
res_fast_chi2 <- vapply(X2, fun.fast.3, skupina=skupina, n=100, numeric(1))
res_fast_chi2 <- unname(res_fast_chi2)

all.equal(res_slow_chi2, res_fast_chi2)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# primerjava hitrosti
microbenchmark(
  fun.start.2(X,m, skupina),
  # apply(X, 2, fun.fast.1),
  # apply(X,2,fun.fast.2),
  #  for(i in 1:m) fun.fast.2(X[,i]),
  # vapply(lapply(seq_len(ncol(X)), function(i) X[,i]), fun.fast.2, numeric(1))
  vapply(X2, fun.fast.3, skupina=skupina, n=100, numeric(1))
  # vapply(as.data.frame(X), fun.fast.2, numeric(1))
)


# profvis({
#   vapply(X2, fun.fast.2, numeric(1))
# 
# })



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# primerjava hitrosti s system.time

# funkcija chi2 (osnovna funkcija)
system.time(for(i in 1:100) {for(i in 1:m) chisq.test(X[ , i], skupina)$stat})
# nova funkcija
system.time(for(i in 1:100) vapply(X2, fun.fast.3, skupina=skupina, n = 100, numeric(1)))




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
chi2_krit <- qchisq(0.95, df=1)

sum(res_fast_chi2 > chi2_krit)/m

