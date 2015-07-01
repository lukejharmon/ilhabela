Anole Solution
========================================================

This is the solution to the anole BiSSe challenge.


```r

library(ape)
library(TreeSim)
```

```
## Loading required package: geiger
```

```r
library(diversitree)
```

```
## Loading required package: deSolve
## Loading required package: subplex
## Loading required package: Rcpp
```

```r
library(laser)
```

```
## Error: there is no package called 'laser'
```

```r

anoleData <- read.csv("~/Documents/teaching/revellClass/2014bogota/continuousModels/anolisDataAppended.csv", 
    row.names = 1)
anoleTree <- read.tree("~/Documents/teaching/revellClass/2014bogota/continuousModels/anolis.phy")

ecomorph <- anoleData[, "ecomorph"]
names(ecomorph) <- rownames(anoleData)

isTG <- as.numeric(ecomorph == "TG")
names(isTG) <- names(ecomorph)

bisseModel <- make.bisse(anoleTree, isTG)
p <- starting.point.bisse(anoleTree)
bisseMLFit <- find.mle(bisseModel, p)

# we can test a constrained model where the character does not affect
# diversification rates

cBisseModel <- constrain(bisseModel, lambda1 ~ lambda0)
cBisseModel <- constrain(cBisseModel, mu1 ~ mu0)

cbMLFit <- find.mle(cBisseModel, p[c(-1, -3)])

# compare models

anova(bisseMLFit, constrained = cbMLFit)
```

```
##             Df lnLik  AIC   ChiSq Pr(>|Chi|)
## full         6   -41 94.1                   
## constrained  4   -41 90.1 0.00108          1
```

```r

prior <- make.prior.exponential(1/(2 * 2))

mcmcRun <- mcmc(bisseModel, bisseMLFit$par, nsteps = 1000, prior = prior, w = 0.1, 
    print.every = 100)
```

```
## 100: {2.5448, 2.7423, 0.1223, 0.7930, 0.2523, 0.0846} -> -53.35378
## 200: {1.8292, 3.3961, 0.0570, 1.1785, 0.0314, 1.2488} -> -54.48636
## 300: {2.1231, 2.7016, 0.0070, 0.8293, 0.0369, 1.3421} -> -54.28697
## 400: {2.8111, 2.0097, 0.7680, 0.0248, 0.0545, 1.4919} -> -59.25650
## 500: {1.9839, 2.5681, 0.0472, 0.0503, 0.1028, 0.1318} -> -52.39312
## 600: {2.6010, 2.7925, 0.4178, 0.8012, 0.1168, 0.1300} -> -54.37971
## 700: {2.1467, 4.1820, 0.0587, 0.3723, 0.3015, 0.3556} -> -57.10820
## 800: {2.1143, 2.2737, 0.0277, 0.9679, 0.0180, 0.7823} -> -55.91828
## 900: {2.4259, 3.0204, 0.3192, 2.1095, 0.0069, 0.6079} -> -58.84477
## 1000: {2.3074, 2.4391, 0.0319, 0.0112, 0.1415, 0.0742} -> -50.98364
```

```r

col <- c("blue", "red")
profiles.plot(mcmcRun[, c("lambda0", "lambda1")], col.line = col, las = 1, legend = "topright")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-11.png) 

```r
profiles.plot(mcmcRun[, c("mu0", "mu1")], col.line = col, las = 1, legend = "topright")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-12.png) 

```r

```

