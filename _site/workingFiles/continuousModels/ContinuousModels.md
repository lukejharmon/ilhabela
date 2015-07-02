Fitting models of continuous character evolution
========================================================

Let's fit some models of continuous character evolution. First, we will learn how to do some tests of "phylogenetic signal," a very common test especially for ecological analyses. Then we will learn how to fit a series of evolutionary models for continuous characters.


```r
library(geiger)
```

```
## Loading required package: ape
```

```r
library(picante)
```

```
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.0-10
## Loading required package: nlme
```

```r
library(phytools)
```

```
## Loading required package: maps
## Loading required package: rgl
```

```
## Warning: failed to assign RegisteredNativeSymbol for getData to getData
## since getData is already defined in the 'phangorn' namespace
```


We will use the same anolis data and phylogenetic tree from previous exercises. If you don't already have them, you can download the files from the following addresses:

<a href="https://drive.google.com/file/d/0B9R4DAZPUvjiV2VhTUxOTlRuQUU/edit?usp=sharing">anolisDataAppended.csv</a> <br>
<a href="https://drive.google.com/file/d/0B9R4DAZPUvjiSkl1aFY2TkNMVFk/edit?usp=sharing">anolis.phy</a>

If you need to, make sure these files are in your working directory and read them in.


```r
anoleData <- read.csv("anolisDataAppended.csv", row.names = 1)
anoleTree <- read.tree("anolis.phy")
```


If you have the data, then the following commands should work:


```r
plot(anoleTree)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```r
name.check(anoleTree, anoleData)
```

```
## [1] "OK"
```


Let's do the two main tests for phylogenetic signal using anole body size. The first test is Blomberg's K, which compares the variance of PICs to what we would espect under a Brownian motion model. K = 1 means that relatives resemble one another as much as we should expect under BM; K < 1 means that there is less "phylogenetic signal" than expected under BM, while K > 1 means that there is more. A significant p-value returned from phylosignal tells you that there is significant phylogenetic signal - that is, close relatives are more similar than random pairs of species. 


```r
anoleSize <- anoleData[, 1]
names(anoleSize) <- rownames(anoleData)
phylosignal(anoleSize, anoleTree)
```

```
##       K PIC.variance.obs PIC.variance.rnd.mean PIC.variance.P
## 1 1.554           0.1389                 0.773          0.001
##   PIC.variance.Z
## 1         -3.913
```

```r
phylosig(anoleTree, anoleSize, method = "K", test = T)
```

```
## $K
## [1] 1.554
## 
## $P
## [1] 0.001
```


Another method for testing phylogenetic signal is Pagel's lambda. Lambda is a tree transformation that stretches tip branches relative to internal branches, making the tree more and more like a complete polytomy. If our estimated lambda = 0, then the traits are inferred to have no phylogenetic signal. Lambda = 1 corresponds to a Brownian motion model; 0 < lambda < 1 is in between.



```r
# First let's look at what lambda does
anoleTreeLambda0 <- rescale(anoleTree, model = "lambda", 0)
anoleTreeLambda5 <- rescale(anoleTree, model = "lambda", 0.5)

par(mfcol = c(1, 3))
plot(anoleTree)
plot(anoleTreeLambda5)
plot(anoleTreeLambda0)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

```r

phylosig(anoleTree, anoleSize, method = "lambda", test = T)
```

```
## $lambda
## [1] 1.017
## 
## $logL
## [1] -3.81
## 
## $logL0
## [1] -60.02
## 
## $P
## [1] 2.893e-26
```

```r

lambdaModel <- fitContinuous(anoleTree, anoleSize, model = "lambda")
```

```
## Loading required package: parallel
```

```
## Warning: Parameter estimates appear at bounds:
## 	lambda
```

```r
brownianModel <- fitContinuous(anoleTree, anoleSize)
nosigModel <- fitContinuous(anoleTreeLambda0, anoleSize)

lambdaModel$opt$aicc
```

```
## [1] 15.65
```

```r
brownianModel$opt$aicc
```

```
## [1] 13.52
```

```r
nosigModel$opt$aicc
```

```
## [1] 124.2
```

```r

# Conclusion: Brownian model is best, no signal model is terrible
```


We can use fitContinuous to fit OU and EB models as well.


```r
brownianModel <- fitContinuous(anoleTree, anoleSize)
OUModel <- fitContinuous(anoleTree, anoleSize, model = "OU")
EBModel <- fitContinuous(anoleTree, anoleSize, model = "EB")

# inspect results
brownianModel
```

```
## GEIGER-fitted comparative model of continuous data
##  fitted 'BM' model parameters:
## 	sigsq = 0.136160
## 	z0 = 4.065918
## 
##  model summary:
## 	log-likelihood = -4.700404
## 	AIC = 13.400807
## 	AICc = 13.524519
## 	free parameters = 2
## 
## Convergence diagnostics:
## 	optimization iterations = 100
## 	failed iterations = 0
## 	frequency of best fit = 1.00
## 
##  object summary:
## 	'lik' -- likelihood function
## 	'bnd' -- bounds for likelihood search
## 	'res' -- optimization iteration summary
## 	'opt' -- maximum likelihood parameter estimates
```

```r
OUModel
```

```
## GEIGER-fitted comparative model of continuous data
##  fitted 'OU' model parameters:
## 	alpha = 0.000000
## 	sigsq = 0.136160
## 	z0 = 4.065918
## 
##  model summary:
## 	log-likelihood = -4.700404
## 	AIC = 15.400807
## 	AICc = 15.650807
## 	free parameters = 3
## 
## Convergence diagnostics:
## 	optimization iterations = 100
## 	failed iterations = 0
## 	frequency of best fit = 0.76
## 
##  object summary:
## 	'lik' -- likelihood function
## 	'bnd' -- bounds for likelihood search
## 	'res' -- optimization iteration summary
## 	'opt' -- maximum likelihood parameter estimates
```

```r
EBModel
```

```
## GEIGER-fitted comparative model of continuous data
##  fitted 'EB' model parameters:
## 	a = -0.736271
## 	sigsq = 0.233528
## 	z0 = 4.066519
## 
##  model summary:
## 	log-likelihood = -4.285970
## 	AIC = 14.571939
## 	AICc = 14.821939
## 	free parameters = 3
## 
## Convergence diagnostics:
## 	optimization iterations = 100
## 	failed iterations = 0
## 	frequency of best fit = 0.41
## 
##  object summary:
## 	'lik' -- likelihood function
## 	'bnd' -- bounds for likelihood search
## 	'res' -- optimization iteration summary
## 	'opt' -- maximum likelihood parameter estimates
```

```r

# calculate AIC weights
bmAICC <- brownianModel$opt$aicc
ouAICC <- OUModel$opt$aicc
ebAICC <- EBModel$opt$aicc

aicc <- c(bmAICC, ouAICC, ebAICC)
aiccD <- aicc - min(aicc)
aw <- exp(-0.5 * aiccD)
aiccW <- aw/sum(aw)
aiccW
```

```
## [1] 0.5353 0.1849 0.2798
```


It is important to realize that measurement error can bias your inferences with fitting these models towards OU. Fortunately, we can easily account for that in fitContinuous.


```r

# We measured 20 anoles per species, and the standard deviation within each
# species was, on average, 0.05
seSize <- 0.05/sqrt(20)

# redo with measurement error
brownianModel <- fitContinuous(anoleTree, anoleSize, SE = seSize)
OUModel <- fitContinuous(anoleTree, anoleSize, model = "OU", SE = seSize)
EBModel <- fitContinuous(anoleTree, anoleSize, model = "EB", SE = seSize)


# calculate AIC weights
bmAICC <- brownianModel$opt$aicc
ouAICC <- OUModel$opt$aicc
ebAICC <- EBModel$opt$aicc

aicc <- c(bmAICC, ouAICC, ebAICC)
aiccD <- aicc - min(aicc)
aw <- exp(-0.5 * aiccD)
aiccW <- aw/sum(aw)
aiccW
```

```
## [1] 0.5346 0.1846 0.2808
```


