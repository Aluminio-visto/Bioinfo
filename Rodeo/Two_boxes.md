---
title: "Rodeo - Plasmid Conjugation"
author: "Jorge RG"
date: "6/7/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Plasmid spreading model

Let's try the 2 boxes basic rodeo model:

```{r}
library(deSolve)
library(rodeo)
data(vars, pars, pros, funs, stoi)

```

## Setting the model

Let's tune the model
```{r}
nBox <- 2    #We've got two connected boxes with media flowing from left to right one.

pars <- rbind(pars, c(name="d", unit="1/hour", description="diffusion parameter"))
pros <- rbind(pros, c(name="diffSub", unit="mg/ml/hour",description="diffusion of substrate", 
                      expression="d * ((left(sub)-sub) + (right(sub)-sub))"))
stoi <- rbind(stoi, c(variable="sub", process="diffSub", expression="1"))

model <- rodeo$new(vars=vars, pars=pars, funs=funs, pros=pros, stoi=stoi, dim=c(nBox))
model$compile(fortran=FALSE)

monod <- function(c, h) { c / (c + h) }      # For limited nutrient availability
rp <- function (x) {rep(x, nBox)}            # For convenient replication

v <- cbind(bac=rp(0.01), sub=rp(0))
model$setVars(v)

p <- cbind(mu=rp(0.8), half=rp(0.1), yield= rp(0.1), vol=c(300, 1000), flow=rp(50), 
           sub_in=rp(1), d=rp(.75))          # Added diffusion (L->R) parameter
model$setPars(p)

```

## Solving and Plotting

We use the simple solver, not FORTRAN, as the model is quite simple:

```{r}
out <- model$dynamics(times=0:120, fortran=FALSE)
```


```{r}
layout(matrix(1:model$lenVars(), nrow=1))
for (vn in model$namesVars()) {
  matplot(out[,"time"], out[,paste(vn, 1:nBox, sep=".")], type="l", xlab="time", ylab=vn, lty=1:nBox, col=1:nBox)
  legend("right", bty="n", lty=1:nBox, col=1:nBox, legend=paste("box",1:nBox))
}

layout(1)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

```{r}

model$plotStoichiometry(box=2, time=0, cex=0.8)

```


