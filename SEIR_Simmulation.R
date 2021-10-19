## SEIR Model

## Overview
## ********
## The SEIR model is a model to recreate the trajectory of a epidemic, in particular following the path of infections for 3 cohorts. 
## The whole population, the 'cautious' 10% of the population defined as 10% of the population with the lowest beta values and a 
## random sample of 0.1% of the population. Defining 4 states, susceptible, exposed, infected and recovered or dead on each day 
## an individual moves from suseptible state to the exposed state a function of their contact rate, the infected populations 
##contact rate and a viral infectivity parameter lamda. Each day an individual can move from being exposed to infected with 
## the disease with probability 1/3 and move from infected to recovered or dead with probability 1/5. Using the SEIR function, 
## alongside an initial run this program aims to plot the variability in these cohorts by running the model 10 separate times
## and plotting the variability on three seperate graphs. 

## SEIR model details
## ******************
## 1. Vectors are initiated whereby 'x' is defined as the storage vector for the whole population each element can either be 
##    0 = Susecptible ## 1 = Exposed ## 2 = Infected ## 3 = Recovered/Dead 
## 2. Vectors are initiated whereby 'S', 'E', 'I', 'R' are vectors showing how many people are in each state each day
## 3. Beta's are assigned to each individual (using rnom) and stored in 'beta' the bottom 10% is identified and stored
## 4. A random saample of 0.1% is also taken from the population for use later (using sample)
## 5. Model starts by using uniform random devidates (runif) to move people from state E -> I and E -> R
## 5. A person, i,  in state S are assigned a probabilty of moving to state E by sum_infectious * lamda * beta(i) and 
##    moved to state S based on random deviates (runif)
## 6. Three count vectors are set up, new_infections, new_infectios_percentile, new_infections_sample and the number of people 
##    moving from E -> I in each group is recorded for that dat
## 7. Model runs for 150 days
## 8. Plots are created using plot, lines, mtext, ablines and legend functions to plot the trajectory of all three samples

seir <-function(n=5500000,ne=10,nt=150) {
  ## SEIR stochastic simulation model.
  ## n = population size; ne = initially exposed; nt = number of days
  ##Sum of betas of indivial in infected state * beta of susepital individual * landa is prob S -> E
  ## 1/3 prob of moving E -> I
  ## 1/5 prob of moving I -> S with the assumption this can happen after 1 day
  ## 0 = Susecptible ## 1 = Exposed ## 2 = Infected ## 3 = Recovered/Dead 
  
  lamda <- 0.4/n  
  x <-rep(0,n) ## initialize to susceptible state
  beta <-rlnorm(n,0,0.5); beta <-beta/mean(beta) ## individual infection rates
  x[1:ne] <- 1 ## create some exposed
  
  new_infections <- new_infections_percentile <- new_infections_sample <- rep(0, nt-1) ## initialize new daily infections count
  ## initialize new daily infections count for bottom 10%
  ## initialize new daily infections count for sample
  
  bottomquantile <- quantile( beta, 0.1 ) ##Beta values in lowest 10%
  low_beta = beta < bottomquantile
  samplepopn <- rep(0, n) ; samplepopn[ sample(1:n, n*0.001) ] <- 1 #Indicator vector identifying 0.1% of popn
  
  for (i in 2:nt) { ## loop over days
    
    prob_exp = sum(beta[x == 2]) * lamda * beta
    u <-runif(n) ## uniform random deviates 
    
    xu = (x == 1 & u < (1/3))
    new_infections[i-1] <- sum(xu) ##New people in infectious state for day i
    new_infections_percentile[i-1] <- sum(xu & low_beta) ##Bottom 10% new people in infectious state for day i
    new_infections_sample[i-1] <- sum(xu & samplepopn == 1) ##Sampled new people in infectious state for day i
    
    x[x == 2 & u < (1/5)] <- 3 ## I -> R with prob 1/5
    x[xu] <- 2 ## E -> I with prob 1/3
    x[x == 0 & u < prob_exp] <- 1  ## S -> E with prob
    
  }
  return(list(new_infections=new_infections, new_infections_percentile=new_infections_percentile, new_infections_sample=new_infections_sample))
} ## seir

## Plot graphs for initial model run 

exec <- seir()

xx <- exec$new_infections*100/5500000
yy <- exec$new_infections_percentile*100/550000
zz <- exec$new_infections_sample*100/5500
day <- c(1:149)

maxxx <- which.max(xx)
maxyy <- which.max(yy)
maxzz <- which.max(zz)

plot(day, xx, pch=19,cex=.5, main="New Infection Trajectories" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)", type = 'l', col = 'brown', ylim = c(0,3))
lines(yy, col = 'green')
lines(zz, col = 'cadetblue')
legend(1, 2, legend=c("Whole population", "Cautious 10% of the popultion", "Random sample of 5,500"),
       col=c("brown", "green", "cadetblue"), lty=1:1, cex=0.5)
abline(v = c(maxxx, maxyy, maxzz) , lty=c(2, 2, 2), lwd=c(0.4, 0.4, 0.4), col=c("grey", "grey"))
mtext(maxxx, side = 1, at = maxxx, cex = 0.7) ; mtext(maxyy, side = 1, at = maxyy, cex = 0.7) ; mtext(maxzz, side = 1, at = maxzz, cex = 0.7) 

##Plot graphs for 10 sample runs

l1 <- list()
l2 <- list()
l3 <- list()

for (i in 1:10) {
  run <- seir()
  l1[i] <- list(run$new_infections*100/5500000)
  l2[i] <- list(run$new_infections_percentile*100/550000)
  l3[i] <- list(run$new_infections_sample*100/5500)
}

a <- pmax(unlist(l1[1]), unlist(l1[2]), unlist(l1[3]), unlist(l1[4]), unlist(l1[5]), unlist(l1[6]), unlist(l1[7]), unlist(l1[8]), unlist(l1[9]), unlist(l1[10]))
b <- pmax(unlist(l2[1]), unlist(l2[2]), unlist(l2[3]), unlist(l2[4]), unlist(l2[5]), unlist(l2[6]), unlist(l2[7]), unlist(l2[8]), unlist(l2[9]), unlist(l2[10]))
c <-  pmax(unlist(l3[1]), unlist(l3[2]), unlist(l3[3]), unlist(l3[4]), unlist(l3[5]), unlist(l3[6]), unlist(l3[7]), unlist(l3[8]), unlist(l3[9]), unlist(l3[10]))

d <- pmin(unlist(l1[1]), unlist(l1[2]), unlist(l1[3]), unlist(l1[4]), unlist(l1[5]), unlist(l1[6]), unlist(l1[7]), unlist(l1[8]), unlist(l1[9]), unlist(l1[10]))
e <- pmin(unlist(l2[1]), unlist(l2[2]), unlist(l2[3]), unlist(l2[4]), unlist(l2[5]), unlist(l2[6]), unlist(l2[7]), unlist(l2[8]), unlist(l2[9]), unlist(l2[10]))
f <-  pmin(unlist(l3[1]), unlist(l3[2]), unlist(l3[3]), unlist(l3[4]), unlist(l3[5]), unlist(l3[6]), unlist(l3[7]), unlist(l3[8]), unlist(l3[9]), unlist(l3[10]))

day <- c(1:149)

plot(day, a, type = "l", ylim = c(0,3), col="grey", pch=19,cex=.5, main="Variability of infection in the whole population" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)")
lines(day, d, type = "l", col="grey")
polygon(c(day, rev(day)), c(a, rev(d)), col = "grey")
plot(day, b, type = "l", ylim = c(0,3), col="grey", pch=19,cex=.5, main="Variability of infection in cautious 10% of the population" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)")
lines(day, e, type = "l", col="grey")
polygon(c(day, rev(day)), c(b, rev(e)), col = "grey")
plot(day, c, type = "l", ylim = c(0,3), col="grey", pch=19,cex=.5, main="Variability of infection in a random sample of 0.1% of population" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)")
lines(day, f, type = "l", col="grey")
polygon(c(day, rev(day)), c(c, rev(f)), col = "grey")


