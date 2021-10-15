## SEIR Model

## Overview
## ********
## The SEIR model is a model to recreate the trajectory of a epidemic, in particular following the path of infections for 3 cohorts. 
## The whole population, the 'cautious' 10% of the population defined as 10% of the population with the lowest beta values and a 
## random sample of 0.1% of the population. Defining 4 states, susceptible, exposed, infected and recovered or dead on each day 
## an individual moves from suseptible state to the exposed state a function of their contact rate, the infected populations 
##contact rate and a viral infectivity parameter lamda. Each day an individual can move from being exposed to infected with 
## the disease with probability 1/3 and move from infected to recovered or dead with probability 1/5

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
## 8. Plots are created using plot, lines, mtext, ablines and legend functions

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
  S <-E <-I <-R <-rep(0,nt) ## set up storage for pop in each state
  S[1] <-n-ne ##Initialize state S
  E[1] <-ne ## initialize state E
  
  new_infections <- rep(0, nt-1) ## initialize new daily infections count
  new_infections_percentile <- rep(0, nt-1) ## initialize new daily infections count for bottom 10%
  new_infections_sample <- rep(0, nt-1) ## initialize new daily infections count for sample
  
  bottomquantile <-quantile(beta, 0.1) ##Beta values in lowest 10%
  samplepopn <- rep(0, n) ; samplepopn[sample(1:n, n*0.001)] <- 1 #Indicator vector identifying 0.1% of popn
  
  for (i in 2:nt) { ## loop over days
    
    sum_infectious = sum(beta[x == 2]) 
    prob_exposed =  sum_infectious * lamda * beta 
    
    u <-runif(n) ## uniform random deviates 
    
    new_infections[i-1] <- sum(x == 1 & u < (1/3)) ##New people in infectious state for day i
    new_infections_percentile[i-1] <- sum(x == 1 & u < (1/3) & beta < bottomquantile) ##Bottom 10% new people in infectious state for day i
    new_infections_sample[i-1] <- sum(x == 1 & u < (1/3) & samplepopn == 1) ##Sampled new people in infectious state for day i
    
    x[x == 2 & u < (1/5)] <- 3 ## I -> R with prob 1/5
    x[x ==1 & u < (1/3)] <- 2 ## E -> I with prob 1/3
    
    x[x ==0 & u < prob_exposed] <- 1  ## S -> E with prob
    
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  list(S=S,E=E,I=I,R=R, new_infections=new_infections, new_infections_percentile=new_infections_percentile, new_infections_sample=new_infections_sample)
} ## seir


exec <- seir()

#Ayush

xx <- exec$new_infections
yy <- exec$new_infections_percentile
zz <- exec$new_infections_sample

trfxx = xx*100/5500000
trfyy = yy*100/550000
trfzz = zz*100/5500
day <- c(1:149)

maxtrfxx <- which.max(trfxx)
maxtrfyy <- which.max(trfyy)
maxtrfzz <- which.max(trfzz)

plot(day, trfxx, pch=19,cex=.5, main="New Infection Trajectories" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)", type = 'l', col = 'brown')
lines(trfyy, col = 'green')
lines(trfzz, col = 'cadetblue')
legend(1, 2, legend=c("Whole population", "Cautious 10% of the popultion", "Random sample of 5,500"),
       col=c("brown", "green", "cadetblue"), lty=1:1, cex=0.5)
abline(v = c(maxtrfxx, maxtrfyy, maxtrfzz) , lty=c(2, 2, 2), lwd=c(0.4, 0.4, 0.4), col=c("grey", "grey"))
mtext(maxtrfxx, side = 1, at = maxtrfxx, cex = 0.7) ; mtext(maxtrfyy, side = 1, at = maxtrfyy, cex = 0.7) ; mtext(maxtrfzz, side = 1, at = maxtrfzz, cex = 0.7) 

#Gowtham

x <- exec$new_infections
y <- exec$new_infections_percentile
z <- exec$new_infections_sample

trfx = (x - mean(x))/(max(x) - min(x))
trfy = (y - mean(y))/(max(y) - min(y))
trfz = (z - mean(z))/(max(z) - min(z))

plot(trfx, type = 'l', col = 'brown')
points(trfy, type = 'l', col = 'red')
lines(trfy, col = 'chocolate1')
points(trfz, type = 'l', col = 'cadetblue1')
lines(trfz, col = 'cadetblue')

legend(1, 0.5, legend=c("Line 1", "Line 2", "Line 3"),
       col=c("brown", "chocolate1", "cadetblue"), lty=1:2, cex=0.8)







