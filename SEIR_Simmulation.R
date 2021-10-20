## SEIR Model

## Overview
## ********
## The SEIR model is a model to recreate the trajectory of a epidemic, in particular following the path of infections for 3 cohorts. 
## The whole population, the 'cautious' 10% of the population defined as 10% of the population with the lowest beta values and a 
## random sample of 0.1% of the population. Defining 4 states, susceptible, exposed, infected and recovered or dead on each day 
## an individual moves from susceptible state to the exposed state a function of their contact rate, the infected populations 
## contact rate and a viral infectivity parameter lamda. Each day an individual can move from being exposed to infected with 
## the disease with probability 1/3 and move from infected to recovered or dead with probability 1/5. Using the SEIR function, 
## alongside an initial run this program aims to plot the variability in these cohorts by running the model 10 separate times
## and plotting the variability on three separate graphs. 

## SEIR model details
## ******************
## 1. Vectors are initiated whereby 'x' is defined as the storage vector for the whole population each element can either be 
##    0 = Susceptible ## 1 = Exposed ## 2 = Infected ## 3 = Recovered/Dead 
## 2. Vectors are initiated whereby 'S', 'E', 'I', 'R' are vectors showing how many people are in each state each day
## 3. Beta's are assigned to each individual (using rnom) and stored in 'beta' the bottom 10% is identified and stored
## 4. A random sample of 0.1% is also taken from the population for use later (using sample)
## 5. Model starts by using uniform random deviates (runif) to move people from state E -> I and E -> R
## 5. A person, i,  in state S are assigned a probability of moving to state E by sum_infectious * lamda * beta(i) and 
##    moved to state S based on random deviates (runif)
## 6. Three count vectors are set up, new_infections, new_infections_percentile, new_infections_sample and the number of people 
##    moving from E -> I in each group is recorded for that dat
## 7. Model runs for 150 days
## 8. Plots are created using plot, lines, mtext, ablines and legend functions to plot the trajectory of all three samples

seir <-function(n=5500000,ne=10,nt=150) {
  
  ## SEIR stochastic simulation model
  ## n = population size; ne = initially exposed; nt = number of days
  ## Sum of betas of individual in infected state * beta of susceptible individual * lamda is prob S -> E
  ## 1/3 prob of moving E -> I
  ## 1/5 prob of moving I -> S with the assumption this can happen after 1 day
  ## 0 = Susceptible ## 1 = Exposed ## 2 = Infected ## 3 = Recovered/Dead 
  
  lamda <- 0.4/n  										## Assuming overall viral infectivity parameter
  x <-rep(0,n) 											## Initializing the vector for persons in susceptible state
  beta <-rlnorm(n,0,0.5); beta <-beta/mean(beta) 		## Assuming probability distribution for transmission rate
  x[1:ne] <- 1 											## Assuming some exposed state persons
  
  ## Initializing new daily infections, the cautious 10% with lowest transmission rate values and the 0.1% random sample of the population
  
  new_infections <- new_infections_percentile <- new_infections_sample <- rep(0, nt-1)
  
  
  bottomquantile <- quantile( beta, 0.1 ) 				## Beta values in lowest 10%
  low_beta = beta < bottomquantile						## low_beta refers to the 10% cautious persons
  samplepopn <- rep(0, n) ; samplepopn[ sample(1:n, n*0.001) ] <- 1 ## Indicator vector identifying 0.1% random sample from the population
  
  ## The infection model
  
  for (i in 2:nt) { 									## Looping over days
    
    prob_exp = sum(beta[x == 2]) * lamda * beta			## Calculating the probability of getting infected after being exposed 
    u <-runif(n) 										## Generating a uniform random deviates vector for model
    
    xu = (x == 1 & u < (1/3))							## Calculating vector xu to store persons after being exposed to infection
    new_infections[i-1] <- sum(xu) 								##New people in infectious state for day i
    new_infections_percentile[i-1] <- sum(xu & low_beta) 		##Bottom 10% new people in infectious state for day i
    new_infections_sample[i-1] <- sum(xu & samplepopn == 1) 	##Sampled new people in infectious state for day i
    
    x[x == 2 & u < (1/5)] <- 3 							## Persons changing states from I -> R with prob 1/5
    x[xu] <- 2 											## Persons changing states from E -> I with prob 1/3 as calculated for xu vector
    x[x == 0 & u < prob_exp] <- 1  						## Persons changing states from S -> E after probability of being exposed
    
  }
  
  ## Returning  list of all required persons/data - new infections per day, infections for the cautious 10% of population and 0.1% random sample from population
  
  return(list(new_infections=new_infections, new_infections_percentile=new_infections_percentile, new_infections_sample=new_infections_sample))
} ## End of SEIR model

## Plotting graphs for the initial model

exec <- seir()

## Normalizing the scales of infections of various cases onto a common scale

xx <- exec$new_infections*100/5500000
yy <- exec$new_infections_percentile*100/550000
zz <- exec$new_infections_sample*100/5500
day <- c(1:149)

## Calculating maximum values(number of days) for each infection peak case

maxxx <- which.max(xx)
maxyy <- which.max(yy)
maxzz <- which.max(zz)

## Plotting various cases on a common scale after normalization

plot(day, xx, pch=19,cex=.5, main="New Infection Trajectories" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)", type = 'l', col = 'brown', ylim = c(0,3))
lines(yy, col = 'green')
lines(zz, col = 'cadetblue')

## Creating a legend for the plot representing the color of each line and corresponding case of the model

legend("topleft", inset = 0.05, legend=c("Whole population", "Cautious 10% of the popultion", "Random sample of 5,500"),
       col=c("brown", "green", "cadetblue"), lty=1:1, cex=0.5)

## Creating a line to highlight the day at which each case peaks in the model

abline(v = c(maxxx, maxyy, maxzz) , lty=c(2, 2, 2), lwd=c(0.4, 0.4, 0.4), col=c("grey", "grey"))
mtext(maxxx, side = 1, at = maxxx, cex = 0.7) ; mtext(maxyy, side = 1, at = maxyy, cex = 0.7) ; mtext(maxzz, side = 1, at = maxzz, cex = 0.7) 

## Plotting graphs replicating 10 model simulations representing variability of the infection in each new infections, 10% cautious people and 0.1% random sample case

## Creating lists to store model outputs for 10 runs

l1 <- list()
l2 <- list()
l3 <- list()
v1 <- rep(0, 10)
v2 <- rep(0, 10)
v3 <- rep(0, 10)

## Looping over to run the model for 10 times with parallel normalization and finding maximum(day at which infection peaks) of each case

for (i in 1:10) {
  run <- seir()
  l1[i] <- list(run$new_infections*100/5500000)
  l2[i] <- list(run$new_infections_percentile*100/550000)
  l3[i] <- list(run$new_infections_sample*100/5500)
  
  v1[i] <- which.max(run$new_infections*100/5500000)
  v2[i] <- which.max(run$new_infections_percentile*100/550000)
  v3[i] <- which.max(run$new_infections_sample*100/550)
}

## Box plot representing the peak values for 10 runs for each of the 3 cases/cohorts

boxplot(v1, v2, v3, main = "Boxplot & Whistlers for the infection peak of 3 cohorts", horizontal = TRUE, col = c("brown","green", "cadetblue"), names = c("Whole population", "Cautious 10%", "Random sample"), at = c(1, 6, 11))


## Below lines of code were written to plot and demonstrate the variability of each case when model is simulated 10 times

## Finding pairwise maximum of each of the simulations to find the peaks of each case/cohort to shade the variability

a <- pmax(unlist(l1[1]), unlist(l1[2]), unlist(l1[3]), unlist(l1[4]), unlist(l1[5]), unlist(l1[6]), unlist(l1[7]), unlist(l1[8]), unlist(l1[9]), unlist(l1[10]))
b <- pmax(unlist(l2[1]), unlist(l2[2]), unlist(l2[3]), unlist(l2[4]), unlist(l2[5]), unlist(l2[6]), unlist(l2[7]), unlist(l2[8]), unlist(l2[9]), unlist(l2[10]))
c <- pmax(unlist(l3[1]), unlist(l3[2]), unlist(l3[3]), unlist(l3[4]), unlist(l3[5]), unlist(l3[6]), unlist(l3[7]), unlist(l3[8]), unlist(l3[9]), unlist(l3[10]))

## Finding pairwise minimum of each of the simulations to find the lowest of each case/cohort to shade the variability

d <- pmin(unlist(l1[1]), unlist(l1[2]), unlist(l1[3]), unlist(l1[4]), unlist(l1[5]), unlist(l1[6]), unlist(l1[7]), unlist(l1[8]), unlist(l1[9]), unlist(l1[10]))
e <- pmin(unlist(l2[1]), unlist(l2[2]), unlist(l2[3]), unlist(l2[4]), unlist(l2[5]), unlist(l2[6]), unlist(l2[7]), unlist(l2[8]), unlist(l2[9]), unlist(l2[10]))
f <- pmin(unlist(l3[1]), unlist(l3[2]), unlist(l3[3]), unlist(l3[4]), unlist(l3[5]), unlist(l3[6]), unlist(l3[7]), unlist(l3[8]), unlist(l3[9]), unlist(l3[10]))

day <- c(1:149)

## Plotting data obtained on simulating the model 10 times and shadowing the variation in each case with 'grey'

## This graph shows variability of infection for the whole population when simulated 10 times

plot(day, a, type = "l", ylim = c(0,3), col="grey", pch=19,cex=.5, main="Variability of infection in the whole population" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)")
lines(day, d, type = "l", col="grey")
polygon(c(day, rev(day)), c(a, rev(d)), col = "grey")

## This graph shows variability of infection for the cautious 10% of the population when simulated 10 times

plot(day, b, type = "l", ylim = c(0,3), col="grey", pch=19,cex=.5, main="Variability of infection in cautious 10% of the population" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)")
lines(day, e, type = "l", col="grey")
polygon(c(day, rev(day)), c(b, rev(e)), col = "grey")

## This graph shows variability of infection for the 0.1% random sample taken from the population

plot(day, c, type = "l", ylim = c(0,3), col="grey", pch=19,cex=.5, main="Variability of infection in a random sample of 0.1% of population" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)")
lines(day, f, type = "l", col="grey")
polygon(c(day, rev(day)), c(c, rev(f)), col = "grey")