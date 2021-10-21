## Work Group 1 - Sanika Baxi - s2159255		Ayush Oza - s2184992		Gowtham Palepu - s2113890
## https://github.com/ayushoza1/SP-Assesment2.git 

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
## and plotting the variability of the runs. 

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
  
  lamda <- 0.4/n ## Assuming overall viral infectivity parameter
  x <-rep(0,n) ## Initializing the vector for persons in susceptible state
  beta <-rlnorm(n,0,0.5); beta <-beta/mean(beta) ## Assuming probability distribution for transmission rate
  x[1:ne] <- 1 ## Assuming some exposed state persons
  
  ## Initializing new daily infections, the cautious 10% with lowest transmission rate values and the 0.1% random sample of the population
  
  new_infections <- new_infections_percentile <- new_infections_sample <- rep(0, nt-1)

  bottomquantile <- quantile( beta, 0.1 ) ## Beta values in lowest 10%
  low_beta = beta < bottomquantile ## low_beta refers to the 10% cautious persons
  samplepopn <- rep(0, n) ; samplepopn[ sample(1:n, n*0.001) ] <- 1 ## Indicator vector identifying 0.1% random sample from the population
  
  ## The infection model
  
  for (i in 2:nt) { ## Looping over days
    
    
    prob_exp = sum(beta[x == 2]) * lamda * beta	## Calculating the probability of getting infected after being exposed 
    u <-runif(n) ## Generating a uniform random deviates vector for model
    
    xu <- (x == 0 & u < prob_exp)	## Calculating vector xu to store persons after being exposed to infection
    x1 <- which(x == 1)
    x2 <- which(x == 2)
    
    new_infections[i-1] <- sum(xu) 	##New people in infectious state for day i
    new_infections_percentile[i-1] <- sum(xu & low_beta) ##Bottom 10% new people in infectious state for day i
    new_infections_sample[i-1] <- sum(xu & samplepopn == 1) ##Sampled new people in infectious state for day i
    
    x[x2][u[x2] < (1/5)] <- 3 ## Persons changing states from I -> R with prob 1/5
    x[x1][u[x1] < (1/3)] <- 2 ## Persons changing states from E -> I with prob 1/3 as calculated for xu vector
    x[xu] <- 1  ## Persons changing states from S -> E after probability of being exposed
    
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

maxxx <- which.max(xx) ##infection peak of whole popn
maxyy <- which.max(yy) ##infection peak of cautious 10%
maxzz <- which.max(zz) ##infection peak of random sample 0.1%

## Plotting various cases on a common scale after normalization

plot(day, xx, pch=19,cex=.5, main="New Infection Trajectories - Initial Run" ,xlab="Day number",ylab="Percentage of cohort newly infected (%)", type = 'l', col = 'brown', ylim = c(0,3.5))
lines(yy, col = 'green')
lines(zz, col = 'cadetblue')

## Creating a legend for the plot representing the color of each line and corresponding case of the model

legend("topleft", inset = 0.05, legend=c("Whole population", "Cautious 10% of the popultion", "Random sample of 5,500"),col=c("brown", "green", "cadetblue"), lty=1:1, cex=0.5)

## Creating a line to highlight the day at which each case peaks in the model

points(maxxx, max(xx), pch=16) ; points(maxyy, max(yy), pch=16) ; points(maxzz, max(zz), pch=16)
abline(v = c(maxxx, maxyy, maxzz) , lty=c(2, 2, 2), lwd=c(0.4, 0.4, 0.4), col=c("grey", "grey"))
text(maxxx, max(xx)-0.2, labels = maxxx, cex = 0.7) ; text(maxyy, max(yy)-0.2, labels = maxyy, cex = 0.7) ; text(maxzz, max(zz)-0.2, labels = maxzz, cex = 0.7)

## 10 sample runs - Overview
## *************************
## Plotting graphs replicating 10 model simulations to represent the variability of the infections in each of the three cohorts. Trajectory
## of the pandemic is plotted for each sample run alongside a stripchart and boxplot. The boxplot shows the variabilty in pandemic peak for the 10 runs 
## , however, it ignores the pairwsie variability or day on day changes in infection peaks. The line plot shows the infection peak in each of the 10 simulations
## The first stripchart shows the day on day variability by standardizing the infection peak movement from one sample run to the next to 10 days movement and 
## showing the corresponding movement in the other cohorts peak. While the second stripshart shows the change in cohorts infection peak from run to run compared
## to the change in infection peak for the whole population.

v1 <- v2 <- v3 <- rep(0, 10) ##initialize vector to store day of infection peak for each cohort for each sample run
d1 <- d2 <- d3 <- rep(0,9) ##initialize vector to store the difference in the infection peak between each run
sd1 <- sd2 <- sd3 <- rep(10,9) ##initialize vector to store the standardized difference in the infection peak between each run coressponding to 10 day movement for whole popn
dd2 <- dd3 <- rep(0,9) ##initialize vector to store the  difference in the infection peak movement compared to the whole populations infection peak moveemnt 

## Looping over to run the model for 10 times with parallel normalization and finding maximum day at which infection peaks of each cohort

for (i in 1:10) {
  run <- seir()
  
  xx <- run$new_infections*100/5500000
  yy <- run$new_infections_percentile*100/550000
  zz <- run$new_infections_sample*100/5500
  
  v1[i] <- which.max(run$new_infections*100/5500000)
  v2[i] <- which.max(run$new_infections_percentile*100/550000)
  v3[i] <- which.max(run$new_infections_sample*100/550)
}

## Calculating the change in the infection peak for each cohort from sample run to the next

for (i in 1:9) { ##Loop over each day to the next
  
  ## If the change in infection peak of whole popn between 2 sample runs is 0 assume 1 day change in infection peak for standardization purposes
  
  if (v1[i+1] == v1[i]) { 
    d1[i] <- 1
  } else {
    d1[i] <- v1[i+1] - v1[i] ## Storage of change in infection peak between two runs
  }
  
  d2[i] <- v2[i+1] - v2[i] ## Storage of change in infection peak for cautious 10% of popn
  d3[i] <- v3[i+1] - v3[i] ## Storage of change in infection peak for sample 0.1% of popn
  
  }

sd2 <- (d2/d1)*10 ##Standardizing the movement in cautious 10% infection peak  corresponding to a 10 day movement in infection peak for the whole population
sd3 <- (d3/d1)*10 ##Standardizing the movement in random 0.1% infection peak  corresponding to a 10 day movement in infection peak for the whole population

dd2 <- d2 - d1 #Difference between move in infection people for whole population compared to cautious 10%
dd3 <- d3 - d1 #Difference between move in infection people for whole population compared to random 0.1%

## Line chart representing the infection peak for each of the 10 sample runs

plot(c(1:10), v1 , pch=19,cex=1, main="Infection peak for sample runs" ,xlab="Run Number",ylab="Infection peak (Day)", type = 'o', col = 'brown', ylim = c(80,105))
lines(v2, col = 'green', pch=19, type = 'o')
lines(v3, col = 'cadetblue', pch=19, type = 'o')
legend("topleft", inset = 0.05, legend=c("Whole population", "Cautious 10% of the popultion", "Random sample of 5,500"),col=c("brown", "green", "cadetblue"), lty=1:1, cex=0.5)

## Box plot representing the peak values for 10 runs for each of the 3 cases/cohorts

boxplot(v1, v2, v3, main = "Boxplot & Whistlers for the infection peak of 3 cohorts", horizontal = TRUE, col = c("brown","green", "cadetblue"), names = c("Whole population", "Cautious 10%", "Random sample"), at = c(1, 6, 11))

## Stripchart representing the movement in infection peaks corresponding to a 10 day move in the infection peak of whole popn

datapeakstd <- list("Whole population"=sd1, "Cautious 10%"=sd2, "Random sample" =sd3) ##list of standardized movements in infection peak for three cohorts

stripchart(datapeakstd, main="Cohort's Infection peak movement Vs \nStandardized 10 day move in population peak ", xlab="Movement in infection peak (days)", ylab="Cohort",method="jitter",col = c("brown","green", "cadetblue"),pch=16)

## Stripchart representing the movement in infection peaks compared to the move in the infection peak of whole popn from sample run to sample run

datapeak <- list("Cautious 10%"=dd2, "Random sample" =dd3) ##list of standardized movements in infection peak for three cohorts

stripchart(datapeak, main="Cohort's Infection peak movement Vs \nmove in population peak ", xlab="Move in infection peak from run to run compared\n to the infection peak movement in whole popn (days)", ylab="Cohort",method="jitter",col = c("brown","green", "cadetblue"),pch=16)


## Conclusion
## The graphs plotted show that infection trajectories simulated using the ZOE data app will have a later peak compared to similation of the whole popluation iven they are cautious 
## (later peaks on line graphs). Infection trajectories using the ZOE app are likely to show lower variability from run to run regarding the infection peak compared to the REACT-2
## trajectories (given the cluster of points around 10 days for cautious 10% on the stripchart compared to the random sample). However, the ZOE app data's simulated movement in infection
## peak will still show variability of +/-2/3 days compared to a simulated model of the whole population (last stripchart and line graph).



