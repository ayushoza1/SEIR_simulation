## 0 = Susecptible ## 1 = Exposed ## 2 = Infected ## 3 = Recovered/Dead 
## x is all of the people 

seir <-function(n=5500000,ne=10,nt=150) {
  ## SEIR stochastic simulation model.
  ## n = population size; ne = initially exposed; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
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







