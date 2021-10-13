## 0 = Susecptible ## 1 = Exposed ## 2 = Infected ## 3 = Recovered/Dead 
## x is all of the people 

seir <-function(n=5500000,ne=10,nt=100) {
  ## SEIR stochastic simulation model.
  ## n = population size; ne = initially exposed; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  lamda <- 0.4/n
  x <-rep(0,n) ## initialize to susceptible state
  beta <-rlnorm(n,0,0.5); beta <-beta/mean(beta) ## individual infection rates
  x[1:ne] <- 1 ## create some exposed
  S <-E <-I <-R <-rep(0,nt) ## set up storage for pop in each state
  S[1] <-n-ne;E[1] <-ne ## initialize
  for (i in 2:nt) { ## loop over days
    

    
    u <-runif(n) ## uniform random deviates
  
    x[x==0 & u < (lamda*(sum(beta[x == 2])*beta[x]))] <- 2  
    x[x==1 & u < 1/3] <- 2 ## E -> I with prob 1/3
    x[x== 2 & u < 1/5] <- 3 ## E -> I with prob 1/5
  
    x[x == 0 & u < beta * I[i - 1]] <-1 ## S -> E with prob beta*I[i-1]
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
} ## seir
  

ne <- 10
nt <- 100
n<-5500000 
x <-rep(0,n) ## initialize to susceptible state
beta <- rlnorm(n,0,0.5); beta <-beta/mean(beta) ## individual infection rates
x[1:ne] <- 1 ## create some exposed
S <-E <-I <-R <-rep(0,nt) ## set up storage for pop in each state
S[1] <-n-ne;E[1] <-ne ## initialize

sum(beta)

f <- c(1,2,3,6)
h <- c(3,4,5,5)

sum(f[h==5])

