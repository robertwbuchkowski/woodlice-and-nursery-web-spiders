# Maybe try running a signal detection model to see whether spiders should move upwards

# also include jumping spiders as IGP with another distribution (probably lower in the canopy)

# maybe also crickets on the ground, but I need a height distribution!

# then use SDT to determine whether isopod incorrect detection changes where the spiders should hunt.

# The two unknown parameters are (1) relative cost of a mistake grap and (2) how good spiders are at identifying isopods as bad prey

# Use Abbott and Sherratt 2013, modify eqn 6 so that the probabilities 
# are our own results of probability of encounter for each species.

PARS <- c(upd = 2.5,
          upu = 0,
          sig.t = 1, 
          sig.e = 5,
          Vca = 7.5,
          Via = -10,
          Vcr = 0,
          Vir = 0,
          Du1 = 5,
          tca1 = 10,
          tia1 = 1,
          tcr1 = 0,
          tir1 = 0,
          # t2 = 1,
          # G2 = 0.5,
          Dd1 = 5,
          # D2 = 1,
          Ploss = 0.25)

outcome <- function(Nin,pars){
  
  with(as.list(c(Nin, pars)),{
    
    N = Nin
    
    Pd = Dd1/(Du1 + Dd1)
    
    sig.p = sqrt(sig.t^2 + (sig.e/sqrt(N))^2)
    
    lambda.star = (sig.p^2)*(log((Vcr - Via)/(Vca-Vir))-
                               log((1 - Pd)/(Pd)))/upd + upd/2
    
    Pia = sum(dnorm(lambda.star:100, mean = upu, sd = sig.p))
    Pcr = 1- Pia
    Pca = sum(dnorm(lambda.star:100, mean = upd, sd = sig.p))
    Pir = 1-Pca
  
    G.N = (1-Pd)*(Pia*Via + Pcr*Vcr) + Pd*(Pca*Vca + Pir*Vir)
    
    t1.N = (1- Pd)*(Pia*tia1 + Pcr*tcr1)+
      (Pd)*(Pca*tca1 + Pir*tir1)
    
    # W.N = ((Du1 + Dd1)*G.N + D2*G2)/(1 + (Du1 + Dd1)*(N+t1.N)+D2*t2)
    
    W.N = G.N*(1-Ploss)^N 
    # We use the function of LOCO, because we don't think the other two are easily 
    # or reasonability applied here 
    # ---loss of future opportunity requires parameterizing it.
    # ---spiders are not likely to face multiple prey at once and 
    # ----predator risk is less important in our empirical case
    
    # browser()
    
    return(W.N)
  })
}

outcome(5, PARS)

plot(sapply(seq(1,40, 1), FUN = outcome, pars = PARS))

optimize(f = outcome,lower = 1, upper = 1000, maximum = T,pars= PARS)

scanp <- function(inPloss){
  PARS <- c(upd = 2.5,
            upu = 0,
            sig.t = 1, 
            sig.e = 5,
            Vca = 7.5,
            Via = -10,
            Vcr = 0,
            Vir = 0,
            Du1 = 5,
            tca1 = 10,
            tia1 = 1,
            tcr1 = 0,
            tir1 = 0,
            Dd1 = 5,
            Ploss = inPloss)
  
  op = optimize(f = outcome,lower = 1, upper = 1000, maximum = T,pars= PARS)
  
  A = outcome(floor(op$maximum), pars=PARS)
  B = outcome(ceiling(op$maximum), pars=PARS)
  
  ppp = ifelse(A>B, 
               c(objective = A, maximum = floor(op$maximum)), 
               c(objective = B, maximum = ceiling(op$maximum)))
  
  return(ppp)
  
  
}

scanp(0.25)

o1 = sapply(seq(0.01, 0.5, 0.01), FUN = scanp)

plot(t(o1))


