# Code to run all three model simulations

require(brms)
require(tidyverse)

# 1. Extract isopod data and build distribution ----

ISOPOD <-
    # behaviour data from 2017
    read_csv("Data/Data_Sheet_Behaviour_FinalFinal.csv") %>% 
      filter(Animal=="ONAS") %>%
      select(-Time, -Hour, -Minute, -Disturbed, -Notes) %>% 
      rename(Species = Animal) %>% 
      left_join( # Add treatment information from the harvest data set
        read_csv("Data/Harvest_SU17_23Sept2017.csv") %>% 
          select(Block, Cage, LTreatment) %>% 
          rename(Treatment = LTreatment)
      ) %>%
  select(-x,-y,-Perch, -Beh, -Observer, -Species) %>%
  mutate(z = replace_na(z,1e6)) %>%
  mutate(z = ifelse(z < 0, 0, z))

IT = ISOPOD$z

IP = ISOPOD$z[ISOPOD$z != 1e6]

hist(IP)

MASS::fitdistr(IP, densfun = "exponential")

ISOPOD %>%
  filter(z < 1e6) %>%
  group_by(Block, Cage) %>%
    summarize(mean(z))

hist(rexp(length(IP), rate = 0.028028028))

# so, isopods are D.I * prop_being_AG * dist_overlap

length(IP)/(length(IT)*6) # probability of being aboveground...convervative
# prop_being_AG = 0.01244444

# this is much lower than expected for P. scaber (Hassall and Tuck 2007),
# where a conservative estimate of 40 -90 (mean = 80) % of the population is sheltering at a given
# time. O. asellus is less heat tolerant, so maybe that is OK.

# 1.1 Load in Jennie's Excell data ----

library(readxl)

# Load 2011 Data

sheets = vector(mode = "list", length = 50)

for(i in 1:50){
  sheets[[i]] = read_xlsx("Data/doi_10/Behavioral observations 2011.xlsx",sheet = (i+1), skip = 3)
}

# View(sheets)

output <- do.call("rbind", sheets)

IDS = excel_sheets("Data/doi_10/Behavioral observations 2011.xlsx")[2:51] %>%
  tibble::enframe(name = NULL) %>%
  separate(value, into=c(NA,"Pred",NA, "Rep"))



output = output %>%
  mutate(Pred = rep(IDS$Pred, each = 34),
         Rep = rep(IDS$Rep, each = 34))

colnames(output) = c("Time", "Pred_x", "Pred_y", "Pred_z", "Pred_Screen", "Pred_Grnd",
                     "Pred_Solidago", "Pred_Grass","GH_x", "GH_y", "GH_z", "GH_Screen", "GH_Grnd",
                     "GH_Solidago", "GH_Grass", "Pred", "Rep")

output2011 = output %>%
  gather(-Time, -Pred, -Rep, key = Var, value = Value) %>%
  separate(Var, into=c("Species", "Var")) %>%
  mutate(Value = as.numeric(Value)) %>% 
  filter(!is.na(Value)) %>%
  left_join(
    data.frame(Value = c(0,1,2,3,1,2,3),
               Var = c(rep("y", 4), rep("x", 3)),
               V2 = c(44,32,19,7,10,20,30))
  ) %>%
  mutate(Year = 2011)

# Load 2012 data

sheets = vector(mode = "list", length = 66)

IDS = excel_sheets("Data/doi_10/Behavioral observations 2012.xlsx")[2:67] %>%
  tibble::enframe(name = NULL) %>%
  separate(value, into=c(NA,"Pred",NA, "Rep"))

for(i in 1:66){
  sheets[[i]] = read_xlsx("Data/doi_10/Behavioral observations 2012.xlsx",sheet = (i+1), skip = 3)
  
  if(i == 6){
    sheets[[i]] = sheets[[i]][,1:15]
  }
  
  sheets[[i]] = sheets[[i]] %>% mutate(Pred = IDS$Pred[i],
                                       Rep = IDS$Rep[i])
  
  colnames(sheets[[i]]) = c("Time", "Pred_x", "Pred_y", "Pred_z", "Pred_Screen", "Pred_Grnd",
                       "Pred_Solidago", "Pred_Grass","GH_x", "GH_y", "GH_z", "GH_Screen", "GH_Grnd",
                       "GH_Solidago", "GH_Grass", "Pred", "Rep")
  
}

output <- do.call("rbind", sheets)

outputf <- output %>%
  gather(-Time, -Pred, -Rep, key = Var, value = Value) %>%
  separate(Var, into=c("Species", "Var")) %>%
  filter(Var %in% c("x", "y", "z")) %>%
  mutate(Value = as.numeric(Value)) %>% 
  filter(!is.na(Value)) %>%
  left_join(
    data.frame(V2 = c(seq(2.5, 75, 2.5),seq(2.5, 50, 2.5),seq(2.5, 50, 2.5)),
               Value = c(seq(1, 30, 1),seq(1, 20, 1),seq(1, 20, 1)),
               Var =c(rep("y", 30),rep("x", 20),rep("z", 20)))
  ) %>%
  mutate(Year = 2012) %>%
  bind_rows(
    output2011
  )

rm(output2011)

outputf %>% filter(Pred == "rim"| Pred == "rima" & Species == "Pred" & Var == "y") %>%
  ggplot(aes(x=V2)) + geom_histogram() + theme_classic()

outputf %>% filter(Pred == "mira" & Species == "Pred" & Var == "y") %>%
  ggplot(aes(x=V2)) + geom_histogram() + theme_classic()

outputf %>% filter(Pred == "rim"| Pred == "rima" & Species == "Pred" & Var == "y") %>%
  group_by(Rep, Year) %>%
  summarize(sd= sd(V2), V2 = mean(V2)) %>%
  filter(!is.na(V2)) %>%
  ggplot(aes(x=V2)) + geom_histogram() + theme_classic()

outputf %>% filter(Pred == "mira" & Species == "Pred" & Var == "y") %>%
  group_by(Rep, Year) %>%
  summarize(sd= sd(V2), V2 = mean(V2)) %>%
  filter(!is.na(V2)) %>%
  ggplot(aes(x=V2)) + geom_histogram() + theme_classic()

# 1.2 Determine the mean and sd heights from combined data ----
# Now I need to fit a random effects model with STAN to get the distributions

# ....1.2.1 CLARUS ----

CLARUS = outputf %>% filter(Pred == "rim"|Pred == "rima" & Species == "Pred" & Var == "y") %>%
  select(Rep, Year, V2) %>%
  filter(!is.na(V2))

make_stancode(
  V2 ~ (1|Year/Rep),
  data=CLARUS
)

fitCLARUS <- brm(
  V2 ~ (1|Rep), chains = 2, cores = 2,
  data=CLARUS,
  control = list(adapt_delta=0.99)
)

plot(fitCLARUS)
summary(fitCLARUS)
pairs(fitCLARUS)

# ....1.2.2 ISOPOD ----

ISOPOD2 = ISOPOD %>%
  filter(z < 1e6) %>%
  mutate(Rep = paste0(Block, Cage)) %>%
  rename(V2 = z) %>%
  select(Rep, V2) %>%
  mutate(V2 = ifelse(V2 == 0, 0.01, V2))

hist(ISOPOD2$V2)

fitISOPOD <- brm(
  V2 ~ (1|Rep), chains = 2, cores = 2,
  data=ISOPOD2,
  family= "gamma",
  control = list(adapt_delta=0.99)
)


plot(fitISOPOD)
summary(fitISOPOD)
pairs(fitISOPOD)

# shape = 1.70
# Intercept (mean) = 2.70
# mean = shape*scale ==> scale = mean/shape = 2.7/1.7

hist(rgamma(1000, shape = 1.7, scale = 1.59))

hist(ISOPOD2$V2)

make_stancode(
  V2 ~ (1|Rep),
  data=ISOPOD2,
  family= "gamma"
)

# .... 1.2.3 MIRA ----
# Now get the P. mira distribution (adding Jennie and our data together)

MIRA = outputf %>% filter(Pred == "mira" & Species == "Pred" & Var == "y") %>%
  select(Rep, V2) %>%
  filter(!is.na(V2)) %>%
  bind_rows(
    behaviour %>% filter(Species == "PIMI" & Temp == "C") %>%
      mutate(Rep = paste0("CS",Block, Year)) %>%
      rename(V2 = z) %>%
      select(Rep, V2) %>%
      filter(!is.na(V2))
  )

mean(MIRA$V2)
sd(MIRA$V2)


fitMIRA <- brm(
  V2 ~ (1|Rep), chains = 2, cores = 2,
  data=MIRA,
  control = list(adapt_delta=0.99)
)

plot(fitMIRA)
summary(fitMIRA)
pairs(fitMIRA)

lme4::lmer(V2~(1|Rep), data = MIRA)

lme4::glmer(V2~(1|Rep), data=ISOPOD2, family = "Gamma")

# .....1.2.4 FEMUR ----
# Now get the grasshopper distribution (adding Jennie and our data together)

FEMUR = outputf %>% filter(Pred == "mira" & Species == "GH" & Var == "y") %>%
  select(Rep, V2) %>%
  filter(!is.na(V2)) %>%
  bind_rows(
    behaviour %>% filter(Species == "MEFE" & Temp == "C") %>%
      mutate(Rep = paste0("CS",Block, Year)) %>%
      rename(V2 = z) %>%
      select(Rep, V2) %>%
      filter(!is.na(V2))
  )

mean(FEMUR$V2)
sd(FEMUR$V2)


fitFEMUR <- brm(
  V2 ~ (1|Rep), chains = 2, cores = 2,
  data=FEMUR,
  control = list(adapt_delta=0.99)
)

plot(fitFEMUR)
summary(fitFEMUR)
pairs(fitMIRA)

summary(lme4::lmer(V2~(1|Rep), data = FEMUR))

# 1.4 Put it all together ----
xs <- seq(0, 100, by = 0.1)
fmira <- dnorm(xs, mean = 41.41, sd = 17.99)
fclarus <- dnorm(xs, mean = 28.95, sd = 19.11)
fisopod <- dgamma(xs, shape = 1.7, scale = 1.59)
ffemur <- dnorm(xs, mean = 54.48, sd = 19.39)

plot(xs~fmira, type = "l", xlim = c(0, max(fmira, fclarus, fisopod)),
     ylab = "Height", xlab = "Frequency")
points(xs~fclarus, type = "l", col ="blue")
points(xs~fisopod, type = "l", col ="orange")
points(xs~ffemur, type = "l", col ="green")

# 1.5 Calculate the probability of encounter based on distributions ----

D1 = data.frame(Species = c("Clarus", "Femur", "Isopod"),
           Htmean = c(28.95, 54.48, NA),
           Htsd = c(19.11, 19.39, NA),
           Htmira = c(41.41, 41.41, 41.41),
           Sdmira = c(17.99, 17.99, 17.99),
           shape = c(NA, NA, 1.7),
           scale = c(NA, NA, 1.59))

D1[,"Overlap2"] = D1[,"Overlap"] = exp((-1*(D1$Htmira - D1$Htmean)^2)/(2*(D1$Sdmira^2 + D1$Htsd^2)))/
  sqrt(2*3.14*(D1$Sdmira^2 + D1$Htsd^2))

for(i in 1:2){
  D1[i,"Overlap2"] = sum(dnorm(seq(0,100,0.1), mean = D1$Htmira[i], sd = D1$Sdmira[i])*
                          dnorm(seq(0,100,0.1), mean = D1$Htmean[i], sd = D1$Htsd[i]))
  
}

D1[3,"Overlap2"] = sum(dnorm(seq(0,100,0.1), mean = D1$Htmira[3], sd = D1$Sdmira[3])*
                         dgamma(seq(0,100,0.1), shape = D1$shape[3], scale = D1$scale[3]))

D1

# Overlap2 is the probability of encounter of each organism.

# Now write the correct version for the function

calc.overlap <- function(htmira){
  D1 = data.frame(Species = c("Clarus", "Femur", "Isopod"),
                  Htmean = c(28.95, 54.48, NA),
                  Htsd = c(19.11, 19.39, NA),
                  Htmira = c(htmira, htmira, htmira),
                  Sdmira = c(17.99, 17.99, 17.99),
                  shape = c(NA, NA, 1.7),
                  scale = c(NA, NA, 1.59),
                  Overlap = c(NA,NA,NA))
  
  for(i in 1:2){
    D1[i,"Overlap"] = sum(dnorm(seq(0,100,0.1), mean = D1$Htmira[i], sd = D1$Sdmira[i])*
                            dnorm(seq(0,100,0.1), mean = D1$Htmean[i], sd = D1$Htsd[i]))
    
  }
  
  D1[3,"Overlap"] = sum(dnorm(seq(0,100,0.1), mean = D1$Htmira[3], sd = D1$Sdmira[3])*
                          dgamma(seq(0,100,0.1), shape = D1$shape[3], scale = D1$scale[3]))
  
  
  return(D1$Overlap)
}

calc.overlap(41.41)


# 1.6 Energy costs (Notes on data sources) ----

# M. Ford 1977 and A. Schmitz 2004 have the metabolic data for other spiders that can be used
# R. Wiegart has the energy content of grasshoppers 5388 cal/g DW
# M Ford 1977 for P. amentata std = 100 x 10^-3 J h-1 and act 450 x 10^-3 J h-1


# Convert Grasshopper body cal to J --> 22543 J/g DW
# Convert to per grasshopper (using weights from Wiegert) --> 
# 22543 J/g x 0.0015 g/# = 33.81 J/#

# Convert grasshopper energy to net gain using active handling time
# Handling time estimated at 20 minutes from Samu (1993). 

# 0.45 - 0.1 J h^-1 = 0.35 J h^-1 extra effort
# 0.35 J h^-1 x 1 h / 60 min = 0.00583333 J min^-1 

# Gain from grasshopper = 33.81 J - (0.0058333 J min^-1 x 20 minutes) = 33.69

# But, not all attacks are successfull, so we need to modify the expected gain:
# start with 25% chance of success.

# Cost of attacking an object:
# 30 seconds to attack and return to resting posture if attack unsuccessful
# cost - 0.00583333 J min^-1 x 0.5 min = 0.002916667 J

# Being higher in the canopy means being hotter and spiders experience increased
# basal respiration rate. Following Bazzaz and Mezga (1973) we estimate that going
# from the ground to 40cm in the canopy increases the temperature by ~5degC.
# For convenience we approximate the increase as 1C per 10cm up to 1-m during the 
# 8 hours of direct sunlight in the middle of the day. 
# We did not observe spiders moving during the hotter parts of the day to lower
# perches, so we don't model an adaptive behavior of shifting height throughout
# the day.

# spider resting metabolism goes up higher in the canopy, so can discount this as well

# 2. Signal detection theory functions ----

source("rosenblatt2019_spiderresp.R")

outcome <- function(VAR,Ploss = 0.25, upd = 2.5,success.rate = 0.25, Woodlice = T, rtPca = F){
  
  with(as.list(c(VAR)),{

    upu = 0
    sig.t = 1
    sig.e = 5
    Vcr = 0
    Vir = 0
    tca1 = 20 # handling time of successful attack
    tia1 = 0.5 # time of incorrect attack or unsuccessfull attack
    tcr1 = 0
    tir1 = 0
    
    TEMP = 25 + exp(htmira)*0.1
    
    TEMP = ifelse(TEMP > 35, 35, TEMP)
    
    stdM = stdE(TEMP)
    actM = actE(TEMP)
    
  
    # Grasshopper energy content times ingestion efficiency and assimilation efficiency
    Vca = (33.81*0.85*0.95 - (actM - stdM)*tca1)*success.rate - # attacked correctly and got it
      (actM - stdM)*tia1*(1-success.rate) # attacked correctly and missed
    Via = -(actM - stdM)*tia1 # attacked incorrectly
    
    N = 1 + exp(Nin)
    
    OVERLAP = calc.overlap(exp(htmira))
    
    Pd = ifelse(OVERLAP[2] > 1e-5, OVERLAP[2]/(OVERLAP[3] + OVERLAP[2]), 1e-5)
    
    sig.p = sqrt(sig.t^2 + (sig.e/sqrt(N))^2)
    
    lambda.star = (sig.p^2)*(log((Vcr - Via)/(Vca-Vir))-
                               log((1 - Pd)/(Pd)))/upd + upd/2
    
    if(lambda.star %in% c(Inf, -Inf)){
      browser()
    }
    
    # browser()
    
    Pia = sum(dnorm(lambda.star:100, mean = upu, sd = sig.p))
    Pcr = 1- Pia
    Pca = sum(dnorm(lambda.star:100, mean = upd, sd = sig.p))
    Pir = 1-Pca
    
    # eatting is so useful the spiders hit everything
    G.N = (1-Pd)*(Pia*Via + Pcr*Vcr) + Pd*(Pca*Vca + Pir*Vir) # in units of Joules
    
    W.N = G.N*(1-Ploss)^N # Still in units of Joules
    
    negcost = -(W.N + (stdM - stdE(25))*8*60)
    
    
    # We use the function of LOCO, because we don't think the other two are easily 
    # or reasonability applied here 
    # ---loss of future opportunity requires parameterizing it.
    # ---spiders are not likely to face multiple prey at once and 
    # ----predator risk is less important in our empirical case
    
    if(rtPca){
      return(Pca)
    }else{
      return(negcost)# return the negative
      }
  })
}


outcome(c(Nin = log(5),htmira = log(41.41)))

op = optim(par = c(c(Nin = log(5),htmira = log(41.41))), fn = outcome)

exp(op$par)

op2 = optim(par = c(c(Nin = log(10),htmira = log(41.41))), fn = outcome, Woodlice = F)

exp(op2$par)

# Explore outcomes over htmira and Nin values

mat1 = log(expand.grid(Nin = seq(1,40,5), htmira = seq(1, 100, 10)))
opt = rep(NA,dim(mat1)[1])
for(i in 1:dim(mat1)[1]){
  opt[i] = outcome(mat1[i,])
}
opt = -1*opt

opt = cbind(exp(mat1), opt)

plot(opt~htmira, data=opt, col = Nin)
plot(opt~Nin, data=opt, col = htmira)

# Scan optimum across 

mat1 = expand.grid(Ploss = seq(0.1,0.9,0.1), upd = seq(2, 20, 2))
opt = matrix(NA,nrow = dim(mat1)[1], ncol = 4)
for(i in 1:dim(mat1)[1]){
  op = optim(par = c(c(Nin = log(5),htmira = log(41.41))), 
             fn = outcome, Ploss = mat1[i,"Ploss"], upd = mat1[i,"upd"])
  
  opt[i,1] = op$value
  opt[i,2:3] = exp(op$par)
  opt[i,4] = outcome(op$par, Ploss = mat1[i,"Ploss"], upd = mat1[i,"upd"],rtPca = T)
  
}

opt = cbind(mat1, opt)
colnames(opt)[3:6] = c("value", "Nin", "htmira", "Pca")

plot(htmira~Ploss, data = subset(opt, htmira < 1000 & upd == 2), type ="l")

png("Plots/FigureS23.png", width = 8, height =5, units = "in", res = 600)

plot(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.2), type ="l", 
     ylim = c(60,100), col = "blue", lwd = 2,
     xlab = "Spider ability to distinguish grasshoppers and woodlices (upd)",
     ylab = "Optimal spider height (cm)")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.4 & Nin > 0.5), type ="p", lwd = 2, col = "blue")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.4), type ="l", lwd = 2, col = "orange")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.4 & Nin > 0.5), type ="p", lwd = 2, col = "orange")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.6), type ="l", lwd = 2, col = "purple")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.6 & Nin > 0.5), type ="p", lwd = 2, col = "purple")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.8), type ="l", lwd = 2, col = "cyan")
points(htmira~upd, data = subset(opt, htmira < 1000 & Ploss == 0.8 & Nin > 0.5), type ="p", lwd = 2, col = "cyan", cex = 2)
legend("topright", legend = seq(0.2, 0.8, by= 0.2),
       col = c("blue", "orange", "purple", "cyan"),
       lwd = 2, lty = 1, title = "Probability of prey loss per chance to attack",
       bty = 'n')
dev.off()



# 3. The net energy gain if spiders always attack ----

alwaysattack <- function(HTMIRA = log(41.41), success.rate = 0.25, 
                         en.per.day = 0.8, # Taken from data in Miller et al. 2014
                         Woodlice = T){
  
  tca1 = 20
  tia1 = 0.5
  
  htmira = exp(HTMIRA)
  
  TEMP = 25 + htmira*0.1
  
  TEMP = ifelse(TEMP > 35, 35, TEMP)
  
  stdM = stdE(TEMP)
  actM = actE(TEMP)
  
  
  # Grasshopper energy content times ingestion efficiency and assimilation efficiency
  Vca = (33.81*0.85*0.95 - (actM - stdM)*tca1)*success.rate - # attacked correctly and got it
    (actM - stdM)*tia1*(1-success.rate) # attacked correctly and missed
  Via = -(actM - stdM)*tia1 # attacked incorrectly
  
  OVERLAP = calc.overlap(htmira)
  
  Pd = ifelse(OVERLAP[2] > 1e-5, OVERLAP[2]/(OVERLAP[3] + OVERLAP[2]), 1e-5)
  
  # eatting is so useful the spiders hit everything
  
  if(Woodlice){
    G.N = en.per.day*(Pd*Vca + (1-Pd)*Via) - (stdM - stdE(25))*8*60
  }else{
    G.N = en.per.day*(Pd*Vca) - (stdM - stdE(25))*8*60
  }
  
  
  return(G.N)
  
}

HT = seq(10,90, 1)

PWL = sapply(log(HT), FUN = alwaysattack)
MWL = sapply(log(HT), FUN = alwaysattack, Woodlice = F)

plot(PWL~HT, type ="n", xlab = "Spider Height (cm)", ylab = "Expected gain (J)")
rect(41.41 - 17.99, 0, 41.41 + 17.99, 7, col = "grey", border = "grey")
abline(v = 41.41, lty = 2)
points(PWL~HT, type ="l", lwd = 3, col = "blue")
points(MWL~HT, type ="l", col = "orange", lty = 3, lwd = 3)
legend("bottomright", legend = c("Woodlice", "No woodlice"), lwd = 3, 
       lty = c(1,2), col = c("blue", "orange"), bty = "n")

# 4. Individual based simulation tracking spider movement ----

#calculate average movement from our data

source("behaviour_March2019.R")

touse = behaviour %>%
  filter(Species == "PIMI") %>%
  mutate(ID = paste0(Year, Block, Cage)) %>%
  select(Time, ID, z)

avg.cage = rep(NA, length(unique(touse$ID)))
freq.move = avg.cage
avg.cage2 = avg.cage

for(i in 1:length(unique(touse$ID))){
  touse2 = subset(touse, ID == unique(touse$ID)[i] & !is.na(z))
  series = abs((touse2$z - lag(touse2$z, k = 1))[-1])
  
  freq.move[i] = length(series[series >0])/length(series)
  
  avg.cage2[i] = mean(series)
  
  if(any(series > 0)){
    avg.cage[i] = mean(series[series > 0])
  }else{
    avg.cage[i] = NA
  }
  
}


rand.mv = 0.1 # fit to overall movement data
dist_mean = mean(avg.cage[!is.na(avg.cage)]) #2.1 in Miller et al. 2014
dist_sd = sd(avg.cage[!is.na(avg.cage)]) #3.2 in Miller et al. 2014

freq.move_mean = mean(freq.move[!is.na(freq.move)])

freq.of.movement.cage = length(avg.cage[!is.na(avg.cage)])/length(avg.cage)


simfunc <- function(steps = 100, reps = 100, Woodlice = T){
  traj = matrix(NA, nrow = steps, ncol = reps)
  
  for(j in 1: reps){
    
    traj[1,j] = runif(1, 0, 100)
    
    for(i in 2:steps){
      OVERLAP = calc.overlap(traj[(i-1),j])
      
      if(runif(1,0,1) <= OVERLAP[2]){
        move = 0
      }else{
        if(runif(1,0,1) <= OVERLAP[3] & Woodlice |
           runif(1,0,1) <= rand.mv){
          move = rnorm(1, mean = dist_mean, sd = dist_sd)*sample(c(-1,1),1)
        }else{
          move = 0
        }
      }
      
      dest = traj[(i-1),j] + move
      
      traj[i,j] = ifelse(dest > 100 | dest < 0, traj[(i-1),j],dest)
      
    }
  }
  return(traj)
}

trajW = simfunc()
traj0 = simfunc(Woodlice = F)

plot(trajW[,1], type ="l", ylim = c(0,100), xlim = c(0, dim(trajW)[2]),
     col = alpha("blue", 0.1),
     xlab = "Time steps (30 minutes)", ylab = "Spider Height")
rect(24, 0, 48, 2, col = "black", border = "black")
for(i in 2:dim(trajW)[1]){
  points(trajW[,i],col = alpha("blue", 0.1), type ="l")
}
for(i in 1:dim(traj0)[1]){
  points(traj0[,i],col = alpha("orange", 0.1), type ="l")
}

polygon(c(seq(1,100,1), rev(seq(1,100,1))),
        c(apply(trajW,1, mean) - apply(trajW,1, sd),
          rev(apply(trajW,1, mean) + apply(trajW,1, sd))),
        border = alpha("blue", 0.2), col = alpha("blue", 0.2))

polygon(c(seq(1,100,1), rev(seq(1,100,1))),
        c(apply(traj0,1, mean) - apply(traj0,1, sd),
          rev(apply(traj0,1, mean) + apply(traj0,1, sd))),
        border = alpha("orange", 0.2), col = alpha("orange", 0.2))

points(apply(trajW,1, mean), type = "l", col = "white", lwd = 4)
points(apply(traj0,1, mean), type = "l", col = "white", lwd = 4)
points(apply(trajW,1, mean), type = "l", col = "blue", lwd = 3)
points(apply(traj0,1, mean), type = "l", col = "orange", lwd = 3)



hist(trajW[1,])
hist(trajW[dim(trajW)[2],])

pathdata = matrix(NA, nrow = dim(trajW)[1], ncol = 4)

for(j in 1: dim(trajW)[1]){
  mv = abs((trajW[,j] - lag(trajW[,j], k = 1))[-1])
  pathdata[j,1] = mean(mv)
  pathdata[j,2] = sd(mv)
  pathdata[j,3] = length(mv[mv >0])/length(mv)
  pathdata[j,4] = ifelse(max(mv[24:48])>0, 1,0)
}

pathdata2 = matrix(NA, nrow = dim(traj0)[1], ncol = 4)

for(j in 1: dim(traj0)[1]){
  mv = abs((traj0[,j] - lag(traj0[,j], k = 1))[-1])
  pathdata2[j,1] = mean(mv)
  pathdata2[j,2] = sd(mv)
  pathdata2[j,3] = length(mv[mv >0])/length(mv)
  pathdata2[j,4] = ifelse(max(mv[24:48])>0, 1,0)
}

check = data.frame(rbind(apply(pathdata,2, mean),
                         apply(pathdata2,2, mean),
      c(mean(avg.cage2[!is.na(avg.cage2)]),sd(avg.cage2[!is.na(avg.cage2)]), freq.move_mean,freq.of.movement.cage)))
colnames(check) = c("Avg_move", "Sd_move", "Prob_move_time_step", "Prob_move_cage")
check[,"Catagory"] = c("SimulationW", "Simulation0", "Expected")
check




write.csv(trajW,"SimRes/trajW_12Dec2019.csv", row.names = F)
write.csv(traj0,"SimRes/traj0_12Dec2019.csv", row.names = F)

# 5. Final plot with all results ----

png("Plots/Figure2.png", width = 10, height =5, units = "in", res = 600)

par(mfrow=c(1,2))

plot(PWL~HT, type ="n", xlab = "Spider Height (cm)", ylab = "Expected gain (J)",
     xlim = c(0,100))
rect(41.41 - 17.99, 0, 41.41 + 17.99, 7, col = "lightgrey", border = "lightgrey")
abline(v = 41.41, lty = 2)
points(PWL~HT, type ="l", lwd = 3, col = "blue")
points(MWL~HT, type ="l", col = "orange", lty = 3, lwd = 3)
legend("bottomright", legend = c("Woodlice", "No woodlice"), lwd = 3, 
       lty = c(1,2), col = c("blue", "orange"), bty = "n")
legend("topleft", legend = "A", bty = "n")


plot(trajW[,1], type ="l", ylim = c(0,100), xlim = c(0, dim(trajW)[2]),
     col = alpha("blue", 0.1),
     xlab = "Time steps (30 minutes)", ylab = "Spider Height")
rect(24, 0, 48, 2, col = "black", border = "black")

for(i in 2:dim(trajW)[1]){
  points(trajW[,i],col = alpha("blue", 0.1), type ="l")
}
for(i in 1:dim(traj0)[1]){
  points(traj0[,i],col = alpha("orange", 0.1), type ="l")
}

polygon(c(seq(1,100,1), rev(seq(1,100,1))),
        c(apply(trajW,1, mean) - apply(trajW,1, sd),
          rev(apply(trajW,1, mean) + apply(trajW,1, sd))),
        border = alpha("blue", 0.2), col = alpha("blue", 0.2))

polygon(c(seq(1,100,1), rev(seq(1,100,1))),
        c(apply(traj0,1, mean) - apply(traj0,1, sd),
          rev(apply(traj0,1, mean) + apply(traj0,1, sd))),
        border = alpha("orange", 0.2), col = alpha("orange", 0.2))

points(apply(trajW,1, mean), type = "l", col = "white", lwd = 4)
points(apply(traj0,1, mean), type = "l", col = "white", lwd = 4)
points(apply(trajW,1, mean), type = "l", col = "blue", lwd = 3)
points(apply(traj0,1, mean), type = "l", col = "orange", lwd = 3)
legend("topleft", legend = "B", bty = "n")

dev.off()
