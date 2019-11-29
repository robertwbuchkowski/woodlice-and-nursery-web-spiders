require(brms)
require(tidyverse)

# Extract isopod data and build distribution ----

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

# Load in Jennie's Excell data ----

library(readxl)

# Load 2011 Data

sheets = vector(mode = "list", length = 50)

for(i in 1:50){
  sheets[[i]] = read_xlsx("Data/doi_10/Behavioral observations 2011.xlsx",sheet = (i+1), skip = 3)
}

View(sheets)

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

# Determine the mean and sd heights from combined data ----
# Now I need to fit a random effects model with STAN to get the distributions

# ....CLARUS ----

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

# ....ISOPOD ----

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

# .... MIRA ----
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

# .....FEMUR ----
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

# Put it all together ----
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

# Calculate the probability of encounter based on distributions ----

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


# .....WRONG VERSION ----
min.f1f2 <- function(x, mu1, mu2, sd1, sd2){
  f1 <- dnorm(x, mean = mu1, sd = sd1)
  f2 <- dnorm(x, mean = mu2, sd = sd2)
  pmin(f1, f2)
}

fitCLARUS
fitMIRA

integrate(min.f1f2, 0, Inf, mu1 = 28.95, mu2 = 41.41, sd1 = 19.11, sd2 = 17.99)


integrate(min.f1f2, -Inf, Inf, mu1 = 28.95, mu2 = 41.41, sd1 = 10, sd2 = 10)


min.f1f3 <- function(x, mu1, shape3, sd1, scale3){
  f1 <- dnorm(x, mean = mu1, sd = sd1)
  f3 <- rgamma(x, shape = shape3, scale = scale3)
  pmin(f1, f3) # NOT CORRECT --> IT IS DOING THE UNION, NOT THE PROBABILITY OF ENCOUNTER!!
}

integrate(min.f1f3, 0, Inf, mu1 = 41.41, shape3 = 1.7, sd1 = 17.99, scale3 = 1.59)
