require(brms)
require(tidyverse)

# Load in data ----
source("behaviour_March2019.R")

includetemp = T
includelowsamplecages = T
shapevec = c("2013" = 3, "2015" = 19, "2017" = 17, "2018" = 15)

head(behaviour)

# clean up data for model ----

behtomod1 <- behaviour %>% select(Time, Cage, Block, Year, Species, z, Treatment, Temp) %>% 
  filter(Treatment %in% c("SIG", "SG")) %>% 
  group_by(Cage, Block, Year, Species, Treatment, Temp) %>%
  summarize(z = mean(z, na.rm=T)) %>% 
  ungroup() %>% spread(key=Species, value=z)

behtomod2 <- behaviour %>% select(Time, Cage, Block, Year, Species, z, Treatment, Temp) %>% 
  filter(Treatment %in% c("SIG", "SG")) %>% 
  group_by(Cage, Block, Year, Species, Treatment, Temp) %>%
  summarize(z = sd(z, na.rm=T)) %>%
  mutate(z =ifelse(z %in% c(0, NA), 0.0001, z)) %>%
  ungroup() %>% mutate(Species = paste0(Species, "sd")) %>% 
  spread(key=Species, value=z)

behtomod <- full_join(behtomod1, behtomod2) %>% 
  filter(MEFEsd > 0.99 |MEFEsd < 0.97)

if(!includelowsamplecages){
  behtomod = behtomod %>% 
    filter(!(Year==2018 & Cage %in% c("17", "21", "22", "37", "38", "42")))
}

if(includetemp){
  hartomod <- harvest %>% select(Block, Year, Treatment, Temp, Mf) %>%
    filter(Treatment %in% c("SIG", "SG")) %>%
    full_join(behtomod) %>% filter(complete.cases(.)) %>% distinct() %>%
    mutate(attack = exp((-1*(PIMI - MEFE)*(PIMI - MEFE))/(2*(PIMIsd*PIMIsd + MEFEsd*MEFEsd)))/sqrt(2*3.14*(PIMIsd*PIMIsd + MEFEsd*MEFEsd)))
}else{
  hartomod <- harvest %>% select(Block, Year, Treatment, Temp, Mf) %>%
    filter(Treatment %in% c("SIG", "SG") & Temp=="C") %>%
    full_join(behtomod) %>% filter(complete.cases(.)) %>% distinct() %>%
    mutate(attack = exp((-1*(PIMI - MEFE)*(PIMI - MEFE))/(2*(PIMIsd*PIMIsd + MEFEsd*MEFEsd)))/sqrt(2*3.14*(PIMIsd*PIMIsd + MEFEsd*MEFEsd)))
}

hartomod2 = hartomod %>% select(Block, Year, Mf, attack)

har2013 <- harvest %>% filter(Year == 2013) %>% 
  select(Treatment, Block, Mf) %>%
  left_join(
    harvest %>% filter(Year == 2013) %>% select(Block) %>% unique() %>%
      mutate(Block2 = seq(1,5))
  ) %>%
  select(-Block) %>% rename(Block = Block2) %>% 
  mutate(Year=2013, Block = paste(Block)) %>% filter(Treatment !="IG") %>%
  mutate(Treatment = ifelse(Treatment=="SG", "CSG", "CSIG"))


# fit model 1 ---- 

# predicts the mean and SD of the height data by treatment
# removing block random effect has no influence on result except reducing R2 a lot!
prior1 <- get_prior(mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
          data=hartomod)

prior1$prior[1] = "normal(0,100)"

make_stancode(
  mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
  data=hartomod, prior=prior1
)

fit1 <- brm(
  mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
  data=hartomod, chains=1, iter= 10000,
  control = list(adapt_delta=0.99), prior=prior1
)

Resid.brm <- residuals(fit1, type='pearson')
Fitted.brm <- fitted(fit1, scale='linear')
ggplot(data=NULL, aes(y=Resid.brm[,,1][,'Estimate'], x=Fitted.brm[,,1][,'Estimate'])) + geom_point() + theme_classic()
ggplot(data=NULL, aes(y=Resid.brm[,,2][,'Estimate'], x=Fitted.brm[,,2][,'Estimate'])) + geom_point() + theme_classic()

# add_loo(fit1)
summary(fit1)

pp_check(fit1, resp="MEFE")
pp_check(fit1, resp="PIMI")
plot(fit1)

# both multi-modal

bayes_R2(fit1)

marginal_effects(fit1, points=T)

# fit model 1.1: Used in the MS! ----

# Including Brandon's prior
# Spider
# mean = -12.7
# sd = 9.4

# Grasshopper
# mean =  -0.13
# sd = 16.1

prior1 <- get_prior(mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
                    data=hartomod)

prior1$prior[1] = "normal(0,100)"
prior1$prior[5] = "normal(-0.13,16.1)" # grasshopper temp
prior1$prior[16] = "normal(-12.7,9.4)" # spider temp

make_stancode(
  mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
  data=hartomod, prior=prior1
)

fit1b <- brm(
  mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
  data=hartomod, chains=1, iter = 10000,
  control = list(adapt_delta=0.99), prior=prior1
)

summary(fit1b)
plot(fit1b)
bayes_R2(fit1b)

# no effect on qualitative outcome!

# Fit model 1.2: removing behavioural cages from 2018 that lost data ----

# Create new dataset without these cages
lostafternoon = c("17", "21",  "22", "37", "38", "42")

temphartomod = hartomod %>% filter(!(Year ==2018 & Cage %in% lostafternoon))

fit1c <- brm(
  mvbind(MEFE, PIMI) ~ Treatment*Temp +(1|Year/Block),
  data=temphartomod, chains=1, iter = 10000,
  control = list(adapt_delta=0.99), prior=prior1
)

summary(fit1c)
plot(fit1c)
bayes_R2(fit1c)

rm(temphartomod, lostafternoon)

# fit model 2----
# predict the relationship between theoretical attack rate and grasshopper survival

#basically attack prediction does not predict the surival well

# test for zero inflation: not zero inflatted!
cnts <- rpois(1000, mean(hartomod2$Mf))
dat.tab <- table(cnts == 0)
dat.tab2 <-table(hartomod2$Mf ==0)
dat.tab/sum(dat.tab) # theoretical
dat.tab2/sum(dat.tab2) # actual
rm(cnts, dat.tab, dat.tab2)

prior2 <- get_prior(
  Mf~attack + (1|Year),
  data=hartomod2,
  family=poisson()
)

prior2$prior[1] <- "normal(0,0.5)"
prior2$prior[3:4] <- "student_t(3, 0, 5)"
prior2

make_stancode(Mf~attack + (1|Year),
              data=hartomod2,
              family=poisson(),
              prior=prior2)

fit2 <- brm(
  Mf~attack + (1|Year),
  data=hartomod2, chains=2, cores=2,
  family=poisson(),
  control = list(adapt_delta=0.99),
  prior=prior2
)

Resid.brm <- residuals(fit2, type='pearson')
Fitted.brm <- fitted(fit2, scale='linear')
ggplot(data=NULL, aes(y=Resid.brm[,'Estimate'], x=Fitted.brm[,'Estimate'])) + 
  geom_point() + theme_classic()

Resid.brm <- residuals(fit2, type='pearson', summary=FALSE)
SSres.brm <- apply(Resid.brm^2,1,sum)
lambda.brm = fitted(fit2, scale='response', summary=FALSE)
YNew.brm <- matrix(rpois(length(lambda.brm), lambda=lambda.brm),nrow=nrow(lambda.brm))

Resid1.brm<-(lambda.brm - YNew.brm)/sqrt(lambda.brm)
SSres.sim.brm<-apply(Resid1.brm^2,1,sum)
mean(SSres.sim.brm>SSres.brm)

# Doesn't work: from http://www.flutterbys.com.au/stats/tut/tut10.6b.html#h4_7 
simRes <- function(lambda, data,n=250, plot=T, family='poisson') {
  require(gap)
  N = nrow(data)
  sim = switch(family,
               'poisson' = matrix(rpois(n*N,apply(lambda,2,mean)),ncol=N, byrow=TRUE)
  )
  a = apply(sim + runif(n,-0.5,0.5),2,ecdf)
  resid<-NULL
  for (i in 1:nrow(data)) resid<-c(resid,a[[i]](data$y[i] + runif(1 ,-0.5,0.5)))
  if (plot==T) {
    par(mfrow=c(1,2))
    gap::qqunif(resid,pch = 2, bty = "n",
                logscale = F, col = "black", cex = 0.6, main = "QQ plot residuals",
                cex.main = 1, las=1)
    plot(resid~apply(lambda,2,mean), xlab='Predicted value', ylab='Standardized residual', las=1)
  }
  resid
}

simRes(lambda.brm, hartomod2, family='poisson')

plot(fit2)
pairs(fit2)
summary(fit2)
bayes_R2(fit2)
marginal_effects(fit2)

# Now estimate the brms model using the example online:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html

# fit model 2.1: without random effects ----
prior2pt0 <- get_prior(
  Mf~attack,
  data=hartomod2,
  family=poisson()
)

prior2pt0$prior[1] = "normal(0,10)"

prior2pt0

make_stancode(
  Mf~attack,
  data=hartomod2,
  family=poisson(),
  prior = prior2pt0
)

fit2pt0 <- brm(
  Mf~attack,
  data=hartomod2,
  family=poisson(),
  prior = prior2pt0,
  control = list(adapt_delta=0.99), chains=1, iter=10000
)

summary(fit2pt0)
plot(fit2pt0)
pairs(fit2pt0)
bayes_R2(fit2pt0)

require(lme4)

model2 = glmer(Mf~attack + (1|Year), data=hartomod2, family=poisson())

summary(model2)

# fit model 2.2 with binomial distribution: Used in MS ----

# Modify the hartomod so it can be a binomial model
hartomod3 = hartomod %>%
  select(Year,Block, Cage, Mf, attack) %>%
  mutate(TMf = ifelse(Year == 2018, 3,2))

prior3pt0 <- get_prior(Mf | trials(TMf) ~ attack  + (1|Year/Block),
                       data = hartomod3,
                       family = binomial())

prior3pt0$prior[1] = "normal(0,10)"

make_stancode(
  Mf | trials(TMf) ~ attack  + (1|Year/Block), 
  data = hartomod3,
  prior = prior3pt0,
  family = binomial()
)

fit3pt0 <- brm(Mf | trials(TMf) ~ attack  + (1|Year/Block), 
               data = hartomod3,
               prior = prior3pt0,
            family = binomial(),
            control = list(adapt_delta=0.99), chains=1, iter=10000)

summary(fit3pt0)
plot(fit3pt0)
bayes_R2(fit3pt0)

# Create SI Figure S2 -----

suppfigdata <- hartomod %>% filter(Temp=="C") %>% group_by(Year, Treatment) %>%
  summarize(sdMf = sd(Mf, na.rm=T),
            Mf = mean(Mf, na.rm=T)) %>%
  bind_rows(
    har2013 %>% 
      mutate(Treatment = ifelse(Treatment=="CSG", "SG", "SIG")) %>%
      select(-Block) %>%
      group_by(Year, Treatment) %>%
      summarize(sdMf = sd(Mf, na.rm=T),
                Mf = mean(Mf, na.rm=T))
  ) %>% 
  gather(-Year, -Treatment, key=Measure, value=Value) %>%
  mutate(TreatM = paste0(Treatment, "_", Measure)) %>%
  select(-Treatment, -Measure) %>%
  spread(key=TreatM, value=Value) %>%
  mutate(Eff = SIG_Mf -SG_Mf,
         sd = sqrt(SG_sdMf^2 + SIG_sdMf^2)) %>% left_join(
  data.frame(Year = c(2013,2015, 2017, 2018),
             SeasonTemp = c(mean(c(74.9,66.5)),
                            mean(c(72.7,72.1)),
                            mean(c(69.9,67.4)),
                            mean(c(73.3,72.9))),
             JulyTemp = c(74.9, 72.7, 69.9, 73.3),
             AugTemp = c(66.5, 72.1, 67.4, 72.9),
             StockDate = c(204, 218,209, 212)
             )
) %>%
  mutate(JulyTemp = (JulyTemp-32)*5/9,
         AugTemp = (AugTemp-32)*5/9)

sp1 <- suppfigdata %>%  
  ggplot(aes(x=StockDate, y=Eff, shape = as.factor(Year))) + 
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_pointrange(aes(ymin = Eff -sd, ymax = Eff + sd),
                  color="grey60") + 
  theme_classic() +
  geom_point(size = 3) + 
  scale_shape_manual(values=shapevec,name="Year") +
  ylab("Effect of woodlice on survival") + 
  xlab("Expt Start (day of the year)")

sp2 <- suppfigdata %>%  
  ggplot(aes(x=JulyTemp, y=Eff, shape = as.factor(Year))) + 
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_pointrange(aes(ymin = Eff -sd, ymax = Eff + sd),
                  color="grey60") + 
  theme_classic() +
  geom_point(size = 3) + 
  scale_shape_manual(values=shapevec,name="Year") +
  ylab("") + 
  scale_x_continuous(limits = c(18,24)) +
  xlab("Average July temperature")

sp3 <- suppfigdata %>%  
  ggplot(aes(x=AugTemp, y=Eff, shape = as.factor(Year))) + 
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_pointrange(aes(ymin = Eff -sd, ymax = Eff + sd),
                  color="grey60") + 
  theme_classic() +
  geom_point(size = 3) + 
  scale_shape_manual(values=shapevec,name="Year") +
  ylab("") + 
  scale_x_continuous(limits = c(18,24)) +
  xlab("Average August temperature")

png("Plots/FigureS2.png", width=8, height=4, units="in", res=600)
cowplot::plot_grid(NULL,
                   cowplot::get_legend(sp1 + theme(legend.position = "top")),
                   NULL,
                   sp1+ theme(legend.position = "none"),
                   sp2+ theme(legend.position = "none"),
                   sp3+ theme(legend.position = "none"), 
                   ncol=3, nrow=2,
                   rel_heights = c(0.1,1),
                   labels=c("", ""," ", "A", "B", "C"))
dev.off()

#stocking dates: 
# 2013 = July 23
# 2015 =  Aug 6
# 2017 = July 28
# 2018 = July 31

hartomod %>% filter(Temp=="C")%>% 
  group_by(Year, Treatment) %>%
  summarize(PIMI = mean(PIMI, na.rm=T)) %>%
  spread(key=Treatment, value=PIMI) %>%
  mutate(Eff = SIG -SG) %>% left_join(
    data.frame(Year = c(2013,2015, 2017, 2018),
               SeasonTemp = c(mean(c(74.9,66.5)),
                              mean(c(70.7,70.7)),
                              mean(c(69.9,67.4)),
                              mean(c(73.3,72.9))),
               JulyTemp = c(74.9, 72.7, 69.9, 73.3),
               AugTemp = c(66.5, 72.1, 67.4, 72.9)
    )
  ) %>%
  ggplot(aes(x=SeasonTemp, y=Eff, shape = as.factor(Year))) + 
  geom_point() + theme_classic() +
  scale_shape_manual(values=shapevec,name="Year")

# Create a graphic showing the P. mira body sizes measured in 2015

png("Plots/FigureS3.png", width=8, height=4, units="in", res=600)
ggpubr::ggarrange(
  read_csv("Data/Cage_data.csv") %>%
    ggplot(aes(x = Pm_body_Sept)) + geom_histogram(bins = 10) + theme_classic() + xlab("Spider body length (cm)"),
  read_csv("Data/Cage_data.csv") %>%
    ggplot(aes(x = Pm_body_Sept, y = MF_Sept)) + geom_point() + theme_classic() + xlab("Spider body length (cm)") + ylab("Grasshoppers surviving in September"),
  labels = "AUTO"
)
dev.off()

# Create Figure 1 graphic ----

hartomod_axis <- hartomod %>% mutate(XAXIS = paste0(Temp,Treatment))

plt_dat <- marginal_effects(fit1b,resp ="PIMI")$`PIMI.PIMI_Treatment:Temp`
plt_dat["XAXIS"] <- paste0(plt_dat$Temp, plt_dat$Treatment)
p1 = hartomod_axis %>% ggplot(aes(x=XAXIS)) + 
  geom_pointrange(data= plt_dat, size=2,
                  aes(y= estimate__, ymin=lower__, ymax=upper__)) +
  geom_jitter(aes(y=PIMI,color=Block, shape=as.factor(Year)),
              size=3, height=0, width=0.2) + theme_classic() +
  ylab(expression(italic(P.~mira)~Height~(cm))) + 
  xlab("Treatment") +
  annotate("text", x = 2.5, y = 95, label = expression(R[Bayes]^2==0.20)) +
  ylim(0,100) +
  scale_shape_manual(values=shapevec,name="Year") +
  scale_x_discrete(breaks=c("CSG","CSIG","WSG", "WSIG"),
                   labels=c("Ambient \n No woodlice",
                            "Ambient \n Woodlice",
                            "Warmed \n No woodlice",
                            "Warmed \n Woodlice"))

plt_dat <- marginal_effects(fit1b,resp ="MEFE")$`MEFE.MEFE_Treatment:Temp`
plt_dat["XAXIS"] <- paste0(plt_dat$Temp, plt_dat$Treatment)
p2 = hartomod_axis %>% ggplot(aes(x=XAXIS)) + 
  geom_pointrange(data= plt_dat, size=2,
                  aes(y= estimate__, ymin=lower__, ymax=upper__)) +
  geom_jitter(aes(y=MEFE,color=Block, shape=as.factor(Year)),
              size=3, height=0, width=0.2) + theme_classic() +
  ylab(expression(italic(M.~femurrubrum)~Height~(cm))) +
  xlab("Treatment") +
  annotate("text", x = 2.5, y = 95, label = expression(R[Bayes]^2==0.11)) +
  ylim(0,100)+
  scale_shape_manual(values=shapevec,name="Year") +
  scale_x_discrete(breaks=c("CSG","CSIG","WSG", "WSIG"),
                   labels=c("Ambient \n No woodlice",
                            "Ambient \n Woodlice",
                            "Warmed \n No woodlice",
                            "Warmed \n Woodlice"))

plt_dat <- marginal_effects(fit3pt0) # added model without random effects
p3 = ggplot() + 
  geom_ribbon(data=plt_dat$attack, aes(x=attack,ymin=lower__, ymax=upper__),alpha=0.1) +
  geom_line(data=plt_dat$attack, aes(x=attack,y=estimate__), size=2) +
  geom_jitter(data= hartomod2, aes(x=attack, y=Mf,color=Block, shape=as.factor(Year)),
              width=0, height=0.1, size=3) + theme_classic() +
  xlab(expression(Predicted~Attack~Rate~(time^-1))) + 
  ylab(expression(italic(M.~femurrubrum)~survival~("#")))+ 
  scale_shape_manual(values=shapevec,name="Year") +
  annotate("text", y = 2.5, x = 0.015, label = expression(R[Bayes]^2==0.24)) + ylim(-0.1, 4)

p4 = ggplot() + 
  geom_violin(data= hartomod_axis, aes(x=XAXIS, y=Mf)) +
  geom_jitter(data= hartomod_axis, aes(x=XAXIS, y=Mf,color=Block, shape=as.factor(Year)),
              width=0.1, height=0.1, size=3) + theme_classic() +
  xlab("Treatment") +
  ylab(expression(italic(M.~femurrubrum)~survival~("#"))) +
  geom_jitter(data = har2013, aes(x=Treatment, y=Mf,
                                  color=Block, shape=as.factor(Year)),
              size=3,width=0.3, height=0.3) + 
  scale_shape_manual(values=shapevec,name="Year") +
  ylim(-0.1,4) +
  scale_x_discrete(breaks=c("CSG","CSIG","WSG", "WSIG"),
                   labels=c("Ambient \n No woodlice",
                            "Ambient \n Woodlice",
                            "Warmed \n No woodlice",
                            "Warmed \n Woodlice"))

png("Plots/Figure1.png", width=8, height=8, units="in", res=600)
cowplot::plot_grid(cowplot::get_legend(p4 + theme(legend.position = "top") +
                                         scale_color_discrete(guide=F)),
                   cowplot::get_legend(p4 + theme(legend.position = "top") +
                                         scale_shape_discrete(guide=F)),
                   p1+ theme(legend.position = "none"),
                   p2+ theme(legend.position = "none"),
                   p3+ theme(legend.position = "none"),
                   p4+ theme(legend.position = "none"),
                   ncol=2, nrow=3,
                   rel_heights = c(0.3,1,1),
                   labels=c("", "", "A", "B", "C", "D"))
dev.off()


# Analyze the wolf spider data -----

read_csv("Data/Behaviour2017.csv") %>%
  filter(Cage %in% c("D", "E")) %>%
  group_by(Animal, Cage) %>%
  summarize(sd = sd(z, na.rm = T), height = mean(z, na.rm = T)) %>%
  filter(Animal != "ONAS")
