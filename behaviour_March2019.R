# New R Script for Behavioural Data on Spider and Grasshopper Interactions

# Load in the libraries
require(tidyverse)

# Load in the available data ----------------------------------------------

# Data from MESc Thesis
MEScdata <- read_csv("Data/MESc_harvestdata.csv")

MEScdata = subset(MEScdata, Order==6 & Treatment=="G"|Order==6 & Treatment=="F"|
                    Order==6 & Treatment=="C")
MEScdata["LTreatment"] = ifelse(MEScdata$Treatment=="C", "SG",ifelse(MEScdata$Treatment=="F", "IG","SIG" ) )

# Harvest data 2015
harvest2015 <- read_csv("Data/Cage_data.csv") %>%
  rename(Mf = MF_Sept, Treatment2 = Treatment) %>% 
  left_join(data.frame(Treatment2 = c("A", "B", "C"),
                       Treatment = c("IG", "SIG", "SG")))

# Harvest data 2017
harvest2017 <- read_csv("Data/Harvest_SU17_23Sept2017.csv")

# Harvest data 2018

harvest2018 <- read_csv("Data/longterm_data_2018.csv") %>% 
  select(Block:Value)

# collect the harvest data

# Block Treatment    Mf  Year Cage  SpHt  MStart

harvest <- MEScdata %>% separate(Plot, into=c("Block", NA), sep=1) %>% 
  select(Block,LTreatment, Mf) %>% rename(Treatment=LTreatment) %>%
  mutate(Year = 2013, Cage = paste(seq(1,15)), SpHt = NA, MStart = 5, Temp="C") %>%
  bind_rows(
    # Add 2015 data
    harvest2015  %>%
      select(Treatment, Cage, Mf) %>% 
      mutate(Cage = paste(Cage), Year=2015, Block = "1", SpHt = NA, MStart=2, Temp="C"),
    # Add 2017 data
    harvest2017 %>% select(-Spider: -Notes, -ugN_gDMES) %>% 
      rename(Mf = Hopper, Treatment = LTreatment) %>%
      mutate(Block = paste(Block),Cage = paste(Cage), 
             Year=2017, SpHt = NA, MStart = 2, Temp="C"),
    # Add 2018 data
    harvest2018 %>% spread(key=Type, value=Value) %>% 
      select(Block, Treatment, GrLive, SpHt) %>%
      separate(Treatment, into=c("Treatment", "Temp"), sep=" ") %>%
      rename(Mf = GrLive) %>%
      mutate(Year = 2018, MStart = 3, Block = paste(Block)) %>%
      filter(Mf != 4) %>% # remove clear typo
      full_join(read_csv("cageID2018.csv") %>% mutate(Block=paste(Block))) %>%
      mutate(Cage= paste(Cage))
    
  ) 

ublock = unique(paste0(harvest$Year, harvest$Block))

harvest <- harvest %>% 
  filter(Treatment %in% c("SIG", "SG", "IG")) %>%
  mutate(ublock = paste0(Year, Block)) %>%
  left_join(
    data.frame(ublock,
               block = paste(seq(1, length(ublock))))
  )

# behaviour data from 2015
behaviour <-read_csv("Data/behaviour_dataR.csv") %>% select(-X1) %>%
  gather(-Time, -Cage, -Doy, -Observer, key=SpMeas, value=Value) %>%
  separate(SpMeas, into=c("Species", "Measure"), sep=1) %>%
  separate(Measure, into=c("ID", "Measure"), sep=1) %>%
  spread(key=Measure, value= Value) %>% 
  mutate(Species = ifelse(Species=="G", "MEFE", "PIMI")) %>%
  rename(Beh = B, Perch = PL) %>%
  mutate(z = as.numeric(V)*5, y = as.numeric(H)*5) %>% #rescale height and width to cm (2" by 2" counting squares)
  select(-V, -H, -ID)  %>% filter(Doy == 217) %>%
  mutate(Block = "1", Year = 2015, x=NA) %>%
  left_join( # get treatment information from the harvest data set
    harvest2015 %>% select(Cage, Treatment)
    ) %>% mutate(Cage = paste(Cage)) %>% filter(!is.na(z)) %>% bind_rows(
      # behaviour data from 2017
      read_csv("Data/Data_Sheet_Behaviour_FinalFinal.csv") %>% 
        filter(Animal %in% c("MEFE", "PIMI") & !is.na(z)) %>%
        select(-Time, -Hour, -Minute, -Disturbed, -Notes) %>% 
        rename(Species = Animal) %>% 
        left_join( # Add treatment information from the harvest data set
          harvest2017 %>% select(Block, Cage, LTreatment) %>% rename(Treatment = LTreatment)
        ) %>%
        mutate(Year = 2017,Doy = 209, Block = paste(Block))
    ) %>% mutate(Temp = "C") %>% bind_rows(
      # add 2018 behaviour data
      read_csv("Data/grasshopper_Stefanie.csv") %>% 
        separate(Time, into=c("Hour", "Minute", "Second"),sep=":") %>%
        mutate(Hour = as.numeric(Hour), Minute= as.numeric(Minute)/60) %>%
        mutate(Time = Hour + Minute) %>% select(-Hour: -Second) %>%
        gather(-Block, -Obs, -Treatment, - Time, -Cage, key=SpM, value=value) %>%
        separate(SpM, into=c("Measure", "ID"), sep=-1) %>% 
        spread(key=Measure, value=value) %>%
        filter(!is.na(Y)) %>%
        rename(Beh= Behavior, Perch = Location, Observer = Obs, Temp = Treatment) %>%
        mutate(Block = paste(Block), Cage= paste(Cage), 
               z = as.numeric(Y), # height
               x = as.numeric(X),
               y = as.numeric(Z)) %>%
        mutate(Species = ifelse(ID=="s", "PIMI", "MEFE"), Doy = 212, Year = 2018) %>%
        select(-X, -Y, -Z, -ID) %>% left_join(
          data.frame(Cage = paste(seq(1,16)),
                     Treatment = c("SIG", "SG", "SIG", "SG",
                                   "SIG", "SG", "SIG", "SG", # IS THIS RIGHT OR A TYPO???
                                   "SIG", "SG", "SIG", "SG",
                                   "SIG", "SG", "SIG", "SG"))
        ) %>%
        mutate(z = z*2.54,
               x = x*2.54,
               y = y*2.54)
      
    ) %>%
  bind_rows(
    read_csv("Data/Short Term Grasshopper Data.csv") %>%
      rename(X_3 = X, Y_3 = Y, Z_3 = Z, Location_3 = Location, 
             Behavior_3=Behavior) %>%
      mutate(AMPM = rep(c(rep("AM", 6), rep("PM", 7)),27)) %>%
      full_join(
        read_csv("Data/Short Term Spider Data.csv") %>%
          rename(X_s = X, Y_s = Y, Z_s = Z, Location_s = Location, 
                 Behavior_s=Behavior) %>%
          mutate(AMPM = rep(c(rep("AM", 6), rep("PM", 7)),27)) %>%
          select(-Notes)
      ) %>% 
      mutate(Treatment = ifelse(Isopod =="Y", "SIG", "SG")) %>%
      left_join(
        data.frame(Cage = seq(1,44),
                   Block = rep(seq(1, 11), each=4))
      ) %>%
      select(-Isopod, -Notes) %>%
      separate(Time, into=c("Hour", "Minute", "Second"),sep=":") %>%
      mutate(Hour = as.numeric(Hour), Minute= as.numeric(Minute)/60) %>%
      mutate(Hour = ifelse(AMPM=="AM", Hour, Hour + 12)) %>%
      mutate(Time = Hour + Minute) %>% select(-Hour: -Second, -AMPM) %>%
      rename(Temp = Temperature) %>%
      gather(-Block,-Cage, -Temp, -Treatment, -Time, key=Measure, value=Value) %>%
      separate(Measure, into=c("Measure", "ID"), sep="_") %>%
      # distinct() %>%
      # filter(!is.na(Value)) %>%
      # mutate(uniqueID = paste0(Block, Cage, Time, ID)) %>%
      spread(key=Measure, value=Value) %>%
      mutate(Block = paste(Block),
             Cage= paste(Cage), 
             z = as.numeric(Y)*2.54, # height
             x = as.numeric(X)*2.54,
             y = as.numeric(Z)*2.54,
             Species = ifelse(ID=="s", "PIMI", "MEFE"),
             Doy = 212, Year = 2018) %>%
      select(-X, -Y, -Z, -ID) %>%
      mutate(Observer = "Unknown") %>%
      rename(Perch = Location, Beh = Behavior)
  )

# The NA warning here is OK, just NA -> NA in transformation

# check for errors due to data loss

behaviour %>% 
  select(Year,Time,Cage, Species, z) %>%
  group_by(Year,Time,Cage, Species) %>%
  summarize(z = sum(z, na.rm=T)) %>%
  ungroup() %>%
  spread(key=Species, value=z) %>%
  filter(MEFE==0 & PIMI==0) %>%
  filter(Time > 7) %>%
  select(Cage, Year) %>%
  table()

# Plot resulting data -----------------------------------------------------

# library(nlme)
# model1 <- lme(Mf ~ Treatment*Year, random=~1|block, data=harvest)
# summary(model1)

# plot(Mf~Year, data=harvest)


# Remove all but the product data frames from workspace -----
rm(list=ls()[!(ls() %in% c("harvest", "behaviour"))])

