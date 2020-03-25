# Analysis of the HOBO data collected in 2018

library(tidyverse)
library(lubridate)

# Extract the HOBO data

d1 <- read.csv("Data/HOBOdata/20.csv", skip = 1)
d1 <- d1[,2:3]
colnames(d1) <- c("Date", "Temp")
d1$Date = mdy_hms(d1$Date)
d1[,"ID"] = "20.csv"

plot(Temp~Date, data= d2, type = "b")

d2 <- read.csv("Data/HOBOdata/10731729.csv", skip = 1)
d2 <- d2[,2:3]
colnames(d2) <- c("Date", "Temp")
d2$Date = mdy_hms(d2$Date)
d2$ID = "10731729.csv"

lf <- list.files("Data/HOBOdata/tempdata080918/")

lf2 <- vector(mode = "list", length = length(lf))

for(i in 1:length(lf)){
  
  da <- read.csv(paste0("Data/HOBOdata/tempdata080918/",lf[i]), skip = 1)
  da <- da[,2:3]
  colnames(da) <- c("Date", "Temp")
  da$Date = mdy_hms(da$Date)
  da[,"ID"] = lf[i]
  
  lf2[[i]] = da
}

View(lf2)

lf2[[24]] = d1
lf2[[25]] = d2

# Add them all together:

LF <- do.call("rbind",lf2)

LF %>% as_tibble() %>% filter(Date > mdy_hms("08/01/17 08:00:00")) %>% ggplot(aes(x=Date, y = Temp, color = ID)) + geom_line()


LF %>% filter(Date < mdy_hms("08/01/17 08:00:00")) %>% dim()
LF %>% filter(Date > mdy_hms("08/01/17 08:00:00")) %>% dim()


# Load in the 44 data

lf44 <- list.files("Data/HOBOdata/44/")

lf442 <- vector(mode = "list", length = length(lf))

for(i in 1:length(lf)){
  
  da <- read.csv(paste0("Data/HOBOdata/44/",lf[i]), skip = 1)
  da <- da[,2:3]
  colnames(da) <- c("Date", "Temp")
  da$Date = mdy_hms(da$Date)
  da[,"ID"] = lf[i]
  
  lf442[[i]] = da
}

LF44 <- do.call("rbind",lf442)

LF44 %>% as_tibble() %>% ggplot(aes(x=Date, y = Temp, color = ID)) + geom_line()

uID = unique(LF44$ID)

for(i in 1:length(uID)){
  print(paste0(uID[i], " is ", all(LF44 %>% filter(ID == uID[i]) == LF %>% filter(ID == uID[i]), na.rm = T)))
}

# LF44 and LF are the same, so only need to use 1
# LF has one extra dataset that LF44 doesn't have and it is 20.csv

# 20.csv is the same as cage20.csv, so use LF44
LF %>% as_tibble() %>% filter(ID %in% c("cage20.csv", "20.csv")) %>% ggplot(aes(x=Date, y = Temp, color = ID, linetype = ID)) + geom_line() + theme_bw()

LF44 = LF44 %>% as_tibble() %>% filter(!is.na(Temp))

write_csv(LF44, "Data/HOBO.csv")

# Starting analyiss with LF44
LF44 = read_csv("Data/HOBO.csv")



LF44 %>% ggplot(aes(x=Date, y = Temp, color = ID)) + geom_line() + theme_bw()

LF44b = LF44[str_detect(LF44$ID, "cage"),]

LF44b = LF44b %>% separate(ID, into=c(NA,"ID", NA), sep=c(4,-4)) %>%
  rename(Cage = ID) %>%
  mutate(Cage = as.numeric(Cage)) %>%
  left_join(
    read_csv("Data/Short Term Spider Data.csv") %>%
      select(Isopod, Temperature, Cage) %>%
      distinct() 
  ) %>% mutate(ID = ifelse(Temperature == "W", "Warmed", "Control"),
               ID2 = ifelse(Isopod == "N", "Woodlouse", "No woodlouse"))

p1 = LF44b %>% ggplot(aes(x=Date, y = Temp, color = ID, linetype = ID2, group = Cage)) + geom_line() + theme_bw() + theme(legend.position = "top") + ylab("Temperature")


p2 = LF44 %>% ggplot(aes(x=Date, y = Temp, color = ID)) + geom_line() + theme_bw() + ylab("Temperature") + theme(legend.position = "top")

png("Plots/HOBO.png", width = 15, height = 8, units = "in", res = 600)
ggpubr::ggarrange(p1,p2)
dev.off()
