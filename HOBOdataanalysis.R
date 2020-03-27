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

mac = data.frame(ID = c("10731719", "10731721","10731724","10731727",
                        "10731723", "10731720", "10731729", "10731729_0", "10731730",
                        "cage20","cage21","cage25","cage26","cage27",
                        "cage29","cage35","cage38","cage44"),
                 Cage = c("25","26","29","38","20","27", "43", "43","44",
                          "20", "21", "25","26", "27", "29", "35", "38", "44"))

LF44 = LF44 %>% separate(ID, into = c("ID", NA), sep = -4) %>%
  left_join(mac) %>%
  filter(!is.na(Cage)) %>% select(-ID)

write_csv(LF44, "Data/HOBO.csv")



# Starting analyiss with LF44
LF44 = read_csv("Data/HOBO.csv") %>% mutate(Cage = as.character(Cage))

LF44 %>% ggplot(aes(x=Date, y = Temp, color = Cage)) + geom_line() + theme_bw()

# Add the treatment identifiers from other database
LF44b = LF44  %>%
  left_join(
    read_csv("Data/Short Term Spider Data.csv") %>%
      select(Isopod, Temperature, Cage) %>%
      distinct() %>%
      mutate(Cage = as.character(Cage))
  ) %>% mutate(ID = ifelse(Temperature == "W", "Warmed", "Control"),
               ID2 = ifelse(Isopod == "N", "Woodlouse", "No woodlouse"))

LF44expt = LF44b %>% filter(Date > mdy_hms("08/09/18 07:00:00") &
                              Date < mdy_hms("08/09/18 19:00:00"))

etavg = LF44b %>% group_by(Cage, Temperature) %>%
  summarize(Temp = mean(Temp)) %>% ungroup() %>% group_by(Temperature) %>%
  summarize(sd = sd(Temp), avg = mean(Temp)) %>%
  mutate(lower = avg - sd, upper= avg + sd) %>% select(-sd)


LF44b = LF44b %>% left_join(etavg)

exptavg = LF44expt %>% group_by(Cage, Temperature) %>%
  summarize(Temp = mean(Temp)) %>% ungroup() %>% group_by(Temperature) %>%
  summarize(sd = sd(Temp), avg = mean(Temp))  %>%
  mutate(lower = avg - sd, upper= avg + sd) %>% select(-sd)

LF44expt = LF44expt %>% left_join(exptavg)

p1 = LF44b %>% 
  ggplot(aes(x=Date, color = ID)) +  
  geom_line(aes(y = avg), size = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) + 
  geom_line(aes(y = Temp, group = Cage, linetype = ID2)) + theme_bw() + theme(legend.position = "top") + ylab("Temperature") + scale_linetype_discrete(name = "") +
  scale_color_discrete(name = "")

p2 = LF44expt %>% 
  ggplot(aes(x=Date, color = ID)) +  
  geom_line(aes(y = avg), size = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) + 
  geom_line(aes(y = Temp, group = Cage, linetype = ID2)) + theme_bw() + theme(legend.position = "top") + ylab("Temperature") + scale_linetype_discrete(name = "") +
  scale_color_discrete(name = "") + xlab("Time")

png("Plots/HOBO.png", width = 15, height = 8, units = "in", res = 600)
ggpubr::ggarrange(p1,p2, common.legend = T, labels= "AUTO")
dev.off()
