# Analysis of the HOBO data collected in 2018

library(tidyverse)
library(lubridate)

# Load in the data
LF44 = read_csv("Data/HOBO.csv") %>% mutate(Cage = as.character(Cage))

# Add the treatment identifiers from other database
LF44b = LF44  %>%
  left_join(
    read_csv("Data/Short Term Spider Data.csv") %>%
      select(Isopod, Temperature, Cage) %>%
      distinct() %>%
      mutate(Cage = as.character(Cage))
  ) %>% mutate(ID = ifelse(Temperature == "W", "Warmed", "Control"),
               ID2 = ifelse(Isopod == "N", "Woodlouse", "No woodlouse"))

# Add averages and sd at Cage level
etavg = LF44b %>% group_by(Cage, Temperature) %>%
  summarize(Temp = mean(Temp)) %>% ungroup() %>% group_by(Temperature) %>%
  summarize(sd = sd(Temp), avg = mean(Temp)) %>%
  mutate(lower = avg - sd, upper= avg + sd) %>% select(-sd)


LF44b = LF44b %>% left_join(etavg)

# Create experimental data for Aug 9th
LF44expt = LF44b %>% filter(Date > mdy_hms("08/09/18 07:00:00") &
                              Date < mdy_hms("08/09/18 19:00:00"))



exptavg = LF44expt %>% group_by(Cage, Temperature) %>%
  summarize(Temp = mean(Temp)) %>% ungroup() %>% group_by(Temperature) %>%
  summarize(sd = sd(Temp), avg = mean(Temp))  %>%
  mutate(lower = avg - sd, upper= avg + sd) %>% select(-sd)

LF44expt = LF44expt %>% left_join(exptavg)

# Plot the results
p1 = LF44b %>% 
  ggplot(aes(x=Date, color = ID)) +  
  geom_line(aes(y = avg), size = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) + 
  geom_line(aes(y = Temp, group = Cage, linetype = ID2)) + theme_bw() + theme(legend.position = "top") + ylab("Temperature") + scale_linetype_discrete(name = "") +
  scale_color_manual(name = "", values = c("blue", "red"))

p2 = LF44expt %>% 
  ggplot(aes(x=Date, color = ID)) +  
  geom_line(aes(y = avg), size = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) + 
  geom_line(aes(y = Temp, group = Cage, linetype = ID2)) + theme_bw() + theme(legend.position = "top") + ylab("Temperature") + scale_linetype_discrete(name = "") +
  scale_color_manual(name = "", values = c("blue", "red")) + xlab("Time")

png("Plots/HOBO.png", width = 15, height = 8, units = "in", res = 600)
ggpubr::ggarrange(p1,p2, common.legend = T, labels= "AUTO")
dev.off()
