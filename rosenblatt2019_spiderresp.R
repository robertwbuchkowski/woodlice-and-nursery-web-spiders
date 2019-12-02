# Analysis of Rosenblatt spider respiration rates

rosenblatt = read_csv("Data/rosenblatt2019.csv") %>%
  filter(!is.na(`20C`)) %>%
  filter(Source == "CT") %>%
  select(-Source) %>%
  mutate(Block = letters[1:15]) %>%
  gather(-Block, -Sex, -Weight, key = Temp, value = Resp) %>%
  mutate(Resp = Resp*Weight) %>%
  select(-Weight, - Sex) %>%
  separate(Temp, into = c("Temp", NA), sep = 2) %>%
  mutate(Temp = as.numeric(Temp)) %>%
  filter(Block != "h")


summary(lme4::lmer(Resp~Temp + (1|Block), data = rosenblatt))
plot(lme4::lmer(Resp~Temp + (1|Block), data = rosenblatt))

# Respiratory quotient of 0.7 from A. Schmitz 2004
# oxycalorific equivalent taken from Ford 1977
# Congruence of the data suggest we can use Ford 1977 activity rates

x = seq(15, 40,1)
y = (0.019891*x - 0.206382)*60*0.7*0.0200832
plot(y~x, type = "l", ylim = c(0,1.5))

act = (30.724*x- 16.572)/1000
std = (9.692*x-5.228)/1000

points(act~x, col = "red", type ="l")
points(std~x, col = "blue", type ="l")

y2 = (act-std) + y # increase in respiration from Ford with different species, then add to rosenblatt data.

points(y2~x, col = "green", type ="l")

y.min = y/60
y2.min = y2/60

stdE = approxfun(x,y.min)
actE = approxfun(x,y2.min)

stdE(20)
actE(20)
