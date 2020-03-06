# Compare the results in 1-D and 3-D calculations

comparedim <- behaviour %>% 
  filter(Treatment %in% c("SIG", "SG")) %>% 
  mutate(x = replace_na(x, 0)) %>% # put zero for all unmeasured x
  gather(-Time:-Perch, -Block, -Year, -Treatment, -Temp,
         key = Dim, value=Height) %>%
  mutate(Spdim = paste0(Species, "_", Dim)) %>%
  select(-Species, -Dim, -Observer, -Beh, -Perch) %>%
  distinct() %>%
  group_by(Cage, Block, Year, Treatment, Temp, Time, Spdim) %>%
  summarise(Height = mean(Height, na.rm=T)) %>%
  spread(key=Spdim, value=Height) %>%
  mutate(Euclidean3D = sqrt(
    (MEFE_x - PIMI_x)^2 +
      (MEFE_y - PIMI_y)^2 +
      (MEFE_z - PIMI_z)^2
  ),
  Euclidean1D = sqrt((MEFE_z - PIMI_z)^2)
  ) %>%
  select(-MEFE_x: -PIMI_z) %>%
  group_by(Cage, Block, Year, Treatment, Temp) %>%
  summarise(Euclidean3D = mean(Euclidean3D, na.rm=T),
            Euclidean1D = mean(Euclidean1D, na.rm=T)) %>%
  filter(Year != 2016)

summary(lm(Euclidean3D~Euclidean1D, data= comparedim))

png("suppfigure1_16May2019.png", width=8, height=8, units="in", res=600)
comparedim %>%
  ggplot(aes(x=Euclidean1D, y=Euclidean3D)) + 
  stat_smooth(method="lm") +
  geom_point(aes(color=Block, shape=as.factor(Year))) + theme_classic() + 
  xlab("Euclidean (1D)") + xlim(0, 70) +
  ylab("Euclidean (3D)") + ylim(0, 70) +
  scale_shape_discrete(name="Year") +
  annotate("text", y = 15, x = 50, label = expression(R[adj]^2==0.72))
dev.off()

comparedim %>% filter(Euclidean3D > 50 & Euclidean1D < 20) %>%
  left_join(behaviour) %>% View()
