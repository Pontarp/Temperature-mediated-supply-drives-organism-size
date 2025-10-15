require(ggplot2)
####Plots for the paper####

finalmatrix2 <- read.csv("BEEP_Poly2024-2.csv")
finalmatrix2$Div <- interaction(finalmatrix2$Div)

#Specialists plot
minplot <- ggplot(data = finalmatrix2, aes(x = tempvect, y = minvect, col = factor(Div))) +
  # geom_point() +
  xlab("Temperature") +
  ylab("Log Body Size") +
  theme_bw() +
  # geom_point(data = finalmatrix2[(simchosen+1):(simchosen+216),], aes(x = tempvect, y = minvect, col = factor(Div), shape=factor(MinSeed))) +
  geom_point(data = finalmatrix2[1:36,], aes(x = tempvect+0.05, y = minvect, col = factor(Div)), shape = 1) +
  geom_point(data = finalmatrix2[613:648,], aes(x = tempvect-0.05, y = minvect, col = factor(Div)), shape = 2) +
  geom_point(data = finalmatrix2[181:216,], aes(x = tempvect, y = minvect, col = factor(Div)), shape = 3) +
  geom_point(data = finalmatrix2[109:144,], aes(x = tempvect-0.05, y = minvect, col = factor(Div)), shape = 4) +
  geom_point(data = finalmatrix2[325:360,], aes(x = tempvect-0.05, y = minvect, col = factor(Div)), shape = 5) +
  geom_point(data = finalmatrix2[1297:1332,], aes(x = tempvect, y = minvect, col = factor(Div)), shape = 6) +
  # geom_point(data = finalmatrix2[216:252,], aes(x = tempvect, y = minvect, col = factor(Div)), shape = 1) +
  # geom_point(data = finalmatrix2[252:288,], aes(x = tempvect, y = minvect, col = factor(Div)), shape = 1) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.25), 
        legend.position = "none") +
  # ggtitle("Specialist Community") +
  scale_color_discrete(name = "Diversity") +
  stat_smooth(method = "lm", se = F, size = 0.5) 
minplot

#Generalists plot
maxplot <- ggplot(data = finalmatrix2, aes(x = tempvect, y = maxvect, col = factor(Div))) + 
  # geom_point() +
  xlab("Temperature") +
  ylab("Log Body Size") +
  theme_bw() +
  geom_point(data = finalmatrix2[1:36,], aes(x = tempvect+0.05, y = maxvect, col = factor(Div)), shape = 1) +
  geom_point(data = finalmatrix2[37:72,], aes(x = tempvect-0.05, y = maxvect, col = factor(Div)), shape = 2) +
  geom_point(data = finalmatrix2[73:108,], aes(x = tempvect, y = maxvect, col = factor(Div)), shape = 3) +
  geom_point(data = finalmatrix2[109:144,], aes(x = tempvect-0.05, y = maxvect, col = factor(Div)), shape = 4) +
  geom_point(data = finalmatrix2[146:180,], aes(x = tempvect-0.1, y = maxvect, col = factor(Div)), shape = 5) +
  geom_point(data = finalmatrix2[181:216,], aes(x = tempvect, y = maxvect, col = factor(Div)), shape = 6) +
  # geom_point(data = finalmatrix2[1:216,], aes(x = tempvect, y = maxvect, col = factor(Div), shape=factor(MaxSeed))) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  # ggtitle("Generalist Community") +
  scale_color_discrete(name = "Diversity") +
  stat_smooth(method = "lm", se = F, size = 0.5) 
maxplot

finalmatrixniche <- read.csv("BEEP_Final-communityniche.csv")
finalmatrixniche$Div <- interaction(finalmatrixniche$Div)

#Specialists plot with full partitioning
nicheplot <- ggplot(data = finalmatrixniche, aes(x = tempvect, y = nichevect, col = Div)) +
  # geom_point() +
  # geom_point(data = finalmatrixniche[1:216,], aes(x = tempvect, y = nichevect, col = factor(Div), shape=factor(Sim))) +
  geom_point(data = finalmatrixniche[1:36,], aes(x = tempvect+0.05, y = nichevect, col = factor(Div)), shape = 1) +
  geom_point(data = finalmatrixniche[253:288,], aes(x = tempvect-0.05, y = nichevect, col = factor(Div)), shape = 2) +
  geom_point(data = finalmatrixniche[73:108,], aes(x = tempvect, y = nichevect, col = factor(Div)), shape = 3) +
  geom_point(data = finalmatrixniche[109:144,], aes(x = tempvect-0.05, y = nichevect, col = factor(Div)), shape = 4) +
  geom_point(data = finalmatrixniche[289:324,], aes(x = tempvect-0.1, y = nichevect, col = factor(Div)), shape = 5) +
  geom_point(data = finalmatrixniche[181:216,], aes(x = tempvect, y = nichevect, col = factor(Div)), shape = 6) +
  xlab("Temperature") +
  ylab("Log Body Size") +
  theme_bw() +
  stat_smooth(method = "lm", se = F, size = 0.5) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  # ggtitle("Specialist Niche-Partition Community") +
  scale_color_discrete(name = "Div") +
  theme(legend.position = "none") 
nicheplot

finalplot <- plot_grid(nicheplot, minplot, maxplot, nrow = 1,
                       label_size = 12, labels = c("A)", "B)", "C)"),
                       label_y = 0.75, label_x = 0.8)

