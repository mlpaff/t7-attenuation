library(tidyverse)
library(cowplot)
library(openxlsx)
library(stringr)

setwd("/Users/Matt/Desktop/Bull_Lab/T7Attenuation/Promoter")

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

fit <- read.xlsx("Fitness.xlsx")

fit$Background <- factor(fit$Background, 
                         levels=c('T7Hi', '10deop', '8st'))

fit$Knockout <- factor(fit$Knockout, 
                       levels=c('wt', 'phi9', 'phi10', 'phi910'), 
                       labels=c(
                         'wt' = 'wt',
                         'phi9' = expression(paste(Delta, phi, '9')),
                         'phi10' = expression(paste(Delta, phi, '10')),
                         'phi910' = expression(paste(Delta, phi, '9/10')))
                       )

fit$Line <- factor(fit$Line,
                   levels=c('T7Hi', '10deop', '910v2', '910v2_L2', '44_910', '8st_9i2', '8st_910'))

back_labs <- c(
  'wt' = 'wt',
  '10deop'=expression(paste('10'['deop'])),
  '8st'=expression(paste('8'[Delta]['stop']))
)
back_colors <- c(
  'T7Hi' = "#E69F00",
  '8st' = "#56B4E9",
  '10deop'='#009E73'
)
back_shapes <- c(
  'T7Hi' = 15,
  '8st' = 16,
  '10deop'= 17
)
line_labs <- c(
  'T7Hi' = 'wt',
  '10deop' = expression(paste('10'['deop'])),
  '910v2' = expression(paste(Delta, phi, '9/10'['wt'])),
  '910v2_L2' = expression(atop(paste(Delta,phi,'9/10'['wt']), '\nL2')),
  '44_910' = expression(paste(Delta, phi, '9/10'['10'[deop]])),
  '8st_9i2' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])), 
  '8st_910' = expression(paste(Delta, phi, '9/10'['8'[Delta]['stop']]))
)

# Initial fitness -- split by knockouts
init_fit_plot <- fit %>% filter(!grepl('evo', Strain), !grepl('L2', Strain)) %>% 
  ggplot(aes(x=Background, y=Fitness, group=Strain)) +
  stat_summary(aes(color=Background, shape=Background), geom="point", fun.y='mean', size=5) + 
  stat_summary(geom="errorbar", fun.data='mean_se', width=0.3) + 
  scale_x_discrete(labels=back_labs) + 
  scale_color_manual(labels=back_labs,
                     values=cbPalette) + 
  scale_shape_manual(labels=back_labs,
                     values=back_shapes) +
  labs(x = 'Genetic Background', y="Fitness (doublings/hour)") +
  scale_y_continuous(limits=c(25, 50), breaks=c(25, 30, 35, 40, 45, 50)) +
  facet_grid(~Knockout, labeller = label_parsed, scales='free_x', space='free_x') + 
  panel_border() + 
  theme(legend.position='right',
        legend.title=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(1.5, 'lines'))

init_fit_plot
save_plot("fit_init.pdf", init_fit_plot, base_width=9)

# Plot fitnesses for evolved lines. 
evo_fit <- fit %>% filter(Line %in% c('910v2', '910v2_L2', '8st_9i2', '8st_910', '44_910') | Strain==c('T7Hi', '11--44')) %>%
  ggplot(aes(x=Line, y=Fitness, group=Strain)) + 
  stat_summary(aes(color=Background, shape=Background), geom='point', fun.data='mean_se', size=5) + 
  stat_summary(geom='errorbar', fun.data='mean_se', width=0.2) +
  geom_segment(aes(x = 3, y = 32.0733, xend = 3, yend = 36.1367), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') +
  geom_segment(aes(x = 3.1, y = 33.41, xend = 3.85, yend = 33.41), arrow = arrow(length=unit(0.35, "cm")), size=1, color='black', linetype=3) +
  geom_segment(aes(x = 4, y = 33.96, xend = 4, yend = 35.1033), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') +
  geom_segment(aes(x = 5, y = 30.473, xend = 5, yend = 32.295), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') +
  geom_segment(aes(x = 6, y = 32.338, xend = 6, yend = 38.0375), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') +
  geom_segment(aes(x = 7, y = 28.88, xend = 7, yend = 35.6275), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') +
  scale_x_discrete(labels=line_labs) +
  scale_color_manual(labels=back_labs,
                     values=cbPalette) + 
  scale_shape_manual(labels=back_labs,
                     values=back_shapes) +
  labs(x = 'Lines', y="Fitness (doublings/hour)") +
  scale_y_continuous(limits=c(25, 50), breaks=c(25, 30, 35, 40, 45, 50)) +
  theme(legend.position='right',
        legend.title=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(1.5, 'lines'),
        axis.text.x = element_text(angle = 45, hjust=1))

evo_fit
save_plot("fit_evo.pdf", evo_fit, base_height=5, base_width=8)

# fit %>% filter(!grepl('evo', Strain)) %>% group_by(Strain) %>% summarize(mean_fit=mean(Fitness))

# Plot fitness relative to WT
fit %>% filter(!grepl('evo', Strain)) %>% group_by(Strain) %>% 
  summarize(Background=Background[1], mean_fit=mean(Fitness)) %>% 
  mutate(rel_fit=mean_fit/max(mean_fit)) %>%
  ggplot(aes(x=Strain, y=rel_fit, fill=Background)) + geom_bar(stat="identity") + 
  scale_x_discrete(labels=line_labs) + 
  scale_color_manual(labels=back_labs,
                     values=back_colors) + 
  scale_shape_manual(labels=back_labs,
                     values=back_shapes) +
  labs(x = 'Genetic Background', y="Relative fitness") +
  panel_border() + 
  theme(legend.position='right',
        legend.title=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(1.5, 'lines'))

# T tests to test for significance in fitness differences
t.test(Fitness~Strain, data=subset(fit, Strain %in% c("T7Hi", "9i2")))
t.test(Fitness~Strain, subset(fit, Strain %in% c("11-44-910v2", "11-44-910v2evo")))
