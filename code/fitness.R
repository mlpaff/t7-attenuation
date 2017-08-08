library(tidyverse)
library(cowplot)
library(openxlsx)
library(stringr)

setwd("/Users/Matt/Desktop/Bull_Lab/T7Attenuation/Promoter")

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
                   levels=c('T7Hi', '10deop', '910v2', '8st_9i2', '8st_910', '44_910'))

back_labs <- c(
  'T7Hi' = 'wt',
  '8st'=expression(paste('8'[Delta]['stop'])),
  '10deop'=expression(paste('10'['deop']))
)
back_colors <- c(
  'T7Hi' = "orange",
  '8st' = "purple",
  '10deop'='lightgreen'
)
back_shapes <- c(
  'T7Hi' = 15,
  '8st' = 16,
  '10deop'= 17
)
line_labs <- c(
  'T7Hi' = 'wt',
  '10deop' = '10deop',
  '910v2' = expression(paste(Delta, phi, '9/10'['wt'])),
  '8st_9i2' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])), 
  '8st_910' = expression(paste(Delta, phi, '9/10'['8'[Delta]['stop']])),
  '44_910' = expression(paste(Delta, phi, '9/10'['10'[deop]]))
)

# Initial fitness -- split by knockouts
init_fit_plot <- fit %>% filter(!grepl('evo', Strain)) %>% 
  ggplot(aes(x=Background, y=Fitness, group=Strain)) +
  stat_summary(aes(color=Background, shape=Background), geom="point", fun.y='mean', size=5) + 
  stat_summary(geom="errorbar", fun.data='mean_se', width=0.25) + 
  scale_x_discrete(labels=back_labs) + 
  scale_color_manual(labels=back_labs,
                     values=back_colors) + 
  scale_shape_manual(labels=back_labs,
                     values=back_shapes) +
  labs(x = 'Genetic Background', y="Fitness (doublings/hour)") +
  scale_y_continuous(limits=c(25, 50), breaks=c(25, 30, 35, 40, 45, 50)) +
  facet_wrap(~Knockout, labeller = label_parsed, scales='free_x', nrow=1) + 
  panel_border() + 
  theme(legend.position='bottom',
        legend.title=element_blank())

init_fit_plot
save_plot("fit_init.pdf", init_fit_plot, base_width=8)

# Plot fitnesses for evolved lines. 
evo_fit <- fit %>% filter(Line %in% c('910v2', '8st_9i2', '8st_910', '44_910') | Strain=='T7Hi') %>%
  ggplot(aes(x=Line, y=Fitness, group=Strain)) + 
  stat_summary(aes(color=Background, shape=Background), geom='point', fun.data='mean_se', size=5) + 
  stat_summary(geom='errorbar', fun.data='mean_se', width=0.2) +
  geom_segment(aes(x = 2, y = 32.0733, xend = 2, yend = 36.1367), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') + 
  geom_segment(aes(x = 3, y = 32.338, xend = 3, yend = 38.0375), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') + 
  geom_segment(aes(x = 4, y = 28.88, xend = 4, yend = 35.6275), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') + 
  geom_segment(aes(x = 5, y = 30.473, xend = 5, yend = 32.295), arrow = arrow(length=unit(0.35, "cm")), size=0.6, color='black') + 
  scale_x_discrete(labels=line_labs) +
  scale_color_manual(values=back_colors) + 
  scale_shape_manual(values=back_shapes) + 
  labs(x = 'Strain', y="Fitness (doublings/hour)") +
  scale_y_continuous(limits=c(25, 50), breaks=c(25, 30, 35, 40, 45, 50)) +
  theme(legend.position='none')
  
evo_fit
save_plot("fit_evo.pdf", evo_fit, base_height=6, base_width=6)
