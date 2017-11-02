library(tidyverse)
library(cowplot)
library(openxlsx)
library(stringr)
library(broom)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

fit <- read.xlsx("../../data/results/Fitness.xlsx")

fit$Background <- factor(fit$Background, 
                         levels=c('T7Hi', '10deop', '8st'))
fit$Line <- factor(fit$Line,
                   levels=c('T7Hi', '10deop', '910v2', '910v2_L2', '44_910', '8st_9i2', '8st_910'))

# Summarize mean and standard error for all strains
fit %>% group_by(Strain) %>% summarize_at(vars(Fitness), funs(mean,sd, se=sd/sqrt(n()))) -> mean_fit
mean_fit

# T tests to test for significance in fitness differences
fit_sig <- function(treatment, control){
  df <- data.frame(comp=paste0(control, '_', treatment),
                   strain=treatment, 
                   tidy(t.test(Fitness~Strain, data=subset(fit, Strain %in% c(control, treatment)))))
}

fit_sig_all <- function(strains, control){
  strains %>% map(fit_sig, control=control) %>%
    bind_rows()
}

strains.list <- as.character(unique(subset(fit, !Strain=="T7Hi")$Strain))
wt_comparison_fit <- fit_sig_all(strains.list, control='T7Hi')

wt_comparison_fit %>% select(strain, estimate, p.value)

deop_comparison_fit <- fit_sig_all(grep('44-', strains.list, value=T), control='11--44')
deop_comparison_fit %>% select(strain, estimate, p.value)




fit$Knockout <- factor(fit$Knockout, 
                       levels=c('wt', 'phi9', 'phi10', 'phi910'), 
                       labels=c(
                         'wt' = 'wt',
                         'phi9' = expression(paste(Delta, phi, '9')),
                         'phi10' = expression(paste(Delta, phi, '10')),
                         'phi910' = expression(paste(Delta, phi, '9/10')))
                       )



back_labs <- c(
  'T7Hi' = 'wt',
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
# line_labs <- c(
#   'T7Hi' = 'wt',
#   '10deop' = expression(paste('10'['deop'])),
#   '910v2' = expression(paste(Delta, phi, '9/10'['wt'])),
#   '910v2_L2' = expression(atop(paste(Delta,phi,'9/10'['wt']), '\nL2')),
#   '44_910' = expression(paste(Delta, phi, '9/10'['10'[deop]])),
#   '8st_9i2' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])), 
#   '8st_910' = expression(paste(Delta, phi, '9/10'['8'[Delta]['stop']]))
# )
line_labs <- c(
  'T7Hi' = 'wt',
  '10deop' = expression(paste('10'['deop'])),
  '910v2' = expression(paste(Delta, phi, '9/10')),
  '910v2_L2' = expression(paste(Delta,phi,'9/10')),
  '44_910' = expression(paste(Delta, phi, '9/10')),
  '8st_9i2' = expression(paste(Delta, phi, '9')), 
  '8st_910' = expression(paste(Delta, phi, '9/10'))
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
save_plot("../../figures/fit_init.pdf", init_fit_plot, base_width=9)

# Plot fitnesses for evolved lines. 
# evo_line <- data.frame(Line = c('T7Hi', '10deop', '44_910', '8st_910', '8st_9i2', '910v2', '910v2_L2'),
#                        evo_line = c('ref', 'ref', 'evo', 'evo', 'evo', 'evo', 'evo'))
# evo_line$evo_line <- factor(evo_line$evo_line,
#                             levels=c("ref", "evo"))

evo_fit <- fit %>% filter(Line %in% c('910v2', '910v2_L2', '8st_9i2', '8st_910', '44_910') | Strain==c('T7Hi', '11--44')) %>%
  ggplot(aes(x=Line, y=Fitness, group=Strain)) + 
  stat_summary(aes(color=Background, shape=Background), geom='point', fun.data='mean_se', size=5) + 
  stat_summary(geom='errorbar', fun.data='mean_se', width=0.2) +
  geom_segment(aes(x=1, xend=1, y=48.5, yend=49)) + geom_segment(aes(x=2, xend=2, y=49, yend=48.5)) +
  geom_segment(aes(x=1, xend=2, y=49, yend=49)) + annotate("text", x = 1.5, y = 50, label = "ancestor", size = 4) +
  geom_curve(aes(x = 2.95, y = 30.7, xend = 2.95, yend = 37.85), arrow = arrow(length=unit(0.3, "cm")), 
             size=0.5, color='grey', curvature=-0.45) +
  geom_segment(aes(x = 3, y = 30.7, xend = 3, yend = 32.61), arrow = arrow(length=unit(0.3, "cm")), size=0.5, color='black') +
  geom_segment(aes(x = 3, y = 34.1, xend = 3, yend = 35), arrow = arrow(length=unit(0.3, "cm")), size=0.5, color='black') +
  geom_segment(aes(x = 4, y = 30.45, xend = 4, yend = 32.6), arrow = arrow(length=unit(0.3, "cm")), size=0.5, color='black') +
  geom_segment(aes(x = 5, y = 31.5, xend = 5, yend = 38.7), arrow = arrow(length=unit(0.3, "cm")), size=0.5, color='black') +
  geom_segment(aes(x = 6, y = 28, xend = 6, yend = 36.8), arrow = arrow(length=unit(0.3, "cm")), size=0.5, color='black') +
  scale_x_discrete(labels=line_labs) +
  scale_color_manual(labels=back_labs,
                     values=cbPalette) + 
  scale_shape_manual(labels=back_labs,
                     values=back_shapes) +
  labs(x='promoter knockout lines', y="Fitness (doublings/hour)") +
  scale_y_continuous(limits=c(25, 50), breaks=c(25, 30, 35, 40, 45, 50)) +
  theme(legend.position='right',
        legend.title=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(1.5, 'lines'))

evo_fit
save_plot("../../figures/fit_evo.pdf", evo_fit, base_height=5, base_width=8)



# Plot fitness relative to WT
# fit %>% filter(!grepl('evo', Strain)) %>% group_by(Strain) %>% 
#   summarize(Background=Background[1], mean_fit=mean(Fitness)) %>% 
#   mutate(rel_fit=mean_fit/max(mean_fit)) %>%
#   ggplot(aes(x=Strain, y=rel_fit, fill=Background)) + geom_bar(stat="identity") + 
#   scale_x_discrete(labels=line_labs) + 
#   scale_color_manual(labels=back_labs,
#                      values=back_colors) + 
#   scale_shape_manual(labels=back_labs,
#                      values=back_shapes) +
#   labs(x = 'Genetic Background', y="Relative fitness") +
#   panel_border() + 
#   theme(legend.position='right',
#         legend.title=element_blank(), 
#         legend.text.align=0,
#         legend.key.size=unit(1.5, 'lines'))


