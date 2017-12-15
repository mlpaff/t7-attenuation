library(tidyverse)
library(cowplot)
library(openxlsx)
library(stringr)
library(broom)
library(colorspace)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

fit <- read.xlsx("../../data/results/Fitness.xlsx")

fit$Background <- factor(fit$Background, 
                         levels=c('T7Hi', '10deop', '8st'))
fit$Line <- factor(fit$Line,
                   levels=c('T7Hi', '10deop', '910v2', '910L2', '44_910', '8st_9i2', '8st_910'))

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

# significance in evolved lines
t.test(Fitness~Strain, data=subset(fit, Strain %in% c('910v2i2', '910L2-evo')))
t.test(Fitness~Strain, data=subset(fit, Strain %in% c('9i2', '9i2evo')))
t.test(Fitness~Strain, data=subset(fit, Strain %in% c('910i3', '910i3evo')))
t.test(Fitness~Strain, data=subset(fit, Strain %in% c('11-44-910v2', '11-44-910v2evo')))


fit$Knockout <- factor(fit$Knockout, 
                       levels=c('wt', 'phi9', 'phi10', 'phi910'), 
                       labels=c(
                         'wt' = 'wt',
                         'phi9' = expression(paste(Delta, phi, '9')),
                         'phi10' = expression(paste(Delta, phi, '10')),
                         'phi910' = expression(paste(Delta, phi, '9/', phi, '10')))
                       )

back_labs <- c(
  'T7Hi' = 'wt',
  '10deop'=expression(paste('10'['deop'])),
  '8st'=expression(paste('8'[Delta]['stop']))
)
back_colors <- c(
  'T7Hi' = "#E69F00",
  '10deop'="#56B4E9",
  '8st' = "#009E73"
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
  '910v2' = expression(paste(Delta, phi, '9/', phi, '10')),
  '910L2' = expression(paste(Delta,phi,'9/', phi, '10-L2')),
  '44_910' = expression(paste(Delta, phi, '9/', phi, '10')),
  '8st_9i2' = expression(paste(Delta, phi, '9')), 
  '8st_910' = expression(paste(Delta, phi, '9/', phi, '10'))
)

# Initial fitness -- split by knockouts
init_fit_plot <- fit %>% filter(!grepl('evo', Strain), !grepl('L2', Strain)) %>% 
  ggplot(aes(x=Background, y=Fitness, group=Strain)) +
  stat_summary(aes(fill=Background), geom="bar", fun.y='mean', size=5) + 
  stat_summary(geom="errorbar", fun.data='mean_se', width=0.3) + 
  scale_x_discrete(labels=back_labs) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,47)) +
  scale_fill_manual(labels=back_labs,
                     values=cbPalette) + 
  # scale_shape_manual(labels=back_labs,
  #                    values=back_shapes) +
  labs(x = 'genetic background', y="fitness (doublings/hour)") +
  #scale_y_continuous(limits=c(25, 50), breaks=c(25, 30, 35, 40, 45, 50)) +
  facet_grid(~Knockout, labeller = label_parsed, scales='free_x', space='free_x') + 
  panel_border() + 
  theme(legend.position='right',
        legend.title=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(1.5, 'lines'))

init_fit_plot
save_plot("../../figures/fit_init.pdf", init_fit_plot, base_width=9)

fit$Strain <- factor(fit$Strain,
                     levels=c('T7Hi', '11--44', '9v2i3', '10i1', '910v2i2', '910L2-init', '910L2-evo', '910v2evo',
                              '11-44-9v2i1', '11-44-10i2', '11-44-910v2', '11-44-910v2evo',
                              '9i2', '9i2evo', '910i3', '910i3evo'))
fit$Passaged <- factor(fit$Passaged, 
                       levels=c("Ancestor", "Evolved"),
                       labels=c("Ancestor strains", "Evolved lines"))
strain_fill <- c(
  '910v2i2' = '#E69F00',
  '910L2-evo' = desaturate("#E69F00", 0.7),
  '11-44-910v2' = '#56B4E9',
  '11-44-910v2evo' = desaturate('#56B4E9', 0.7),
  '9i2' = '#009E73',
  '9i2evo' = desaturate("#009E73", 0.7),
  '910i3' = '#009E73',
  '910i3evo' = desaturate("#009E73", 0.7))
line_labs2 <- c(
  '910v2' = expression(paste(Delta, phi, '9/', phi, '10'['wt'])),
  '910L2' = expression(paste(Delta,phi,'9/', phi, '10'['wt'],'-L2')),
  '44_910' = expression(paste(Delta, phi, '9/', phi, '10'['10'['deop']])),
  '8st_9i2' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])), 
  '8st_910' = expression(paste(Delta, phi, '9/', phi, '10'['8'[Delta]['stop']]))
)

evo_fit <- fit %>% filter(Line %in% c('910v2', '910L2', '8st_9i2', '8st_910', '44_910') | Strain==c('T7Hi', '11--44'), 
                          !Strain %in% c("910L2-init", "910v2evo", "T7Hi", "11--44")) %>%
  ggplot(aes(x=Line, y=Fitness, group=Strain)) +
  geom_segment(aes(x=-Inf, xend=Inf, y=45.08333, yend=45.08333), color='#E69F00') +
  annotate("text", x=1, y=46.08333, label="wt ancestor", size=3.5) + 
  stat_summary(aes(fill=Strain), geom='bar', position='dodge', fun.data='mean_se') +
  stat_summary(geom='errorbar', position=position_dodge(0.9), fun.data='mean_se', width=0.25) +
  scale_x_discrete(labels=line_labs2) +
  #scale_color_manual(values=back_colors) +
  scale_fill_manual(values=strain_fill) +
  labs(x='promoter knockout lines', y="fitness (doublings/hour)") +
  scale_y_continuous(expand=c(0,0), limits=c(0, 50)) +
  theme(legend.position='none',
        legend.title=element_blank(),
        legend.text.align=0,
        legend.key.size=unit(1.5, 'lines'))

evo_fit
save_plot("../../figures/fit_evo.pdf", evo_fit, base_height=5, base_width=5)


# calculate relative fitness compared to WT
fit %>% filter(!Strain=='910L2-init') %>% group_by(Strain) %>%
  summarize(Background=Background[1], mean_fit=mean(Fitness)) %>%
  mutate(rel_fit=mean_fit/max(mean_fit))


