# Statistical analysis of hisat2-aligned tpm data
library(tidyverse)
library(cowplot)
library(broom)

rna <- read.csv("../../data/results/counts_rna_abundance.csv")

rna$gene <- factor(rna$gene, levels=unique(rna$gene), ordered=TRUE)


# pairwise t-test comparing wt to all other strains (correction for multiple testing)
pvals <- function(a, b){
  rna %>% filter(strain %in% c(a, b)) %>%
    spread(strain, tpm) %>% group_by(gene) %>%
    nest() %>%
    mutate(t_test = map(data, ~ tidy(t.test(.[[a]], .[[b]], data = ., paired = F)))) %>%
    select(-data) %>%
    unnest() %>%
    select(gene, p.value, estimate) %>%
    mutate(control=a, strain=b, padj = p.adjust(p.value, method = "fdr")) %>%
    #filter(padj < 0.05) %>%
    arrange(padj) %>%
    mutate_at(vars(p.value, estimate, padj), signif, digits = 2) %>%
    select(control, strain, gene, estimate, p.value, padj)
}
get_all_pvals <- function(refer, strain.list){
  strain.list %>% map(pvals, a=refer) %>%
    bind_rows()
}
strains <- as.character(unique(subset(rna, !strain=="T7Hi")$strain))


# consolidate all significantly diff expressed genes for pairwise comparisons with wild-type
vs_wt_pvals <- get_all_pvals(ref="T7Hi", strain.list=strains)
# t-test comparing knockouts in 10deop background to the 10deop ancestor
deop_pvals <- get_all_pvals(ref="11-44", strain.list=c("11-44-phi9v2", "11-44-phi10"))
# compare evolved to initial strains for wt-910 and 8st-910 lines
evo_pvals <- bind_rows(get_all_pvals(ref="phi910v2", strain.list="910L2evo"), 
                       get_all_pvals(ref="8st-910", strain.list="8st-910evo"))

# combine all comparisons into a single df
all_pvals <- bind_rows(vs_wt_pvals, deop_pvals, evo_pvals)

# add in information on treatment, background, promoter knockout type to data frame. 
labels <-data.frame(
  strain=c('wt', 'phi9v2', 'phi910v2','phi10', '11-44', '11-44-phi9v2', '11-44-phi10', '8st-9', '8st-910', '8st-910evo', '910L2evo'),
  treatment=c('wt', 'wt-9', 'wt-910', 'wt-10', 'deop', 'deop-9', 'deop-10', '8st-9', '8st-910', '8st-910evo', '910L2evo'),
  background=c('wt', 'wt', 'wt', 'wt', '10deop', '10deop', '10deop', '8st', '8st', '8st', 'wt'),
  knockout=c('wt', 'phi9', 'phi910', 'phi10', 'wt', 'phi9', 'phi10', 'phi9', 'phi910', 'phi910', 'phi910')
)
all_pvals <- inner_join(labels, all_pvals) 

write.csv(all_pvals, "../../data/results/pairwise_t_vs_wt.csv", row.names = F)
write.csv(evo_pvals, "../../data/results/pairwise_t_evo.csv", row.names=F)

all_pvals %>% filter(gene %in% c("8", "9", "10A", "11", "12", "13"))
all_pvals %>% filter(strain=="11-44")

