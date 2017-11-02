# Statistical analysis of hisat2-aligned tpm data
library(tidyverse)
library(cowplot)
library(broom)

rna <- read.csv("../../data/results/hisat2_rna_abundance.csv")

rna$gene <- factor(rna$gene, levels=unique(rna$gene), ordered=TRUE)

# Normalize the raw RNA counts, first to the sum of raw counts for all genes up to and including gene 7.7, then by TPM
rna %>% group_by(rep, strain) %>% 
  mutate(nfactor=sum(counts[gene<='7.7'])) %>%
  mutate(normcounts=counts/nfactor) %>%
  mutate(rpk=normcounts/((stop-start)/1000)) %>%
  mutate(rpm=sum(rpk)) %>%
  mutate(tpm=rpk/rpm) -> counts

counts$gene <- factor(counts$gene, levels=rna$gene[1:60])

# anova analysis, assuming intactions between gene, background, and knockout
counts %>% select(strain, background, knockout, gene, rep, tpm) %>% mutate(gene1=paste('gene', gene)) %>%
  filter(!strain %in% c("910L2evo", "8st-910evo")) -> df

fit <- lm(tpm ~ gene1*background*knockout, data=df)

glance(fit)
anova(fit)
summary(fit)


counts %>% filter(strain %in% c("T7Hi", "phi9v2")) %>%
  spread(strain, tpm) %>% group_by(gene) %>%
  nest() %>%
  #mutate(t_test=map(.x ~tidy(t.test(.$phi910v2, .$`910L2evo`, data=.x, paired = F)))) %>%
  mutate(t_test = map(data, ~ tidy(t.test(.$T7Hi, .$phi9v2, data = ., paired = F)))) %>%
  select(-data) %>%
  unnest() %>%
  select(gene, p.value, estimate) %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) %>%
  filter(p.value_adj < 0.05) %>%
  arrange(p.value_adj) %>%
  mutate_at(vars(p.value, estimate, p.value_adj), signif, digits = 2) %>%
  select(gene, estimate, p.value, p.value_adj)

# pairwise t-test comparing wt to all other strains (correction for multiple testing)
pvals <- function(a, b){
  counts %>% filter(strain %in% c(a, b)) %>%
    spread(strain, tpm) %>% group_by(gene) %>%
    nest() %>%
    #mutate(t_test=map(.x ~tidy(t.test(.$phi910v2, .$`910L2evo`, data=.x, paired = F)))) %>%
    mutate(t_test = map(data, ~ tidy(t.test(.[[a]], .[[b]], data = ., paired = F)))) %>%
    select(-data) %>%
    unnest() %>%
    select(gene, p.value, estimate) %>%
    mutate(control=a, strain=b, padj = p.adjust(p.value, method = "fdr")) %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    mutate_at(vars(p.value, estimate, padj), signif, digits = 2) %>%
    select(control, strain, gene, estimate, p.value, padj)
}
get_all_pvals <- function(refer, strain.list){
  strain.list %>% map(pvals, a=refer) %>%
    bind_rows()
}
strains <- as.character(unique(subset(counts, !strain=="T7Hi")$strain))


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

all.pvals %>% filter(gene %in% c("8", "9", "10A", "11", "12", "13"))
all.pvals %>% filter(strain=="11-44")




t.test(tpm~strain, subset(df, strain %in% c("T7Hi", "11-44") & gene=='10A'))

get_p_vals <- function(treat, cont, gene.name, dfa){
  stat <- t.test(tpm~strain, subset(dfa, strain %in% c(cont, treat) & gene==gene.name))
  sigs <- data.frame(comparison=paste0(cont, '_', treat),
                     strain=treat,
                     gene=gene.name,
                     mean.diff=as.numeric(stat$estimate[1]-stat$estimate[2]),
                     p.value=stat$p.value)
  return(sigs)
}
# Feed in list of strains to compare to control, returns p-values for a specific gene
get_pval_for_gene <- function(strains, cont, gene.name, dfa){
  strains %>% map(get_p_vals, cont=cont, gene.name=gene.name, dfa=dfa) %>%
    bind_rows()
}

get_all_pvals <- function(genes, cont, strains, dfa){
  genes %>% map(get_pval_for_gene, cont=cont, strains=strains, dfa=dfa) %>%
    bind_rows()
}

# Function that takes a list of genes and performs a t-test for each gene between 2 strains
gene.list <- as.character(unique(df$gene))
strains.list <- as.character(unique(subset(df, !strain=="T7Hi")$strain))

all.pvals <- get_all_pvals(genes=gene.list, cont="T7Hi", strains=strains.list, dfa=df)

all.pvals %>% filter(strain=="phi9v2") %>% mutate(padj=p.adjust(p.value, method='fdr')) %>%
  filter(padj<0.05) %>% arrange(padj)


# t test comparing initial and final population for evolved lines
counts %>% select(strain, background, knockout, gene, rep, tpm) %>% mutate(gene1=paste('gene', gene)) %>%
  filter(rep != '2', strain %in% c("wt", "phi910v2", "8st-910",  "910L2evo", "8st-910evo")) -> evo_counts
wt_910_vs_evo <- get_all_pvals(gene.list, control="phi910v2", strains="910L2evo", df=evo_counts)
st_910_vs_evo <- get_all_pvals(gene.list, control="8st-910", strains="8st-910evo", df=evo_counts)

wt_910_vs_evo <- wt_910_vs_evo %>% mutate(p.adj=p.adjust(p.value, method='fdr')) %>%
  filter(p.adj < 0.1) %>%
  arrange(p.adj) %>%
  mutate_at(vars(p.value, mean.diff, p.adj), signif, digits = 2) %>%
  select(gene, mean.diff, p.value, p.adj)

 
