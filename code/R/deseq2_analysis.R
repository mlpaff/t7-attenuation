# Run differential expression analysis on the raw counts data generated from the hisat2 alignments using DEseq2
library(tidyverse)
library(DESeq2)
library(cowplot)
library(stringr)

counts <- read.csv("../../data/results/hisat2_rna_abundance.csv", header=TRUE)
t7_genes <- read.csv("../../data/t7_genes.csv", header=TRUE)

# Create the counts matrix with condition, gene, counts
# replicate 2 has issues with incorrect labelling so I am excluding it from the analysis until the samples are re-sequenced
counts %>% mutate(condition=paste0(strain, '-', Sample)) %>% select(condition, gene, counts) %>%
  spread(condition, counts) -> mtx_counts

# Create matrix using the counts data for input into DESeq2
mtx_counts %>% remove_rownames() %>% column_to_rownames(var="gene") -> mtx_counts

# Set conditions table with condition and corresponding sample
counts %>% mutate(condition=paste0(strain, '-', Sample)) %>% 
  select(condition, gene, treatment) %>% spread(gene, treatment) -> meta

coldata <- data.frame(
  row.names = colnames(mtx_counts),
  condition = meta$`0.3`
)

# Check that columns of the count matrix and the rows of the column data are in the same order
#all(rownames(coldata) %in% colnames(mtx_counts))
#all(rownames(coldata) == colnames(mtx_counts))

# Construct the DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = mtx_counts,
                              colData = coldata,
                              design = ~ condition)

# Pre-filter low-count reads
dds <- dds[ rowSums(counts(dds)) > 1, ]
# Set the 'wt' strain as the reference level condition -- the point that you want to make the comparisons against
#dds$condition <- relevel(dds$condition, ref="wt")
# remove levels that do not have samples when working with subsets of the data
#dds$condition <- droplevels(dds$condition)

# Run differential expression analysis
dds <- DESeq(dds)

# Function extracts results for single pair-wise comparison between "condition" and wt T7 -- returns data frame of results
get_results <- function(con){
  # given a condition, return results comparing to wt
  res <- results(dds, contrast=c("condition", con, "wt"), alpha=0.05)
  # Re-order the data frame in order of ascending adjusted p-value
  #resOrdered <- res[order(res$padj),]
  #df <- data.frame(resOrdered, gene=row.names(resOrdered), treatment=con)
  df <- data.frame(res, gene=row.names(res), treatment=con) %>% arrange(padj)
}

# Consolidate all results into a single data frame of differential expression for all conditions compared to wt
get_all_results <- function(list_strains){
  # Create a data frame with results from all pairwise comparisons of diff expression from each promoter knockout strain
  list_strains %>% map(get_results) %>%
    bind_rows()
}

# Extract the conditions (T7 strains) that will be compared to wt T7
treats <- str_extract(paste(unique(coldata$condition)), '.*-.*|deop')
treats <- treats[!is.na(treats)]

results_data_frame <- get_all_results(treats)

# add in information on strain, treatment, background, promoter knockout type to data frame. 
labels <-data.frame(
  strain=c('wt', 'phi9v2', 'phi910v2','phi10', '11-44', '11-44-phi9v2', '11-44-phi10', '8st-9', '8st-910', '8st-910evo', '910L2evo'),
  treatment=c('wt', 'wt-9', 'wt-910', 'wt-10', 'deop', 'deop-9', 'deop-10', '8st-9', '8st-910', '8st-910evo', '910L2evo'),
  background=c('wt', 'wt', 'wt', 'wt', '10deop', '10deop', '10deop', '8st', '8st', '8st', 'wt'),
  knockout=c('wt', 'phi9', 'phi910', 'phi10', 'wt', 'phi9', 'phi10', 'phi9', 'phi910', 'phi910', 'phi910')
)

results_data_frame <- inner_join(labels, results_data_frame) %>% filter(padj<=0.05)


# Set stars to indicate significance level from adjusted p-values
results_data_frame$star <- ""
results_data_frame$star[results_data_frame$padj <= .05]  <- "*"
results_data_frame$star[results_data_frame$padj <= .01]  <- "**"
results_data_frame$star[results_data_frame$padj <= .001] <- "***"

# Output differential expression comparing different knockout strains to wt (T7Hi)
write.csv(results_data_frame, "../../data/results/deseq2_results_vs_wt.csv", row.names=F)

results_data_frame %>% filter(knockout=='phi9')

# Extract results for comparisons of differential expression between evolved and initial strains
# Results for wt-910 vs wt-910evo
res_910 <- results(dds, contrast=c("condition", "910L2evo", "wt-910"), alpha=0.05)
df_910 <- data.frame(res_910, gene=row.names(res_910), treatment="910L2evo")
# Results for 8st-910 vs 8st-910evo
res_8st <- results(dds, contrast=c("condition", "8st-910evo", "8st-910"), alpha=0.05)
df_8st <- data.frame(res_8st, gene=row.names(res_8st), treatment="8st-910evo")

# Combine data frames for evolution results and add labels
df_evo <- inner_join(labels, bind_rows(df_910, df_8st))

df_evo$star <- ""
df_evo$star[df_evo$padj <= .05]  <- "*"
df_evo$star[df_evo$padj <= .01]  <- "**"
df_evo$star[df_evo$padj <= .001] <- "***"

# Save comparison of differential expression between initial and evolved strains 
write.csv(df_evo, "../../data/results/deseq2_results_evolved.csv")

# Return normalized counts data for all strains 
# Function for extracting plotCounts data from dds
get_gene_abundance <- function(geneName){
  gene_ab <- plotCounts(dds, gene=geneName, intgroup="condition", returnData=TRUE) %>% 
    mutate(gene1=paste("gene", geneName), gene=geneName) %>%
    left_join(labels, by=c("condition"="treatment")) 
    #left_join(results_data_frame %>% filter(gene==geneName))
}

get_all_abundance <- function(genes_list){
  genes_list %>% map(get_gene_abundance) %>%
    bind_rows()
}

t7_abundance <- get_all_abundance(as.character(t7_genes$genes))
t7_abundance$gene1 <- factor(t7_abundance$gene1, levels=unique(t7_abundance$gene1), ordered=TRUE)

write.csv(t7_abundance, "../../data/results/tot_t7_abundance.csv")
