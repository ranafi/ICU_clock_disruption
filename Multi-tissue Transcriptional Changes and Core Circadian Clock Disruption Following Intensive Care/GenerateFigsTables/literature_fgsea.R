#This file generates results from table4, where we see if our differentially
#expressed IC genes, ranked by significance, are enriched in existing gene sets
#associated with sleep loss and fasting (2 common biproducts of IC). The final
#table stored in variable 'out.' Notice in final paper table we do not consider
#skeletal muscle tissue, this is because other publications suggest the post-mortem
#interval seriously changes gene expression in this tissue (and some GTEx subjects
#have PMI > a few hours).

library(fgsea)
library(data.table)
library(ggplot2)
library(stringr)
library(fgsea)  #bioconductor install this if not already done
library(rstudioapi)
path = dirname(dirname(rstudioapi::getActiveDocumentContext()$path)) #path of parent folder
setwd(paste0(path, "/GenerateFigsTables/compare_to_existing_literature/"))

pathways <- gmtPathways("literature_gene_set.gmt")

rnk_files = c("AdiposeSubcutaneousFDR_preranked.rnk")

run_fgsea = function(file, gsea_param = 1){
  ranks <- read.table(file, header=T, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$rnk, ranks$...1)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 1,
                    eps = 0.0,
                    gseaParam = gsea_param,
                    maxSize  = 1000)
}

fgseaRes = lapply(rnk_files, run_fgsea, gsea_param=1)
names(fgseaRes) = c("Adipose-Subcutaneous_FDR_ranked")

pvals = sapply(fgseaRes, `[[`, 2)
colnames(pvals) = c("Adipose-Subcutaneous_FDR_ranked_pvals")

out = round(pvals, 3)
rownames(out) = c("up_and_down_genes_assoc_w_fasting", "up_genes_assoc_w_fasting", "down_genes_assoc_w_fasting", "up_and_down_genes_assoc_w_sleep_loss_in_musc_skel", "up_and_down_genes_assoc_w_sleep_loss_in_adipose_sub", "up_genes_assoc_w_sleep_loss_in_adipose_sub", "down_genes_assoc_w_sleep_loss_in_adipose_sub", "up_genes_assoc_w_sleep_loss_in_musc_skel","down_genes_assoc_w_sleep_loss_in_musc_skel" )

#####if you want to plot a particular fgsea plot do so with the following code #####
plot_pathway = function(file, pathway, gsea_param = 1){
  name = names(pathway)
  ranks <- read.table(file, header=T, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$rnk, ranks$...1)
  plotEnrichment(unlist(pathway), gseaParam = gsea_param, ranks) + labs(title=paste(name, str_remove(file, '.rnk')))
  
}
plot_pathway(rnk_files[1], pathways[2])
