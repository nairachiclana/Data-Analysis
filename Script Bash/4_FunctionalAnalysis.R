options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
type.expression = args[1]
lfc = as.numeric(args[2])
pval = as.numeric(args[3])


library(tibble)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyr)
library(rlang)
library(clusterProfiler)
library(curl)
library(dplyr)
library(ggplot2)


load(file="savedData/filtered_df.RData")
load(file="savedData/DEA.RData")

#OUR GENES
genes<-names(filtered.df)[11:ncol(filtered.df)]
df.exp.genes<-diff.exp.df[(diff.exp.df$gene.name %in% genes) , ]
names(df.exp.genes) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "hgnc_symbol")

#ENRICHMENT ANALYSIS

orgdb = "org.Hs.eg.db"
biomart_dataset = "hsapiens_gene_ensembl"
keggname = "hsa"
mart=useMart(biomart = "ensembl", dataset = biomart_dataset)

entrezsymbol=getBM(attributes = c("entrezgene", "hgnc_symbol"), mart=mart)
entrezsymbol$entrezgene = as.character(entrezsymbol$entrezgene)

summarize_cp = function(df, comparison) {
  summaries = data.frame()
  for (ont in names(df)) {
    ontsum = summary(df[[ont]])
    ontsum$ont = ont
    summaries = rbind(summaries, ontsum)
  }
  summaries$comparison = comparison
  return(summaries)
}

enrich_cp = function(df, comparison, type, lfc, pval) {
  df = df %>% data.frame()  %>% left_join(entrezsymbol, by = "hgnc_symbol") %>% filter(!is.na(entrezgene))
  if(type=="all")
    df<-df %>% filter(abs(logFC) > lfc & adj.P.Val < pval) # lfc and pval threshold defined above in the volcano plot
  if(type=="over")
    df<-df %>% filter(logFC > lfc  & adj.P.Val < pval)
  if(type=="under")
    df<-df %>% filter(logFC < -lfc & adj.P.Val < pval)
  genes = df$entrezgene
  mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  cc = enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  bp = enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
  all = list(mf = mf, cc = cc, bp = bp, kg = kg)
  all[["summary"]] = summarize_cp(all, comparison)
  return(all)
}


enrich_rs<-enrich_cp(df.exp.genes, "TNBC/LumA", type=type.expression, lfc, pval)
ggsave(dotplot(enrich_rs$kg, x="count", showCategory=10))

