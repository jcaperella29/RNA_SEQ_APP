# get_annotation.R
args <- commandArgs(trailingOnly = TRUE)
ids <- unlist(strsplit(args[1], ","))

library(biomaRt)

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = ids,
  mart = mart
)

write.csv(annot, file = "annotation_output.csv", row.names = FALSE)
