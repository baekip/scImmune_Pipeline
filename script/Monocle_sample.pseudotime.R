#/TBI/Tools/R/current/bin/R
#load(file = "~/Projects/datasets/pbmc3k_final.Rda")
###Preparation of Seurat Process###
suppressMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
    library(Cairo)
    library(monocle)
})

args = commandArgs(TRUE)
options (bitmapType='cairo')

###analysis argument
sample.id = args[1]
input.path = args[2]
output.path = args[3]



###preparation process
basic.path = paste(c(output.path,"/Pseudotime/"), collapse="")

###load OC.Rda (Ordering Cells)
sample.subset.OC = readRDS (file = paste(c(input.path,"/Rdata/",sample.id,".subset.OC.Rda"), collapse=""))

###Finding Genes that 
my_pseudotime_de <- differentialGeneTest(sample.subset.OC,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)")
my_pseudotime_de %>% arrange(qval) %>% head()

###save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$gene_short_name


my_pseudotime_gene <- my_pseudotime_gene$gene_short_name

###draw significant gene plot 
###!!Error
my_gene <- row.names (subset(fData(sample.subset.OC),
                      my_pseudotime_gene))
CairoPNG (file = paste(c(basic.path, "plot_genes_in_pseudotime.png"), collapse = ""), width = 800, height = 600)
plot_genes_in_pseudotime (sample.subset.OC[my_genes], color_by = "res.0.6")
dev.off()

###Analyzing Branches in Single-Cell Trajectories
BEAM_res <- BEAM (sample.subset.OC, branch_point = 1)
BEAM_res_order <- BEAM_res [order(BEAM_res$qval),]
BEAM_res_order <- BEAM_res [,c("gene_short_name", "pval", "qval")]

###for a second (abt 10min)
CairoPNG (file = paste (c(basic.path, "plot_genes_branched_heatmap.png"), collapse = ""), width = 800, height = 600)
my_branched_heatmap <- plot_genes_branched_heatmap (sample.subset.OC[row.names(subset(BEAM_res_order, qval < 1e-15)),], branch_point = 1, num_clusters = 4, use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()

###Save RDS file
saveRDS (my_pseudotime_de, file = paste(c(output.path,"/Rdata/",sample.id,".my_pseudotime_de.Rda"), collapse=""))
saveRDS (BEAM_res, file = paste(c(output.path,"/Rdata/",sample.id,".BEAM_res.Rda"), collapse=""))
saveRDS (my_branched_heatmap, file = paste(c(output.path,"/Rdata/",sample.id,".BEAM_res.Rda"), collapse=""))


