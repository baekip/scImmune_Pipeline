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
basic.path = paste(c(output.path,"/Trajectory/"), collapse="")

###load Cluster.Rda
sample.cluster = readRDS (file = paste(c(input.path,"/Rdata/",sample.id,".Cluster.Rda"), collapse=""))

###format monocle data
sample.import <- importCDS (sample.cluster)
my_feat <- fData(sample.import)
names(my_feat) <- c('gene_short_name')
sample.cds <- newCellDataSet (exprs(sample.import),
                              phenoData = new("AnnotatedDataFrame", data = pData(sample.import)),
                              featureData = new("AnnotatedDataFrame", data = my_feat),
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
sample.cds <- estimateSizeFactors(sample.cds)
sample.cds <- estimateDispersions(sample.cds)

###variance estimation steps
sample.cds.DG <- detectGenes (sample.cds, min_expr = 0.1)
disp_table <- dispersionTable (sample.cds)
table (disp_table$mean_expression >= 0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
unsup_clustering_genes

###Marks genes for clustering
sample.cds.DG.OF <- setOrderingFilter (sample.cds.DG, unsup_clustering_genes$gene_id)
CairoPNG (file = paste(c(basic.path, "plot_ordering_gene.png"), collapse=""),
    width=800, height=600)
plot_ordering_genes(sample.cds.DG.OF)
dev.off()

CairoPNG (file = paste(c(basic.path, "plot_pc_variance_explained.png"), collapse=""),
    width=800, height=600)
plot_pc_variance_explained(sample.cds.DG.OF, return_all = FALSE)
dev.off()

###Constructing Single Cell Trajectories
expressed_genes <- row.names (subset(fData(sample.cds.DG.OF),
                              num_cells_expressed >=10))
expressed_genes
sample.subset <- sample.cds [expressed_genes]
sample.subset <- detectGenes (sample.subset, min_expr=0.1)

fData(sample.subset)$use_for_ordering <- fData(sample.subset)$num_cells_expressed > 0.05*ncol(sample.subset)
table(fData(sample.subset)$use_for_ordering)

###trajectory analysis (for a second time, abt 10min)
clustering_DEG_genes <- differentialGeneTest(sample.subset,
                                             fullModelFormulaStr = '~res.0.6')
my_ordering_genes <- row.names (clustering_DEG_genes)[order(clustering_DEG_genes$qval)[1:1000]]
sample.subset.OF <- setOrderingFilter (sample.subset, ordering_genes = my_ordering_genes)
sample.subset.RD <- reduceDimension (sample.subset.OF, method = 'DDRTree')

###the warnings were for use of deprecated code
sample.subset.OC <- orderCells (sample.subset.RD)

#GM_state <- function(cds){
#    if (length(unique(pData(cds)$State)) > 1){
#        T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
#        return(as.numeric(names(T0_counts)[which
#              (T0_counts == max(T0_counts))]))
#    } else {
#      return (1)
#    }
#}

###for identifying the State which contains most of the cells from time zero
#sample.subset.OC <- orderCells(sample.subset.OC,
#                    root_state = GM_state(sample.subset.OC))

CairoPNG (file = paste(c(basic.path, "plot_cell_trajectory_orig.png"), collapse=""),
plot_cell_trajectory(sample.subset.OC, color_by = "res.0.6")
dev.off()

CairoPNG (file = paste(c(basic.path, "plot_cell_trajectory_pseudotime.png"), collapse=""),
plot_cell_trajectory(sample.subset.OC, color_by = "Pseudotime")
dev.off()

###Save RDS file
saveRDS (clustering_DEG_genes, file = paste(c(output.path,"/Rdata/",sample.id,".Cluster.Rda"), collapse=""))
saveRDS (sample.subset.OC, file = paste(c(output.path,"/Rdata/",sample.id,".subset.OC.Rda"), collapse=""))


