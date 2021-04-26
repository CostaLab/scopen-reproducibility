library(ggplot2)
library(cicero)

df <- read.table("GM12878_imputed.txt", header = TRUE)
input_cds <- make_atac_cds(df, binarize = FALSE)

set.seed(2018)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none",
                             check_duplicates = FALSE)

tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

data("human.hg19.genome")
conns <- run_cicero(cicero_cds, human.hg19.genome) # Takes a few minutes to run

#conns_sub <- subset(conns, coaccess > 0)
write.table(conns, file = "GM12878_Cicero.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)

conns <- read.table("GM12878_Cicero.txt", header = TRUE)
gene_annotation <- read.table("hg19.annotation.bed", header = TRUE, sep = "\t")
data(gene_annotation_sample)

gene_annotation <- subset(gene_annotation, feature %in% unique(gene_annotation_sample$feature))

gene_annotation$strand <- as.character(gene_annotation$strand)
gene_annotation$feature <- as.character(gene_annotation$feature)
gene_annotation$gene <- as.character(gene_annotation$gene)
gene_annotation$transcript <- as.character(gene_annotation$transcript)
gene_annotation$symbol <- as.character(gene_annotation$symbol)

pdf("GM12878.pdf", height = 4, width = 4)
plot_connections(conns, "chr19", 42340000, 42562000,
                 #gene_model = gene_annotation,
                 coaccess_cutoff = 0.0, 
                 connection_width = 1.0, 
                 collapseTranscripts = "longest",
                 #viewpoint = "chr19_42379190_42380000",
                 include_axis_track = FALSE)
dev.off()
