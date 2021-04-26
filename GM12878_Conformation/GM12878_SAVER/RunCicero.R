library(ggplot2)
library(cicero)

df <- read.table("GM12878_imputed.txt", header = TRUE)
input_cds <- make_atac_cds(df, binarize = FALSE)

set.seed(2018)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "PCA")
input_cds <- reduce_dimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none", 
                             preprocess_method = "PCA")

tsne_coords <- reducedDims(input_cds)$tSNE
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
