library(ggplot2)
library(cicero)

df <- read.table("GSM2970932_sciATAC_GM12878_counts.txt", header = TRUE)
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

write.table(conns, file = "./GM12878.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)

gene_annotation <- read.table("gencode.v19.txt", header = TRUE)
gene_annotation$strand <- as.character(gene_annotation$strand)
gene_annotation$feature <- as.character(gene_annotation$feature)
gene_annotation$gene <- as.character(gene_annotation$gene)
gene_annotation$transcript <- as.character(gene_annotation$transcript)
gene_annotation$symbol <- as.character(gene_annotation$symbol)

pdf("GM12878.pdf", height = 4, width = 6)
plot_connections(conns, "chr19", 42348190, 42571957,
                 gene_model = gene_annotation,
                 coaccess_cutoff = 0.0, 
                 connection_width = .5, 
                 collapseTranscripts = "gene",
                 viewpoint = "chr19_42379190_42380000",
                 include_axis_track = FALSE)
dev.off()