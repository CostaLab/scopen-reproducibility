library(ggplot2)
library(cicero)
library(optparse)

option_list = list(
    make_option("--input", type="character", default=NULL, help="input feature file"),
    make_option("--output", type="character", default=NULL, help="output file name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

df <- read.table(opt$input, header = TRUE)
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

write.table(conns, file = opt$output, 
            sep = "\t", quote = FALSE, row.names = FALSE)
