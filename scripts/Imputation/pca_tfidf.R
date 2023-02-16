library(optparse)
library(methods)
library(stringr)
library(missMDA)
library(Signac)
library(Seurat)

option_list = list(
  make_option("--input", type="character", default=NULL, help="input feature file"),
  make_option("--output", type="character", default=NULL, help="output file name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input
x <- read.table(input_file, header = TRUE)

## run TF-IDF transformation
chrom_assay <- CreateChromatinAssay(
  counts = as.matrix(x),
  sep = c("_", "_"),
  genome = 'hg19'
)

obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

obj <- RunTFIDF(obj)

x <- obj@assays$peaks@data
x[x==0] <- NA

## Imputation
x.comp <- imputePCA(as.data.frame(x), ncp = 30, seed = 2020, scale = FALSE,
                    init = 2)

x_complete <- x.comp$completeObs
colnames(x_complete) <- colnames(x)
rownames(x_complete) <- rownames(x)

write.table(x_complete, file = opt$output, 
            quote = FALSE, sep = "\t")
