library(optparse)
library(cisTopic)
library(methods)
library(stringr)
suppressMessages(library(patchwork))
suppressMessages(library(glue))
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(cowplot))
suppressMessages(library(genomation))
suppressMessages(library(dplyr))

option_list = list(
  make_option("--input", type="character", default=NULL, help="input feature file"),
  make_option("--output_dir", type="character", default=NULL, help="output file name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input

x <- read.table(input_file, header = TRUE, check.names = FALSE)
region <- as.data.frame(t(as.data.frame(str_split(rownames(x), "_"))))
colnames(region) <- c("chrom", "p1", "p2")
rownames(x) <- paste0(region$chrom, ":", region$p1, "-", region$p2)
cisTopicObject <- createcisTopicObject(count.matrix = as.matrix(x))

cisTopicObject <- runWarpLDAModels(cisTopicObject, 
                                   topic = c(2:30), 
                                   seed = 987, 
                                   nCores = 4, 
                                   addModels=FALSE,
                                   iterations = 500)

cisTopicObject <- selectModel(cisTopicObject)

t1 <- as.data.frame(cisTopicObject@selected.model$document_expects)
t2 <- as.data.frame(t(cisTopicObject@selected.model$topics))

# normalization
t1 <- t1 / rowSums(t1)
t2 <- t2 / colSums(t2)

colnames(t1) <- colnames(x)

x_complete <- as.matrix(t2) %*% as.matrix(t1)

colnames(x_complete) <- colnames(x)
rownames(x_complete) <- rownames(x)

write.table(x_complete, file = sprintf("%s/cisTopic.txt", opt$output_dir), 
            quote = FALSE, sep = "\t")
write.table(t1, file = sprintf("%s/documents.txt", opt$output_dir), 
            quote = FALSE, sep = "\t")
