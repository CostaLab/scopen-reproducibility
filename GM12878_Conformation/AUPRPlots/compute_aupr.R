library(GenomicRanges)
library(cicero)
library(PRROC)
library(ggplot2)
library(cowplot)
library(gridExtra)

df_hic <- read.table("../ConnPlots/GM12878_HiC.txt", header = TRUE)
df_chia_pet <- read.table("../ConnPlots/GM12878_CHIAPET.txt", header = TRUE)

get_pr <- function(df_true, df_pred){
    df_pred$coaccess[is.na(df_pred$coaccess)] <- 0
    df_pred$isLinked <- compare_connections(df_pred, df_true, 
                                            maxgap = 1000)
    
    df_class0 <- subset(df_pred, isLinked == TRUE)
    df_class1 <- subset(df_pred, isLinked == FALSE)
    
    pr <- pr.curve(scores.class0 = df_class0$coaccess,
                   scores.class1 = df_class1$coaccess,
                   curve = TRUE)
    
    pr$backgroup_precision <- nrow(df_class0) / nrow(df_pred)
    
    return(pr)
}

for (method in c("Raw", "scOpen", "MAGIC", "cisTopic", "DCA", "PCA", "SAVER", "SCALE", "scBFA", "scImpute")) {
    if(method == "Raw"){
        df <- read.table("../GM12878_Cicero/GM12878.txt", header = TRUE)    
    } else{
        df <- read.table(glue::glue("../GM12878_{method}/GM12878_Cicero.txt"), header = TRUE)
    }
    
    pr_hic <- get_pr(df_true = df_hic, df_pred = df)
    pr_chiapet <- get_pr(df_true = df_chia_pet, df_pred = df)
    
    saveRDS(pr_hic, glue::glue("{method}_hic.Rds"))
    saveRDS(pr_chiapet, glue::glue("{method}_chiapet.Rds"))
    
    pdf(glue::glue("{method}_AUPR.pdf"))
    plot(pr_hic)
    plot(pr_chiapet)
    dev.off()
}
