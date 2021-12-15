library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)

df <- read.table("../Diff/All_statistics.txt", header = TRUE)

df$CD_PC <- df$Protection_Score_CD_PC + df$TC_CD_PC
df$CNT <- df$Protection_Score_CNT + df$TC_CNT
df$DCT <- df$Protection_Score_DCT + df$TC_DCT
df$DL_TAL <- df$Protection_Score_DL_TAL + df$TC_DL_TAL
df$EC <- df$Protection_Score_EC + df$TC_EC
df$Fib <- df$Protection_Score_Fib + df$TC_Fib
df$IC <- df$Protection_Score_IC + df$TC_IC
df$Injured_PT <- df$Protection_Score_Injured_PT + df$TC_Injured_PT
df$Lymphoid <- df$Protection_Score_Lymphoid + df$TC_Lymphoid
df$Mac <- df$Protection_Score_Mac + df$TC_Mac
df$Pod <- df$Protection_Score_Pod + df$TC_Pod
df$PT_S1 <- df$Protection_Score_PT_S1 + df$TC_PT_S1
df$PT_S2 <- df$Protection_Score_PT_S2 + df$TC_PT_S2
df$PT_S3 <- df$Protection_Score_PT_S3 + df$TC_PT_S3
df$TAL <- df$Protection_Score_TAL + df$TC_TAL


sel_tfs = grep("var", df$Motif,invert=TRUE)

df=df[sel_tfs, ]

# df <- subset(df, select = c("Motif", "Num",
#                            "CD_PC", "CNT", "DCT",
#                            "DL_TAL", "EC", "Fib",
#                            "IC", "Injured_PT", "Lymphoid",
#                            "Mac", "Pod", "PT_S1",
#                            "PT_S2", "PT_S3", "TAL"))

df <- subset(df, select = c("Motif", "Num",
                            "Fib", "Injured_PT", "Lymphoid",
                            "Mac", "PT_S1",
                            "PT_S2", "PT_S3"))


df$Mean <- apply(df[3:ncol(df)], 1, mean)
df$Var <- apply(df[3:ncol(df)], 1, sd)


write.table(df, file = "TF.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# 


rownames(df) <- df$Motif
df$Motif <- NULL

df_plot <- subset(df, Num > 1000 & Var > 0.05)
df_plot$Mean <- NULL
df_plot$Var <- NULL
df_plot$Num <- NULL

df_plot_scale <- t(apply(df_plot, 1, scale))
colnames(df_plot_scale) <- colnames(df_plot)
rownames(df_plot_scale) <- rownames(df_plot)

library(dendsort)
col_dend = dendsort(hclust(dist(t(df_plot_scale))))

p <- Heatmap(as.matrix(df_plot_scale),
             name = "TF",
             cluster_columns = col_dend,
             cluster_rows = TRUE,
             show_row_names = FALSE,
             column_names_gp = gpar(fontsize = 10),
             show_column_dend = TRUE,
             show_row_dend = FALSE,
             clustering_method_rows = "ward.D2",
             col = rev(brewer.pal(n = 11, name = "RdBu")))

pdf("heatmap_without_name.pdf", width = 3, height = 6)
draw(p)
dev.off()


p <- Heatmap(as.matrix(df_plot_scale),
             name = "TF",
             cluster_columns = col_dend,
             cluster_rows = TRUE,
             show_row_names = TRUE,
             row_names_gp = gpar(fontsize = 4),
             column_names_gp = gpar(fontsize = 10),
             show_column_dend = TRUE,
             show_row_dend = FALSE,
             clustering_method_rows = "ward.D2",
             col = rev(brewer.pal(n = 11, name = "RdBu")))


pdf("heatmap_with_name.pdf", width = 4, height = 30)
draw(p)
dev.off()
