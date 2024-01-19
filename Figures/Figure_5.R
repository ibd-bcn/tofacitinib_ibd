options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(readxl)

source('Figures/functions_plots.R')
# Figure 5E -------------------------

## Scatter plots of our DE data and __ paper
setwd('Analysis/')
de_data <- readRDS('data/01_DE/REPASO/new_complete.RDS')

# Obtaining data of upp-regulated genes and down-regulated ones from both Responders (R) and Non-responders (NR) from the M2 population
de_data2 <- de_data[de_data$cluster == 'M2' &
                      de_data$annotation == 'annotation_intermediate' &
                      de_data$comp %in% c('w0R_vs_POSTR', 'w0NR_vs_POSTNR', 'RPOST_NRPOST') &
                      de_data$sign %in% c('UPP', 'UP' ,'DWW','DW'), c("p_val","avg_log2FC",  "sign", "comp", "gene")]

# Responders
df_w0R_vs_POSTR_UPP <- subset(de_data2, comp == "w0R_vs_POSTR" & sign == "UPP" | comp == "w0R_vs_POSTR" & sign == "UP" , select = c("avg_log2FC", "p_val", "gene"))
filtered_df_w0R_vs_POSTR_UPP <- subset(df_w0R_vs_POSTR_UPP, p_val < 0.05)
filtered_df_w0R_vs_POSTR_UPP$condition <- rep("responder_UPP", nrow(filtered_df_w0R_vs_POSTR_UPP))

df_w0R_vs_POSTR_DWW <- subset(de_data2, comp == "w0R_vs_POSTR" & sign == "DWW" | comp == "w0R_vs_POSTR" & sign == "DW", select = c("avg_log2FC", "p_val", "gene"))
filtered_df_w0R_vs_POSTR_DWW <- subset(df_w0R_vs_POSTR_DWW, p_val < 0.05)
filtered_df_w0R_vs_POSTR_DWW$condition <- rep("respoonder_DWW", nrow(filtered_df_w0R_vs_POSTR_DWW))


dataframe_R <- rbind(filtered_df_w0R_vs_POSTR_UPP, filtered_df_w0R_vs_POSTR_DWW) # Final dataframe Responders

# Non-responders
df_w0NR_vs_POSTNR_UPP <- subset(de_data2, comp == "w0NR_vs_POSTNR" & sign == "UPP" | comp == "w0NR_vs_POSTNR" & sign == "UP" , select = c("avg_log2FC", "p_val", "gene"))
filtered_df_w0NR_vs_POSTNR_UPP <- subset(df_w0NR_vs_POSTNR_UPP, p_val< 0.05)
filtered_df_w0NR_vs_POSTNR_UPP$condition <- rep("non_responder_UPP", nrow(filtered_df_w0NR_vs_POSTNR_UPP))

df_w0NR_vs_POSTNR_DWW <- subset(de_data2, comp == "w0NR_vs_POSTNR" & sign == "DWW" | comp == "w0NR_vs_POSTNR" & sign == "DW", select = c("avg_log2FC", "p_val", "gene"))
filtered_df_w0NR_vs_POSTNR_DWW <- subset(df_w0NR_vs_POSTNR_DWW, p_val < 0.05)
filtered_df_w0NR_vs_POSTNR_DWW$condition <- rep("non_responder_DWW", nrow(filtered_df_w0NR_vs_POSTNR_DWW))

dataframe_NR <- rbind(filtered_df_w0NR_vs_POSTNR_UPP, filtered_df_w0NR_vs_POSTNR_DWW) # Final dataframe Non-responders

# Data from other paper
res_DE <- read_csv("~/Azu/rnaseq2/res_DE.csv") # DEGs nota para angela de subirlo al git
res_DE <- na.omit(res_DE)

res_DE_upp <- res_DE[res_DE$padj < 0.05 & res_DE$log2FoldChange >=1,]
res_DE_dww <- res_DE[res_DE$padj < 0.05 & res_DE$log2FoldChange <= -1,]
res_DE_upp$condition <- rep("paper_upp",nrow(res_DE_upp))
res_DE_dww$condition <- rep("paper_dww",nrow(res_DE_dww))
resfinal <- rbind(res_DE_upp, res_DE_dww)



# Tidying the data

names(dataframe_R)[names(dataframe_R) == "avg_log2FC"] <- "log2FoldChange"
resfinal <- resfinal[order(names(resfinal))]
names(resfinal)[names(resfinal) == "...1"] <- "gene"
names(resfinal)[names(resfinal) == "pvalue"] <- "p_val"
dataframe_R <- dataframe_R[order(names(dataframe_R))]

names(dataframe_NR)[names(dataframe_NR) == "avg_log2FC"] <- "log2FoldChange"
dataframe_NR <- dataframe_NR[order(names(dataframe_NR))]



resfinal$lfcSE <- NULL
resfinal$baseMean <- NULL
resfinal$padj <- NULL
resfinal$stat <- NULL

responders <- rbind(resfinal, dataframe_R)

nonresponders <- rbind(resfinal, dataframe_NR)

responders <- na.omit(responders)
responders$p_val <- NULL
responders2 <- pivot_wider(responders, names_from = "condition", values_from = "log2FoldChange")


nonresponders <- na.omit(nonresponders)
nonresponders$p_val <- NULL
nonresponders2 <- pivot_wider(nonresponders, names_from = "condition", values_from = "log2FoldChange")

filtered_p3 <- subset(nonresponders2, !is.na(paper_upp) & !is.na(non_responder_UPP))
filtered_p4 <- subset(nonresponders2, !is.na(paper_dww) & !is.na(non_responder_DWW))

filtered_p3$paper_dww<- NULL
filtered_p3$non_responder_DWW <- NULL
filtered_p4$paper_upp <- NULL
filtered_p4$non_responder_UPP <- NULL
final <- merge(filtered_p3, filtered_p4, by.y = "gene", all.x = TRUE, all.y  = TRUE)

# Plot

plot <- ggplot(data = final) +
  geom_point(aes(x = paper_upp, y = non_responder_UPP), color = "red", size = 1) +
  geom_point(aes(x = paper_dww, y = non_responder_DWW), color = "blue", size = 1) +
  geom_vline(xintercept = 0, colour = "gray") +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_label_repel(aes(label = gene, x = paper_upp, y = non_responder_UPP),
                   box.padding = 0.5, size = 2.5, color = "black",
                   nudge_y = 0.01, segment.size = 0.1,
                   arrow = arrow(length = unit(0.005, "npc"), type = "closed", ends = "last"), max.overlaps = Inf) +
  geom_label_repel(data = final, aes(label = gene, x = paper_dww, y = non_responder_DWW),
                   box.padding = 0.5, size = 2.5, color = "black",
                   nudge_y = 0.01, segment.size = 0.1,
                   arrow = arrow(length = unit(0.005, "npc"), type = "closed", ends = "last"), max.overlaps =  Inf) +

  labs(x = "log2FoldChange", y = "avg_log2FoldChange") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6)
  )


print(plot)
save_sizes(plot =plot , filename = '5E', device = 'jpeg')
save_sizes(plot = plot, filename = '5E', device = 'tiff')
save_sizes(plot = plot, filename = '5E', device = 'svg')
save_sizes(plot = plot, filename = '5E', device = 'pdf')




labs(title = "Responders", x = "log2FoldChange", y = "avg_log2FoldChange") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),       # Increase axis text size
    axis.title = element_text(size = 16),      # Increase axis label size
    plot.title = element_text(size = 18)       # Increase title size
  )
# Figure 5B --------------------------------------------------------------------

#
# Heatmap DMSO
#

DMSO <- read_excel("Figures/extra_data/230916 Base de datos macrofagos Madrid FC 3.xlsx",
                   sheet = "FC Stimuli-Control")

#Heatmap
tiff(filename = "Figures/output/DMSO_heatmap_larga.tiff", width=6,height=10, units = "in",res = 1200)
heatmap_plot("DMSO", DMSO)

jpeg(filename = "Figures/output/DMSO_heatmap_larga.jpeg",width=6,height=10,units = "in",res = 1200)
heatmap_plot("DMSO", DMSO)

svg(filename = "Figures/output/DMSO_heatmap_larga.svg",width=3,height=5)
heatmap_plot("DMSO", DMSO)

pdf(file= "Figures/output/DMSO_heatmap_larga.pdf",width=6,height=10)
heatmap_plot("DMSO", DMSO)


#Legend
tiff(filename = "Figures/output/DMSO_legend_cuadrada.tiff",width=3,height=3,units = "in",res = 1200)
legend <- heatmap_lgd("DMSO")

jpeg(filename = "Figures/output/DMSO_legend_cuadrada.jpeg",width=3,height=3,units = "in",res = 1200)
legend <- heatmap_lgd("DMSO")

svg(filename = "Figures/output/DMSO_legend_cuadrada.svg",width=3,height=3)
legend <- heatmap_lgd("DMSO")

pdf(file = "Figures/output/DMSO_legend_cuadrada.pdf",width=3,height=3)
legend <- heatmap_lgd("DMSO")

# Figure 5D --------------------------------------------------------------------

#
# Heatmap TOFA
#

TOFA <- read_excel("Figures/extra_data/230916 Base de datos macrofagos Madrid FC 3.xlsx",
                   sheet = "FC Tofa stimuli-stimuli")

#Heatmap
tiff(filename = "Figures/output/TOFA_heatmap_larga.tiff", width=6,height=10, units = "in",res = 1200)
heatmap_plot("TOFA", TOFA)

jpeg(filename = "Figures/output/TOFA_heatmap_larga.jpeg",width=6,height=10,units = "in",res = 1200)
heatmap_plot("TOFA", TOFA)

svg(filename = "Figures/output/TOFA_heatmap_larga.svg",width=3,height=5)
heatmap_plot("TOFA", TOFA)

pdf(file= "Figures/output/TOFA_heatmap_larga.pdf",width=6,height=10)
heatmap_plot("TOFA", TOFA)

#Legend
tiff(filename = "Figures/output/TOFA_legend_cuadrada.tiff",width=3,height=3,units = "in",res = 1200)
legend <- heatmap_lgd("TOFA")

jpeg(filename = "Figures/output/TOFA_legend_cuadrada.jpeg",width=3,height=3,units = "in",res = 1200)
legend <- heatmap_lgd("TOFA")

svg(filename = "Figures/output/TOFA_legend_cuadrada.svg",width=3,height=3)
legend <- heatmap_lgd("TOFA")

pdf(file = "Figures/output/TOFA_legend_cuadrada.pdf",width=3,height=3)
legend <- heatmap_lgd("TOFA")

