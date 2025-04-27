options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(readxl)
library(tidyr)
library(ggalluvial)
library(Seurat)
library(plyr)
library(dplyr)
library(pheatmap)
library(viridis)

source('Figures/functions_plots.R')

# Figure 5A Heatmap anti IL10

library(readxl)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
options(bitmapType='cairo')
# LEE LA BASE DE DATOS QUE TE DE Y LA HOJA QUE TOQUE
DMSO <-  read_excel("~/Elisa/il10_elisa.xlsx",
                    sheet = "FC Macs  antiIL-10")
DMSO <- DMSO %>%
  dplyr::mutate(across(where(is.numeric), log2))

colnames(DMSO)[2] <- "Sample"
colnames(DMSO)[3] <- "Condition"
DMSO <- DMSO[,2:ncol(DMSO)]

# Identify numeric columns in the elisa dataframe
numeric_cols <- sapply(DMSO, is.numeric)

#LPS anti
LPS_anti <- subset(DMSO, Condition == "LPS+anti-IL10")
LPS_anti <- data.frame(LPS_anti, row.names = NULL)

#TNF anti
TNF_anti <- subset(DMSO, Condition == "TNF+anti-IL10")
TNF_anti <- data.frame(TNF_anti, row.names = NULL)

#IFN anti
IFN_anti <- subset(DMSO, Condition == "IFN+anti-IL10")
IFN_anti <- data.frame(IFN_anti, row.names = NULL)


results <- LPS_anti %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS_anti[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS_anti[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNF_anti %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNF_anti[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNF_anti[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}


results <- IFN_anti %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFN_anti[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFN_anti[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}


heatmap_tt <- rbind(LPS_anti, TNF_anti, IFN_anti)
heatmap_tt <- heatmap_tt %>% relocate(Sample)

openxlsx::write.xlsx(heatmap_tt, "~/stats_il10new.xlsx", rowNames=T)




col_fun = colorRamp2(c(-3, 0, 5), c("green", "black", "red"))


row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <-
  c(
    "TNF",
    "SOCS3",
    "MX1",
    "CXCL1",
    "CXCL10",
    "CXCL8",
    "IL1B",
    "IL23A",
    "IL6",
    "INHBA",
    "CLEC5A",
    "IL10"
  )
numeric_cols <- sapply(DMSO, is.numeric)

DMSO$Condition <- as.character(DMSO$Condition)

lps_anti <- subset(DMSO, Condition == "LPS+anti-IL10")
lps_anti <- data.frame(lps_anti, row.names = NULL)
rownames(lps_anti) <- c("LPS+anti-IL10_1", "LPS+anti-IL10_2", "LPS+anti-IL10_3")
lps_anti <- lps_anti[, 3:length(colnames(lps_anti))]
lps_anti <- colMeans(lps_anti)

tnf_anti <- subset(DMSO, Condition == "TNF+anti-IL10")
tnf_anti <- data.frame(tnf_anti, row.names = NULL)
rownames(tnf_anti) <- c("TNF+anti-IL10_1", "TNF+anti-IL10_2", "TNF+anti-IL10_3")
tnf_anti <- tnf_anti[, 3:length(colnames(tnf_anti))]
tnf_anti <- colMeans(tnf_anti)

ifn_anti <- subset(DMSO, Condition == "IFN+anti-IL10")
ifn_anti <- data.frame(ifn_anti, row.names = NULL)
rownames(ifn_anti) <- c("IFN+anti-IL10_1", "IFN+anti-IL10_2", "IFN+anti-IL10_3")
ifn_anti <- ifn_anti[, 3:length(colnames(ifn_anti))]
ifn_anti <- colMeans(ifn_anti)


t_DMSO <- t(data.frame( LPS_anti = lps_anti, TNF_anti = tnf_anti, IFN_anti = ifn_anti))

t_DMSO <- t_DMSO[, gene_vector]
col_names <- colnames(t_DMSO)

valid_genes <- intersect(gene_vector, colnames(t_DMSO))
t_DMSO <- t_DMSO[, valid_genes, drop = FALSE]
rownames(t_DMSO) <- gsub("LPS_anti", "LPS+anti-IL10", rownames(t_DMSO))
rownames(t_DMSO) <- gsub("TNF_anti", "TNF+anti-IL10", rownames(t_DMSO))
rownames(t_DMSO) <- gsub("IFN_anti", "IFN+anti-IL10", rownames(t_DMSO))

#Statistics
macros_dmso_stats <-
  as.data.frame(read_excel("~/stats_il10new.xlsx"))
genes_adjusted <- paste(valid_genes, "_p_adjusted", sep = "")
macros_dmso_stats <-
  macros_dmso_stats[, c("Condition", genes_adjusted)]
macros_dmso_stats <- unique(macros_dmso_stats)
rownames(macros_dmso_stats) <- macros_dmso_stats$Condition
macros_dmso_stats <-
  macros_dmso_stats[c( "LPS+anti-IL10", "TNF+anti-IL10", "IFN+anti-IL10"),]
macros_dmso_stats <- macros_dmso_stats[, 2:ncol(macros_dmso_stats)]
sig_mat <-
  ifelse(macros_dmso_stats < 0.05,
         ifelse(
           macros_dmso_stats < 0.005,
           ifelse(
             macros_dmso_stats < 0.001,
             ifelse(macros_dmso_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")

png(
  "Figures/output/il10.png",
  width = 10,
  height = 4,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_DMSO,
  na_col = "white",
  name = "Legend",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  row_split = as.factor(c("IFN+anti-IL10","LPS+anti-IL10",  "TNF+anti-IL10")),
  row_title = NULL,
  row_names_gp = gpar(fontsize = 25),
  column_names_gp =  gpar(fontsize = 25, fontface = "italic"),
  column_names_rot = 60,
  border_gp = gpar(col = "white", lwd = 2),
  heatmap_legend_param = list(
    at = c(-7, 0, 7),  # Force the legend range
    color_bar = "continuous"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] == "1") {
      grid.text("*", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "2") {
      grid.text("**", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "3") {
      grid.text("***", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
    if (sig_mat[i, j] == "4") {
      grid.text("****", x  , y  , gp = gpar(fontsize = 40, col = "white"))
    }
  }
)




draw(heatmap)
dev.off()


# Figure 5B --------------------------------------------------------------------

## Scatter plots of our DE data and Cuevas et al.
setwd('Analysis/')
de_data <- readRDS('data/01_DE/REPASO/new_complete.RDS')

# Obtaining data of upp-regulated genes and down-regulated ones from Non-responders (NR) for the M2 population
de_data2 <- de_data[de_data$cluster == 'M2' &
                      de_data$annotation == 'annotation_intermediate' &
                      de_data$comp %in% 'w0NR_vs_POSTNR' &
                      de_data$sign %in% c('UPP', 'UP' , 'DWW', 'DW'), c("p_val", "avg_log2FC",  "sign", "comp", "gene")]

# Non-responders
df_w0NR_vs_POSTNR_UPP <-
  subset(
    de_data2,
    comp == "w0NR_vs_POSTNR" &
      sign == "UPP" |
      comp == "w0NR_vs_POSTNR" &
      sign == "UP" ,
    select = c("avg_log2FC", "p_val", "gene")
  )
filtered_df_w0NR_vs_POSTNR_UPP <-
  subset(df_w0NR_vs_POSTNR_UPP, p_val < 0.05)
filtered_df_w0NR_vs_POSTNR_UPP$condition <-
  rep("non_responder_UPP", nrow(filtered_df_w0NR_vs_POSTNR_UPP))

df_w0NR_vs_POSTNR_DWW <-
  subset(
    de_data2,
    comp == "w0NR_vs_POSTNR" &
      sign == "DWW" |
      comp == "w0NR_vs_POSTNR" &
      sign == "DW",
    select = c("avg_log2FC", "p_val", "gene")
  )
filtered_df_w0NR_vs_POSTNR_DWW <-
  subset(df_w0NR_vs_POSTNR_DWW, p_val < 0.05)
filtered_df_w0NR_vs_POSTNR_DWW$condition <-
  rep("non_responder_DWW", nrow(filtered_df_w0NR_vs_POSTNR_DWW))

dataframe_NR <-
  rbind(filtered_df_w0NR_vs_POSTNR_UPP,
        filtered_df_w0NR_vs_POSTNR_DWW) # Final dataframe Non-responders

# DEGs data from Cuevas et al
res_DE <-
  read_csv("Figures/extra_data/res_DE.csv") # DEGs nota para angela de subirlo al git
res_DE <- na.omit(res_DE)

res_DE_upp <-
  res_DE[res_DE$padj < 0.05 & res_DE$log2FoldChange >= 1, ]
res_DE_dww <-
  res_DE[res_DE$padj < 0.05 & res_DE$log2FoldChange <= -1, ]
res_DE_upp$condition <- rep("paper_upp", nrow(res_DE_upp))
res_DE_dww$condition <- rep("paper_dww", nrow(res_DE_dww))
resfinal <- rbind(res_DE_upp, res_DE_dww)



# Tidying the data for matching


resfinal <- resfinal[order(names(resfinal))]
names(resfinal)[names(resfinal) == "...1"] <- "gene"
names(resfinal)[names(resfinal) == "pvalue"] <- "p_val"
resfinal$lfcSE <- NULL
resfinal$baseMean <- NULL
resfinal$padj <- NULL
resfinal$stat <- NULL


names(dataframe_NR)[names(dataframe_NR) == "avg_log2FC"] <-
  "log2FoldChange"
dataframe_NR <- dataframe_NR[order(names(dataframe_NR))]

nonresponders <- rbind(resfinal, dataframe_NR)


nonresponders <- na.omit(nonresponders)
nonresponders$p_val <- NULL
nonresponders2 <-
  pivot_wider(nonresponders, names_from = "condition", values_from = "log2FoldChange")

filtered_up <-
  subset(nonresponders2,!is.na(paper_upp) &
           !is.na(non_responder_UPP))
filtered_dw <-
  subset(nonresponders2,!is.na(paper_dww) &
           !is.na(non_responder_DWW))


filtered_up$paper_dww <- NULL
filtered_up$non_responder_DWW <- NULL
filtered_dw$paper_upp <- NULL
filtered_dw$non_responder_UPP <- NULL
final <-
  merge(
    filtered_up,
    filtered_dw,
    by.y = "gene",
    all.x = TRUE,
    all.y  = TRUE
  )

# Plot

plot <- ggplot(data = final) +
  geom_point(aes(x = paper_upp, y = non_responder_UPP),
             color = "red",
             size = 1) +
  geom_point(aes(x = paper_dww, y = non_responder_DWW),
             color = "blue",
             size = 1) +
  geom_vline(xintercept = 0, colour = "gray") +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_label_repel(
    aes(label = gene, x = paper_upp, y = non_responder_UPP),
    box.padding = 0.5,
    size = 2.5,
    color = "black",
    nudge_y = 0.01,
    segment.size = 0.1,
    max.overlaps = Inf
  ) +

  geom_label_repel(
    data = final,
    aes(label = gene, x = paper_dww, y = non_responder_DWW),
    box.padding = 0.5,
    size = 2.5,
    color = "black",
    nudge_y = 0.01,
    segment.size = 0.1,
    max.overlaps = Inf
  ) +

  labs(x = "log2FoldChange", y = "avg log2FoldChange") +
  theme_bw() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6)
  )


print(plot)
save_sizes(plot = plot ,
           filename = '5E',
           device = 'jpeg')
save_sizes(plot = plot,
           filename = '5E',
           device = 'tiff')
save_sizes(plot = plot,
           filename = '5E',
           device = 'svg')
save_sizes(plot = plot,
           filename = '5E',
           device = 'pdf')

# Figure 5C --------------------------------------------------------------------
cellinfo <- read.csv('Figures/extra_data/cellinfo.csv')

cellinfo2 <- pivot_longer(cellinfo[cellinfo$TX == 'PRE', ],
                          cols = 1:3,
                          values_to = 'poblacio',
                          names_to = 'xaxis')
cellinfo2$xaxis <-
  factor(cellinfo2$xaxis,
         levels = c('RECEPTOR_R', 'SENDER', 'RECEPTOR_NR'))
cellinfo2 <- cellinfo2[!is.na(cellinfo2$poblacio), ]
cellinfo2$PROP[cellinfo2$xaxis == 'SENDER' &
                 cellinfo2$MEAN == 0 & cellinfo2$poblacio != 'M2'] <- 0

pre <- ggplot(
  cellinfo2,
  aes(
    x = xaxis,
    stratum = poblacio,
    alluvium = INDIVIDUAL,
    y = PROP,
    fill = MEAN,
    label = poblacio
  )
) +
  scale_x_discrete(expand = c(-.1,-.1)) +
  geom_flow() +
  geom_stratum(
    alpha = 1,
    width = 0.5,
    fill = c(
      'gray',
      # M2
      "#F8C05F",
      # IDA macrophages
      "#A28ED5",
      # DN
      "#FEF1CA",
      # cy my
      "gray",
      #M2
      "#6784A6",
      # M1
      "#EB7799",
      # Infl mono
      "#ABC57B",
      # Infl fibro
      "#F8C05F",
      # IDA macro
      "#A28ED5",
      # DN
      'gray',
      #M2
      "#F8C05F",
      # IDA macrophages
      "#A28ED5",
      # DN
      "#FEF1CA" # cy my
    )
  ) +
  scale_colour_gradient(
    low = "#EEFCFB00",
    high = "#21A59F",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_void() +
  labs(title = 'PRE tx')

cellinfo_post <- read.csv('Figures/extra_data/cellinfo2.csv')
cellinfo_post$xaxis <-
  factor(cellinfo_post$xaxis,
         levels = c('RECEPTOR_R', 'SENDER', 'RECEPTOR_NR'))


save_sizes(plot = pre ,
           filename = '5E',
           device = 'jpeg')
save_sizes(plot = pre,
           filename = '5E',
           device = 'tiff')
save_sizes(plot = pre,
           filename = '5E',
           device = 'svg')
save_sizes(plot = pre,
           filename = '5E',
           device = 'pdf')

# Figure 5D --------------------------------------------------------------------
#Read Single-Cell with cyt annotation (go to cytosig analysis)
todas <- readRDS("/path/to/single_cell_object")
todas <- todas[, todas$annotation_cyt != "NA"]
#Split by cond
todas_R <- todas[, todas$week == 'W0' & todas$response == "R"]
todas_NR <- todas[, todas$week == 'W0' & todas$response == "NR"]

todas_R@active.ident <- as.factor(todas$annotation_cyt)
a <-
  DotPlot(todas_R, features = c("IL10", "IL10RA", "IL10RB"))$data
a$gene <- paste(a$features.plot, "_R", sep = "")

todas_NR@active.ident <- as.factor(todas$annotation_cyt)
a1 <-
  DotPlot(todas_NR, features = c("IL10", "IL10RA", "IL10RB"))$data
a1$gene <- paste(a1$features.plot, "_NR", sep = "")

combined_a <- rbind(a, a1)

combined_a <- combined_a[combined_a$id != "Eosinophils",]

combined_a$id <-
  factor(
    combined_a$id,
    levels = c(
      "Colonocytes",
      "Goblet",
      "Tuft",
      "Enteroendocrines",
      "Stem",
      "Undifferentiated_epithelium",
      "Cycling_TA",
      "Macrophage",
      "Inflammatory_monocytes",
      "Neutrophils",
      "Mast",
      "DCs",
      "S1",
      "S2",
      "S3",
      "Inflammatory_fibroblasts",
      "Myofibroblasts",
      "Endothelium",
      "Perycites",
      "Glia",
      "CD4",
      "CD8",
      "CD4_CD8_IFIT3",
      "Thf",
      "Th17",
      "Tregs",
      "ILC",
      "DN_EOMES",
      "Naive_T_cells",
      "B_cell",
      "Cycling_B_cell",
      "PC_IgA",
      "PC_IgG",
      "Plasmablast_IgA",
      "Plasmablast_IgG"
    )
  )

png(
  "Figures/output/dotplot.png",
  width = 10,
  height = 4,
  units = "in",
  res = 600
)
ggplot(data = combined_a, aes(
  x = id,
  y = gene,
  color = avg.exp.scaled,
  size = pct.exp
)) +
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()
  ) +
  scale_color_gradientn(
    colours = viridis::viridis(20),
    limits = c(-0.93, 2.5),
    oob = scales::squish,
    name = 'avg.exp'
  )
dev.off()

# Figure 5E --------------------------------------------------------------------
#currently using todas_R and todas_NR from Figure 4D

cells_R <- colnames(todas_R)
cells_NR <- colnames(todas_NR)

cell_annotations_R <-
  todas_R@meta.data[, c("cell_name", "annotation_cyt")]
data_R <-
  as.data.frame(t(as.data.frame(todas_R[c("IL10", "IL10RA", "IL10RB"),]@assays$RNA@data)))
data_R$cell_name <- rownames(data_R)
data_R$cell_type <-
  mapvalues(
    x = data_R$cell_name ,
    from = cell_annotations_R$cell_name ,
    to = cell_annotations_R$annotation_cyt
  )
data_R$condition <- "R"

cell_annotations_NR <-
  todas_NR@meta.data[, c("cell_name", "annotation_cyt")]
data_NR <-
  as.data.frame(t(as.data.frame(todas_NR[c("IL10", "IL10RA", "IL10RB"),]@assays$RNA@data)))
data_NR$cell_name <- rownames(data_NR)
data_NR$cell_type <-
  mapvalues(
    x = data_NR$cell_name ,
    from = cell_annotations_NR$cell_name ,
    to = cell_annotations_NR$annotation_cyt
  )
data_NR$condition <- "NR"

data <- rbind(data_R, data_NR)

#Do it for each condition c("IL10","IL10RA","IL10RB")
results <- data %>%
  group_by(cell_type) %>%
  do({
    df <- .
    if (length(unique(df$condition)) == 2) {
      wilcox_result <-
        wilcox.test(IL10RB ~ condition, data = df, exact = FALSE)
      tibble(p_value = wilcox_result$p.value)
    } else {
      tibble(p_value = NA)
    }
  }) %>%
  ungroup()

q_values_adjusted <-
  p.adjust(na.omit(results$p_value), method = "bonferroni")
q_values_full <- rep(NA_real_, nrow(results))
q_values_full[!is.na(results$p_value)] <- q_values_adjusted
results$q_value <- q_values_full
results <- results[results$q_value < 0.05,]

#Cytosig part of the plot
todas_w0_R_vs_w0_NR_annotation_cyt <-
  read_excel("Analysis/Cytosig/comparisons_results/todas_W0_NR_vs_W0_R_annotation_cyt.xlsx")
R_vs_NR <-
  todas_w0_R_vs_w0_NR_annotation_cyt[todas_w0_R_vs_w0_NR_annotation_cyt$cytokine == "IL10",]
R_vs_NR <- R_vs_NR[!is.na(R_vs_NR$annotation), ]

#Eliminate EO
R_vs_NR <- R_vs_NR[R_vs_NR$annotation != "Eosinophils" ,]

#Grab significant ones
significant <- R_vs_NR[R_vs_NR$q_value < 0.05,]$annotation

df_hm <- data.frame(
  conditon = c(rep("IL10_R", times = 35), rep("IL10_NR", times = 35)),
  cyt_val = c(R_vs_NR$median_todas_W0_R, R_vs_NR$median_todas_W0_NR),
  annotation = c(R_vs_NR$annotation, R_vs_NR$annotation)
)
df_hm$value <- NA

new_row <- data.frame(
  conditon = "",
  cyt_val = c(10, 10),
  # Adjust this if you have specific values for this row
  annotation = c("Macrophage", "PC_IgA"),
  value = c("***", "*")
)

# Combine df_hm with the new row
df_hm <- rbind(df_hm, new_row)

df_hm$conditon <-
  factor(df_hm$conditon, levels = c("IL10_NR", "IL10_R", ""))
df_hm$annotation <-
  factor(
    df_hm$annotation,
    levels = c(
      "Colonocytes",
      "Goblet",
      "Tuft",
      "Enteroendocrines",
      "Stem",
      "Undifferentiated_epithelium",
      "Cycling_TA",
      "Macrophage",
      "Inflammatory_monocytes",
      "Neutrophils",
      "Mast",
      "DCs",
      "S1",
      "S2",
      "S3",
      "Inflammatory_fibroblasts",
      "Myofibroblasts",
      "Endothelium",
      "Perycites",
      "Glia",
      "CD4",
      "CD8",
      "CD4_CD8_IFIT3",
      "Thf",
      "Th17",
      "Tregs",
      "ILC",
      "DN_EOMES",
      "Naive_T_cells",
      "B_cell",
      "Cycling_B_cell",
      "PC_IgA",
      "PC_IgG",
      "Plasmablast_IgA",
      "Plasmablast_IgG"
    )
  )


png(
  "Figures/output/heatmap.png",
  width = 10,
  height = 3,
  units = "in",
  res = 600
)

ggplot(df_hm, aes(annotation, conditon, fill = cyt_val)) +
  geom_tile() + cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.key.size = unit(0.4, "cm"),
    axis.title.x = element_blank()
  ) +
  ylab('') +
  theme(axis.ticks = element_blank()) + scale_fill_viridis(
    option = "D",
    limits = c(-1.75, 5.72),
    na.value = "white"
  ) +
  labs(fill = NULL) + geom_text(aes(label = value), color = "black", size = 4)

dev.off()

