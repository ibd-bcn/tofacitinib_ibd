#Libraries
library(Seurat)
library(DropletUtils)
library(plyr)
library(readxl)
library(dplyr)
library(tidyr)
library(future.apply)
library(tibble)
library(readr)

# Supplementary 12 ----------------------------------------
#Object TODAS
todas <- readRDS("~/Thommasetal.RDS")

#First part: Prepare INPUT cytosig and RUN; -----------------------------------

#w0_R
todas_w0_R <- todas[, todas$Remission_status == "Remission"]
#Obtain COUNT MATRIX
counts <- todas_w0_R@assays$RNA@counts
# Calculate TPM
TPM <- counts / Matrix::colSums(counts) * 1e6
# Perform the log2(TPM/10 + 1) transformation
log_TPM <- log2(TPM / 10 + 1)
# You can add this as a new assay to your Seurat object
todas_w0_R[["log_TPM"]] <- CreateAssayObject(counts = log_TPM)
#Write 10x Counts
print("TODAS_w0_R: Writting 10x Counts")
write10xCounts(x = todas_w0_R@assays$log_TPM@data,
               path = "~/Elisa/adalimumab/todas_w0_R",
               version = "3")
#Run Cytosig
print("TODAS: Run Cytosig")
system(
  "taskset -c 0-20 ~/CytoSig/CytoSig/CytoSig_run.py -i ~/Elisa/adalimumab/todas_w0_R/ -o Analysis/Cytosig/results/todas_w0_R.xlsx -r 1000 -a 10000  -e 1  -s 1 -c 1000 -z 0.95"
)



#w0_NR
todas_w0_NR <-todas[, todas$Remission_status == "Non_Remission"]
#Obtain COUNT MATRIX
counts <- todas_w0_NR@assays$RNA@counts
# Calculate TPM
TPM <- counts / Matrix::colSums(counts) * 1e6
# Perform the log2(TPM/10 + 1) transformation
log_TPM <- log2(TPM / 10 + 1)
# You can add this as a new assay to your Seurat object
todas_w0_NR[["log_TPM"]] <- CreateAssayObject(counts = log_TPM)
#Write 10x Counts
print("TODAS_w0_NR: Writting 10x Counts")
write10xCounts(x = todas_w0_NR@assays$log_TPM@data,
               path = "~/Elisa/adalimumab/todas_w0_NR",
               version = "3")
#Run Cytosig
print("TODAS_w0_NR Run Cytosig")
system(
  "taskset -c 0-10 ~/CytoSig/CytoSig_run.py -i ~/Elisa/adalimumab/todas_w0_NR/ -o ~/Elisa/adalimumab/todas_w0_NR.xlsx -r 1000 -a 10000  -e 1  -s 1 -c 1000 -z 0.95"
)

#Second part:Convert Cytosig output TO EXCELL-----------------------------------
print("TODAS: Converting OUTPUT")
system("python3 convertoexcel.py")

#Third part:Calculate statistical analysis -------------------------------------
print("TODAS: Statistical analysis")
anotacio_cytosig <-
  read_csv("Analysis/Cytosig/anot.csv")

#Modify TODAS and ADD column for annotation_cyt // cell names // condition
todas$cell_name <- colnames(todas)

# todas$annotation_cyt <- mapvalues(
#   x =  todas$annotation_intermediate,
#   from =  anotacio_cytosig$intermediate,
#   to = anotacio_cytosig$intermediate_Cytosig
# )
#
# todas$condition <- paste(todas$week_3, "_", todas$response, sep = "")

#Specify Conditions/Combinations/Annotation
conditions <- c(unique(todas$Remission_status))
combinations <- combn(conditions, 2)
annotation <-
  c("subset",
    "annotation")

#Read cytosig data
cyt_data1 <-
  read_csv("~/Elisa/adalimumab/todas_w0_R_converted.csv")
colnames(cyt_data1)[which(names(cyt_data1) == "...1")] <- "cytokine"

cyt_data2 <-
  read_csv("~/Elisa/adalimumab/todas_w0_NR_converted.csv")
cyt_data2 <- cyt_data2[, 2:ncol(cyt_data2)]

cyt_data <- cbind(cyt_data1,cyt_data2)

cyt_data_long <-
  pivot_longer(
    cyt_data,
    cols = -cytokine,
    names_to = "cell_name",
    values_to = "z_score"
  )

#Run the statistical analysis CONDITIONS
for (com in 1:6) {
  for (anot in annotation) {
    #Conditions
    cond1 <- combinations[1, com]
    cond2 <- combinations[2, com]

    # Prepare cell annotations
    cell_annotations <- todas@meta.data %>%
      dplyr::select(cell_name = cell_name, annotation = anot)

    #Select cells from your condition
    cells_cond1 <-
      todas@meta.data[todas@meta.data$Remission_status == cond1, ]$cell_name
    cells_cond2 <-
      todas@meta.data[todas@meta.data$Remission_status == cond2, ]$cell_name

    #Create group1_long & group2_long
    group1_long <-
      cyt_data_long[cyt_data_long$cell_name %in% cells_cond1, ]
    group1_long$group <- "Group 1"
    group2_long <-
      cyt_data_long[cyt_data_long$cell_name %in% cells_cond2, ]
    group2_long$group <- "Group 2"

    # Combine and annotate data
    combined_data <- bind_rows(group1_long, group2_long) %>%
      left_join(cell_annotations, by = "cell_name")

    # Perform Wilcoxon rank-sum tests and calculate medians
    results <- combined_data %>%
      group_by(cytokine, annotation) %>%
      do({
        data <- .
        if (length(unique(data$group)) == 2) {
          wilcox_result <-
            wilcox.test(z_score ~ group, data = data, exact = FALSE)
          tibble(
            p_value = wilcox_result$p.value,
            median_group1 = median(data$z_score[data$group == "Group 1"]),
            median_group2 = median(data$z_score[data$group == "Group 2"])
          )
        } else {
          tibble(
            p_value = NA_real_,
            median_group1 = NA_real_,
            median_group2 = NA_real_
          )
        }
      }) %>%
      ungroup()

    # Adjust p-values for multiple testing (Bonferroni correction)
    q_values_adjusted <-
      p.adjust(na.omit(results$p_value), method = "bonferroni")

    # Prepare a vector for q-values matching the original data frame
    q_values_full <- rep(NA_real_, nrow(results))
    q_values_full[!is.na(results$p_value)] <- q_values_adjusted

    # Assign adjusted q-values to the results
    results$q_value <- q_values_full

    # Save the results
    colnames(results)[which(names(results) == "median_group1")]  <-
      paste0("median_todas_", cond1, sep = "")
    colnames(results)[which(names(results) == "median_group2")]  <-
      paste0("median_todas_", cond2, sep = "")
    openxlsx::write.xlsx(
      x = results,
      file = paste0(
        "~/Elisa/adalimumab/",
        cond1,
        "_vs_",
        cond2,
        "_",
        anot,
        ".xlsx",
        sep = ""
      )
    )

  }




}



todas_R@active.ident <- as.factor(todas$annotation)
a <-
  DotPlot(todas_R, features = c("IL10", "IL10RA", "IL10RB"))$data
a$gene <- paste(a$features.plot, "_R", sep = "")

todas_NR@active.ident <- as.factor(todas$annotation)
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


ggplot(data = combined_a, aes(
  x = id,
  y = gene,
  color = avg.exp.scaled,
  size = pct.exp
)) +
  geom_point() + cowplot::theme_cowplot() +
  # theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  # theme(
  #   axis.ticks = element_blank(),
  #   axis.text.x = element_blank(),
  #   axis.title = element_blank()
  #
  scale_color_gradientn(
    colours = viridis::viridis(20),
    limits = c(-0.93, 2.5),
    oob = scales::squish,
    name = 'avg.exp'
  )




cells_R <- colnames(todas_R)
cells_NR <- colnames(todas_NR)

cell_annotations_R <-
  todas_R@meta.data[, c("cell_name", "annotation")]
data_R <-
  as.data.frame(t(as.data.frame(todas_R[c("IL10", "IL10RA", "IL10RB"),]@assays$RNA@data)))
data_R$cell_name <- rownames(data_R)
data_R$cell_type <-
  mapvalues(
    x = data_R$cell_name ,
    from = cell_annotations_R$cell_name ,
    to = cell_annotations_R$annotation
  )
data_R$condition <- "R"

cell_annotations_NR <-
  todas_NR@meta.data[, c("cell_name", "annotation")]
data_NR <-
  as.data.frame(t(as.data.frame(todas_NR[c("IL10", "IL10RA", "IL10RB"),]@assays$RNA@data)))
data_NR$cell_name <- rownames(data_NR)
data_NR$cell_type <-
  mapvalues(
    x = data_NR$cell_name ,
    from = cell_annotations_NR$cell_name ,
    to = cell_annotations_NR$annotation
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


todas_w0_R_vs_w0_NR_annotation_cyt <-
  read_excel("adalimumab/Non_Remission_vs_Remission_annotation.xlsx")
R_vs_NR <-
  todas_w0_R_vs_w0_NR_annotation_cyt[todas_w0_R_vs_w0_NR_annotation_cyt$cytokine == "IL10",]
R_vs_NR <- R_vs_NR[!is.na(R_vs_NR$annotation), ]

significant <- R_vs_NR[R_vs_NR$q_value < 0.05,]$annotation





df_hm <- data.frame(
  conditon = c(rep("IL10_R", times = 30), rep("IL10_NR", times = 30)),
  cyt_val = c(R_vs_NR$median_todas_Remission, R_vs_NR$median_todas_Non_Remission),
  annotation = c(R_vs_NR$annotation, R_vs_NR$annotation)
)
df_hm$value <- NA

df_hm$conditon <-
  factor(df_hm$conditon, levels = c("IL10_NR", "IL10_R", ""))
new_row <- data.frame(
  conditon = "",
  cyt_val = c(10, 10),

  annotation = c("Macrophages", "Glia", "Goblet", "Stem"),
  value = c("***", "*", "*", "***")
)

df_hm$annotation <- as.character(df_hm$annotation)
df_hm <- rbind(df_hm, new_row)
df_hm$annotation <-
  factor(
    df_hm$annotation,
    levels = c(
      "Colonocytes",
      "Goblet",
      "Tuft",
      "Stem",
      "Enteroendocrines",
      "Cycling TA",
      "Macrophages",
      "Inf monocytes",
      "Mast",
      "DCs",
      "Cycling",
      "S1",
      "S2",
      "S3",
      "Myofibroblasts",
      "Endothelium",
      "Perycites",
      "Glia",
      "Naive T cells",
      "CD8",
      "Th17",
      "Tregs",
      "ILCs",
      "DN EOMES",
      "NKs",
      "B cell",
      "Cycling B cells",
      "PC IgA",
      "PC IgG",
      "Plasmablasts"
    )
  )

df_hm$annotation
png(
  "~/Elisa/adalimumab/heatmap2.png",
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
    limits = c(-4, 6),
    na.value = "white"
  ) +
  labs(fill = NULL) + geom_text(aes(label = value), color = "black", size = 4)

dev.off()
