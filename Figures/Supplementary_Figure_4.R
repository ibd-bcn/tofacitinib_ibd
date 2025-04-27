options(stringsAsFactors = FALSE,bitmapType = "cairo")
library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)

# Supplementary figure 4A: UMAP -----------------
colors <- c("Plasma cell" = "#DB324D",
            "T cell" = "#911EB4",
            "Epithelium" = "#66999B",
            "Macrophages" = "#FF69B4",
            "Mast cells" = "#A65628",
            "Inflammatory monocytes" = "#FF7F00",
            "DCs" = "#666666",
            "Cycling" = "#DCBEFF",
            "Fibroblasts" = "#44BBA4",
            "Endothelium" = "#029386",
            "Glia" = "#6495ED")

cell_type_dictionary <- list(
  "PC IgG"= "Plasma cell",
  "PC IgA" = "Plasma cell",
  "Plasmablasts" = "Plasma cell",
  "Cycling B cells" = "Plasma cell",
  "B cell" = "Plasma cell",
  "Naive T cells" = "T cell",
  "Tregs" = "T cell",
  "Th17" = "T cell",
  "CD8" = "T cell",
  "ILCs" = "T cell",
  "DN EOMES" = "T cell",
  "NKs" = "T cell",
  "Colonocytes" = "Epithelium",
  "Tuft" = "Epithelium",
  "Stem" = "Epithelium",
  "Cycling TA" = "Epithelium",
  "Enteroendocrines" = "Epithelium",
  "Goblet" = "Epithelium",
  "Macrophages" = "Macrophages",
  "Mast" = "Mast cells",
  "Inf monocytes" = "Inflammatory monocytes",
  "DCs" = "DCs",
  "Cycling" = "Cycling",
  "S3" = "Fibroblasts",
  "S2" = "Fibroblasts",
  "Perycites" = "Endothelium",
  "Myofibroblasts" = "Fibroblasts",
  "S1" = "Fibroblasts",
  "Endothelium" = "Endothelium",
  "Glia" = "Glia"

)

todas <- readRDS("~/thommas_et_al.RDS")
todas$reduced_anot <- todas$annotation

for (pattern in names(cell_type_dictionary)) {
  todas$reduced_anot <- gsub(pattern, cell_type_dictionary[[pattern]], todas$reduced_anot)
}



supfig4a <- DimPlot(todas, group.by = 'reduced_anot') +scale_color_manual(values = colors) +
  theme_umap()

save_sizes(plot = supfig4a, filename = 'supfig4a', device = 'tiff')
save_sizes(plot = supfig4a, filename = 'supfig4a', device = 'svg')
save_sizes(plot = supfig4a, filename = 'supfig4a', device = 'jpeg')
save_sizes(plot = supfig4a, filename = 'supfig4a', device = 'pdf')

supfig4a_legend <- supfig4a +
  theme(legend.position = 'right')

save_sizes(plot = supfig4a_legend, filename = 'supfig4a_legend', device = 'tiff')
save_sizes(plot = supfig4a_legend, filename = 'supfig4a_legend', device = 'svg')
save_sizes(plot = supfig4a_legend, filename = 'supfig4a_legend', device = 'jpeg')
save_sizes(plot = supfig4a_legend, filename = 'supfig4a_legend', device = 'pdf')


# Supplementary figure 4B: Progeny -----------------

prog_sub <-
  progeny(
    todas,
    scale = FALSE,
    organism = "Human",
    top = 500,
    perm = 1,
    return_assay = TRUE
  )


DefaultAssay(object = prog_sub) <- "progeny"

prog_sub <- ScaleData(prog_sub)
prog_sub@assays$progeny@data <- prog_sub@assays$progeny@scale.data
prog_sub$condition <- paste0(prog_sub$sample_id,"_", prog_sub$Remission_status)

progeny_scores_df <-
  as.data.frame(t(GetAssayData(
    prog_sub, layer = "data",
    assay = "progeny"
  ))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity,-Cell)


todas$condition <- paste0(todas$sample_id,"_", todas$Remission_status)
meta <- todas@meta.data
meta$Cell <- rownames(meta)

#New annotation
cell_type_dictionary <- list(
  "PC IgG"= "Plasma cell",
  "PC IgA" = "Plasma cell",
  "Plasmablasts" = "Plasma cell",
  "Cycling B cells" = "Plasma cell",
  "B cell" = "Plasma cell",
  "Naive T cells" = "T cell",
  "Tregs" = "T cell",
  "Th17" = "T cell",
  "CD8" = "T cell",
  "ILCs" = "T cell",
  "DN EOMES" = "T cell",
  "NKs" = "T cell",
  "Colonocytes" = "Epithelium",
  "Tuft" = "Epithelium",
  "Stem" = "Epithelium",
  "Cycling TA" = "Epithelium",
  "Enteroendocrines" = "Epithelium",
  "Goblet" = "Epithelium",
  "Macrophages" = "Macrophages",
  "Mast" = "Mast cells",
  "Inf monocytes" = "Inflammatory monocytes",
  "DCs" = "DCs",
  "Cycling" = "Cycling",
  "S3" = "Fibroblasts",
  "S2" = "Fibroblasts",
  "Perycites" = "Endothelium",
  "Myofibroblasts" = "Fibroblasts",
  "S1" = "Fibroblasts",
  "Endothelium" = "Endothelium",
  "Glia" = "Glia"

)

meta$reduced_anot <-
  cell_type_dictionary[meta$annotation]

#Merge metadata and progenyscore by cell
merged_df <- merge(progeny_scores_df, meta, by = "Cell")

#Split data per conditions
todas_w0_NR <-
  todas[, todas$Remission_status == "Non_Remission"]
todas_w0_R <- todas[,  todas$Remission_status == "Remission"]



#Todas_Pre_R
pscores_w0_R <- merged_df[merged_df$Cell %in% colnames(todas_w0_R), ]
summarized_pscores_w0_R <- pscores_w0_R %>%
  group_by(Pathway, reduced_anot) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_pscores_w0_R <- summarized_pscores_w0_R %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
summarized_pscores_w0_R$row_names <-
  rownames(summarized_pscores_w0_R)
write_csv(
  summarized_pscores_w0_R,
  paste0(
    "~/todas_reduced_anot_w0_R.csv",
    sep = ""
  )
)

#Todas_Pre_NR
pscores_w0_NR <-
  merged_df[merged_df$Cell %in% colnames(todas_w0_NR), ]
summarized_pscores_w0_NR <- pscores_w0_NR %>%
  group_by(Pathway, reduced_anot) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_pscores_w0_NR <- summarized_pscores_w0_NR %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
summarized_pscores_w0_NR$row_names <-
  rownames(summarized_pscores_w0_NR)
write_csv(
  summarized_pscores_w0_NR,
  paste0(
    "~/todas_reduced_anot_w0_NR.csv",
    sep = ""
  )
)


conditions <- c(unique(todas$condition))
combinations <- combn(conditions, 2)
annotation <- c("reduced_anot")

for (com in 1:6) {
  for (anot in annotation) {
    #Conditions
    cond1 <- combinations[1, com]
    cond2 <- combinations[2, com]
    dd <-
      merged_df[merged_df$condition %in% c(cond1, cond2), c("Cell", "Pathway", "Activity", anot, "condition")]
    colnames(dd) <-
      c("Cell", "Pathway", "Activity", "annotation", "condition")
    dd$group <- ifelse(dd$condition == cond1, "Group 1", "Group 2")

    # Perform Wilcoxon rank-sum tests and calculate medians
    results <- dd %>%
      group_by(Pathway, annotation) %>%
      do({
        data <- .
        if (length(unique(data$group)) == 2) {
          wilcox_result <-
            wilcox.test(Activity ~ group, data = data, exact = FALSE)
          tibble(
            p_value = wilcox_result$p.value,
            median_group1 = median(data$Activity[data$group == "Group 1"]),
            median_group2 = median(data$Activity[data$group == "Group 2"])
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


    q_values_adjusted <-
      p.adjust(na.omit(results$p_value), method = "bonferroni")


    q_values_full <- rep(NA_real_, nrow(results))
    q_values_full[!is.na(results$p_value)] <- q_values_adjusted


    results$q_value <- q_values_full

    colnames(results)[which(names(results) == "median_group1")]  <-
      paste0("median_todas_", cond1, sep = "")
    colnames(results)[which(names(results) == "median_group2")]  <-
      paste0("median_todas_", cond2, sep = "")

    openxlsx::write.xlsx(
      x = results,
      file = paste0(
        "~/progeny/",
        cond1,
        "_vs_",
        cond2,
        "_",
        gsub("annotation_", "", anot),
        ".xlsx",
        sep = ""
      )
    )

  }

}



heatmap_progeny <- function(stimuli, anot = "reduced_anot") {
  #Open dfs
  todas_intermediate_w0_NR <-
    read_csv(paste0("~/todas_reduced_anot_w0_NR.csv", sep = ""))
  cn_w0_NR <- todas_intermediate_w0_NR$row_names
  jak_w0_NR <- todas_intermediate_w0_NR[[stimuli]]

  todas_intermediate_w0_R <-
    read_csv(paste0("~/todas_reduced_anot_w0_R.csv", sep = ""))
  cn_w0_R <- todas_intermediate_w0_R$row_names
  jak_w0_R <- todas_intermediate_w0_R[[stimuli]]


  #Modify df
  df_w0_NR <- data.frame(cn = cn_w0_NR, jak_w0_NR = jak_w0_NR)
  df_w0_R <- data.frame(cn = cn_w0_R, jak_w0_R = jak_w0_R)


  combined_df <- df_w0_NR %>%
    full_join(df_w0_R, by = "cn")


  colnames(combined_df) <- c("cn", "w0_NR", "w0_R")
  cells <- combined_df$cn
  rownames(combined_df) <- combined_df$cn
  combined_df <- combined_df[, -c(1)]
  data_t <- t(combined_df)
  colnames(data_t) <- cells
  data_t <- data_t[c("w0_NR", "w0_R"), ]

  available_cells <- colnames(data_t)

  desired_cells <- c(
    "T cell", "Plasma cell", "B cell", "Macrophages",
    "Neutrophils", "Mast cells", "Inflammatory monocytes", "DCs",
    "Eosinophils", "Fibroblasts", "Endothelium", "Glia", "Epithelium"
  )

  filtered_cells <- desired_cells[desired_cells %in% available_cells]

  data <- data_t[, filtered_cells, drop = FALSE]

  paletteLength = 100
  myColor <-
    colorRamp2(range(na.omit(data)),
               hcl_palette = "Reds",
               reverse = TRUE)
  data <- data[c("w0_NR", "w0_R"), ]
  colnames(data) <-
    gsub(pattern = "Inflammatory monocytes",
         replacement = "Inf mono",
         x = colnames(data))

  #Heatmap
  progeny_hmap <- Heatmap(
    data,
    name = "PROGENy (500)",
    col = myColor,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = NA),
    row_title = NULL,
    column_title = NULL,
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 18)
  )

  # Draw the heatmap
  draw(
    progeny_hmap,
    heatmap_legend_side = "right",
    annotation_legend_side = "left"
  )
}
heatmap_progeny(stimuli = "TNFa") #change per condition

cond1 <- "Remission"
cond2 <- "Non_Remission"


dd <- merged_df[merged_df$Remission_status %in% c(cond1, cond2), c("Cell", "Pathway", "Activity", anot, "Remission_status")]

print(dim(dd))

colnames(dd) <- c("Cell", "Pathway", "Activity", "annotation", "Remission_status")
dd$group <- ifelse(dd$Remission_status == cond1, "Group 1", "Group 2")


results <- dd %>%
  group_by(Pathway, annotation) %>%
  do({
    data <- .
    if (length(unique(data$group)) == 2) {
      wilcox_result <- wilcox.test(Activity ~ group, data = data, exact = FALSE)
      tibble(
        p_value = wilcox_result$p.value,
        median_group1 = median(data$Activity[data$group == "Group 1"]),
        median_group2 = median(data$Activity[data$group == "Group 2"])
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


q_values_adjusted <- p.adjust(na.omit(results$p_value), method = "bonferroni")
q_values_full <- rep(NA_real_, nrow(results))
q_values_full[!is.na(results$p_value)] <- q_values_adjusted
results$q_value <- q_values_full

colnames(results)[which(names(results) == "median_group1")] <- paste0("median_todas_", cond1)
colnames(results)[which(names(results) == "median_group2")] <- paste0("median_todas_", cond2)

openxlsx::write.xlsx(x = results, file = paste0("~/progeny/", cond1, "_vs_", cond2, "_", gsub("annotation_", "", anot), ".xlsx"))


