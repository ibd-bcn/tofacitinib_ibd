#Libraries
library(progeny)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(dplyr)

#Object TODAS
todas <- readRDS("/path/to/single_cell_object")
todas@meta.data[todas@meta.data$annotation_intermediate == "HS"  &
                  todas@meta.data$subset == "myeloids",]$annotation_intermediate <-
  "HS_Mac"
todas@meta.data[todas@meta.data$annotation_intermediate == "HS"  &
                  todas@meta.data$subset == "tcells",]$annotation_intermediate <-
  "HS_Tcells"

#First part: Run PROGENy; ------------------------------------------------------
prog_sub <-
  progeny(
    todas,
    scale = FALSE,
    organism = "Human",
    top = 500,
    perm = 1,
    return_assay = TRUE
  )

#Second part: Analyse PROGENy; -------------------------------------------------
# Change assay
DefaultAssay(object = prog_sub) <- "progeny"
# Scale the data
prog_sub <- ScaleData(prog_sub)
prog_sub@assays$progeny@data <- prog_sub@assays$progeny@scale.data
prog_sub$condition <-
  paste(prog_sub$week_3, "_", prog_sub$response, sep = "")

#Convert into DF
progeny_scores_df <-
  as.data.frame(t(GetAssayData(
    prog_sub, slot = "data",
    assay = "progeny"
  ))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity,-Cell)

#Add new annotation of cells
todas$condition <- paste(todas$week_3, "_", todas$response, sep = "")
meta <- todas@meta.data
meta$Cell <- rownames(meta)

#New annotation
cell_type_dictionary <- list(
  "Plasmablast_IGKC" = "Plasma cell",
  "PC_IGLL5" = "Plasma cell",
  "PC_IER" = "Plasma cell",
  "PC_PSAT1" = "Plasma cell",
  "Plasmablast_IgG" = "Plasma cell",
  "Plasmablast_IgA" = "Plasma cell",
  "PC_IgG" = "Plasma cell",
  "PC_heat_shock" = "Plasma cell",
  "PC_IgA" = "Plasma cell",
  "PC_IGLV6-57" = "Plasma cell",
  "PC_IFIT1" = "Plasma cell",
  "B_cell" = "B cell",
  "Cycling_B_cell" = "B cell",
  "CD4" = "T cell",
  "Th17" = "T cell",
  "HS_Tcells" = "T cell",
  "HS_Mac" = "Macrophages",
  "Tregs" = "T cell",
  "CD8" = "T cell",
  "Ribhi T cells" = "T cell",
  "NK" = "T cell",
  "MT hi IER" = "T cell",
  "Thf" = "T cell",
  "Naïve T cells" = "T cell",
  "CD4 CD8 IFIT3" = "T cell",
  "ILC3" = "T cell",
  "DN EOMES" = "T cell",
  "Tuft" = "Epithelium",
  "Colonocytes" = "Epithelium",
  "Secretory progenitor" = "Epithelium",
  "Epithlium Rib hi" = "Epithelium",
  "IER Epithelium" = "Epithelium",
  "Cycling TA" = "Epithelium",
  "Undifferentiated epithelium" = "Epithelium",
  "Stem" = "Epithelium",
  "Goblet" = "Epithelium",
  "Enteroendocrines" = "Epithelium",
  "Mature goblet" = "Epithelium",
  "APOA4" = "Epithelium",
  "M2" = "Macrophages",
  "M0" = "Macrophages",
  "Macrophage NRG1" = "Macrophages",
  "Neutrophils" = "Neutrophils",
  "M1" = "Macrophages",
  "Rib hi myeloids" = "Macrophages",
  "Mast" = "Mast cells",
  "Inflammatory monocytes" = "Inflammatory monocytes",
  "DCs" = "DCs",
  "Eosinophils" = "Eosinophils",
  "pDC" = "DCs",
  "Cycling myeloids" = "Macrophages",
  "S3" = "Fibroblasts",
  "Inflammatory_fibroblasts" = "Fibroblasts",
  "S2" = "Fibroblasts",
  "MT fibroblasts" = "Fibroblasts",
  "IER_fibroblasts" = "Fibroblasts",
  "Perycites" = "Endothelium",
  "Myofibroblasts" = "Fibroblasts",
  "S1" = "Fibroblasts",
  "Endothelium" = "Endothelium",
  "Glia" = "Glia",
  "Cycling_fibroblasts" = "Fibroblasts",
  "Cycling_T_cell" = "T cell",
  "Cycling_PC" = "Plasma cell",
  "Cycling_Myeloid" = "Macrophages"
)

meta$reduced_anot <-
  cell_type_dictionary[meta$annotation_intermediate]

#Merge metadata and progenyscore by cell
merged_df <- merge(progeny_scores_df, meta, by = "Cell")

#Split data per conditions
todas_w8_NR <-
  todas[, todas$week_3 == "POST" & todas$response == "NR"]
todas_w8_R <-
  todas[, todas$week_3 == "POST" & todas$response == "R"]
todas_w0_NR <-
  todas[, todas$week_3 == "W0" & todas$response == "NR"]
todas_w0_R <- todas[, todas$week_3 == "W0" & todas$response == "R"]

#Todas_Post_NR
pscores_w8_NR <-
  merged_df[merged_df$Cell %in% colnames(todas_w8_NR), ]
summarized_pscores_w8_NR <- pscores_w8_NR %>%
  group_by(Pathway, reduced_anot) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_pscores_w8_NR <- summarized_pscores_w8_NR %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
summarized_pscores_w8_NR$row_names <-
  rownames(summarized_pscores_w8_NR)
write_csv(
  summarized_pscores_w8_NR,
  paste0(
    "~/tofacitinib_ibd/Analysis/PROGENy/data/todas_reduced_anot_w8_NR.csv",
    sep = ""
  )
)

#Todas_Post_R
pscores_w8_R <- merged_df[merged_df$Cell %in% colnames(todas_w8_R), ]
summarized_pscores_w8_R <- pscores_w8_R %>%
  group_by(Pathway, reduced_anot) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_pscores_w8_R <- summarized_pscores_w8_R %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
summarized_pscores_w8_R$row_names <-
  rownames(summarized_pscores_w8_R)
write_csv(
  summarized_pscores_w8_R,
  paste0(
    "~/tofacitinib_ibd/Analysis/PROGENy/data/todas_reduced_anot_w8_R.csv",
    sep = ""
  )
)

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
    "~/tofacitinib_ibd/Analysis/PROGENy/data/todas_reduced_anot_w0_R.csv",
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
    "~/tofacitinib_ibd/Analysis/PROGENy/data/todas_reduced_anot_w0_NR.csv",
    sep = ""
  )
)

#Third part: PROGENy Statistics; ----------------------------------------------
#Specify Conditions/Combinations/Annotation
conditions <- c(unique(todas$condition))
combinations <- combn(conditions, 2)
annotation <- c("reduced_anot")

#Run the statistical analysis CONDITIONS
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
        "~/tofacitinib_ibd/Analysis/PROGENy/data/todas_",
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

