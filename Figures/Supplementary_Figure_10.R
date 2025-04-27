library(readxl)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
options(bitmapType='cairo')


# Supplementary figure 10 A ------------------------

library(readxl)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
options(bitmapType='cairo')
# LEE LA BASE DE DATOS QUE TE DE Y LA HOJA QUE TOQUE
DMSO <-
  read_excel("~/Elisa/250212 Base de datos UC Macs JCC Revision Angela Heatmap 2.xlsx",
             sheet = "FC Macs inhibitors")



DMSO <- DMSO %>%
  dplyr::mutate(across(where(is.numeric), log2))

colnames(DMSO)[2] <- "Sample"
colnames(DMSO)[3] <- "Condition"
DMSO <- DMSO[DMSO$Condition != "M-DMSO-IL4", ]
DMSO <- DMSO[,2:ncol(DMSO)]

# Identify numeric columns in the elisa dataframe
numeric_cols <- sapply(DMSO, is.numeric)

#LPS tofa
LPS_tofa <- subset(DMSO, Condition == "LPS+TOFA")
LPS_tofa <- data.frame(LPS_tofa, row.names = NULL)

# filgo
LPS_filgo <- subset(DMSO, Condition == "LPS+FILGO")
LPS_filgo <- data.frame(LPS_filgo, row.names = NULL)

# upa
lps_upa <- subset(DMSO, Condition == "LPS+UPA")
lps_upa <- data.frame(lps_upa, row.names = NULL)

results <- LPS_tofa %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS_tofa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS_tofa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- lps_upa %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  lps_upa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  lps_upa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS_filgo %>%
  select(-Condition, -Sample) %>%
  dplyr::summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS_filgo[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS_filgo[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}


heatmap_tt <- rbind(LPS_tofa, lps_upa, LPS_filgo)

heatmap_tt <- heatmap_tt %>% relocate(Sample)


openxlsx::write.xlsx(heatmap_tt, "~/stats_inhibitors_new.xlsx", rowNames=T)




col_fun = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))




# colnames(DMSO)[1] <- "Condition"
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
    "IL10",
    "CLEC5A"
  )
numeric_cols <- sapply(DMSO, is.numeric)

#TOFA
DMSO$Condition <- as.character(DMSO$Condition)
TOFA <- subset(DMSO, Condition == "LPS+TOFA")
TOFA <- data.frame(TOFA, row.names = NULL)
rownames(TOFA) <- c("LPS+TOFA_1", "LPS+TOFA_2", "LPS+TOFA_3", "LPS+TOFA_4", "LPS+TOFA_5", "LPS+TOFA_6")
TOFA <- TOFA[, 4:length(colnames(TOFA))]
TOFA <- colMeans(TOFA)
#FILGO
FILGO <- subset(DMSO, Condition == "LPS+FILGO")
FILGO <- data.frame(FILGO, row.names = NULL)
rownames(FILGO) <- c("LPS+FILGO_1", "LPS+FILGO_2", "LPS+FILGO_3", "LPS+FILGO_4", "LPS+FILGO_5", "LPS+FILGO_6")
FILGO <- FILGO[, 3:length(colnames(FILGO))]
FILGO <- colMeans(FILGO)
#UPA
UPA <- subset(DMSO, Condition == "LPS+UPA")
UPA <- data.frame(UPA, row.names = NULL)
rownames(UPA) <- c("LPS+UPA", "LPS+UPA_2", "LPS+UPA_3", "LPS+UPA_4", "LPS+UPA_5", "LPS+UPA_6")
UPA <- UPA[, 3:length(colnames(UPA))]
UPA <- colMeans(UPA)

t_DMSO <- t(data.frame(
  FILGO = FILGO,
  UPA = UPA))

t_DMSO <- t_DMSO[, gene_vector]
col_names <- colnames(t_DMSO)

valid_genes <- intersect(gene_vector, colnames(t_DMSO))
t_DMSO <- t_DMSO[, valid_genes, drop = FALSE]
## SOCS3 CXCL10 CXCL1 ILB IL23A IL6 IL10



#Statistics
macros_dmso_stats <-
  as.data.frame(read_excel("~/stats_inhibitors_new.xlsx"))
genes_adjusted <- paste(valid_genes, "_p_adjusted", sep = "")
macros_dmso_stats <-
  macros_dmso_stats[, c("Condition", genes_adjusted)]
macros_dmso_stats <- unique(macros_dmso_stats)
rownames(macros_dmso_stats) <- macros_dmso_stats$Condition
macros_dmso_stats <-
  macros_dmso_stats[c( "LPS+FILGO", "LPS+UPA"),]
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
  "~/inhibitors_new.png",
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
  row_split = as.factor(c("FILGO", "UPA")),
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



# Supplementary figure 10 B ------------------------------------


df <-read_excel("Elisa/250212 Base de datos UC Macs JCC Revision Angela Heatmap 4.xlsx", sheet = "FC Fibs JAKin")

col_fun = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))

colnames(df)[1] <- "Sample"
colnames(df)[2] <- "Condition"
df <- df %>%
  dplyr::mutate(across(where(is.numeric), log2))

# colnames(df)[1] <- "Condition"
row_title <- gpar(fontsize = 14)
col_names <- gpar(fontsize = 14)
gene_vector <-
  c(
    "SOCS3", "CXCL10", "CXCL1", "IL6", "IL1B", "TNF", "CHI3L1", "WNT5A"

  )
numeric_cols <- sapply(df, is.numeric)


#FILGO
FILGO <- subset(df, Condition == "LPS+FILGO")
FILGO <- data.frame(FILGO, row.names = NULL)
rownames(FILGO) <- c("LPS+FILGO_1", "LPS+FILGO_2", "LPS+FILGO_3", "LPS+FILGO_4", "LPS+FILGO_5","LPS+FILGO_6", "LPS+FILGO_7")
FILGO <- FILGO[, 3:length(colnames(FILGO))]
FILGO <- colMeans(FILGO)
#UPA
UPA <- subset(df, Condition == "LPS+UPA")
UPA <- data.frame(UPA, row.names = NULL)
rownames(UPA) <- c("LPS+UPA", "LPS+UPA_2", "LPS+UPA_3", "LPS+UPA_4", "LPS+UPA_5", "LPS+UPA_6", "LPS+UPA_7")
UPA <- UPA[, 3:length(colnames(UPA))]
UPA <- colMeans(UPA)

t_df <- t(data.frame(
  FILGO = FILGO,
  UPA = UPA))

t_df <- t_df[, gene_vector]
col_names <- colnames(t_df)

valid_genes <- intersect(gene_vector, colnames(t_df))
t_df <- t_df[, valid_genes, drop = FALSE]
## SOCS3 CXCL10 CXCL1 ILB IL23A IL6 IL10



#Statistics
fibs_df_stats <-
  as.data.frame(read_excel("~/stats_inhibitors_fibs_0319.xlsx"))
genes_adjusted <- paste(valid_genes, "_p_adjusted", sep = "")
fibs_df_stats <- fibs_df_stats[, c("Condition", genes_adjusted)]
fibs_df_stats <- unique(fibs_df_stats)
rownames(fibs_df_stats) <- fibs_df_stats$Condition
fibs_df_stats <-
  fibs_df_stats[c( "LPS+FILGO", "LPS+UPA"),]
fibs_df_stats <- fibs_df_stats[, 2:ncol(fibs_df_stats)]
sig_mat <-
  ifelse(fibs_df_stats < 0.05,
         ifelse(
           fibs_df_stats < 0.005,
           ifelse(
             fibs_df_stats < 0.001,
             ifelse(fibs_df_stats < 0.0001, "4", "3"),
             "2"
           ),
           "1"
         ),
         "")
# col_fun = colorRamp2(c(-7, 0, 7), c("green", "black", "red"))
## Obtain database
png(
  "~/fibs_inhibitors_new.png",
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
heatmap <- Heatmap(
  t_df,
  na_col = "white",
  name = "Legend",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  row_split = as.factor(c( "FILGO", "UPA")),
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
