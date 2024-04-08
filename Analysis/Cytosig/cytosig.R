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

#Object TODAS
todas <- readRDS("/path/to/single_cell_object")

#First part: Prepare INPUT cytosig and RUN; ------------------------------------

#w0_R
todas_w0_R <- todas[,todas$week_3 == "W0" & todas$response == "R" ]
#Obtain COUNT MATRIX
counts <- todas_w0_R@assays$RNA@counts
# Calculate TPM
TPM <- counts / Matrix::colSums(counts) * 1e6
# Perform the log2(TPM/10 + 1) transformation
log_TPM <- log2(TPM/10 + 1)
# You can add this as a new assay to your Seurat object
todas_w0_R[["log_TPM"]] <- CreateAssayObject(counts = log_TPM)
#Write 10x Counts
print("TODAS_w0_R: Writting 10x Counts")
write10xCounts(x = todas_w0_R@assays$log_TPM@data, path = "/tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w0_R", version="3")
#Run Cytosig
print("TODAS: Run Cytosig")
system("taskset -c 0-20 CytoSig/CytoSig/CytoSig_run.py -i /tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w0_R/ -o /tofacitinib_ibd/Analysis/Cytosig/results/todas_w0_R.xlsx -r 1000 -a 10000  -e 1  -s 1 -c 1000 -z 0.95")

#w0_NR
todas_w0_NR <- todas[,todas$week_3 == "W0" & todas$response == "NR" ]
#Obtain COUNT MATRIX
counts <- todas_w0_NR@assays$RNA@counts
# Calculate TPM
TPM <- counts / Matrix::colSums(counts) * 1e6
# Perform the log2(TPM/10 + 1) transformation
log_TPM <- log2(TPM/10 + 1)
# You can add this as a new assay to your Seurat object
todas_w0_NR[["log_TPM"]] <- CreateAssayObject(counts = log_TPM)
#Write 10x Counts
print("TODAS_w0_NR: Writting 10x Counts")
write10xCounts(x = todas_w0_NR@assays$log_TPM@data, path = "/tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w0_NR", version="3")
#Run Cytosig
print("TODAS_w0_NR Run Cytosig")
system("taskset -c 0-20 CytoSig/CytoSig/CytoSig_run.py -i /tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w0_NR/ -o /tofacitinib_ibd/Analysis/Cytosig/results/todas_w0_NR.xlsx -r 1000 -a 10000  -e 1  -s 1 -c 1000 -z 0.95")


#w8_R
todas_w8_R <- todas[,todas$week_3 == "POST" & todas$response == "R" ]
#Obtain COUNT MATRIX
counts <- todas_w8_R@assays$RNA@counts
# Calculate TPM
TPM <- counts / Matrix::colSums(counts) * 1e6
# Perform the log2(TPM/10 + 1) transformation
log_TPM <- log2(TPM/10 + 1)
# You can add this as a new assay to your Seurat object
todas_w8_R[["log_TPM"]] <- CreateAssayObject(counts = log_TPM)
#Write 10x Counts
print("TODAS_w8_R: Writting 10x Counts")
write10xCounts(x = todas_w8_R@assays$log_TPM@data, path = "/tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w8_R", version="3")
#Run Cytosig
print("TODAS_w8_R Run Cytosig")
system("taskset -c 0-20 CytoSig/CytoSig/CytoSig_run.py -i /tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w8_R/ -o /tofacitinib_ibd/Analysis/Cytosig/results/todas_w8_R.xlsx -r 1000 -a 10000  -e 1  -s 1 -c 1000 -z 0.95")

#w8_NR
todas_w8_NR <- todas[,todas$week_3 == "POST" & todas$response == "NR" ]
#Obtain COUNT MATRIX
counts <- todas_w8_NR@assays$RNA@counts
# Calculate TPM
TPM <- counts / Matrix::colSums(counts) * 1e6
# Perform the log2(TPM/10 + 1) transformation
log_TPM <- log2(TPM/10 + 1)
# You can add this as a new assay to your Seurat object
todas_w8_NR[["log_TPM"]] <- CreateAssayObject(counts = log_TPM)
#Write 10x Counts
print("TODAS_w8_NR: Writting 10x Counts")
write10xCounts(x = todas_w8_NR@assays$log_TPM@data, path = "/tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w8_NR", version="3")
#Run Cytosig
print("TODAS_w8_NR Run Cytosig")
system("taskset -c 0-20 CytoSig/CytoSig/CytoSig_run.py -i /tofacitinib_ibd/Analysis/Cytosig/10xcounts/todas_w8_NR/ -o /tofacitinib_ibd/Analysis/Cytosig/results/todas_w8_NR.xlsx -r 1000 -a 10000  -e 1  -s 1 -c 1000 -z 0.95")


#Second part:Convert Cytosig output TO EXCELL-----------------------------------
print("TODAS: Converting OUTPUT")
system("python3 convertoexcel.py")

#Third part:Calculate statistical analysis -------------------------------------
print("TODAS: Statistical analysis")
anotacio_cytosig <- read_csv("/tofacitinib_ibd/Analysis/Cytosig/anot.csv")

#Modify TODAS and ADD column for annotation_cyt // cell names // condition
todas$cell_name <- colnames(todas)
todas$annotation_intermediate <- gsub(" ","_", todas$annotation_intermediate)

todas$annotation_cyt <- mapvalues(
  x =  todas$annotation_intermediate,
  from =  anotacio_cytosig$intermediate,
  to = anotacio_cytosig$intermediate_Cytosig
)

todas$condition <- paste(todas$week_3,"_",todas$response,sep = "")

#Specify Conditions/Combinations/Annotation
conditions <- c(unique(todas$condition))
combinations <- combn(conditions, 2)
annotation <- c("subset","annotation_intermediate","annotation_refined","annotation_cyt")

#Read cytosig data
cyt_data1 <- read_csv("/tofacitinib_ibd/Analysis/Cytosig/results/todas_w8_R_converted.csv")
colnames(cyt_data1)[which(names(cyt_data1) == "...1")] <- "cytokine"

cyt_data2 <- read_csv("/tofacitinib_ibd/Analysis/Cytosig/results/todas_w8_NR_converted.csv")
cyt_data2 <- cyt_data2[,2:ncol(cyt_data2)]

cyt_data3 <- read_csv("/tofacitinib_ibd/Analysis/Cytosig/results/todas_w0_R_converted.csv")
cyt_data3 <- cyt_data3[,2:ncol(cyt_data3)]

cyt_data4 <- read_csv("/tofacitinib_ibd/Analysis/Cytosig/results/todas_w0_NR_converted.csv")
cyt_data4 <- cyt_data4[,2:ncol(cyt_data4)]

cyt_data <- cbind(cyt_data1,cyt_data2,cyt_data3,cyt_data4)

cyt_data_long <- pivot_longer(cyt_data, cols = -cytokine, names_to = "cell_name", values_to = "z_score")

#Run the statistical analysis CONDITIONS
for(com in 1:6){
  for(anot in annotation){
    #Conditions
    cond1 <- combinations[1,com]
    cond2 <- combinations[2,com]

    # Prepare cell annotations
    cell_annotations <- todas@meta.data %>%
      dplyr::select(cell_name = cell_name, annotation = anot)

    #Select cells from your condition
    cells_cond1 <- todas@meta.data[todas@meta.data$condition == cond1,]$cell_name
    cells_cond2 <- todas@meta.data[todas@meta.data$condition == cond2,]$cell_name

    #Create group1_long & group2_long
    group1_long <- cyt_data_long[cyt_data_long$cell_name %in% cells_cond1,]
    group1_long$group <- "Group 1"
    group2_long <- cyt_data_long[cyt_data_long$cell_name %in% cells_cond2,]
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
          wilcox_result <- wilcox.test(z_score ~ group, data = data, exact = FALSE)
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
    q_values_adjusted <- p.adjust(na.omit(results$p_value), method = "bonferroni")

    # Prepare a vector for q-values matching the original data frame
    q_values_full <- rep(NA_real_, nrow(results))
    q_values_full[!is.na(results$p_value)] <- q_values_adjusted

    # Assign adjusted q-values to the results
    results$q_value <- q_values_full

    # Save the results
    colnames(results)[which(names(results) == "median_group1")]  <- paste0("median_todas_",cond1,sep = "")
    colnames(results)[which(names(results) == "median_group2")]  <- paste0("median_todas_",cond2,sep = "")
    openxlsx::write.xlsx(x = results, file = paste0("/tofacitinib_ibd/Analysis/Cytosig/comparisons_results/todas_",cond1,"_vs_",cond2,"_",anot,".xlsx",sep = ""))

  }




}

