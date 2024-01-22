#
# Libraries
#
message('Loading libraries')

library(DropletUtils)
library(SummarizedExperiment)
library(Seurat)
library(plyr)
library(ggplot2)
library(beepr)
library(celda)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scDblFinder)
library(viridis)
library(MASS)
library(harmony)
library(readr)
library(readxl)
library(openxlsx)
library(dplyr)
library(WriteXLS)

#
# Functions
#

#Obtain cell_phone data function
obtain_cell_phone <- function(directory, pvalues, significant_means) {
  message(directory)
  message("Working on that directory")
  #Create empty dataframe for cell_phone
  columns <- c(
    "Cluster_1",
    "Cluster_2",
    "Subset_1",
    "Subset_2",
    "id_cp_interaction",
    "interacting_pair",
    "partner_a",
    "partner_b",
    "gene_a",
    "gene_b",
    "secreted",
    "receptor_a",
    "receptor_b",
    "annotation_strategy",
    "is_integrin",
    "variable",
    "mean",
    "pvalue",
    "log10p"
  )
  cell_phone <- data.frame(matrix(vector(), ncol = 19))
  colnames(cell_phone) <- columns

  #Anotation file
  anotation <- read_excel("TOFA_cells_annotation.xlsx")
  anotation$cluster <- gsub(" ", "_", anotation$cluster)

  anotation[61, 2] <- "Naive_T_cells"
  anotation[38, 2] <- "PC_IGLV6_57"
  anotation[1, 2] <- "Cycling_B_cell_cycling"
  anotation[20, 2] <- "HS_myeloids"
  anotation[31, 2] <- "Cycling_B_cell_plasmas"
  anotation[58, 2] <- "HS_tcells"

  #Extract all the IDs interactions
  id_interactions <- significant_means$id_cp_interaction

  #Extract cell interactions and informative columns
  extract_cell_interactions <- colnames(significant_means)
  information_columns <- extract_cell_interactions[1:12]
  extract_cell_interactions <-
    extract_cell_interactions[13:length(extract_cell_interactions)]

  #Function to obtain the cell_phone
  y <- 0
  #Iterate over the cellxcell
  for (x in 1:length(extract_cell_interactions)) {
    #Obtain dataframe with a wanted information
    information <-
      significant_means[1:12] #Grab the 12 first columns where all the information is found
    interested_column <-
      significant_means[, 12 + x] #Grab the column of the interaction that you are interested with
    name_of_column <-
      extract_cell_interactions[x] #Extract the name of the interaction
    information[[name_of_column]] <-
      interested_column #Add to the information dataframe, the column of the interested interaction
    inf_df <-
      information[complete.cases(information[[name_of_column]]),] #Delete all the NA rows in the column of the interested interaction
    #Number of rows
    row_number <-
      nrow(inf_df) #If the number of rows is 0, then go next, because no interaction found
    if (row_number == 0) {
      next
    }
    #Iterate over all the interactions inside the cellxcell
    for (p in 1:length(row_number)) {
      y <- y + 1
      word_cell <-
        as.list(strsplit(extract_cell_interactions[x], ".", fixed = TRUE))[[1]]
      #give variables a name to add to the dataframe
      cell_type1 <- word_cell[1]
      cell_type2 <- word_cell[2]
      cell_phone[y, 1] <- cell_type1
      cell_phone[y, 2] <- cell_type2
      if (cell_type1 == "Naïve_T_cells") {
        cell_type1 <- "Naive_T_cells"

      }
      if (cell_type2 == "Naïve_T_cells") {
        cell_type2 <- "Naive_T_cells"

      }
      #subset_type
      subset1 <- anotation[anotation$cluster == cell_type1, ]
      subset1 <- subset1$subset
      subset2 <- anotation[anotation$cluster == cell_type2, ]
      subset2 <- subset2$subset
      cell_phone[y, 3] <- subset1
      cell_phone[y, 4] <- subset2
      #id_cp_interactions
      cell_phone[y, 5] <- inf_df[p, 1]

      #interacting_pairs
      cell_phone[y, 6] <- inf_df[p, 2]

      #partner_a
      cell_phone[y, 7] <- inf_df[p, 3]

      #parnet_b
      cell_phone[y, 8] <- inf_df[p, 4]

      #gene_a
      cell_phone[y, 9] <- inf_df[p, 5]

      #gene_b
      cell_phone[y, 10] <- inf_df[p, 6]

      #secreted
      cell_phone[y, 11] <- inf_df[p, 7]

      #receptor_a
      cell_phone[y, 12] <- inf_df[p, 8]

      #receptor_b
      cell_phone[y, 13] <- inf_df[p, 9]

      #anotation_strategy
      cell_phone[y, 14] <- inf_df[p, 10]

      #is_integrin
      cell_phone[y, 15] <- inf_df[p, 11]

      #variable
      cell_phone[y, 16] <- extract_cell_interactions[x]

      #mean
      cell_phone[y, 17] <- inf_df[p, extract_cell_interactions[x]]

      #pvalue
      pvalue_interaction <-
        pvalues[pvalues$id_cp_interaction == inf_df[p, 1], extract_cell_interactions[x]]
      cell_phone[y, 18] <- pvalue_interaction

      #log10p
      cell_phone[y, 19] <- log10(cell_phone[y, 18])
    }

  }

  return(cell_phone)
}

#Fuse cellphone output into one dataframe
fuse_all_cellphones <-
  function(cell_phone_group1,
           cell_phone_group2,
           cell_phone_group3,
           cell_phone_group4) {
    dataframe <- data.frame(matrix(ncol = 28, nrow = 1))
    colnames_df <- colnames(cell_phone_group4)
    colnames_df <-
      c(
        colnames_df[1:16],
        "mean1",
        "mean2",
        "mean3",
        "mean4",
        "pvalue1",
        "pvalue2",
        "pvalue3",
        "pvalue4",
        "log10p1",
        "log10p2",
        "log10p3",
        "log10p4"
      )
    colnames(dataframe) <- colnames_df
    #Add group 1
    for (x in 1:nrow(cell_phone_group1)) {
      dataframe[x, 1] <- cell_phone_group1[x, 1]
      dataframe[x, 2] <- cell_phone_group1[x, 2]
      dataframe[x, 3] <- cell_phone_group1[x, 3]
      dataframe[x, 4] <- cell_phone_group1[x, 4]
      dataframe[x, 5] <- cell_phone_group1[x, 5]
      dataframe[x, 6] <- cell_phone_group1[x, 6]
      dataframe[x, 7] <- cell_phone_group1[x, 7]
      dataframe[x, 8] <- cell_phone_group1[x, 8]
      dataframe[x, 9] <- cell_phone_group1[x, 9]
      dataframe[x, 10] <- cell_phone_group1[x, 10]
      dataframe[x, 11] <- cell_phone_group1[x, 11]
      dataframe[x, 12] <- cell_phone_group1[x, 12]
      dataframe[x, 13] <- cell_phone_group1[x, 13]
      dataframe[x, 14] <- cell_phone_group1[x, 14]
      dataframe[x, 15] <- cell_phone_group1[x, 15]
      dataframe[x, 16] <- cell_phone_group1[x, 16]
      dataframe[x, 17] <- cell_phone_group1[x, 17]
      dataframe[x, 21] <- cell_phone_group1[x, 18]
      dataframe[x, 25] <- cell_phone_group1[x, 19]
    }

    dataframe_names <-
      c("cell_phone_group2",
        "cell_phone_group3",
        "cell_phone_group4")

    #Add group 2
    for (name in dataframe_names) {
      message(name)
      p <- get(name)

      for (x in 1:nrow(p)) {
        message(name)
        message(x)
        coincidence1 <- p[x, 1]
        coincidence2 <- p[x, 2]
        coincidence3 <- p[x, 3]
        coincidence4 <- p[x, 4]
        coincidence5 <- p[x, 5]
        filtered_df <-
          subset(
            dataframe,
            dataframe$Cluster_1 == coincidence1 &
              dataframe$Cluster_2 == coincidence2 &
              dataframe$Subset_1 == coincidence3 &
              dataframe$Subset_2 == coincidence4 &
              dataframe$id_cp_interaction == coincidence5
          )
        number_rows <- as.numeric(nrow(filtered_df))
        if (number_rows != 0) {
          condition <-
            dataframe$Cluster_1 == coincidence1 &
            dataframe$Cluster_2 == coincidence2 &
            dataframe$Subset_1 == coincidence3 &
            dataframe$Subset_2 == coincidence4 &
            dataframe$id_cp_interaction == coincidence5
          if (name == "cell_phone_group2") {
            dataframe$mean2[condition] <- p[x, 17]
            dataframe$pvalue2[condition] <- p[x, 18]
            dataframe$log10p2[condition] <- p[x, 19]
          }
          if (name == "cell_phone_group3") {
            dataframe$mean3[condition] <- p[x, 17]
            dataframe$pvalue3[condition] <- p[x, 18]
            dataframe$log10p3[condition] <- p[x, 19]
          }
          if (name == "cell_phone_group4") {
            dataframe$mean4[condition] <- p[x, 17]
            dataframe$pvalue4[condition] <- p[x, 18]
            dataframe$log10p4[condition] <- p[x, 19]
          }

        } else{
          dataframe[nrow(dataframe) + 1, 1] <- p[x, 1]
          dataframe[nrow(dataframe), 2] <- p[x, 2]
          dataframe[nrow(dataframe), 3] <- p[x, 3]
          dataframe[nrow(dataframe), 4] <- p[x, 4]
          dataframe[nrow(dataframe), 5] <- p[x, 5]
          dataframe[nrow(dataframe), 6] <- p[x, 6]
          dataframe[nrow(dataframe), 7] <- p[x, 7]
          dataframe[nrow(dataframe), 8] <- p[x, 8]
          dataframe[nrow(dataframe), 9] <- p[x, 9]
          dataframe[nrow(dataframe), 10] <- p[x, 10]
          dataframe[nrow(dataframe), 11] <- p[x, 11]
          dataframe[nrow(dataframe), 12] <- p[x, 12]
          dataframe[nrow(dataframe), 13] <- p[x, 13]
          dataframe[nrow(dataframe), 14] <- p[x, 14]
          dataframe[nrow(dataframe), 15] <- p[x, 15]
          dataframe[nrow(dataframe), 16] <- p[x, 16]
          if (name == "cell_phone_group2") {
            dataframe[nrow(dataframe), 18] <- p[x, 17]
            dataframe[nrow(dataframe), 22] <- p[x, 18]
            dataframe[nrow(dataframe), 26] <- p[x, 19]
          }
          if (name == "cell_phone_group3") {
            dataframe[nrow(dataframe), 19] <- p[x, 17]
            dataframe[nrow(dataframe), 23] <- p[x, 18]
            dataframe[nrow(dataframe), 27] <- p[x, 19]
          }
          if (name == "cell_phone_group4") {
            dataframe[nrow(dataframe), 20] <- p[x, 17]
            dataframe[nrow(dataframe), 24] <- p[x, 18]
            dataframe[nrow(dataframe), 28] <- p[x, 19]
          }

        }
      }
    }
    return(dataframe)
  }


#
#Cellphone input preparation
#

#Create directories
dir.create('~/TOFA/plots/cellphone_output')
dir.create('~/TOFA/plots/cellphone_output/raw')

#Obtain genes
dummydata <-
  read10xCounts('/home/acorraliza/TOFA_data/output_cellranger/TOF-005-W0/')
geneinfo <- rowData(dummydata)

#Read the subset objects
epi <-
  readRDS('Analysis/data/00_annotation_process/00_anotadas/epi.RDS')
stroma <-
  readRDS('Analysis/data/00_annotation_process/00_anotadas/stroma.RDS')
plasmas <-
  readRDS('Analysis/data/00_annotation_process/00_anotadas/plasmas.RDS')
tcells <-
  readRDS('Analysis/data/00_annotation_process/00_anotadas/tcells.RDS')
cycling <-
  readRDS('Analysis/data/00_annotation_process/00_anotadas/cycling.RDS')
myeloids <-
  readRDS('Analysis/data/00_annotation_process/00_anotadas/myeloids.RDS')

#Modify certain cell names to avoid overlapping
myeloids$annotation_intermediate <-
  gsub("HS", "HS_myeloids", myeloids$annotation_intermediate)
tcells$annotation_intermediate <-
  gsub("HS", "HS_tcells", tcells$annotation_intermediate)
tcells$annotation_intermediate <-
  gsub("Naïve T cells",
       "Naive_T_cells",
       tcells$annotation_intermediate)
cycling$annotation_intermediate <-
  gsub("Cycling_B_cell",
       "Cycling_B_cell_cycling",
       cycling$annotation_intermediate)
plasmas$annotation_intermediate <-
  gsub("Cycling_B_cell",
       "Cycling_B_cell_plasmas",
       plasmas$annotation_intermediate)

#List of interactions we want to check
thelist <- list(
  #Myeloids
  c('myeloids', 'myeloids'),
  c('myeloids', 'plasmas'),
  c('myeloids', 'stroma'),
  c('myeloids', 'tcells'),
  c('myeloids', 'cycling'),
  c('myeloids', 'epi'),
  #Tcells
  c('tcells', 'tcells'),
  c('tcells', 'epi'),
  c('tcells', 'stroma'),
  c('tcells', 'plasmas'),
  c('tcells', 'cycling'),
  #Epi
  c('epi', 'epi'),
  c('epi', 'plasmas'),
  c('epi', 'stroma'),
  c('epi', 'cycling'),
  #Plasma
  c('plasmas', 'plasmas'),
  c('plasmas', 'stroma'),
  c('plasmas', 'cycling'),
  #Cycling
  c('cycling', 'cycling'),
  c('cycling', 'stroma'),
  #Stroma
  c('stroma', 'stroma')
)

#List of the objects
thedata <- list(tcells,
                plasmas,
                epi,
                stroma,
                myeloids,
                cycling)

names(thedata) <- c('tcells',
                    'plasmas',
                    'epi',
                    'stroma',
                    'myeloids',
                    'cycling')



#Obtain folders with data
for (i in 1:length(thelist)) {
  print(thelist[[i]])
  new <- merge(thedata[thelist[[i]][1]][[1]],
               thedata[thelist[[i]][2]][[1]])
  comparison <-
    paste(thelist[[i]][1], thelist[[i]][2], sep = '_and_')
  #Rename metadata
  new <- RenameCells(new, new.names = gsub(' ', '_', colnames(new)))
  new <- RenameCells(new, new.names = gsub('-', '_', colnames(new)))
  #Set directory for input cellphone
  setwd("~/TOFA/plots/cellphone_output/raw")
  dir.create(comparison)
  setwd(comparison)

  #Create cellphonedb_meta.txt
  print('normalizing_counts')
  count_raw <- GetAssayData(object = new, slot = "counts")
  count_norm <- apply(count_raw, 2, function(x)
    (x / sum(x)) * 10000)
  count_norm <- as.data.frame(count_norm)
  count_norm$Gene <-
    suppressWarnings(mapvalues(rownames(count_norm), geneinfo$Symbol, geneinfo$ID))
  count_norm <-
    count_norm[, c(ncol(count_norm), 1:(ncol(count_norm) - 1))]
  write.table(
    count_norm,
    'cellphonedb_count.txt',
    sep = '\t',
    quote = F,
    row.names = F
  )

  #Create cellphonedb_meta.txt
  print('generating meta file')
  meta_data <- cbind(rownames(new@meta.data),
                     new@meta.data[, 'annotation_intermediate', drop = F])
  colnames(meta_data) <- c('Cell',    'cell_type')
  meta_data$cell_type <- gsub(" ", "_", meta_data$cell_type)
  meta_data$cell_type <- gsub("-", "_", meta_data$cell_type)
  write.table(
    meta_data,
    'cellphonedb_meta.txt',
    sep = '\t',
    quote = F,
    row.names = F
  )


}

#Data frame of patients
sample_type <- data.frame(
  PATIENT = c(
    "TOF_005",
    "TOF_005",
    "TOF_005",
    "TOF_009",
    "TOF_009",
    "TOF_010",
    "TOF_011",
    "TOF_011",
    "TOF_012",
    "TOF_012",
    "TOF_012",
    "TOF_013",
    "TOF_013",
    "TOF_015",
    "TOF_015",
    "TOF_016",
    "TOF_016",
    "TOF_019",
    "TOF_019",
    "TOF_019",
    "TOF_022",
    "TOF_023",
    "TOF_023"
  ),
  TIMEPOINT = c(
    "W0",
    "W8",
    "W24",
    "W0",
    "w48",
    "W0",
    "W0",
    "W8",
    "W0",
    "W8",
    "W16",
    "W0",
    "W8",
    "W0",
    "W8",
    "W0",
    "W8",
    "W0",
    "W8",
    "W16",
    "W0",
    "W0",
    "W8"
  ),
  week_2 = c(
    "W0",
    "W8",
    "W16",
    "W0",
    "W16",
    "W0",
    "W0",
    "W8",
    "W0",
    "W8",
    "W16",
    "W0",
    "W8",
    "W0",
    "W8",
    "W0",
    "W8",
    "W0",
    "W8",
    "W16",
    "W0",
    "W0",
    "W8"
  ),
  SAMPLE = c(
    "TOF-005-W0",
    "TOF-005-W8",
    "TOF-005-W24",
    "TOF-009-W0",
    "TOF-009-W48",
    "TOF-010-W0",
    "TOF-011-W0",
    "TOF-011-W8",
    "TOF-012-W0",
    "TOF-012-W8",
    "TOF-012-W16",
    "TOF-013-W0",
    "TOF-013-W8",
    "TOF-015-W0",
    "TOF-015-W8",
    "TOF-016-W0",
    "TOF-016-W8",
    "TOF-019-W0",
    "TOF-019-W8",
    "TOF-019-W16",
    "TOF-022-W0",
    "TOF-023-W0",
    "TOF-023-W8"
  ),
  Response = c(
    "NR",
    "NR",
    "NR",
    "R",
    "R",
    "NR",
    "R",
    "R",
    "R",
    "R",
    "R",
    "NR",
    "NR",
    "R",
    "R",
    "NR",
    "NR",
    "NR",
    "NR",
    "NR",
    "NR",
    "R",
    "R"
  )
)

sample_type[] <- lapply(sample_type, function(x)
  gsub("-", "_", x))

#Create and separate for groups
sample_type$group <- NA
for (x in 1:nrow(sample_type)) {
  if (sample_type[x, 2] == "W0") {
    if (sample_type[x, 5] == "NR") {
      sample_type[x, 6] <- 2
    } else{
      sample_type[x, 6] <- 1
    }
  } else{
    if (sample_type[x, 5] == "NR") {
      sample_type[x, 6] <- 4
    } else{
      sample_type[x, 6] <- 3
    }
  }
}

#Input
dir.create('~/TOFA/plots/cellphone_output/input/')

#Obtain final input for cellphone run
for (x in 1:length(thelist)) {
  print(thelist[x])
  #Obtain cells and paths
  cell1 <- thelist[[x]][1]
  cell2 <- thelist[[x]][2]
  path1 <-
    paste(
      "~/TOFA/plots/cellphone_output/raw/",
      cell1,
      "_and_",
      cell2,
      "/cellphonedb_count.txt",
      sep = ""
    )
  path2 <-
    paste(
      "~/TOFA/plots/cellphone_output/raw/",
      cell1,
      "_and_",
      cell2,
      "/cellphonedb_meta.txt",
      sep = ""
    )


  #Extract files
  db_count <- read.delim(path1, header = TRUE)
  db_meta <- read.delim(path2, header = TRUE)

  #Create directories
  dir.create(
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group1",
      sep = ""
    )
  )
  dir.create(
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group2",
      sep = ""
    )
  )
  dir.create(
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group3",
      sep = ""
    )
  )
  dir.create(
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group4",
      sep = ""
    )
  )

  #Start separating files
  group1_sample <- sample_type[sample_type$group == 1, ]$SAMPLE
  group2_sample <- sample_type[sample_type$group == 2, ]$SAMPLE
  group3_sample <- sample_type[sample_type$group == 3, ]$SAMPLE
  group4_sample <- sample_type[sample_type$group == 4, ]$SAMPLE

  db_count_colnames <- colnames(db_count)
  #Dataframe group_1
  finding1 <-
    paste(
      group1_sample[1],
      "|",
      group1_sample[2],
      "|",
      group1_sample[3],
      "|",
      group1_sample[4],
      "|",
      group1_sample[5],
      "|Gene",
      sep = ""
    )
  db_count_colnames_group1 <-
    db_count_colnames[grepl(finding1, db_count_colnames)]
  selected_db_count1 <- db_count[db_count_colnames_group1]
  selected_db_meta1 <-
    db_meta[db_meta$Cell %in%  db_count_colnames_group1[2:length(db_count_colnames_group1)], ]
  write.table(
    selected_db_meta1,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group1/cellphonedb_meta.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count1,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group1/cellphonedb_count.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  #Dataframe group_2
  finding2 <-
    paste(
      group2_sample[1],
      "|",
      group2_sample[2],
      "|",
      group2_sample[3],
      "|",
      group2_sample[4],
      "|",
      group2_sample[5],
      "|",
      group2_sample[6],
      "|Gene",
      sep = ""
    )
  db_count_colnames_group2 <-
    db_count_colnames[grepl(finding2, db_count_colnames)]
  selected_db_count2 <- db_count[db_count_colnames_group2]
  selected_db_meta2 <-
    db_meta[db_meta$Cell %in%  db_count_colnames_group2[2:length(db_count_colnames_group2)], ]
  write.table(
    selected_db_meta2,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group2/cellphonedb_meta.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count2,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group2/cellphonedb_count.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  #Dataframe group_3
  finding3 <-
    paste(
      group3_sample[1],
      "|",
      group3_sample[2],
      "|",
      group3_sample[3],
      "|",
      group3_sample[4],
      "|",
      group3_sample[5],
      "|",
      group3_sample[6],
      "|Gene",
      sep = ""
    )
  db_count_colnames_group3 <-
    db_count_colnames[grepl(finding3, db_count_colnames)]
  selected_db_count3 <- db_count[db_count_colnames_group3]
  selected_db_meta3 <-
    db_meta[db_meta$Cell %in%  db_count_colnames_group3[2:length(db_count_colnames_group3)], ]
  write.table(
    selected_db_meta3,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group3/cellphonedb_meta.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count3,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group3/cellphonedb_count.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  #Dataframe group_4
  finding4 <-
    paste(
      group4_sample[1],
      "|",
      group4_sample[2],
      "|",
      group4_sample[3],
      "|",
      group4_sample[4],
      "|",
      group4_sample[5],
      "|",
      group4_sample[6],
      "|Gene",
      sep = ""
    )
  db_count_colnames_group4 <-
    db_count_colnames[grepl(finding4, db_count_colnames)]
  selected_db_count4 <- db_count[db_count_colnames_group4]
  selected_db_meta4 <-
    db_meta[db_meta$Cell %in%  db_count_colnames_group4[2:length(db_count_colnames_group4)], ]
  write.table(
    selected_db_meta4,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group4/cellphonedb_meta.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count4,
    paste(
      '~/TOFA/plots/cellphone_output/input/',
      cell1,
      "_",
      cell2,
      "_group4/cellphonedb_count.txt",
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

}

#
# Run bash script
#

system("bash cellphone.sh")

#
# Modify input
#

# Get a list of all directories within the parent directory
directory_groups <- "~/TOFA/Cell_phone/output"
directories <- list.dirs(directory_groups, recursive = FALSE)

#Create the 4 cellphones
cell_phone_group1 <- NA
cell_phone_group2 <- NA
cell_phone_group3 <- NA
cell_phone_group4 <- NA

# Iterate over the directories
for (dir in directories) {
  message(dir)
  message("Putting together. ")
  #Obtain directory for function
  stripped_dir <-
    sub("^/home/mmoro/TOFA/plots/cellphone_output/output/",
        "",
        dir)
  split_string <- strsplit(stripped_dir, "_")[[1]]
  first_element <- split_string[1]
  second_element <- split_string[2]
  directory <- c(first_element, second_element)
  #Obtain 4 excells for the different groups
  if (grepl("group1", dir)) {
    #Group1
    pvalues_dir <- paste(dir, "/out/pvalues.txt", sep = "")
    pvalues <- read.delim(pvalues_dir, header = TRUE)
    significant_means_dir <-
      paste(dir, "/out/significant_means.txt", sep = "")
    significant_means <-
      read.delim(significant_means_dir, header = TRUE)
    cb_group1 <-
      obtain_cell_phone(directory, pvalues, significant_means)
    cell_phone_group1 <- rbind(cell_phone_group1, cb_group1)
  }
  if (grepl("group2", dir)) {
    #Group2
    pvalues_dir <- paste(dir, "/out/pvalues.txt", sep = "")
    pvalues <- read.delim(pvalues_dir, header = TRUE)
    significant_means_dir <-
      paste(dir, "/out/significant_means.txt", sep = "")
    significant_means <-
      read.delim(significant_means_dir, header = TRUE)
    cb_group2 <-
      obtain_cell_phone(directory, pvalues, significant_means)
    cell_phone_group2 <- rbind(cell_phone_group2, cb_group2)
  }
  if (grepl("group3", dir)) {
    #Group3
    pvalues_dir <- paste(dir, "/out/pvalues.txt", sep = "")
    pvalues <- read.delim(pvalues_dir, header = TRUE)
    significant_means_dir <-
      paste(dir, "/out/significant_means.txt", sep = "")
    significant_means <-
      read.delim(significant_means_dir, header = TRUE)
    cb_group3 <-
      obtain_cell_phone(directory, pvalues, significant_means)
    cell_phone_group3 <- rbind(cell_phone_group3, cb_group3)
  }
  if (grepl("group4", dir)) {
    #Group4
    pvalues_dir <- paste(dir, "/out/pvalues.txt", sep = "")
    pvalues <- read.delim(pvalues_dir, header = TRUE)
    significant_means_dir <-
      paste(dir, "/out/significant_means.txt", sep = "")
    significant_means <-
      read.delim(significant_means_dir, header = TRUE)
    cb_group4 <-
      obtain_cell_phone(directory, pvalues, significant_means)
    cell_phone_group4 <- rbind(cell_phone_group4, cb_group4)
  }

}

dir.create('/home/mmoro/TOFA/plots/cellphone_output/results')

#Save the obtained dataframes

#Cellphone group 1 dataframe
df_filtered <- cell_phone_group1[order(cell_phone_group1$pvalue),]
cell_phone_group1 <-
  df_filtered[!duplicated(df_filtered[, c("Cluster_1",
                                          "Cluster_2",
                                          "id_cp_interaction",
                                          "Subset_1",
                                          "Subset_2")]),]
cell_phone_group1 <- cell_phone_group1[-nrow(cell_phone_group1), ]

file_path1 <-
  "~/TOFA/plots/cellphone_output/results/cell_phone1.xlsx"

# Export the dataframe to Excel
write.xlsx(cell_phone_group1, file_path1, rowNames = FALSE)

#Cellphone group 2 dataframe
df_filtered <- cell_phone_group2[order(cell_phone_group2$pvalue),]
cell_phone_group2 <-
  df_filtered[!duplicated(df_filtered[, c("Cluster_1",
                                          "Cluster_2",
                                          "id_cp_interaction",
                                          "Subset_1",
                                          "Subset_2")]),]
cell_phone_group2 <- cell_phone_group2[-nrow(cell_phone_group2), ]

file_path2 <-
  "~/TOFA/plots/cellphone_output/results/cell_phone2.xlsx"

# Export the dataframe to Excel
write.xlsx(cell_phone_group2, file_path2, rowNames = FALSE)

#Cellphone group 3 dataframe
df_filtered <- cell_phone_group3[order(cell_phone_group3$pvalue),]
cell_phone_group3 <-
  df_filtered[!duplicated(df_filtered[, c("Cluster_1",
                                          "Cluster_2",
                                          "id_cp_interaction",
                                          "Subset_1",
                                          "Subset_2")]),]
cell_phone_group3 <- cell_phone_group3[-nrow(cell_phone_group3), ]

file_path3 <-
  "~/TOFA/plots/cellphone_output/results/cell_phone3.xlsx"

# Export the dataframe to Excel
write.xlsx(cell_phone_group3, file_path3, rowNames = FALSE)

#Cellphone group 4 dataframe
df_filtered <- cell_phone_group4[order(cell_phone_group4$pvalue),]
cell_phone_group4 <-
  df_filtered[!duplicated(df_filtered[, c("Cluster_1",
                                          "Cluster_2",
                                          "id_cp_interaction",
                                          "Subset_1",
                                          "Subset_2")]),]
cell_phone_group4 <- cell_phone_group4[-nrow(cell_phone_group4), ]

file_path4 <-
  "~/TOFA/plots/cellphone_output/results/cell_phone4.xlsx"

# Export the dataframe to Excel
write.xlsx(cell_phone_group4, file_path4, rowNames = FALSE)

#Calling function
final_cellphone <-
  fuse_all_cellphones(cell_phone_group1,
                      cell_phone_group2,
                      cell_phone_group3,
                      cell_phone_group4)

#Path
file_path_final <-
  "~/TOFA/plots/cellphone_output/results/final_cell_phone.xlsx"

# Export the dataframe to Excel
write.xlsx(final_cellphone, file_path_final, rowNames = FALSE)
