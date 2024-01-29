#
#Libraries----------------------------------------------------------------------
#

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

#
## Cellphone input generation---------------------------------------------------
#

#Cell ranger output
dummydata <-
  read10xCounts('output_cellranger/TOF-005-W0/')
geneinfo <- rowData(dummydata)

#Annotated subsets
epi <-
  readRDS(
    'epi.RDS'
  )
stroma <-
  readRDS(
    'stroma.RDS'
  )
plasmas <-
  readRDS(
    'plasmas.RDS'
  )
tcells <-
  readRDS(
    'tcells.RDS'
  )
cycling <-
  readRDS(
    'cycling.RDS'
  )
myeloids <-
  readRDS(
    'myeloids.RDS'
  )

#Modify names to avoid overlapping
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

#All combination of subsets
thelist <- list(
  #Myeloids
  c('myeloids', 'myeloids'),
  c('myeloids', 'plasmas'),
  c('myeloids', 'stroma'),
  c('myeloids', 'tcells'),
  c('cycling', 'myeloids'),
  c('myeloids', 'epi'),
  #Tcells
  c('tcells', 'tcells'),
  c('tcells', 'epi'),
  c('tcells', 'stroma'),
  c('tcells', 'plasmas'),
  c('cycling', 'tcells'),
  #Epi
  c('epi', 'epi'),
  c('epi', 'plasmas'),
  c('epi', 'stroma'),
  c('epi', 'cycling'),
  #Plasma
  c('plasmas', 'plasmas'),
  c('stroma', 'plasmas'),
  c('cycling', 'plasmas'),
  #Cycling
  c('cycling', 'cycling'),
  c('cycling', 'stroma'),
  #Stroma
  c('stroma', 'stroma')
)

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

#Separate samples
group1_sample <- sample_type[sample_type$group == 1,]$SAMPLE
group2_sample <- sample_type[sample_type$group == 2,]$SAMPLE
group3_sample <- sample_type[sample_type$group == 3,]$SAMPLE
group4_sample <- sample_type[sample_type$group == 4,]$SAMPLE


#List of seurat objects
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

#Create directory
dir.create('output/')
setwd('output/')


#Obtain input
for (i in 1:length(thelist)) {
  print(thelist[[i]])
  new <- merge(thedata[thelist[[i]][1]][[1]],
               thedata[thelist[[i]][2]][[1]])
  comparison <-
    paste(thelist[[i]][1], thelist[[i]][2], sep = '_and_')

  #Rename metadata
  new <- RenameCells(new, new.names = gsub(' ', '_', colnames(new)))
  new <- RenameCells(new, new.names = gsub('-', '_', colnames(new)))

  #Cells we are working with
  cell1 <- thelist[[i]][1]
  cell2 <- thelist[[i]][2]

  #Create directories for comparison
  dir.create(paste(cell1, "_", cell2, "_group1", sep = ""))
  dir.create(paste(cell1, "_", cell2, "_group2", sep = ""))
  dir.create(paste(cell1, "_", cell2, "_group3", sep = ""))
  dir.create(paste(cell1, "_", cell2, "_group4", sep = ""))

  #Create cellphonedb_meta.txt
  print('normalizing_counts')
  count_raw <- new@assays$RNA@counts
  count_norm <- apply(count_raw, 2, function(x)
    (x / sum(x)) * 10000)
  count_norm <- as.data.frame(count_norm)
  count_norm$Gene <-
    suppressWarnings(mapvalues(
      rownames(count_norm),
      geneinfo$Symbol,
      geneinfo$ID,
      warn_missing = F
    ))
  count_norm <-
    count_norm[, c(ncol(count_norm), 1:(ncol(count_norm) - 1))]

  #Create cellphonedb_meta.txt
  print('generating meta')
  meta_data <- cbind(rownames(new@meta.data),
                     new@meta.data[, 'annotation_intermediate', drop = F])
  colnames(meta_data) <- c('Cell',    'cell_type')
  meta_data$cell_type <- gsub(" ", "_", meta_data$cell_type)
  meta_data$cell_type <- gsub("-", "_", meta_data$cell_type)

  #Save per groups
  db_count_colnames <- colnames(count_norm)
  #input group_1
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
  selected_db_count1 <- count_norm[db_count_colnames_group1]
  selected_db_meta1 <-
    meta_data[meta_data$Cell %in%  db_count_colnames_group1[2:length(db_count_colnames_group1)],]
  write.table(
    selected_db_meta1,
    paste(cell1, "_", cell2, "_group1/cellphonedb_meta.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count1,
    paste(cell1, "_", cell2, "_group1/cellphonedb_count.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  #input group_2
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
  selected_db_count2 <- count_norm[db_count_colnames_group2]
  selected_db_meta2 <-
    meta_data[meta_data$Cell %in%  db_count_colnames_group2[2:length(db_count_colnames_group2)],]
  write.table(
    selected_db_meta2,
    paste(cell1, "_", cell2, "_group2/cellphonedb_meta.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count2,
    paste(cell1, "_", cell2, "_group2/cellphonedb_count.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  #input group_3
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
  selected_db_count3 <- count_norm[db_count_colnames_group3]
  selected_db_meta3 <-
    meta_data[meta_data$Cell %in%  db_count_colnames_group3[2:length(db_count_colnames_group3)],]
  write.table(
    selected_db_meta3,
    paste(cell1, "_", cell2, "_group3/cellphonedb_meta.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count3,
    paste(cell1, "_", cell2, "_group3/cellphonedb_count.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  #input group_4
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
  selected_db_count4 <- count_norm[db_count_colnames_group4]
  selected_db_meta4 <-
    meta_data[meta_data$Cell %in%  db_count_colnames_group4[2:length(db_count_colnames_group4)],]
  write.table(
    selected_db_meta4,
    paste(cell1, "_", cell2, "_group4/cellphonedb_meta.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    selected_db_count4,
    paste(cell1, "_", cell2, "_group4/cellphonedb_count.txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


}

#
## Cellphone run----------------------------------------------------------------
#

#Run bash script
system("bash run_cellphone.sh")

