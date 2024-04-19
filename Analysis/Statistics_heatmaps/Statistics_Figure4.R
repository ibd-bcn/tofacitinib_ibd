library(readxl)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(dplyr)
options(bitmapType='cairo')

#Macrophages DMSO---------------------------------------------------------------
#Read Database
DMSO <- read_excel("Figures/extra_data/230916 Base de datos macrofagos Madrid FC 3.xlsx",
                   sheet = "FC Stimuli-Control")

DMSO <- DMSO %>%
  mutate(across(where(is.numeric), log2))

colnames(DMSO)[2] <- "Sample"
colnames(DMSO)[3] <- "Condition"
DMSO <- DMSO[DMSO$Condition != "M-DMSO-IL4", ]
DMSO <- DMSO[,2:ncol(DMSO)]

numeric_cols <- sapply(DMSO, is.numeric)
#LPS
LPS <- subset(DMSO, Condition == "M-DMSO-LPS")
LPS <- data.frame(LPS, row.names = NULL)
#TNFa
TNFa <- subset(DMSO, Condition == "M-DMSO-TNF?")
TNFa <- data.frame(TNFa, row.names = NULL)
#IFNg
IFNG <- subset(DMSO, Condition == "M-DMSO-IFN?")
IFNG <- data.frame(IFNG, row.names = NULL)

#Statistical analysis
results <- IFNG %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFNG[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFNG[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNFa %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNFa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNFa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

#Final matrix
heatmap_tt <- rbind(LPS, IFNG, TNFa)

#Reorganize columns
heatmap_tt <- heatmap_tt %>% relocate(Sample)

#Save
openxlsx::write.xlsx(heatmap_tt, "Figures/extra_data/Heatmap/macros_dmso_stats.xlsx", rowNames=T)


# Macrophages TOFA--------------------------------------------------------------
TOFA <- read_excel("230916 Base de datos macrofagos Madrid FC 3.xlsx",
                   sheet = "FC Tofa stimuli-stimuli")

TOFA <- TOFA %>%
  mutate(across(where(is.numeric), log2))

colnames(TOFA)[2] <- "Sample"
colnames(TOFA)[3] <- "Condition"
TOFA <- TOFA[TOFA$Condition != "M-TOFA-IL4", ]

TOFA <- TOFA[,2:ncol(TOFA)]
numeric_cols <- sapply(TOFA, is.numeric)

#LPS
LPS <- subset(TOFA, Condition == "M-TOFA-LPS")
LPS <- data.frame(LPS, row.names = NULL)

#TNFa
TNFa <- subset(TOFA, Condition == "M-TOFA-TNF?")
TNFa <- data.frame(TNFa, row.names = NULL)

#IFNg
IFNG <- subset(TOFA, Condition == "M-TOFA-IFN?")
IFNG <- data.frame(IFNG, row.names = NULL)

#Statistics
results <- IFNG %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFNG[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFNG[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNFa %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNFa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNFa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

#Final Matrix
heatmap_tt <- rbind(LPS, IFNG, TNFa)

#Reorganize columns
heatmap_tt <- heatmap_tt %>% relocate(Sample)

#Save
openxlsx::write.xlsx(heatmap_tt, "Figures/extra_data/macros_tofa_stats.xlsx", rowNames=T)

#DMSO Fibros--------------------------------------------------------------------
DMSO <- read_excel("Figures/extra_data/Base de datos fibroblasts qPCR 2.xlsx",
                   sheet = "FC DMSO")


DMSO <- DMSO %>%
  mutate(across(where(is.numeric), log2))

colnames(DMSO)[2] <- "Sample"
colnames(DMSO)[1] <- "Condition"

DMSO <- DMSO[DMSO$Sample != "FIB1",]


#LPS
LPS <- subset(DMSO, Condition == "LPS")
LPS <- data.frame(LPS, row.names = NULL)
#TNFa
TNFa <- subset(DMSO, Condition == "TNF")
TNFa <- data.frame(TNFa, row.names = NULL)
#IFNg
IFNG <- subset(DMSO, Condition == "IFN")
IFNG <- data.frame(IFNG, row.names = NULL)

#Statistics

results <- IFNG %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFNG[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFNG[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNFa %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNFa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNFa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

#Final Matrix
heatmap_tt <- rbind(LPS, IFNG, TNFa)
#Reorganize columns
heatmap_tt <- heatmap_tt %>% relocate(Sample)
#Save
openxlsx::write.xlsx(heatmap_tt, "Figures/extra_data//fibros_dmso_stats.xlsx", rowNames=T)


#TOFA Fibros--------------------------------------------------------------------
TOFA <- read_excel("Figures/extra_data/Base de datos fibroblasts qPCR 2.xlsx",
                   sheet = "FC TOFA")

TOFA <- TOFA %>%
  mutate(across(where(is.numeric), log2))
colnames(TOFA)[2] <- "Sample"
colnames(TOFA)[1] <- "Condition"
TOFA <- TOFA[TOFA$Sample != "FIB1",]

#LPS
LPS <- subset(TOFA, Condition == "LPS+TOFA")
LPS <- data.frame(LPS, row.names = NULL)
#TNFa
TNFa <- subset(TOFA, Condition == "TNF+TOFA")
TNFa <- data.frame(TNFa, row.names = NULL)
#IFNg
IFNG <- subset(TOFA, Condition == "IFN+TOFA")
IFNG <- data.frame(IFNG, row.names = NULL)

#Statistics
results <- IFNG %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  IFNG[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  IFNG[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- LPS %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  LPS[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  LPS[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

results <- TNFa %>%
  select(-Condition, -Sample) %>%
  summarise(across(everything(), ~list(t.test(.x, mu = 0)$p.value))) %>%
  mutate(across(everything(), ~unlist(.x)))

all_p_values <- unlist(results)
adjusted_p_values <- p.adjust(all_p_values, method = "fdr")

for (gene in names(results)) {
  TNFa[[paste(gene, "p_value", sep = "_")]] <- results[[gene]]
  TNFa[[paste(gene, "p_adjusted", sep = "_")]] <- adjusted_p_values[names(adjusted_p_values) == gene]
}

#Final Matrix
heatmap_tt <- rbind(LPS, IFNG, TNFa)
#Reorganize columns
heatmap_tt <- heatmap_tt %>% relocate(Sample)
#Save
openxlsx::write.xlsx(heatmap_tt, "Figures/extra_data/fibros_tofa_stats.xlsx", rowNames=T)
