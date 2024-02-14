library(dplyr)
library(ggplot2)
library(readr)


## blood RNA --------------------------------------------------------------------
blood_rna <- read_excel("Analysis/data/Blood RNA supplementary figure 1 B 1.xlsx",
                        na = "X")

blood_rna_responders <- blood_rna %>%
  dplyr::filter(Response == "R") %>%
  as.data.frame()

resultats <- NULL
for(gene in c("SOCS1",
              "SOCS3",
              "IRF1",
              "MS4A1",
              "CD3E",
              "PROK2",
              "FCG3RA")) {
  print(shapiro.test(blood_rna_responders %>% dplyr::pull(gene)))
  kw <- kruskal.test(blood_rna_responders[,gene] ~ blood_rna_responders$`Time point`)
  res <- data.frame(
    response = 'R',
    gene = gene,
    comp = 'kruskall',
    p.val = kw$p.value)

  resultats <- rbind(resultats, res)

  plot <- ggpubr::ggboxplot(blood_rna_responders, x = 'Time point', y = gene, add = 'jitter') +
    ggplot2::labs(title = paste('Responders', gene)) +
    ggpubr::stat_compare_means()
  print(plot)
}

blood_rna_nonresponders <- blood_rna %>%
  dplyr::filter(Response == "NR") %>%
  as.data.frame()

# resultats <- NULL
for(gene in c("SOCS1",
              "SOCS3",
              "IRF1",
              "MS4A1",
              "CD3E",
              "PROK2",
              "FCG3RA")) {
  print(shapiro.test(blood_rna_nonresponders %>% dplyr::pull(gene)))
  kw <- kruskal.test(blood_rna_nonresponders[,gene] ~ blood_rna_nonresponders$`Time point`)
  res <- data.frame(
    response = 'NR',
    gene = gene,
    comp = 'kruskall',
    p.val = kw$p.value)

  resultats <- rbind(resultats, res)

  plot <- ggpubr::ggboxplot(blood_rna_nonresponders, x = 'Time point', y = gene, add = 'jitter') +
    ggplot2::labs(title = paste('No Responders', gene)) +
    ggpubr::stat_compare_means()
  print(plot)
}

resultats$padj_fdr <- p.adjust(resultats$p.val)

## blood populations --------------------------------------------------------------------

blood_cells <- read_excel("Analysis/data/Analíticas supplemmentary figure 1 A 1.xlsx",
                          na = "X")

resultats <- NULL
cells_pop <- colnames(blood_cells)[4:ncol(blood_cells)]

blood_cells_responders <- blood_cells %>%
  dplyr::filter(RESPONSE == "R") %>%
  as.data.frame()

for(cells in cells_pop) {
  print(shapiro.test(blood_cells_responders %>% dplyr::pull(cells)))
  kw <- kruskal.test(blood_cells_responders[,cells] ~ blood_cells_responders$`WEEK 2`)
  res <- data.frame(
    response = 'R',
    gene = cells,
    comp = 'kruskall',
    p.val = kw$p.value)

  resultats <- rbind(resultats, res)

  plot <- ggpubr::ggboxplot(blood_cells_responders, x = 'WEEK 2', y = formulaic::add.backtick(cells), add = 'jitter') +
    ggplot2::labs(title = paste('Responders', cells)) +
    ggpubr::stat_compare_means()
  print(plot)
}

blood_cells_nonresponders <- blood_cells %>%
  dplyr::filter(RESPONSE == "NR") %>%
  as.data.frame()

for(cells in cells_pop) {
  print(shapiro.test(blood_cells_nonresponders %>% dplyr::pull(cells)))
  kw <- kruskal.test(blood_cells_nonresponders[,cells] ~ blood_cells_nonresponders$`WEEK 2`)
  res <- data.frame(
    response = 'NR',
    gene = cells,
    comp = 'kruskall',
    p.val = kw$p.value)

  resultats <- rbind(resultats, res)

  plot <- ggpubr::ggboxplot(blood_cells_nonresponders, x = 'WEEK 2', y = formulaic::add.backtick(cells), add = 'jitter' ) +
    ggplot2::labs(title = paste('No Responders', cells)) +
    ggpubr::stat_compare_means()
  print(plot)
}

resultats$padj_fdr <- p.adjust(resultats$p.val)
