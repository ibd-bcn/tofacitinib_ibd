


# 4C Volcano plots stromal cells

de_data <- readRDS('/home/acorraliza/TOFA_data/20220222_TOFAS_23/01_DE/REPASO/new_complete.RDS')

# Volcano plots  myeloid cells

# S1 responders
cluster <- "S1"
comp <- "w0R_vs_POSTR"
filtered_genes <- c("ABCA8", "ADAMDEC1", "FOS", "GSN", "FN1",
                    "SELENOP", "EGR1", "ADIRF", "COL15A1", "AC007952.4", "CHI3L1",
                    "IL13RA2", "CXCL6", "OSMR", "IFITM3", "INHBA",
                    "GBP1", "NRG1", "WNT5A", "MMP3", "CXCL1", "CCL19", "FTH1", "F3", "IGFBP7", "CXCL14", "PDLIM4")

fig4c_S1 <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_S1)

save_sizes(plot =fig4c_S1 , filename = 'fig4c_S1R', device = 'jpeg')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1R', device = 'tiff')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1R', device = 'svg')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1R', device = 'pdf')

# S1 non-responders
cluster <- "S1"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("CHI3L1", "POSTN", "FTH1", "RPS4Y1", "PTMA",
                    "PLCG2", "TIMP1", "PPIB", "RPL8", "RARRES2", "TPM2",
                    "MT-CO3", "MT-ND1", "MT-ND3", "RPL23A", "MTRNR2L8",
                    "IFI27", "CD63", "XIST", "RPL21", "RPL9")

fig4c_S1NR <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_S1NR)

save_sizes(plot =fig4c_S1NR , filename = 'fig4c_S1NR', device = 'jpeg')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'tiff')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'svg')
save_sizes(plot = fig4c_S1NR, filename = 'fig4c_S1NR', device = 'pdf')





# S2 Responders

cluster <- "S2"
comp <- "w0R_vs_POSTR"
filtered_genes <- c("CHI3L1", "OSMR", "CXCL1", "MX1", "C1R",
                    "PDLIM4", "CST3", "IL15RA", "VIM", "COL6A2", "TIMP1",
                    "IGFBP7", "IFI6", "CFD", "HMGB1", "MALAT1",
                    "EGR1", "RPL39", "FOS", "TXNIP", "H3F3B", "FOSB", "IGFBP3", "RPL34")

fig4c_S2R <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_S2R)

save_sizes(plot =fig4c_S2R , filename = 'fig4c_S2R', device = 'jpeg')
save_sizes(plot = fig4c_S2R, filename = 'fig4c_S2R', device = 'tiff')
save_sizes(plot = fig4c_S2R, filename = 'fig4c_S2R', device = 'svg')
save_sizes(plot = fig4c_S2R, filename = 'fig4c_S2R', device = 'pdf')

#S2 non-responders

cluster <- "S2"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("PLCG2", "FTH1", "TMEM176A", "RPS8", "PFN1",
                    "RPL11", "VIM", "CXCL14", "IFITM3", "TPI1", "GBP3",
                    "GBP1", "MT-ND3", "MT-CO3", "MT-ND1", "MT-ATP6",
                    "MT-CO2", "MT-CYB", "MTRNR2L8", "MT-CO1", "MT-ND4")

fig4c_S2NR <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_S2NR)

save_sizes(plot =fig4c_S2NR , filename = 'fig4c_S2NR', device = 'jpeg')
save_sizes(plot = fig4c_S2NR, filename = 'fig4c_S2NR', device = 'tiff')
save_sizes(plot = fig4c_S2NR, filename = 'fig4c_S2NR', device = 'svg')
save_sizes(plot = fig4c_S2NR, filename = 'fig4c_S2NR', device = 'pdf')

# S3 Responders

cluster <- "S3"
comp <- "w0R_vs_POSTR"
filtered_genes <- c("ADAMDEC1", "ABCA8", "CXCL12", "TXNIP", "SELENOP",
                    "IGFBP3", "OGN", "FXYD1", "RPL30", "RPS28", "RPS21",
                    "CHI3L1", "CCL19", "CXCL12", "CXCL6", "CXCL2",
                    "CXCL14", "CXCL8", "FAP", "CLU", "TIMP1", "CXCL1", "CHRDL2", "RARRES1", "FTH1", "IGFBP4")

fig4c_S3R <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_S3R)

save_sizes(plot =fig4c_S3R , filename = 'fig4c_S3R', device = 'jpeg')
save_sizes(plot = fig4c_S3R, filename = 'fig4c_S3R', device = 'tiff')
save_sizes(plot = fig4c_S3R, filename = 'fig4c_S3R', device = 'svg')
save_sizes(plot = fig4c_S3R, filename = 'fig4c_S3R', device = 'pdf')

#S3 non-responders

cluster <- "S3"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("CXCL14", "PLCG2", "RPS4Y1", "PTMA", "TMEM176A",
                    "SERF2", "RARRES2", "RPL8", "RPL12", "TMEM14C", "LGALS3",
                    "CD99", "SOCS3", "MT-CO3", "MT-CO2", "MT-ATP6",
                    "MT-CO2", "MT-CYB", "MTRNR2L8", "MT-CO1", "MT-CO1", "TPT1", "MT-ATP6", "MT-CYB", "RPL9",
                    "RPL21", "CD63", "EMILIN1", "RPS27")

fig4c_S3NR <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_S3NR)

save_sizes(plot =fig4c_S3NR , filename = 'fig4c_S3NR', device = 'jpeg')
save_sizes(plot = fig4c_S3NR, filename = 'fig4c_S3NR', device = 'tiff')
save_sizes(plot = fig4c_S3NR, filename = 'fig4c_S3NR', device = 'svg')
save_sizes(plot = fig4c_S3NR, filename = 'fig4c_S3NR', device = 'pdf')


# Inflammatory fibroblasts responders

cluster <- "Inflammatory_fibroblasts"
comp <- "w0R_vs_POSTR"
filtered_genes <- c("CHI3L1", "PDPN", "PLAU", "ACTG1", "RRBP1",
                    "IGFBP4", "BGN", "ADM", "COL6A2", "FBLN2", "HLA-C",
                    "TMEM165", "YIPF2")

fig4c_IFR <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_IFR)

save_sizes(plot =fig4c_IFR , filename = 'fig4c_IFR', device = 'jpeg')
save_sizes(plot = fig4c_IFR, filename = 'fig4c_IFR', device = 'tiff')
save_sizes(plot = fig4c_IFR, filename = 'fig4c_IFR', device = 'svg')
save_sizes(plot = fig4c_IFR, filename = 'fig4c_IFR', device = 'pdf')

#Inflammatory_fibroblasts non-responders

cluster <- "Inflammatory_fibroblasts"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("PLCG2", "RPS26", "TPI1", "RPL19", "RPL11",
                    "FTH1", "HLA-B", "MT2A", "S100A11", "FAU", "TNFRSF12A",
                    "IL32", "TIMP1", "MT-ND1", "MT-CO3", "RPL10A",
                    "TPT1", "MT-ND4", "MT-ATP6", "RPS16", "SOCS3", "CXCL13", "IER2")
fig4c_IFNR <- volcano_plot(cluster, comp, filtered_genes)
print(fig4c_IFNR)

save_sizes(plot =fig4c_IFNR , filename = 'fig4c_IFNR', device = 'jpeg')
save_sizes(plot = fig4c_IFNR, filename = 'fig4c_IFNR', device = 'tiff')
save_sizes(plot = fig4c_IFNR, filename = 'fig4c_IFNR', device = 'svg')
save_sizes(plot = fig4c_IFNR, filename = 'fig4c_IFNR', device = 'pdf')
