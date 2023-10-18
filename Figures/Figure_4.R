


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

save_sizes(plot =fig4c_S1 , filename = 'fig4c_S1', device = 'jpeg')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1', device = 'tiff')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1', device = 'svg')
save_sizes(plot = fig4c_S1, filename = 'fig4c_S1', device = 'pdf')

# M2 non-responders
cluster <- "M2"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("MMP9", "INHBA", "CD300E", "SPP1", "IDO1",
                    "MMP12", "CLEC5A", "PLCG2", "ITGAX", "PLAUR", "IL1RN",
                    "CXCL3", "VIM", "IL7R", "TMSB4X", "TPT1",
                    "RPL9", "MS4A6A", "SELENOP", "FUCA1", "CYBA")

fig3b_M2NR <- volcano_plot(cluster, comp, filtered_genes)
print(fig3b_M2NR)

save_sizes(plot =fig3b_M2NR , filename = 'fig3b_M2NR', device = 'jpeg')
save_sizes(plot = fig3b_M2NR, filename = 'fig3b_M2NR', device = 'tiff')
save_sizes(plot = fig3b_M2NR, filename = 'fig3b_M2NR', device = 'svg')
save_sizes(plot = fig3b_M2NR, filename = 'fig3b_M2NR', device = 'pdf')



# M1 non-responders
cluster <- "M1"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("PLCG2", "MIF", "FTH1", "ENO1", "TPI1",
                    "MT-ATP8", "P4HA1", "GAPDH", "SPP1", "MT-CO3", "MT-CYB",
                    "MT-CO2", "MT-ATP6", "MT-CO1", "MT-ND1", "MT-ND3",
                    "IFI6", "OAS1", "S100A9", "MX1")

fig3b_M1NR <- volcano_plot(cluster, comp, filtered_genes)
print(fig3b_M1NR)

save_sizes(plot =fig3b_M1NR , filename = 'fig3b_M1NR', device = 'jpeg')
save_sizes(plot = fig3b_M1NR, filename = 'fig3b_M1NR', device = 'tiff')
save_sizes(plot = fig3b_M1NR, filename = 'fig3b_M1NR', device = 'svg')
save_sizes(plot = fig3b_M1NR, filename = 'fig3b_M1NR', device = 'pdf')


# M0 Responders

cluster <- "M0"
comp <- "w0R_vs_POSTR"
filtered_genes <- c("S100A9", "CD14", "MMP12", "CXCL9", "IFI27",
                    "CCL18", "C1QA", "C1QB", "STAT1", "ISG15", "FCGR2A",
                    "FABP1", "LYZ", "IGF1", "RLPSSS", "MARCKS",
                    "HLA-DQA2", "RPL39", "RPS3A", "RPS28", "RPL34", "TMSB4X")

fig3b_M0R <- volcano_plot(cluster, comp, filtered_genes)
print(fig3b_M0R)

save_sizes(plot =fig3b_M0R , filename = 'fig3b_M0R', device = 'jpeg')
save_sizes(plot = fig3b_M0R, filename = 'fig3b_M0R', device = 'tiff')
save_sizes(plot = fig3b_M0R, filename = 'fig3b_M0R', device = 'svg')
save_sizes(plot = fig3b_M0R, filename = 'fig3b_M0R', device = 'pdf')

#M0 non-responders

cluster <- "M0"
comp <- "w0NR_vs_POSTNR"
filtered_genes <- c("MMP9", "MALAT1", "NEAT1", "C15orf48", "FABP5",
                    "S100A8", "IFITM2", "CYBA", "MTRNR2L8", "MTRNR2L12", "RPS2")

fig3b_M0NR <- volcano_plot(cluster, comp, filtered_genes)
print(fig3b_M0NR)

save_sizes(plot =fig3b_M0NR , filename = 'fig3b_M0NR', device = 'jpeg')
save_sizes(plot = fig3b_M0NR, filename = 'fig3b_M0NR', device = 'tiff')
save_sizes(plot = fig3b_M0NR, filename = 'fig3b_M0NR', device = 'svg')
save_sizes(plot = fig3b_M0NR, filename = 'fig3b_M0NR', device = 'pdf')
