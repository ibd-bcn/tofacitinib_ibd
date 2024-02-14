data_by_patient <-data.frame(
          stringsAsFactors = FALSE,
                                     Patient = c("227",
                                                 "7903","8507","9247","9691",
                                                 "11174","11831","12135",
                                                 "12378","12479","13087","13134",
                                                 "13384","13566","TOF_001",
                                                 "TOF_002","TOF_003","TOF_005",
                                                 "TOF_006","TOF_007","TOF_009",
                                                 "TOF_010","TOF_011",
                                                 "TOF_012","TOF_013","TOF_015",
                                                 "TOF_016","TOF_019","TOF_022",
                                                 "TOF_023","TOF_024"),
                                      Cohort = c("LV","LV",
                                                 "LV","LV","LV","LV","LV",
                                                 "LV","LV","LV","LV","LV",
                                                 "LV","LV","BCN","BCN","BCN",
                                                 "BCN","BCN","BCN","BCN",
                                                 "BCN","BCN","BCN","BCN",
                                                 "BCN","BCN","BCN","BCN","BCN",
                                                 "BCN"),
                               Milagro.score = c("R","R",
                                                 "NR","NR","R","NR","R",
                                                 "NR","R","R","NR","R","NR",
                                                 "NR","NR","NR","R","NR",
                                                 "NR","R","R","NR","R","R",
                                                 "NR","R","NR","NR","NR","R",
                                                 "NR"),
                                         Sex = c("F","M",
                                                 "M","F","F","M","M","M",
                                                 "F","M","M","M","F","M",
                                                 "F","F","F","M","F","M","M",
                                                 "F","M","M","M","F","F",
                                                 "M","M","F","F"),
                                         Age = c(61.21013005,43.01711157,62.83915127,
                                                 52.96098563,42.38193018,
                                                 45.7412731,44.36413415,18.11635866,
                                                 28.05475702,67.63312799,
                                                 73.62902122,45.42094456,
                                                 48.47091034,48.56125941,31,44,43,54,
                                                 23,47,33,49,40,46,32,32,
                                                 48,23,41,22,55),
                            Disease.duration = c(17.44832307,16.92265572,10.02874743,
                                                 12.24914442,7.244353183,
                                                 27.56741958,3.789185489,3.403148528,
                                                 16.38603696,25.66187543,
                                                 12.77481177,3.496235455,
                                                 1.979466119,2.412046543,9,5,4,6,2,
                                                 23,8,14,4,2,5,8,28,7,
                                                 14,10,20),
          UC.extension = c("Left-side colitis",
                           "Left-side colitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Proctitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Pancolitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Left-side colitis",
                           "Pancolitis",
                           "Proctitis",
                           "Pancolitis",
                           "Pancolitis",
                           "Proctitis",
                           "Left-side colitis",
                           "Proctitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Pancolitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Pancolitis",
                           "Proctitis",
                           "Pancolitis",
                           "Left-side colitis",
                           "Left-side colitis",
                           "Proctitis"),
                    Total.biologicals.before = c(3L,2L,2L,
                                                 2L,2L,3L,2L,3L,2L,3L,
                                                 3L,3L,2L,3L,3L,2L,2L,2L,
                                                 1L,2L,0L,2L,4L,2L,1L,2L,
                                                 1L,2L,3L,2L,1L)
                  )

data_by_sample <- data.frame(
                     stringsAsFactors = FALSE,
                          check.names = FALSE,
                                                Patient = c("227","227",
                                                            "7903","7903","8507",
                                                            "8507","9247",
                                                            "9247","9691","9691",
                                                            "11174","11174",
                                                            "11831","11831",
                                                            "12135","12135",
                                                            "12378","12378","12479",
                                                            "12479","13087",
                                                            "13087","13134",
                                                            "13134","13384",
                                                            "13384","13566",
                                                            "13566","TOF_001",
                                                            "TOF_001","TOF_002",
                                                            "TOF_002","TOF_003",
                                                            "TOF_003","TOF_005",
                                                            "TOF_005","TOF_006",
                                                            "TOF_006",
                                                            "TOF_007","TOF_007",
                                                            "TOF_009","TOF_009",
                                                            "TOF_010","TOF_010",
                                                            "TOF_011","TOF_011",
                                                            "TOF_012",
                                                            "TOF_012","TOF_013",
                                                            "TOF_013","TOF_015",
                                                            "TOF_015","TOF_016",
                                                            "TOF_016","TOF_019",
                                                            "TOF_019","TOF_022",
                                                            "TOF_022",
                                                            "TOF_023","TOF_023",
                                                            "TOF_024","TOF_024"),
                                                 Cohort = c("LV","LV","LV",
                                                            "LV","LV","LV",
                                                            "LV","LV","LV",
                                                            "LV","LV","LV","LV",
                                                            "LV","LV","LV",
                                                            "LV","LV","LV",
                                                            "LV","LV","LV","LV",
                                                            "LV","LV","LV",
                                                            "LV","LV","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN","BCN","BCN",
                                                            "BCN"),
                                          Milagro.score = c("R","R","R",
                                                            "R","NR","NR","NR",
                                                            "NR","R","R",
                                                            "NR","NR","R","R",
                                                            "NR","NR","R","R",
                                                            "R","R","NR",
                                                            "NR","R","R","NR",
                                                            "NR","NR","NR",
                                                            "NR","NR","NR","NR",
                                                            "R","R","NR",
                                                            "NR","NR","NR","R",
                                                            "R","R","R","NR",
                                                            "NR","R","R","R",
                                                            "R","NR","NR",
                                                            "R","R","NR","NR",
                                                            "NR","NR","NR",
                                                            "NR","R","R","NR",
                                                            "NR"),
                                              Timepoint = c("Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx","Post-tx",
                                                            "Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx","Pre-tx",
                                                            "Post-tx",
                                                            "Pre-tx","Post-tx"),
                                                    Sex = c("F","F","M",
                                                            "M","M","M","F",
                                                            "F","F","F","M",
                                                            "M","M","M","M",
                                                            "M","F","F","M",
                                                            "M","M","M","M",
                                                            "M","F","F","M",
                                                            "M","F","F","F",
                                                            "F","F","F","M",
                                                            "M","F","F","M",
                                                            "M","M","M","F",
                                                            "F","M","M","M",
                                                            "M","M","M","F",
                                                            "F","F","F","M",
                                                            "M","M","M","F",
                                                            "F","F","F"),
                                                    Age = c(61.21013005,
                                                            61.21013005,43.01711157,
                                                            43.01711157,
                                                            62.83915127,62.83915127,
                                                            52.96098563,
                                                            52.96098563,42.38193018,
                                                            42.38193018,
                                                            45.7412731,45.7412731,
                                                            44.36413415,
                                                            44.36413415,18.11635866,
                                                            18.11635866,28.05475702,
                                                            28.05475702,
                                                            67.63312799,67.63312799,
                                                            73.62902122,
                                                            73.62902122,45.42094456,
                                                            45.42094456,
                                                            48.47091034,48.47091034,
                                                            48.56125941,
                                                            48.56125941,31,31,44,44,
                                                            43,43,54,54,23,
                                                            23,47,47,33,33,
                                                            49,49,40,40,46,
                                                            46,32,32,32,32,
                                                            48,48,23,23,41,
                                                            41,22,22,55,55),
                                       Disease.duration = c(17.44832307,
                                                            17.44832307,16.92265572,
                                                            16.92265572,
                                                            10.02874743,10.02874743,
                                                            12.24914442,
                                                            12.24914442,7.244353183,
                                                            7.244353183,
                                                            27.56741958,27.56741958,
                                                            3.789185489,
                                                            3.789185489,3.403148528,
                                                            3.403148528,
                                                            16.38603696,16.38603696,
                                                            25.66187543,
                                                            25.66187543,12.77481177,
                                                            12.77481177,3.496235455,
                                                            3.496235455,
                                                            1.979466119,1.979466119,
                                                            2.412046543,
                                                            2.412046543,9,9,5,5,4,
                                                            4,6,6,2,2,23,
                                                            23,8,8,14,14,4,
                                                            4,2,2,5,5,8,
                                                            8,28,28,7,7,14,
                                                            14,10,10,20,20),
                     UC.extension = c("Left-side colitis",NA,
                                      "Left-side colitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Proctitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Pancolitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Left-side colitis",
                                      NA,"Pancolitis",
                                      NA,"Proctitis",NA,
                                      "Pancolitis",NA,
                                      "Pancolitis",NA,
                                      "Proctitis",NA,
                                      "Left-side colitis",NA,
                                      "Proctitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Pancolitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Pancolitis",NA,
                                      "Proctitis",NA,
                                      "Pancolitis",NA,
                                      "Left-side colitis",NA,
                                      "Left-side colitis",NA,
                                      "Proctitis"),
                    Current.treatment.without.biologics = c(NA,NA,
                                                            "Pentasa 4g",NA,NA,NA,
                                                            "Colitofalk 2g",NA,
                                                            "Pentasa 4g",NA,NA,
                                                            NA,"Pentasa 4g",NA,
                                                            NA,NA,NA,NA,NA,
                                                            NA,
                                                            "Claversal 3g systemic",NA,NA,
                                                            NA,NA,NA,
                                                            "Pyridoxine + Isoniazide (positive quantiferon without any signs of active TB)",NA,
                                                            "Azathioprine  + Prednisone",NA,
                                                            "Azathioprine",NA,NA,NA,
                                                            "Topic mesalazine + Topic budesonida",NA,
                                                            "Beclometasond",NA,
                                                            "Mesalazine",NA,
                                                            "Mesalazine + Beclometasond",NA,
                                                            "Mesalazine",NA,
                                                            "Beclometasond",NA,NA,NA,
                                                            NA,NA,
                                                            "Azathioprine",NA,
                                                            "Azathioprine",NA,
                                                            "Azathioprine",NA,
                                                            "Beclometasond",NA,
                                                            "Azathioprine",NA,NA,NA),
                               Total.biologicals.before = c(3L,3L,2L,2L,
                                                            2L,2L,2L,2L,2L,
                                                            2L,3L,3L,2L,2L,
                                                            3L,3L,2L,2L,3L,
                                                            3L,3L,3L,3L,3L,
                                                            2L,2L,3L,3L,3L,
                                                            3L,2L,2L,2L,2L,
                                                            2L,2L,1L,1L,2L,
                                                            2L,0L,0L,2L,2L,
                                                            4L,4L,2L,2L,1L,
                                                            1L,2L,2L,1L,1L,
                                                            2L,2L,3L,3L,2L,
                                                            2L,1L,1L),
                                            `CRP.mg/dl` = c(0.85,0.26,0.09,
                                                            0.03,3.8,0.07,
                                                            1.1,1.14,0.38,0.04,
                                                            1.7,0.44,0.36,
                                                            0.4,0.36,2.76,0.06,
                                                            0.03,0.72,0.04,
                                                            0.26,0.28,0.72,
                                                            0.27,0.09,0.57,1.25,
                                                            3.91,0.4,0.4,1.8,
                                                            4.1,0.4,0.7,0.4,
                                                            0.68,1.5,1.54,1,
                                                            0.4,0.4,0.4,0.4,
                                                            1.27,0.47,0.4,
                                                            0.41,0.4,0.58,1.15,
                                                            0.4,0.4,0.4,0.4,
                                                            0.95,3.5,0.4,
                                                            0.45,0.4,0.4,0.4,
                                                            0.66),
                                        Endoscopic.Mayo = c("2","3","3",
                                                            "1","3","3","3",
                                                            "3","3","3","3",
                                                            "3","3","1","3",
                                                            "3","3","1","3",
                                                            "1","3","3","3",
                                                            "0","3","3","3",
                                                            "3","3","3","3",
                                                            NA,"3","3","3",
                                                            "3","2","3","3",
                                                            "3","3","0","3",
                                                            NA,"3","3","3",
                                                            "1","3","2","3",
                                                            "3","3","3","3",
                                                            "3","3",NA,
                                                            "3","0","3","3"),
                                             Total.Mayo = c("9","5","10",
                                                            "2","11","6","12",
                                                            "11","12","3",
                                                            "11","11","10","3",
                                                            "11","12","8",
                                                            "3","9","1","8",
                                                            "10","11","2","10",
                                                            "12","12","12",
                                                            "5","7","10",NA,
                                                            "4","3","9","9",
                                                            "5","7","6","4",
                                                            "7","0","8",
                                                            NA,"8","8","10",
                                                            "2","10","9","10",
                                                            "5","9","9","6",
                                                            "8","9",NA,
                                                            "8","0","11","9")
                  )

treatments <- data.frame(
  stringsAsFactors = FALSE,
         Treatment = c("Aminosalicylates",
                       "Corticosteroids","Immunomodulator","Monotherapy",
                       "Combination therapy"),
               BCN = c(3L, 6L, 6L, 10L, 3L),
                LV = c(5L, 0L, 0L, 5L, 0L)
)

# sex and response -----------------------------------------
fi_sex <- table(data_by_patient$Milagro.score, data_by_patient$Sex)
fisher.test(fi_sex)
chisq.test(fi_sex)

# sex and leuven -----------------------------------------
fi_sex <- table(data_by_patient$Cohort, data_by_patient$Sex)
fisher.test(fi_sex)
chisq.test(fi_sex)

# age and response -----------------------------------------

shapiro.test(data_by_patient$Age)
t.test(data_by_patient$Age ~ data_by_patient$Milagro.score)
wilcox.test(data_by_patient$Age ~ data_by_patient$Milagro.score)

# age and cohort -----------------------------------------

shapiro.test(data_by_patient$Age[data_by_patient$Cohort == 'LV'])
shapiro.test(data_by_patient$Age[data_by_patient$Cohort == 'BCN'])
t.test(data_by_patient$Age ~ data_by_patient$Cohort)
wilcox.test(data_by_patient$Age ~ data_by_patient$Cohort)

# disease duration and response -----------------------------------------

shapiro.test(data_by_patient$Disease.duration)
wilcox.test(data_by_patient$Disease.duration ~ data_by_patient$Milagro.score)

# disease duration and cohort -----------------------------------------

shapiro.test(data_by_patient$Disease.duration[data_by_patient$Cohort == 'LV'])
shapiro.test(data_by_patient$Disease.duration[data_by_patient$Cohort == 'BCN'])

wilcox.test(data_by_patient$Disease.duration ~ data_by_patient$Cohort)

# UC extension and leuven -----------------------------------------
fi_ext <- table(data_by_patient$Cohort, data_by_patient$UC.extension)
fisher.test(fi_ext)
chisq.test(fi_ext)


# prev drugs and response -----------------------------------------

shapiro.test(data_by_patient$Total.biologicals.before)
wilcox.test(data_by_patient$Total.biologicals.before ~ data_by_patient$Milagro.score)
.
# prev drugs and cohort -----------------------------------------

shapiro.test(data_by_sample_w0$Total.biologicals.before)
wilcox.test(data_by_sample_w0$Total.biologicals.before ~ data_by_sample_w0$Cohort)

# crp by group ------------------------------------------------------------------
library(ggplot2)
ggplot(data_by_sample, aes(x = `CRP.mg/dl`, fill = Timepoint, color= Milagro.score)) +
  geom_histogram()

data_by_sample$myvar <- paste(data_by_sample$Milagro.score, data_by_sample$Timepoint, sep='_')
shapiro.test(data_by_sample$`CRP.mg/dl`) # no normal
kruskal.test(data_by_sample$`CRP.mg/dl`~data_by_sample$myvar)

library(FSA)
dunnTest(data_by_sample$`CRP.mg/dl`~data_by_sample$myvar)



wilcox.test(x = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "R_Pre-tx"],
            y = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "R_Post-tx"], paired = TRUE, alternative = "two.sided")

wilcox.test(x = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "NR_Pre-tx"],
            y = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "NR_Post-tx"], paired = TRUE, alternative = "two.sided")


# crp by cohort ------------------------------------------------------------------
library(ggplot2)
ggplot(data_by_sample_w0, aes(x = `CRP.mg/dl`, fill = Cohort, color= Cohort)) +
  geom_histogram()

data_by_sample$myvar <- paste(data_by_sample$Cohort, data_by_sample$Timepoint, sep='_')
shapiro.test(data_by_sample$`CRP.mg/dl`) # no normal
kruskal.test(data_by_sample$`CRP.mg/dl`~data_by_sample$myvar)

library(FSA)
dunnTest(data_by_sample$`CRP.mg/dl`~data_by_sample$myvar)



wilcox.test(x = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "LV_Pre-tx"],
            y = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "BCN_Pre-tx"],
            paired = FALSE)

wilcox.test(x = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "LV_Post-tx"],
            y = data_by_sample$`CRP.mg/dl`[data_by_sample$myvar == "BCN_Post-tx"],
            paired = FALSE)



# endoscopic mayo by group ------------------------------------------------------------------
library(ggplot2)
data_by_sample$Endoscopic.Mayo <- as.numeric(data_by_sample$Endoscopic.Mayo)
ggplot(data_by_sample, aes(x = Endoscopic.Mayo, fill = Timepoint, color= Milagro.score)) +
  geom_histogram()

data_by_sample$myvar <- paste(data_by_sample$Milagro.score, data_by_sample$Timepoint, sep='_')
shapiro.test(data_by_sample$Endoscopic.Mayo) # no normal
kruskal.test(data_by_sample$Endoscopic.Mayo~data_by_sample$myvar)

library(FSA)
dunnTest(data_by_sample$Endoscopic.Mayo~data_by_sample$myvar)

wilcox.test(x = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "R_Pre-tx"],
            y = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "R_Post-tx"], paired = TRUE, alternative = "two.sided")

wilcox.test(x = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "NR_Pre-tx"],
            y = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "NR_Post-tx"], paired = TRUE, alternative = "two.sided")

# Endoscopic.Mayo by cohort ------------------------------------------------------------------
library(ggplot2)
ggplot(data_by_sample_w0, aes(x = as.numeric(Endoscopic.Mayo), fill = Cohort, color= Cohort)) +
  geom_histogram()

data_by_sample$myvar <- paste(data_by_sample$Cohort, data_by_sample$Timepoint, sep='_')

data_by_sample$myvar <- paste(data_by_sample$Milagro.score, data_by_sample$Cohort, data_by_sample$Timepoint, sep='_')
shapiro.test(data_by_sample$Endoscopic.Mayo) # no normal
kruskal.test(data_by_sample$Endoscopic.Mayo~data_by_sample$myvar)

library(FSA)
dunnTest(data_by_sample$Endoscopic.Mayo~data_by_sample$myvar)


wilcox.test(x = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "LV_Pre-tx"],
            y = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "BCN_Pre-tx"],
            paired = FALSE)

wilcox.test(x = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "LV_Post-tx"],
            y = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "BCN_Post-tx"],
            paired = FALSE)

wilcox.test(x = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "NR_LV_Pre-tx"],
            y = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "NR_BCN_Pre-tx"],
            paired = FALSE)

wilcox.test(x = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "NR_LV_Post-tx"],
            y = data_by_sample$Endoscopic.Mayo[data_by_sample$myvar == "NR_BCN_Post-tx"],
            paired = FALSE)

# Total.Mayo by cohort ------------------------------------------------------------------
library(ggplot2)
ggplot(data_by_sample_w0[data_by_sample_w0$Milagro.score == 'NR',], aes(y = as.numeric(Total.Mayo), x = Cohort, color= Cohort)) +
  geom_boxplot() +
  labs(title= 'NR pre')+
  ggpubr::stat_compare_means(method='wilcox')

ggplot(data_by_sample[data_by_sample$Milagro.score == 'R' &
                        data_by_sample$Timepoint == 'Pre-tx',],
       aes(y = as.numeric(Total.Mayo), x = Cohort, color= Cohort)) +
  geom_boxplot() +
  labs(title= 'R pre')+
  ggpubr::stat_compare_means(method='wilcox')

ggplot(data_by_sample[data_by_sample$Milagro.score == 'NR' &
                        data_by_sample$Timepoint == 'Post-tx',],
       aes(y = as.numeric(Total.Mayo), x = Cohort, color= Cohort)) +
  geom_boxplot() +
  labs(title= 'NR post')+
  ggpubr::stat_compare_means(method='wilcox')

ggplot(data_by_sample[data_by_sample$Milagro.score == 'R' &
                        data_by_sample$Timepoint == 'Post-tx',],
       aes(y = as.numeric(Total.Mayo), x = Cohort, color= Cohort)) +
  geom_boxplot() +
  labs(title= 'R post') +
  ggpubr::stat_compare_means(method='wilcox')

data_by_sample$Total.Mayo <- as.numeric(data_by_sample$Total.Mayo)

data_by_sample$myvar <- paste(data_by_sample$Cohort, data_by_sample$Timepoint, sep='_')
shapiro.test(data_by_sample$Total.Mayo) # no normal
kruskal.test(data_by_sample$Total.Mayo~data_by_sample$myvar)

library(FSA)
dunnTest(data_by_sample$Total.Mayo~data_by_sample$myvar)


wilcox.test(x = data_by_sample$Total.Mayo[data_by_sample$myvar == "LV_Pre-tx"],
            y = data_by_sample$Total.Mayo[data_by_sample$myvar == "BCN_Pre-tx"],
            paired = FALSE)

wilcox.test(x = data_by_sample$Total.Mayo[data_by_sample$myvar == "LV_Post-tx"],
            y = data_by_sample$Total.Mayo[data_by_sample$myvar == "BCN_Post-tx"],
            paired = FALSE)



# total mayo by group ------------------------------------------------------------------
library(ggplot2)
data_by_sample$Total.Mayo <- as.numeric(data_by_sample$Total.Mayo)
ggplot(data_by_sample, aes(x = Total.Mayo, fill = Timepoint, color= Milagro.score)) +
  geom_histogram()

# data_by_sample$myvar <- paste(data_by_sample$Milagro.score, data_by_sample$Timepoint, sep='_')
shapiro.test(data_by_sample$Total.Mayo) # no normal
kruskal.test(data_by_sample$Total.Mayo~data_by_sample$myvar)

library(FSA)
dunnTest(data_by_sample$Total.Mayo~data_by_sample$myvar)



wilcox.test(x = data_by_sample$Total.Mayo[data_by_sample$myvar == "R_Pre-tx"],
            y = data_by_sample$Total.Mayo[data_by_sample$myvar == "R_Post-tx"], paired = TRUE, alternative = "two.sided")

wilcox.test(x = data_by_sample$Total.Mayo[data_by_sample$myvar == "NR_Pre-tx"],
            y = data_by_sample$Total.Mayo[data_by_sample$myvar == "NR_Post-tx"], paired = TRUE, alternative = "two.sided")

