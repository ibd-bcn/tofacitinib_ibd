library(readr)

corr_tofa <- read_delim("Analysis/data/corr_tofa_sangre.csv",
                        delim = ";", escape_double = FALSE,
                        locale = locale(decimal_mark = ",",
                                        grouping_mark = "."),
                        na = c("ND","n.d"), trim_ws = TRUE)

# socs 1 -------------------------------------------------------------------------

y <-  corr_tofa$SOCS1[!is.na(corr_tofa$`Conc tofa (ng/ml)`) & !is.na(corr_tofa$SOCS1)]
x <-  corr_tofa$`Conc tofa (ng/ml)`[!is.na(corr_tofa$`Conc tofa (ng/ml)`) & !is.na(corr_tofa$SOCS1)]

cx <- x - mean(x)

model <- lm( y ~ log(x))

#view the output of the model
summary(model)


# socs 3 -------------------------------------------------------------------------

y <-  corr_tofa$SOCS3[!is.na(corr_tofa$`Conc tofa (ng/ml)`) & !is.na(corr_tofa$SOCS3)]
x <-  corr_tofa$`Conc tofa (ng/ml)`[!is.na(corr_tofa$`Conc tofa (ng/ml)`) & !is.na(corr_tofa$SOCS3)]

model <- lm(y ~ log(x))

#view the output of the model
summary(model)


# IRF1 -------------------------------------------------------------------------

y <-  corr_tofa$IRF1[!is.na(corr_tofa$`Conc tofa (ng/ml)`) & !is.na(corr_tofa$IRF1)]
x <-  corr_tofa$`Conc tofa (ng/ml)`[!is.na(corr_tofa$`Conc tofa (ng/ml)`) & !is.na(corr_tofa$IRF1)]

model <- lm(y ~ log(x))

#view the output of the model
summary(model)
