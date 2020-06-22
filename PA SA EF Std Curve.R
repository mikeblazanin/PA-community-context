library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.csv("2020-02-17 PA SA EF Std Curve Tidy.csv")

data$OD600_plt_adj <- data$OD600_plt-data$OD600_plt[data$Bacteria == "Blank"]
data$OD600_spec_adj <- data$OD600_spec-data$OD600_spec[data$Bacteria == "Blank"]

pa_std_curve_plt <- lm(cfu.mL ~ OD600_plt_adj + 0,
                   data[data$Bacteria == "P. aeruginosa" &
                          data$OD600_plt_adj < 1, ])
sa_std_curve_plt <- lm(cfu.mL ~ OD600_plt_adj + 0,
                   data[data$Bacteria == "S. aureus" &
                          data$OD600_plt_adj < 1, ])
ef_std_curve_plt <- lm(cfu.mL ~ OD600_plt_adj + 0,
                   data[data$Bacteria == "E. faecalis", ])

pa_std_curve_spec <- lm(cfu.mL ~ OD600_spec_adj + 0,
                       data[data$Bacteria == "P. aeruginosa" &
                              data$OD600_spec_adj < 3, ])
sa_std_curve_spec <- lm(cfu.mL ~ OD600_spec_adj + 0,
                       data[data$Bacteria == "S. aureus" &
                              data$OD600_spec_adj < 2, ])
ef_std_curve_spec <- lm(cfu.mL ~ OD600_spec_adj + 0,
                       data[data$Bacteria == "E. faecalis", ])

pa_std_curve_plt
sa_std_curve_plt
ef_std_curve_plt
pa_std_curve_spec
sa_std_curve_spec
ef_std_curve_spec

plot_lines <- data.frame(intercept = 0,
                         Bacteria = c("P. aeruginosa",
                                      "S. aureus",
                                      "E. faecalis"),
                         spec_slope = c(pa_std_curve_spec$coefficients,
                                        sa_std_curve_spec$coefficients,
                                        ef_std_curve_spec$coefficients),
                         plt_slope = c(pa_std_curve_plt$coefficients,
                                       sa_std_curve_plt$coefficients,
                                       ef_std_curve_plt$coefficients))

ggplot(data = data[data$Bacteria != "Blank", ], 
       aes(x = OD600_plt_adj, y = cfu.mL)) +
  geom_point() +
  geom_abline(data = plot_lines, aes(intercept = intercept, slope = plt_slope)) +
  facet_wrap(~Bacteria, scales = "free") +
  NULL

ggplot(data = data[data$Bacteria != "Blank", ], 
       aes(x = OD600_spec_adj, y = cfu.mL)) +
  geom_point() +
  geom_abline(data = plot_lines, aes(intercept = intercept, slope = spec_slope)) +
  facet_wrap(~Bacteria, scales = "free")

##2020-06 standard curve ----
od_2 <- read.csv("2020-06 PA SA EF Std Curve ODs.csv",
                 stringsAsFactors = FALSE)
cfu_2 <- read.csv("2020-06 PA SA EF Std Curve Counts.csv",
                  stringsAsFactors = FALSE)

#Reshape od
od_2 <- pivot_longer(od_2, cols = starts_with("X"),
                     names_to = "dilution", values_to = "OD600")
od_2$od_dilution <- sapply(od_2$dilution, function(x) {substr(x, 2, nchar(x))})

#run calculations
od_2$OD600_adj <- od_2$OD600 - od_2$Blank_avg
cfu_2$cfu_ml <- cfu_2$cfu/cfu_2$Dilution
cfu_2 <- group_by(cfu_2, Bacteria, Rep)
cfu_2_sum <- summarize(cfu_2,
                       cfu_ml_avg = mean(cfu_ml, na.rm = TRUE))

#Join cfu & od data
std_curve2 <- left_join(od_2, cfu_2_sum, by = c("Bacteria", "Rep"))
#Calculate cfu/ml in dilution series
std_curve2$cfu_ml_dil <- as.numeric(std_curve2$od_dilution)*std_curve2$cfu_ml_avg
#Drop unneeded cols & rows
std_curve2 <- select(std_curve2, Machine, Bacteria, Rep, 
                     OD600_adj, cfu_ml_dil)
std_curve2 <- std_curve2[complete.cases(std_curve2), ]

#Run linear models
lm2_plt <- lm(data = std_curve2[std_curve2$Machine == "Plate" &
                                  std_curve2$OD600_adj <= 1, ],
              cfu_ml_dil ~ OD600_adj + OD600_adj:Bacteria + 0)
lm2_spec <- lm(data = std_curve2[std_curve2$Machine == "Spec" &
                                   std_curve2$OD600_adj <= 1, ],
               cfu_ml_dil ~ OD600_adj + OD600_adj:Bacteria + 0)

#Calc lines
std_curve2$pred_cfu <- NA
std_curve2$pred_cfu[std_curve2$Machine == "Plate"] <-
  predict(lm2_plt, std_curve2)[std_curve2$Machine == "Plate"]
std_curve2$pred_cfu[std_curve2$Machine == "Spec"] <-
  predict(lm2_spec, std_curve2)[std_curve2$Machine == "Spec"]

#Plot lines
ggplot(data = std_curve2, aes(x = cfu_ml_dil, y = OD600_adj)) +
  geom_point() +
  facet_grid(Machine ~ Bacteria, scales = "free") +
  geom_line(aes(x = pred_cfu))

ggplot(data = std_curve2, aes(x = OD600_adj, y = cfu_ml_dil)) +
  geom_point() +
  facet_grid(Bacteria ~ Machine, scales = "free") +
  geom_line(aes(y = pred_cfu))

#Print models
lm2_plt
lm2_spec
