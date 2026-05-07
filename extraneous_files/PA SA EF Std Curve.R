stop("all paths need to be corrected to './extraneous_files/'")

library(ggplot2)
library(tidyr)
library(dplyr)

##2020-02-17 std curve ----
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
cfu_2 <- group_by(cfu_2, Date, Bacteria, Rep)
cfu_2_sum <- summarize(cfu_2,
                       cfu_ml_avg = mean(cfu_ml, na.rm = TRUE))

#Join cfu & od data
std_curve2 <- left_join(od_2, cfu_2_sum, 
                        by = c("Date", "Bacteria", "Rep"))
#Calculate cfu/ml in dilution series
std_curve2$cfu_ml_dil <- as.numeric(std_curve2$od_dilution)*std_curve2$cfu_ml_avg
#Drop unneeded cols & rows
std_curve2 <- select(std_curve2, Date, Machine, Bacteria, Rep, 
                     OD600_adj, cfu_ml_dil)
std_curve2 <- std_curve2[complete.cases(std_curve2), ]

#Run linear models
lm2 <- lm(data = std_curve2[(std_curve2$Machine == "Plate" &
                              std_curve2$OD600_adj <= 0.5) |
                              (std_curve2$Machine == "Spec" &
                              std_curve2$OD600_adj <= 1), ],
          cfu_ml_dil ~ OD600_adj:Bacteria:Machine + 0)

#Calc lines
std_curve2$pred_cfu <- predict(lm2, std_curve2)

#Plot data & lines in linear
ggplot(data = std_curve2, 
       aes(x = OD600_adj, y = cfu_ml_dil)) +
  geom_point(aes(color = Date)) +
  facet_grid(Bacteria ~ Machine, scales = "free") +
  geom_line(aes(y = pred_cfu)) +
  # scale_y_continuous(trans = "log10") +
  # scale_x_continuous(trans = "log10") +
  NULL

#Plot data & lines in log-log space
ggplot(data = std_curve2, 
       aes(x = OD600_adj, y = cfu_ml_dil)) +
  geom_point(aes(color = Date)) +
  facet_grid(Bacteria ~ Machine, scales = "free") +
  geom_line(aes(y = pred_cfu)) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  NULL

#Print model information
lm2

# 2020-07 Standard Curve (static, Bolton Broth) ----
std_curve3_od <- read.csv("2020-07-12 PA SA EF Std Curve Tidy.csv",
                       stringsAsFactors = FALSE)
std_curve3_cfu <- read.csv("2020-07-12 PA SA EF Std Curve CFU.csv",
                           stringsAsFactors = FALSE)

std_curve3_od <- group_by(std_curve3_od,
                       Machine, Bacteria, Rep, Dilution)
std_curve3_od_sum <- summarize(std_curve3_od,
                            "OD600" = mean(OD600))
std_curve3_od_sum$OD600_adj <- NA
std_curve3_od_sum$OD600_adj[std_curve3_od_sum$Machine == "Plate"] <-
  std_curve3_od_sum$OD600[std_curve3_od_sum$Machine == "Plate"] -
  std_curve3_od_sum$OD600[std_curve3_od_sum$Machine == "Plate" & 
                         std_curve3_od_sum$Bacteria == "Blank"]
std_curve3_od_sum$OD600_adj[std_curve3_od_sum$Machine == "Spec"] <-
  std_curve3_od_sum$OD600[std_curve3_od_sum$Machine == "Spec"] -
  std_curve3_od_sum$OD600[std_curve3_od_sum$Machine == "Spec" & 
                         std_curve3_od_sum$Bacteria == "Blank"]
std_curve3_od_sum <- std_curve3_od_sum[std_curve3_od_sum$Bacteria != "Blank",]

std_curve3_cfu$cfu_ml <- std_curve3_cfu$cfu/std_curve3_cfu$Dilution
  
std_curve3 <- left_join(std_curve3_od_sum, std_curve3_cfu,
                        by = c("Bacteria", "Rep"))
std_curve3$cfu_ml_dil <- std_curve3$cfu_ml*std_curve3$Dilution.x

lm3 <- lm(data = std_curve3, 
                         cfu_ml_dil ~ OD600_adj:Bacteria:Machine + 0)
std_curve3$pred_cfu <- predict(lm3, std_curve3)

ggplot(data = std_curve3, 
       aes(x = OD600_adj, y = cfu_ml_dil, color = Rep)) +
  geom_point() +
  facet_grid(Bacteria ~ Machine) +
  geom_line(aes(y = pred_cfu))

ggplot(data = std_curve3, 
       aes(x = OD600_adj, y = cfu_ml_dil, color = Rep)) +
  geom_point() +
  facet_grid(Bacteria ~ Machine) +
  geom_line(aes(y = pred_cfu)) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10")

lm3


