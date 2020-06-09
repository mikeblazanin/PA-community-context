library(ggplot2)
library(tidyr)

setwd("C:/users/mikeb/Google Drive/Research Projects/11 P aeruginosa Community/")
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

plot(data[data$Bacteria == "P. aeruginosa", "cfu.mL"],
     data[data$Bacteria == "P. aeruginosa", "OD600_adj"])

plot(data[data$Bacteria == "P. aeruginosa", "OD600_adj"],
     data[data$Bacteria == "P. aeruginosa", "cfu.mL"])
abline(pa_std_curve)
