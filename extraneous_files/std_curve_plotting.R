stop("all paths need to be corrected to './extraneous_files/'")

library(ggplot2)

std_cv_data <- read.csv("2020-02-18_std_curve.csv",
                        stringsAsFactors = F)
colnames(std_cv_data)[4:11] <- 2:9

std_cv_data$cfu_ml <- NA

for (row in 1:nrow(std_cv_data)) {
  for (col in 11:4) {
    if (!is.na(std_cv_data[row, col]) &
        !is.na(as.numeric(std_cv_data[row, col]))) {
      std_cv_data$cfu_ml[row] <- as.numeric(std_cv_data[row, col])*
        10**(1+as.numeric(colnames(std_cv_data)[col]))
    }
  }
}

std_cv_data$rel_cfu_ml <- NA
for (bact in unique(std_cv_data$Bacteria)) {
  myrows <- std_cv_data$Bacteria == bact
  std_cv_data$rel_cfu_ml[myrows] <- std_cv_data$cfu_ml[myrows]/
    std_cv_data$cfu_ml[myrows & std_cv_data$Media.Type == "Rich"]
}

tiff("std_curve.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data = std_cv_data,
       aes(x = Bacteria, y = cfu_ml, color = Media.Type)) +
  geom_point(size = 3, position = position_dodge(width = 0.1), alpha = 0.7) +
  theme_bw() +
  labs(y = "Overnight Culture Density (cfu/mL)") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("black", "gray60")) +
  NULL
dev.off()
