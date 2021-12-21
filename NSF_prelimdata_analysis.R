library(ggplot2)

#Read data
H6growth <- xlsx::read.xlsx("./Katie_data/2021.12.13-H6DataForMike.xlsx",
                            sheetIndex = 1, startRow = 2)
H6InSeq <- xlsx::read.xlsx("./Katie_data/2021.12.13-H6DataForMike.xlsx",
                           sheetIndex = 2, colIndex = 1:6)
H6EOP <- xlsx::read.xlsx("./Katie_data/2021.12.13-H6DataForMike.xlsx",
                         sheetIndex = 3, startRow = 2, colIndex = 1:3)

#Single step phage growth curve
H6growth <- tidyr::pivot_longer(H6growth, cols = !Time..min.,
                                names_sep = "_", names_to = c("Treat", "Rep"),
                                values_to = "pfu_ml")
H6growth <- H6growth[complete.cases(H6growth), ]
lm(log10(pfu_ml) ~ Time..min., H6growth[H6growth$Treat %in% c("Filtered") & 
                                          H6growth$Time..min. %in% c(0, 10), ])

png("H6_singlestep_growth.png", width = 4, height = 4, units = "in", res = 150)
ggplot(data = H6growth[H6growth$Treat %in% c("Control.Filtered", "Filtered"), ], 
       aes(x = Time..min., y = pfu_ml,
                            color = Treat, group = paste(Rep, Treat))) +
  geom_line(lwd = 2, alpha = 0.75) +
  scale_y_continuous(trans = "log10", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Time (min)", y = "Phage particle\nconcentration (pfu/mL)") +
  scale_color_discrete(breaks = c("Filtered", "Control.Filtered"),
    labels = c("Bacteria\nAdded", "No-bacteria\nControl"),
    name = "Treatment") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.2),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.background = element_rect(color = "black")) +
  NULL
dev.off()

#Inseq results for phage receptor
H6InSeq$plot_color <- ifelse(H6InSeq$Pval < 0.05, "sig", "nonsig")
H6InSeq$plot_color[H6InSeq$Pval < 0.05 & 
                     grepl("pil*", H6InSeq$Gene.Name)] <- "pilus"
png("H6_InSeq.png", width = 4, height = 4, units = "in", res = 150)
ggplot(data = H6InSeq,
       aes(x = Input.Average..1, y = Output.Average..1,
           color = plot_color)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_x_continuous(trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Gene input frequency", y = "Gene output frequency") +
  scale_color_manual(breaks = c("pilus", "sig", "nonsig"),
                     values = c("red", "black", "gray"),
                     labels = c("Pilus significant hit", 
                                "Non-pilus\nsignificant hit",
                                "Non-significant hit")) +
  theme_bw() +
  theme(legend.position = c(0.275, 0.85),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(color = "black"),
        legend.margin = margin(t = -0.05, l = 0.05, 
                               r = 0.05, b = 0.05, unit = "in"),
        plot.margin = margin(t = 5.5, l = 5.5, r = 10, b = 5.5, unit = "pt")) +
  NULL
dev.off()

#EOP results
H6EOP_sum <- dplyr::summarize(dplyr::group_by(H6EOP, Bacteria),
                   mean_titer = mean(H6.Titer))
png("H6_EOP.png", width = 2, height = 4, units = "in", res = 150)
ggplot(data = H6EOP_sum[H6EOP_sum$Bacteria %in% c("PA14", "PA14^pilA"), ], 
                    aes(x = Bacteria, y = mean_titer+1)) +
  geom_col() +
  scale_y_continuous(trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +  labs(y = "Phage plaquing (pfu/mL)") +
  scale_x_discrete(breaks = c("PA14", "PA14^pilA"),
                   labels = c("PA14", "PA14 \u0394pilA")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1)) +
  NULL
dev.off()
