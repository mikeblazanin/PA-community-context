library(gcplyr)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

if("growth_curves" %in% list.files()) {setwd("./growth_curves")}

#Tidy ----
dat <- read_wides(files = list.files(pattern = "Data.csv"),
                  startrow = 41, startcol = "B")
dat <- lapply(dat, FUN = function(x) {select(x, !3)})
dat <- trans_wide_to_tidy(dat, id_cols = c("file", "Time"))
dat <- merge_dfs(dat, collapse = TRUE)
dat$Time <- lubridate::time_length(lubridate::hms(dat$Time))
dat <- mutate(dat, 
              date = str_split_i(file, " ", 1),
              Time_hr = Time/60/60)

#Designs ----
designs <- read_blocks(list.files(pattern = "Design.csv"))
designs <- trans_block_to_wide(designs)
designs <- trans_wide_to_tidy(designs, id_cols = "block_name", 
                              values_to_numeric = FALSE)
designs <- separate_tidy(
  designs, "Measurements", sep = ", ",
  into = c("Bact_pop", "Bact_timepoint", "Bact_treat", "Bact_isol"))

designs <- mutate(designs,
                  date = str_split_i(block_name, " ", 1),
                  Bact_pop = gsub("^P", "", Bact_pop),
                  Bact_timepoint = gsub("^T", "", Bact_timepoint),
                  Bact_isol = gsub("^#", "", Bact_isol),
                  Bact_treat_C = factor(
                    ifelse(Bact_treat %in% c("NCP", "NN"), "NC",
                           ifelse(Bact_treat %in% c("CNP", "CP"), "C",
                                  Bact_treat)),
                    levels = c("Anc", "NC", "C", "Blank")),
                  Bact_treat_P = factor(
                    ifelse(Bact_treat %in% c("NN", "CNP"), "NP",
                           ifelse(Bact_treat %in% c("NCP", "CP"), "P",
                                  Bact_treat)),
                    levels = c("Anc", "NP", "P", "Blank")),
                  Bact_treat = factor(
                    Bact_treat,
                    levels = c("Anc", "NN", "CNP", "NCP", "CP", "Blank")))
                               
dat <- merge_dfs(dat, designs, drop = TRUE)

#Exclude contaminated pops
dat <- filter(dat,
              !((Bact_pop == 9 & Bact_treat == "NN" & Bact_timepoint >= 4) |
                  (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

#Subtract blank ----
dat <- mutate(group_by(dat, file),
              meas_blanked = Measurements - 
                min(Measurements[Bact_treat == "Blank"], na.rm = TRUE))

ggplot(data = dat,
       aes(x = Time_hr, y = Measurements, group = Well)) +
  geom_line() +
  facet_wrap(~file)
ggplot(data = dat,
       aes(x = Time_hr, y = meas_blanked, group = Well)) +
  geom_line() +
  facet_wrap(~date)

#Filter data, excluding any bad wells
dat <- filter(dat,
              Bact_treat != "Blank",
              #Bad wells
              !(date == "2024-01-26" & Bact_pop == 6 &
                  Bact_treat == "CNP" & Bact_timepoint == 8 &
                  Bact_isol == 1))

#Smooth ----
dat <- mutate(group_by(dat, date, Well),
              sm_meas = smooth_data(x = Time_hr, y = meas_blanked,
                                    sm_method = "moving-median",
                                    window_width_n = 3),
              sm_meas = smooth_data(x = Time_hr, y = sm_meas,
                                     sm_method = "moving-average",
                                     window_width_n = 5))

ggplot(data = dat,
       aes(x = Time_hr, y = sm_meas, group = Well)) +
  geom_line() +
  facet_wrap(~date)

#Deriv ----
dat <- mutate(group_by(dat, date, Well),
              deriv = calc_deriv(y = sm_meas, x = Time_hr, window_width_n = 5),
              deriv_percap = calc_deriv(
                y = sm_meas, x = Time_hr, window_width_n = 5,
                percapita = TRUE, trans_y = 'log', blank = 0))

#Assess noise
for (mydate in unique(dat$date)) {
  png(paste0("./plots/sm_deriv_derivpercap_", mydate, ".png"),
      width = 4, height = 10, units = "in", res = 150)
  print(
    cowplot::plot_grid(ncol = 1,
      ggplot(data = filter(dat, date == mydate),
             aes(x = Time_hr, y = sm_meas, color = Bact_isol)) +
        geom_line(aes(group = Well), alpha = 0.6) +
        facet_grid(Bact_timepoint ~ Bact_treat),
      ggplot(data = filter(dat, date == mydate),
             aes(x = Time_hr, y = deriv, color = Bact_isol)) +
        geom_line(aes(group = Well), alpha = 0.6) +
        facet_grid(Bact_timepoint ~ Bact_treat),
      ggplot(data = filter(dat, date == mydate),
             aes(x = Time_hr, y = deriv_percap, color = Bact_isol)) +
        geom_line(aes(group = Well), alpha = 0.6) +
        facet_grid(Bact_timepoint ~ Bact_treat)
    ))
  dev.off()
}

#Summarize ----
dat_sum <- summarize(
  group_by(filter(dat, Time_hr < 24), 
           date, Well, Bact_treat, Bact_treat_C, Bact_treat_P,
           Bact_timepoint, Bact_pop, Bact_isol),
  max_percap = max_gc(deriv_percap),
  max_dens = max_gc(sm_meas[Time_hr < 23]),
  lag = lag_time(x = Time_hr, y = sm_meas, deriv = deriv_percap,
                 y0 = min_gc(sm_meas)),
  auc = auc(x = Time_hr, y = sm_meas),
  max_percap_time = Time_hr[which_max_gc(deriv_percap)],
  max_percap_dens = sm_meas[which_max_gc(deriv_percap)],
  min_dens = min_gc(sm_meas))
  
#Plot ----
png("../plots/max_dens.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = dat_sum,
       aes(x = Bact_timepoint, y = max_dens,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ Bact_pop) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Maximum Density (OD600)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

png("../plots/max_percap.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = dat_sum,
       aes(x = Bact_timepoint, y = max_percap,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ Bact_pop) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Maximum Growth Rate (/hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

ggplot(data = dat_sum,
       aes(x = Bact_timepoint, y = max_percap_time, color = Bact_treat)) +
  geom_point(size =3, alpha = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(~ Bact_pop) +
  theme_bw()
       
png("../plots/auc.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = dat_sum,
       aes(x = Bact_timepoint, y = auc,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ Bact_pop) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Area Under the Curve\n(OD600 hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

for(pop in unique(dat$Bact_pop)) {
  png(paste0("lagcheck_pop", pop, ".png"), 
             width = 7, height = 12, units = "in", res = 300)
  print(
    ggplot(data = filter(dat, Bact_pop == pop),
           aes(x = Time_hr, y = log(sm_meas))) +
      geom_point(size = 0.5, alpha = 0.8) +
      theme_bw() +
      facet_wrap(~ Bact_treat * Bact_timepoint * Bact_isol) +
      geom_hline(data = filter(dat_sum, Bact_pop == pop),
                 aes(yintercept = log(min_dens)), 
                 color = "red", linewidth = 0.5, alpha = 0.3) +
      geom_vline(data = filter(dat_sum, Bact_pop == pop),
                 aes(xintercept = lag), 
                 color = "red", linewidth = 0.5, alpha = 0.3) +
      geom_abline(data = filter(dat_sum, Bact_pop == pop),
                  color = "red", alpha = 0.3,
                  aes(slope = max_percap, 
                      intercept = log(max_percap_dens) - max_percap*max_percap_time)) +
      geom_point(data = filter(dat_sum, Bact_pop == pop), color = "red", size = 0.5,
                 aes(x = max_percap_time, y = log(max_percap_dens)))
  )
  dev.off()
}

png("../plots/lag_time.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = dat_sum,
       aes(x = Bact_timepoint, y = lag,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ Bact_pop) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Lag Time (hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

write.csv(dat_sum, "../growth_curve_metrics.csv", row.names = FALSE)
