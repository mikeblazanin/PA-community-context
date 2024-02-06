library(gcplyr)
library(ggplot2)
library(dplyr)
library(stringr)

setwd("./growth_curves")

#Tidy ----
dat <- read_wides(files = list.files(pattern = "Data.csv"),
                  startrow = 41, startcol = "B")
dat <- lapply(dat, FUN = function(x) {select(x, !3)})
dat <- trans_wide_to_tidy(dat, id_cols = c("file", "Time"))
dat <- merge_dfs(dat, collapse = TRUE)
dat$Time <- lubridate::time_length(lubridate::hms(dat$Time))
dat2 <- mutate(read_tidys("2024-02-02 William Growth Curve Pop 7 Tidy.csv"),
               Measurements = as.numeric(Measurements),
               Time = as.numeric(Time))
dat <- merge_dfs(dat, dat2)
dat <- mutate(dat, date = str_split_i(file, " ", 1))

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
                    levels = c("Anc", "NC", "C")),
                  Bact_treat_P = factor(
                    ifelse(Bact_treat %in% c("NN", "CNP"), "NP",
                           ifelse(Bact_treat %in% c("NCP", "CP"), "P",
                                  Bact_treat)),
                    levels = c("Anc", "NP", "P")),
                  Bact_treat = factor(Bact_treat,
                                      levels = c("Anc", "NN", "CNP", "NCP", "CP")))
                               
dat <- merge_dfs(dat, designs, drop = TRUE)

#Subtract blank ----
dat <- mutate(group_by(dat, file),
              meas_blanked = Measurements - mean(Measurements[Treat == "Blank"]))
dat <- filter(dat, Treat != "Blank")

#Smooth ----
dat <- mutate(group_by(dat, file, Well),
              sm_meas = smooth_data(x = Time, y = meas_blanked,
                                    sm_method = "moving-median",
                                    window_width_n = 3),
              sm_meas = smooth_data(x = Time, y = sm_meas,
                                     sm_method = "moving-average",
                                     window_width_n = 5))

ggplot(data = dat, aes(x = Time, y = sm_meas,
                       color = Treat, group = paste(file, Well))) +
  geom_line() +
  facet_grid(Pop~Timepoint)

#Deriv ----
dat <- mutate(group_by(dat, file, Well),
              deriv = calc_deriv(y = sm_meas, x = Time, window_width_n = 7),
              deriv_percap = calc_deriv(y = Measurements, x = Time,
                                        window_width_n = 7,
                                        percapita = TRUE, trans_y = 'log',
                                        blank = 0))

#Assess noise
ggplot(data = dat, aes(x = Time, y = deriv,
                       color = Pop, group = paste(file, Well))) +
  geom_line() +
  facet_wrap(~ Well)

ggplot(data = dat, aes(x = Time, y = deriv_percap,
                       color = Pop, group = paste(file, Well))) +
  geom_line() +
  facet_wrap(~ Well)

#Summarize ----
dat_sum <- summarize(
  group_by(dat, file, Well, Pop, Treat, bact_treat, phage_treat, Timepoint, Isol),
  max_percap = max_gc(deriv_percap),
  max_dens = max_gc(sm_meas),
  lag = lag_time(x = Time, y = sm_meas, deriv = deriv_percap),
  auc = auc(x = Time, y = sm_meas))
  
ggplot(data = dat_sum,
       aes(x = Timepoint, y = max_dens, color = Treat)) +
  geom_point(size =3, alpha = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(~ Pop) +
  theme_bw()

ggplot(data = dat_sum,
       aes(x = Timepoint, y = max_percap, color = Treat)) +
  geom_point(size =3, alpha = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(~ Pop) +
  theme_bw()
       
ggplot(data = dat_sum,
       aes(x = Timepoint, y = lag, color = Treat)) +
  geom_point(size =3, alpha = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(~ Pop) +
  theme_bw()       

ggplot(data = dat_sum,
       aes(x = Timepoint, y = auc, color = Treat)) +
  geom_point(size =3, alpha = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(~ Pop) +
  theme_bw()       
       