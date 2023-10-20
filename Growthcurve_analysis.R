library(gcplyr)
library(ggplot2)
library(dplyr)

setwd("./growth_curves")

#Tidy ----
dat <- read_wides(files = c("William 2023-08-12 PA Community, Pop 3.csv",
                            "William An 2023-10-15 Growth Curve Pop 4.csv"),
                  startrow = 41, startcol = "B")

dat <- trans_wide_to_tidy(dat, id_cols = c("file", "Time", "T° 600"))
dat <- lapply(dat, FUN = function(x) {
  colnames(x) <- gsub("T° 600", "Temp", colnames(x)); return(x)})
dat[[1]]$Pop <- "3"
dat[[2]]$Pop <- "4"

dat <- merge_dfs(dat, collapse = TRUE)

#Designs ----
trt_designs <- make_design(
  nrows = 8, ncol = 12, output_format = "tidy",
  Timepoint = mdp(c(0, 1, 2, 4, 8),
                  2:5, 2:10,
                  "110022222222333333334444444455555555",
                  FALSE),
  Treat = mdp(c("Anc", "NN", "NCP", "CNP", "CP"),
              2:5, 2:10,
              "110022334455223344552233445522334455",
              FALSE),
  Timepoint = mdp("Blank", 7, 11),
  Treat = mdp("Blank", 7, 11))

trt_isols3 <- make_design(
  nrows = 8, ncol = 12, output_format = "tidy",
  Pop = mdp(3, 1:8, 1:12),
  Isol = mdp(1:6, 2:5, 2:10,
             "120012341212121212121212121212121212",
             FALSE),
  Isol = mdp("Blank", 7, 11))
trt_isols4 <- make_design(
  nrows = 8, ncol = 12, output_format = "tidy",
  Pop = mdp(4, 1:8, 1:12),
  Isol = mdp(1:6, 2:5, 2:10,
             "120012121213121212121213121212151212",
             FALSE),
  Isol = mdp("Blank", 7, 11))

trt_isols <- do.call(rbind, list(trt_isols3, trt_isols4))

dat <- merge_dfs(dat, trt_designs)
dat <- merge_dfs(dat, trt_isols, drop = TRUE)


#Time ----
dat$Time <- lubridate::time_length(lubridate::hms(dat$Time))

#Treatments split ----
dat <- mutate(dat,
              bact_treat = 
                ifelse(Treat == "Anc", "Anc",
                       ifelse(Treat %in% c("NCP", "NN"), "NC", "C")),
              phage_treat =
                ifelse(Treat == "Anc", "Anc",
                       ifelse(Treat %in% c("CNP", "NN"), "NP", "P")))

#Smooth ----
dat <- mutate(group_by(dat, file),
              meas_blanked = Measurements - mean(Measurements[Treat == "Blank"]))
dat <- filter(dat, Treat != "Blank")
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
  lag = lag_time(x = Time, y = sm_meas, deriv = deriv_percap))
  
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
       