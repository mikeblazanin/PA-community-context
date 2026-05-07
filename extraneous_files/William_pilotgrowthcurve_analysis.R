library(gcplyr)
library(ggplot2)
library(dplyr)

setwd("./growth_curves")

dat <- read_wides(files = c("William growth curves pilot 2023-02-08.csv",
                            "2023-02-18 William growth curve.csv"),
                  startrow = c(41, 35), endrow = c(128, 180),
                  startcol = c("B", "A"))

dat <- trans_wide_to_tidy(dat, id_cols = c("file", "Time", "T 600"))

dat <- merge_dfs(dat, collapse = TRUE)
dat2 <- filter(dat, file == "2023-02-18 William growth curve")
dat <- filter(dat, file == "William growth curves pilot 2023-02-08")

designs <- make_design(nrows = 8, ncols = 12,
                       output_format = "tidy",
                       Bacteria = list(c("P6 T0", "CP P3", "CNP P7", "NCP P5", "NN P10"),
                            2:7, 2:11,
                            rep(1:5, each = 12),
                            FALSE),
                       Bacteria = list("Control", 7, 11, "1"),
                       start_od = list(0.1*c(1, 1/2, 1/4, 1/8, 1/16),
                                       2:7, 2:11,
                                       "112233445500",
                                       FALSE),
                       start_od = list(0, 7, 11, "1"))

dat <- merge_dfs(dat, designs)

design2 <- make_design(nrows = 8, ncols = 12,
                       output_format = "tidy",
                       timepoint = make_designpattern(c(1, 2, 4),
                                                      3:5, 2:11,
                                                      "123", FALSE),
                       timepoint = make_designpattern(0, 2, 2:4, "1"),
                       timepoint = list(c(1, 2, 4), "G", 3:10,
                                        "12300123"),
                       pop = make_designpattern(c(3, 4, 5), 3:5, 2:11, 
                                                "1231231231"),
                       pop = list(4, "G", 3:6, "1"),
                       pop = list(5, "G", 8:11, "1"),
                       treat = list(c("CP", "CNP", "NCP", "NN"), c("C", "D", "E"),
                                    2:11, "1112223334"),
                       treat = list("Anc", "B", 2:4, "1"),
                       treat = list("NN", "G", 3:10, "11100111"),
                       treat = list("Ctrl", "B", 11, "1"))

dat2 <- merge_dfs(dat2, design2)
dat2 <- filter(dat2, treat != "NA")


dat <- dplyr::filter(dat, start_od != "NA")

dat$Time <- lubridate::time_length(lubridate::hms(dat$Time))
                            

ggplot(data = dat, aes(x = Time, y = Measurements, 
                       color = start_od, group = Well)) +
  geom_line() +
  facet_wrap(~ Bacteria) +
  scale_y_log10()

min(filter(dat, Bacteria == "Control")$Measurements)

dat <- mutate(group_by(dat, Well, Bacteria, start_od),
              deriv = calc_deriv(y = Measurements, x = Time, x_scale = 3600,
                                 window_width_n = 9),
              deriv_percap = calc_deriv(y = Measurements, x = Time, 
                                        x_scale = 3600, window_width_n = 9,
                                        percapita = TRUE, trans_y = 'log',
                                        blank = 0.119))

ggplot(data = dat, aes(x = Time, y = deriv, 
                       color = start_od, group = Well)) +
  geom_line() +
  facet_wrap(~ Bacteria)

ggplot(data = dat, aes(x = Time, y = deriv_percap, 
                       color = start_od, group = Well)) +
  geom_line() +
  facet_wrap(~ Bacteria)

dat_sum <- summarize(group_by(dat, Well, Bacteria, start_od),
                     max_percap = max(deriv_percap, na.rm = TRUE),
                     max_dens = max(Measurements, na.rm = TRUE),
                     lag = lag_time(rate = max_percap,
                                    y0 = min(Measurements),
                                    y1 = Measurements[which.max(deriv_percap)],
                                    x1 = Time[which.max(deriv_percap)]))

ggplot(dat_sum, aes(x = Bacteria, y = max_percap, color = start_od)) +
  geom_point(size = 2, position = position_dodge(width = 0.4))

ggplot(dat_sum, aes(x = Bacteria, y = max_dens, color = start_od)) +
  geom_point(size = 2, position = position_dodge(width = 0.4))
                     
ggplot(dat_sum, aes(x = Bacteria, y = lag/3600, color = start_od)) +
  geom_point(size = 2, position = position_dodge(width = 0.4))

# dat2 analysis ----
dat2$Time <- lubridate::time_length(lubridate::hms(dat2$Time))
dat2 <- 
  mutate(group_by(dat2, Well, timepoint, pop, treat),
         deriv = calc_deriv(y = Measurements, x = Time, x_scale = 3600,
                            window_width_n = 9),
         deriv_percap = 
           calc_deriv(y = Measurements, x = Time, 
                      x_scale = 3600, window_width_n = 9,
                      percapita = TRUE,
                      blank = 
                        min(filter(dat2, treat == "Ctrl")$Measurements, 
                            na.rm = TRUE)))

ggplot(data = dat2, aes(x = Time, y = Measurements)) +
  geom_line() +
  facet_wrap(~Well)

ggplot(data = dat2, aes(x = Time, y = deriv)) +
  geom_line() +
  facet_wrap(~Well)

ggplot(data = dat2, aes(x = Time, y = deriv_percap)) +
  geom_line() +
  facet_wrap(~Well)

dat2_sum <- summarize(group_by(dat2, Well, timepoint, pop, treat),
                     max_percap = max(deriv_percap, na.rm = TRUE),
                     max_dens = max(Measurements, na.rm = TRUE),
                     lag = lag_time(x = Time, y = Measurements,
                                    deriv = deriv_percap))

ggplot(filter(dat2_sum, treat != "Ctrl"), aes(x = pop, y = lag, color = pop)) +
  geom_point(size = 3) +
  facet_grid(timepoint ~ treat)

ggplot(dat2_sum, aes(x = pop, y = lag, color = pop)) +
  geom_point(size = 2) +
  facet_grid(treat ~ timepoint)

ggplot(data = filter(dat2, treat != "Ctrl"), 
       aes(x = Time, y = deriv_percap, color = pop)) +
  geom_line(lwd = 1.5, alpha = 0.8) +
  facet_grid(timepoint ~ treat)
