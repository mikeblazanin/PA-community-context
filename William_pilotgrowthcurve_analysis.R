library(gcplyr)
library(ggplot2)
library(dplyr)

lag_time <- function(rate, y0, x1, y1, trans_y = "log") {
  if(!canbe.numeric(rate)) {stop("rate must be numeric")}
  if(!canbe.numeric(x1)) {stop("x1 must be numeric")}
  if(!canbe.numeric(y1)) {stop("y1 must be numeric")}
  if(!canbe.numeric(y0)) {stop("y0 must be numeric")}
  rate <- as.numeric(rate)
  y0 <- as.numeric(y0)
  x1 <- as.numeric(x1)
  y1 <- as.numeric(y1)
  if(length(unique(sapply(list(rate, y0, x1, y1), length))) != 1) {
    stop("rate, y0, x1, y1 must all be the same length")}
  
  if(!trans_y %in% c("linear", "log")) {
    stop("trans_y must be one of c('linear', 'log')")}
  if(trans_y == "log") {y0 <- log(y0); y1 <- log(y1)}
  
  return((y0-y1)/rate + x1)
}

canbe.numeric <- function(x) {
  if(is.numeric(x) | 
     tryCatch(
       {
         as.numeric(x)
         #This only gets executed if the line above doesn't warn or error
         TRUE 
       }, 
       warning = function(w) {FALSE},
       error = function(e) {stop(e)}
     )) {
    return(TRUE) 
  } else {
    return(FALSE)
  }
}

dat <- read_wides(files = "William growth curves pilot 2023-02-08.csv",
                  startrow = 41, endrow = 128,
                  startcol = "B")

dat <- trans_wide_to_tidy(dat, id_cols = c("file", "Time", "T 600"))

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
