library(ggplot2)
library(dplyr)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Read data
H6growth <- xlsx::read.xlsx("./Katie_data/2021.12.13-H6DataForMike.xlsx",
                            sheetIndex = 1, startRow = 2)
H6InSeq <- xlsx::read.xlsx("./Katie_data/2021.12.13-H6DataForMike.xlsx",
                           sheetIndex = 2, colIndex = 1:6)
H6EOP <- xlsx::read.xlsx("./Katie_data/2021.12.13-H6DataForMike.xlsx",
                         sheetIndex = 3, startRow = 2, colIndex = 1:3)

#Single step phage growth curve ----
H6growth <- tidyr::pivot_longer(H6growth, cols = !Time..min.,
                                names_sep = "_", names_to = c("Treat", "Rep"),
                                values_to = "pfu_ml")
H6growth <- H6growth[complete.cases(H6growth), ]
adsorb_lm <- lm(log10(pfu_ml) ~ Time..min., 
                H6growth[H6growth$Treat %in% c("Filtered") & 
                           H6growth$Time..min. %in% c(0, 10), ])
adsorb_lm

#adsorbtion rate:
2.3/(4.7*10**7) * adsorb_lm$coefficients[2]/log10(exp(1))

H6growth$predict <- NA
H6growth$predict[which(H6growth$Treat == "Filtered")[1:6]] <- 10**predict(adsorb_lm)

png("H6_singlestep_growth.png", width = 4, height = 4, units = "in", res = 150)
ggplot(data = H6growth[H6growth$Treat %in% c("Control.Filtered", "Filtered"), ], 
       aes(x = Time..min., y = pfu_ml,
                            color = Treat, group = paste(Rep, Treat))) +
  geom_line(lwd = 2, alpha = 0.75) +
  #geom_line(aes(y = predict), color = "black") +
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

#Inseq results for phage receptor ----
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
                     labels = c("Sig. hit in\npilus gene", 
                                "Sig. hit in\nnon-pilus gene",
                                "Non-sig. hit")) +
  theme_bw() +
  theme(legend.position = c(0.225, 0.85),
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

#EOP results ----
H6EOP_sum <- dplyr::summarize(dplyr::group_by(H6EOP, Bacteria),
                   mean_titer = mean(H6.Titer),
                   std_err = sd(H6.Titer)/sqrt(n()),
                   n = n())
png("H6_EOP.png", width = 2, height = 4, units = "in", res = 150)
ggplot(data = H6EOP_sum[H6EOP_sum$Bacteria %in% c("PA14", "PA14^pilA"), ], 
                    aes(x = Bacteria, y = mean_titer+1)) +
  geom_col() +
  geom_hline(yintercept = 1, lty = 2) +
  geom_text(
    data = data.frame(Bacteria = "PA14^pilA", mean_titer = 5, label = "bd"),
            aes(label = label), size = 6) +
  geom_errorbar(aes(ymin = mean_titer - std_err,
                    ymax = mean_titer + std_err),
                width = 0.4) +
  scale_y_continuous(trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +  
  labs(y = "Phage plaquing (pfu/mL)") +
  scale_x_discrete(breaks = c("PA14", "PA14^pilA"),
                   labels = c("PA14", "PA14 \u0394pilA")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1)) +
  NULL
dev.off()

#pyocyanin data (visually estimated from Gail's plots since
# I don't have the data itself)
pyo_data_mock <- data.frame(Treat = c("Ctrl", "H6_heat_inac", "H6"),
                            pyocyanin = c(1, 0.98, 0.1),
                            pyocyanin_err = c(NA, 0.02, 0.015))
png("H6_pyocyanin.png", width = 3, height = 4, units = "in", res = 150)
ggplot(data = pyo_data_mock,
       aes(x = Treat, y = pyocyanin)) +
  scale_x_discrete(limits = c("Ctrl", "H6_heat_inac", "H6"),
                   labels = c("Control\n(no phage)", "Heat-\ninactivated\nH6",
                              "H6")) +
  geom_col() +
  geom_errorbar(aes(x = Treat, ymin = pyocyanin - pyocyanin_err,
                    ymax = pyocyanin + pyocyanin_err),
                width = 0.2) +
  geom_hline(yintercept = 1, lty = 2) +
  labs(y = "Pyocyanin Production\nRelative to Control") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1)) +
  NULL
dev.off()

#Hypothetical growth curve data ----
#Define function that calculates derivatives
calc_deriv <- function(density, percapita = FALSE,
                       subset_by = NULL, time = NULL,
                       time_normalize = NULL) {
  #Note! density values must be sorted sequentially into their unique sets already
  
  #Provided a vector of density values, this function returns (by default) the
  # difference between sequential values
  #if percapita = TRUE, the differences of density are divided by density
  #if subset_by is provided, it should be a vector (same length as density),
  # the unique values of which will separate calculations
  #if time_normalize is specified, time should be provided as a simple 
  # numeric (e.g. number of seconds) in some unit
  #Then the difference will be normalized for the time_normalize value
  #(e.g. if time is provided in seconds and the difference per hour is wanted,
  # time_normalize should = 3600)
  
  #Check inputs
  if (!is.null(time) & !is.numeric(time)) {
    stop("time is not numeric")
  }
  if (!is.null(time_normalize)) {
    if (!is.numeric(time_normalize)) {
      stop("time_normalize is not numeric")
    } else if (is.null(time)) {
      stop("time_normalize is specified, but time is not provided")
    }
  }
  
  #Calc derivative
  ans <- c(density[2:length(density)]-density[1:(length(density)-1)])
  #Percapita (if specified)
  if (percapita) {
    ans <- ans/density[1:(length(density)-1)]
  }
  #Time normalize (if specified)
  if (!is.null(time_normalize)) {
    ans <- ans/
      (c(time[2:length(time)]-time[1:(length(time)-1)])/time_normalize)
  }
  #Subset by (if specified)
  if (!is.null(subset_by)) {
    ans[subset_by[2:length(subset_by)] != subset_by[1:(length(subset_by)-1)]] <- NA
  }
  return(c(ans, NA))
}

baranyi_func2 <- function(r, k, v, q0, d0, t_vals) {
  #Copied from Ram et al 2019
  #where m = r
  if (anyNA(c(r, k, v, q0, d0, t_vals))) {return(NA)}
  if (q0 < 0) {q0 <- 0}
  t_vals_hrs <- t_vals/3600
  a <- t_vals_hrs + 1/r*log((exp(-r*t_vals_hrs)+q0)/(1+q0))
  d <- k/(1-(1-((k/d0)**v))*exp(-r*v*a))**(1/v)
  return(d)
}

baranyi_fit2_err <- function(params, t_vals, dens_vals) {
  #params <- c("logk" = ..., "logd0" = ..., "r" = ..., "v" = ...,
  #            "q0" = ...)
  #   (m = r for this one)
  pred_vals <- baranyi_func2(r = params["r"],
                             k = 10**params["logk"],
                             v = params["v"],
                             q0 = params["q0"],
                             d0 = 10**params["logd0"],
                             t_vals = t_vals)
  pred_vals[pred_vals < 0] <- 0
  err <- sum((log10(pred_vals) - log10(dens_vals))**2)
  if (is.infinite(err) | is.na(err)) {return(2*10**300)} else {return(err)}
}

#Simulate data
set.seed(1)
gc_data <- data.frame(time = 3600*seq(from = 0, to = 12, by = 1/2))
gc_data$dens = baranyi_func2(r = 1, k = 10**9, v = 1, q0 = 0.05, d0 = 10**7,
                             t_vals = gc_data$time) +
  10**7 * runif(nrow(gc_data), min = - 0.8, max = 3)

#Fit data
gc_fit <- optim(par = c("logk" = log10(10**9),
                       "logd0" = log10(10**7),
                       "r" = 1,
                       "v" = 1, 
                       "q0" = 0.05),
      fn = baranyi_fit2_err,
      dens_vals = gc_data$dens,
      t_vals = gc_data$time,
      method = "BFGS")

#Calculate lag time
t_vals <- seq(from = 0, to = max(gc_data$time), by = 1)
pred_vals <- baranyi_func2(r = gc_fit$par["r"],
                           k = 10**gc_fit$par["logk"],
                           v = gc_fit$par["v"],
                           q0 = gc_fit$par["q0"],
                           d0 = 10**gc_fit$par["logd0"],
                           t_vals = t_vals)
pred_deriv <- calc_deriv(density = pred_vals)
exp_grow_index <- which.max(pred_deriv)
lagtime <- 
  (((10**gc_fit$par["logd0"]) - pred_vals[exp_grow_index])/pred_deriv[exp_grow_index]) + 
  t_vals[exp_grow_index]

#Plot
png("gcdata.png", width = 4, height = 4, units = "in", res = 150)
ggplot(data = gc_data, aes(x = time, y = dens)) +
  geom_point(color = "black", size = 3) +
  scale_y_continuous(breaks = c(0, 10**8, 10**9),
                     labels = c(0, parse(text = "10^8"), parse(text = "10^9"))) +
  scale_x_continuous(breaks = 3600*seq(0, 12, 3), labels = seq(0, 12, 3)) +
  geom_line(data = data.frame(time = t_vals, dens = pred_vals), 
            lwd = 2, alpha = 0.5, color = my_cols[6]) +
  #geom_vline(xintercept = lagtime) +
  geom_hline(yintercept = 10**gc_fit$par["logd0"], lty = 2) +
  geom_abline(slope = pred_deriv[exp_grow_index],
              intercept = pred_vals[exp_grow_index] - 
                pred_deriv[exp_grow_index]*t_vals[exp_grow_index],
              color = my_cols[5], lwd = 1.5, alpha = 0.5) +
  # geom_point(data = data.frame(time = lagtime, dens = 10**gc_fit$par["logd0"]),
  #            color = my_cols[5], size = 3, pch = 8, stroke = 2, alpha = 0.5) +
  labs(x = "Time (hr)", y = "Density (cfu/mL)") +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 16)) +
  NULL
dev.off()
  
#Genome sequencing results ----
resis_muts <- read.table("./Ben_data/rTIVPH6_only_relevantSigHits.txt",
                         sep = "\t", header = TRUE)
geneprod_lookup <- read.csv("./Ben_data/gene.product_lookup.csv",
                            na.strings = c(""))
#This info on geneprod taken from the Pseudomonas database
#at pseudomonas.com, tho I'm sure other sites have the same info
genome_length <- 6537648

#note that this pivot loses info bc there are some rows with the same
#gene.product name (but obvi dift mutations)
resis_muts_lng <- 
  tidyr::pivot_longer(resis_muts,
                      cols = -gene.product,
                      names_to = c("strain", "mutant", NA),
                      names_sep = "_")
resis_muts_lng$strain <- gsub(resis_muts_lng$strain,
                              pattern = "X", replacement = "")
resis_muts_lng <- full_join(resis_muts_lng, geneprod_lookup)
resis_muts_lng <- tidyr::separate(resis_muts_lng, col = Coordinates,
                  into = c("coord_start", "coord_end", "coord_strand"),
                  sep = "[\\. \\(\\)]+",
                  remove = FALSE)
                  
resis_muts_lng$qntl <-
  ecdf(resis_muts_lng$value)(resis_muts_lng$value)

resis_muts_lng$coord_middle <-
  (as.numeric(resis_muts_lng$coord_start) + as.numeric(resis_muts_lng$coord_end))/2

resis_muts_sum <-
  summarize(group_by(resis_muts_lng, 
                     strain, gene.product, Locus.Tag, Feature.Type, 
                     coord_start, coord_middle, coord_end, 
                     coord_strand, Gene.Name, Product.Description),
            val_mean = mean(value),
            val_logmean = 10**mean(log10(value)),
            val_sd = sd(value),
            val_logsd = 10**sd(log10(value)),
            frac_signif = sum(value < 0.05)/n(),
            n = n())

resis_muts_sum_sum <-
  summarize(group_by(resis_muts_sum, 
                     gene.product, Locus.Tag, Feature.Type, 
                     coord_start, coord_middle, coord_end, 
                     coord_strand, Gene.Name, Product.Description),
            frac_signif_avg = mean(frac_signif),
            frac_signif_se = sd(frac_signif)/sqrt(n()))

#Plots
ggplot(data = resis_muts_lng,
       aes(x = gene.product, y = paste(strain, mutant), 
           color = -qntl)) +
  geom_point() +
  theme(axis.text = element_blank())

ggplot(data = resis_muts_lng,
       aes(x = gene.product, y = qntl, color = strain)) +
  geom_point(alpha = 0.5) +
  theme(axis.text.x = element_blank())

ggplot(data = resis_muts_lng,
       aes(x = gene.product, y = value, color = strain)) +
  geom_point(alpha = 0.5) +
  theme(axis.text.x = element_blank())

ggplot(data = resis_muts_lng[resis_muts_lng$value < 0.05, ],
       aes(x = gene.product, y = value, color = strain)) +
  geom_point(alpha = 0.5) +
  theme(axis.text.x = element_blank())

ggplot(data = resis_muts_sum,
       aes(x = gene.product, y = frac_signif)) +
  geom_point(alpha = 0.5) +
  theme(axis.text.x = element_blank())

png("./Ben_data/frac_signif.png", height = 6, width = 3, units = "in", res = 300)
ggplot(data = resis_muts_sum_sum,
       aes(x = reorder(gene.product, -frac_signif_avg), y = frac_signif_avg)) +
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymin = frac_signif_avg - 1.96*frac_signif_se,
                    ymax = frac_signif_avg + 1.96*frac_signif_se)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0,
                                   size = 5))
dev.off()

ggplot(data = resis_muts_sum_sum,
       aes(x = coord_middle/genome_length)) +
  #coord_polar(theta = "x") +
  geom_dotplot(binwidth = 1/65) +
  theme_bw()

#To make polar dotplot (binned dots), need to do it manually
#So define the bins, then run through & add rows to a data.frame
# for each dot we need (iterating height so dots stack)
nbins <- 120
bins <- data.frame(
  bin_starts <- genome_length * seq(from = 0, to = 1-1/nbins, length.out = nbins),
  bin_cntrs <- 
  genome_length*(seq(from = 0, to = 1-1/nbins, length.out = nbins) + 1/(2*nbins)),
  bin_ends <- genome_length * seq(from = 1/nbins, to = 1, length.out = nbins))
y_start <- 1
y_incr <- 0.04
resis_muts_ss_dotplot <- 
  data.frame(bin_cntr =  0, height = 0, frac_signif_avg = 1)
for (i in 1:nrow(bins)) {
  rows = which(!is.na(resis_muts_sum_sum$coord_middle) &
                 resis_muts_sum_sum$coord_middle >= bins$bin_starts[i] &
                 resis_muts_sum_sum$coord_middle < bins$bin_ends[i])
  rows <- rows[order(-resis_muts_sum_sum$frac_signif_avg[rows])]
  if(length(rows) > 0) {
    heights <- seq(from = y_start, by = y_incr, length.out = length(rows))
    for (j in 1:length(rows)) {
      row <- rows[j]
      resis_muts_ss_dotplot <- 
        rbind(resis_muts_ss_dotplot,
              data.frame(bin_cntr = bins$bin_cntrs[i],
                         height = heights[j],
                         frac_signif_avg = resis_muts_sum_sum$frac_signif_avg[row]))
    }
  }
}

png("./Ben_data/mutation_locations.png", width = 5, height = 5, 
    units = "in", res = 300)
ggplot(data = resis_muts_ss_dotplot,
       aes(x = bin_cntr, y = height, color = -frac_signif_avg)) +
  coord_polar(theta = "x") +
  #scale_color_manual(values = c("black", "white")) +
  geom_point(size = 1.5) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  guides(color = FALSE) +
  NULL
dev.off()

png("./Ben_data/mutation_locations_forscale.png", width = 7, height = 5, 
    units = "in", res = 300)
ggplot(data = resis_muts_ss_dotplot,
       aes(x = bin_cntr, y = height, color = frac_signif_avg)) +
  coord_polar(theta = "x") +
  scale_color_gradient(low = "#56B1F7", high = "#132B43",
                       name = "Frequency of\nResistance\nMutation") +
  geom_point(size = 1.5) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  NULL
dev.off()
