#TODO:

## Import libraries ----
library(deSolve)
library(ggplot2)
library(dplyr)
library(gcplyr)
library(tidyr)
library(scales)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

## Define model ----
baranyi_phage <- function(t, state, parameters) {
  state[state < 10**-3] <- 0

  with(as.list(c(state, parameters)), {

    #Note that lag = 1/qr * ln(1 + 1/q0)
    #          q0 = 1/(e^(qr*lag)-1)
    #calc_q0 <- function(qr, lag) {return(1/(exp(qr*lag)-1))}
    
    # Implicit resource (with recycling)
    M = 1 - (S1+I1)/k1 - (S2+I2)/k2
    #N <- k - S1*c1 - I1*c1 - S2*c2 - I2*c2
    
    # Strain 1
    if(lag1 == 0) {
      adj1 <- 1
    } else {
      q01 <- 1/(exp(qr1*lag1)-1)
      q1 <- q01 * exp(qr1 * t)
      adj1 <- q1 / (1 + q1)
    }
    
    inf1 <- a1 * S1 * P
    #S1grow <- r1 * adj1 * S1 * (N / (k * c1))
    S1grow <- r1 * adj1 * S1 * M
    
    dS1 <- S1grow - inf1
    dI1 <- inf1 - I1 / tau1
    
    # Strain 2
    if(lag2 == 0) {
      adj2 <- 1
    } else {
      q02 <- 1/(exp(qr2*lag2)-1)
      q2 <- q02 * exp(qr2 * t)
      adj2 <- q2 / (1 + q2)
    }
    
    inf2 <- a2 * S2 * P
    #S2grow <- r2 * adj2 * S2 * (N / (k * c2))
    S2grow <- r2 * adj2 * S2 * M
    
    dS2 <- S2grow - inf2
    dI2 <- inf2 - I2 / tau2
    
    # Phage dynamics
    dP <- 
      (beta1 * (I1 / tau1)
       + beta2 * (I2 / tau2)
       - inf1 - inf2
       - P * y)
    
    #Note: y ~ 0.01/hr
    #  10.1371/journal.pbio.0040193 
    #  10.1093/evolut/qpag022
    #  10.1371/journal.ppat.1000253
    
    list(c(dS1, dI1, dS2, dI2, dP))
  })
}

#Simple example ----
if(F) {
  parameters <- c(
    # growth
    r1 = 0.75, r2 = 0.75,
    k1 = 1e11, k2 = 1e9,
    
    #lag
    lag1 = 0.9, lag2 = 1.5,
    qr1 = 0.75, qr2 = 0.75,
    
    # infection
    a1 = 1.4e-8,
    a2 = 1.4e-8,
    
    # lysis
    tau1 = 1,
    tau2 = 1,
    
    # burst sizes
    beta1 = 11,
    beta2 = 11,
    
    #phage decay rate
    y = 0.01
  )
  
  state <- c(
    S1 = 1e6,
    I1 = 0,
    S2 = 1,
    I2 = 0,
    P  = 0
  )
  
  times = seq(0, 48, by = 0.01)
  
  out <- ode(y = state, times = times,
             func = baranyi_phage, parms = parameters)
  
  outdf <- tidyr::pivot_longer(as.data.frame(out), cols = -c(time), 
                               names_to = "Pop", values_to = "Density")
  
  ggplot(data = outdf,
         aes(x = time, y = Density, color = Pop)) +
    geom_line(lwd = 1.5, alpha = 0.5) +
    scale_y_log10() +
    theme_bw() +
    scale_color_manual(values = my_cols)
}

#Define simulations function ----
run_sims_ode <- function(par_list,
                         init_state,
                         times,
                         ode_func = baranyi_phage,
                         combinatorial = TRUE,
                         n_passages = 1,
                         dilution_rate = 1,
                         print_info = TRUE) {
  
  # helper: safe ODE wrapper
  myTryCatch <- function(expr) {
    warn <- err <- NULL
    value <- withCallingHandlers(
      tryCatch(expr, error = function(e) {
        err <<- e
        NULL
      }),
      warning = function(w) {
        warn <<- w
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, warning = warn, error = err)
  }
  
  # build parameter grid
  if (combinatorial) {
    param_grid <- expand.grid(par_list, stringsAsFactors = FALSE)
  } else {
    param_grid <- as.data.frame(
      lapply(par_list, function(x) rep_len(x, max(lengths(par_list))))
    )
  }
  
  n_sims <- nrow(param_grid)
  
  if (print_info) {
    message(n_sims, " simulations will be run")
  }
  
  # storage
  results <- vector("list", n_sims)
  failures <- vector("list", 0)
  
  # main loop over parameter sets
  for (i in seq_len(n_sims)) {
    
    pars <- as.list(param_grid[i, ])
    
    current_state <- init_state
    sim_list <- vector("list", n_passages)
    
    time_offset <- 0
    
    for (p in seq_len(n_passages)) {
      
      sim <- myTryCatch(
        ode(y = current_state,
            times = times,
            func = ode_func,
            parms = pars)
      )
      
      if (!is.null(sim$error)) {
        failures[[length(failures) + 1]] <- cbind(run = i, passage = p, param_grid[i, ])
        break
      }
      
      df <- as.data.frame(sim$value)
      
      #adjust time so passages are continuous
      df$time <- df$time + time_offset
      
      #track metadata
      df$run <- i
      df$passage <- p
      
      for (nm in names(pars)) {
        df[[nm]] <- pars[[nm]]
      }
      
      #store
      sim_list[[p]] <- df
      
      #prepare next passage
      final_state <- as.numeric(df[nrow(df), names(current_state)])
      names(final_state) <- names(current_state)
      
      # apply dilution
      current_state <- final_state * dilution_rate
      
      # update time offset
      time_offset <- max(df$time)
      
      #avoid duplicate timepoints
      # if (p > 1) {
      #   sim_list[[p]] <- sim_list[[p]][-1, ]
      # }
    }
    
    # combine passages for this parameter set
    results[[i]] <- do.call(rbind, sim_list)
    
    if (print_info && i %% max(1, floor(n_sims / 10)) == 0) {
      message(round(100 * i / n_sims), "% done")
    }
  }
  
  # combine all simulations
  ybig <- do.call(rbind, results)
  
  yfail <- if (length(failures) > 0) {
    do.call(rbind, failures)
  } else {
    NULL
  }
  
  list(
    ybig = ybig,
    yfail = yfail
  )
}

# Run 1 ----
run1 <- run_sims_ode(
  par_list = list(
    r1 = 0.8,
    r2 = 0.8,
    qr1 = 0.8,
    qr2 = 0.8,
    lag1 = 0.9,
    lag2 = 1.5,
    #lag2 = seq(0.9, 1.5, length.out = 2),
    #Estimated adsorption rate is 8.4e-7/hr
    a1 = signif(10**seq(-10, -12, by = -0.1), 3),
    a2 = 0, #1e-17,
    tau1 = 1,
    tau2 = 1,
    beta1 = 11,
    beta2 = 11,
    k1 = 1.4e10,
    k2 = 1e10,
    #k2 = seq(1, 1.4, length.out = 2)*1e10,
    y = 0),
  init_state <- c(
    S1 = 1e6,
    I1 = 0,
    S2 = 1,
    I2 = 0,
    P  = 1e6),
  times = seq(0, 24, by = 2),
  n_passages = 50, dilution_rate = 0.01
)

r1_out <- tidyr::pivot_longer(run1[[1]], 
                              cols = c(S1, I1, S2, I2, P), 
                              names_to = "Pop", values_to = "Density")

r1_sum <- summarize(group_by(r1_out, across(-c(time, Density, passage))),
                    final_dens = Density[time == max(time)])
                    

ggplot(data = filter(r1_out, Pop %in% c("S1", "S2", "P")),
       aes(x = time, y = log10(Density), color = Pop)) +
  geom_line(lwd = 0.9, alpha = 0.8) +
  coord_cartesian(ylim = c(-1, NA)) +
  #scale_y_log10(limits = c(0.1, NA)) +
  theme_bw() +
  scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
                     values = my_cols[c(1, 4, 5, 2, 7)]) +
  facet_wrap(. ~ a1)

# ggplot(data = r1_out,
#        aes(x = time, y = Density, color = Pop)) +
#   geom_line(lwd = 1.5, alpha = 0.8) +
#   coord_cartesian(ylim = c(0, NA)) +
#   theme_bw() +
#   scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
#                      values = my_cols[c(1, 4, 5, 2, 7)]) +
#   facet_grid(k2 ~ lag2)

ggplot(data = mutate(filter(r1_sum, Pop %in% c("S1", "S2", "P")),
                     final_dens = ifelse(final_dens < 0.01, 0, final_dens)),
       aes(x = log10(a1), y = final_dens, color = paste(lag2, k2))) +
  geom_line() +
  facet_grid(~ Pop) +
  scale_y_log10(limits = c(0.01, NA))


run1 <- run_sims_ode(
  par_list = list(
    r1 = 0.8,
    r2 = 0.8,
    qr1 = 0.8,
    qr2 = 0.8,
    lag1 = 0.9,
    lag2 = 1.5,
    #lag2 = seq(0.9, 1.5, length.out = 2),
    #Estimated adsorption rate is 8.4e-7/hr
    a1 = 0, #signif(10**seq(-10, -12, by = -0.1), 3),
    a2 = 0, #1e-17,
    tau1 = 1,
    tau2 = 1,
    beta1 = 11,
    beta2 = 11,
    k1 = 1.4e10, 
    k2 = 1e9,
    #k2 = seq(1, 1.4, length.out = 2)*1e10,
    y = 0),
  init_state <- c(
    S1 = 1,
    I1 = 0,
    S2 = 1e6,
    I2 = 0,
    P  = 0),
  times = seq(0, 24, by = 2),
  n_passages = 50, dilution_rate = 0.01
)

r1_out <- tidyr::pivot_longer(run1[[1]], 
                              cols = c(S1, I1, S2, I2, P), 
                              names_to = "Pop", values_to = "Density")

r1_sum <- summarize(group_by(r1_out, across(-c(time, Density, passage))),
                    final_dens = Density[time == max(time)])


ggplot(data = r1_out,
       aes(x = time, y = log10(Density), color = Pop)) +
  geom_line(lwd = 1.5, alpha = 0.8) +
  coord_cartesian(ylim = c(-1, NA)) +
  #scale_y_log10(limits = c(0.1, NA)) +
  theme_bw() +
  scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
                     values = my_cols[c(1, 4, 5, 2, 7)]) +
  facet_wrap(. ~ a1)

#Run 2 ----
run2 <- run_sims_ode(
  par_list = list(
    r1 = 0.8,
    r2 = 0.8,
    qr1 = 0.8,
    qr2 = 0.8,
    lag1 = 0.9,
    lag2 = 0.9*seq(1, 1.8, by = 0.1),
    #Estimated adsorption rate is 8.4e-7/hr
    a1 = signif(10**seq(-10.9, -11.7, by = -0.025), 5),
    a2 = 0, #1e-17,
    tau1 = 1,
    tau2 = 1,
    beta1 = 11,
    beta2 = 11,
    k1 = 1.4e10,
    k2 = 1.4*seq(0.5, 1, by = 0.1)*1e10,
    y = 0.01),
  init_state <- c(
    S1 = 1e6,
    I1 = 0,
    S2 = 1,
    I2 = 0,
    P  = 1e6),
  times = seq(0, 24, by = 3),
  n_passages = 50, dilution_rate = 0.01
)

r2_out <- tidyr::pivot_longer(run2[[1]], 
                              cols = c(S1, I1, S2, I2, P), 
                              names_to = "Pop", values_to = "Density")
r2_out <- mutate(ungroup(r2_out),
                 Dens1 = ifelse(Density < 1, 1, Density))

r2_sum <- summarize(group_by(r2_out, across(-c(time, Density, Dens1, passage))),
                    final_dens = Dens1[time == max(time)],
                    avgdens_lastpassage = mean(Dens1[passage == 50]))

r2_sum_sum <- 
  pivot_wider(r2_sum,
              names_from = Pop,
              values_from = c(final_dens, avgdens_lastpassage))
r2_sum_sum <- mutate(ungroup(r2_sum_sum),
                     S1_S2_ratio =
                       avgdens_lastpassage_S1/avgdens_lastpassage_S2,
                     S1_freq = avgdens_lastpassage_S1/
                       (avgdens_lastpassage_S1 + avgdens_lastpassage_S2),
                     S2_freq = avgdens_lastpassage_S2/
                       (avgdens_lastpassage_S1 + avgdens_lastpassage_S2))
                     
ggplot(data = filter(r2_out, lag2 == 1.5, k2 == 1e10),
       aes(x = time, y = log10(Density), color = Pop)) +
  geom_line(lwd = 1.5, alpha = 0.8) +
  coord_cartesian(ylim = c(-1, NA)) +
  #scale_y_log10(limits = c(0.1, NA)) +
  theme_bw() +
  scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
                     values = my_cols[c(1, 4, 5, 2, 7)]) +
  facet_wrap(. ~ a1)

ggplot(data = filter(r2_out, 
                     Pop == "P", 
                     k2 %in% c(min(r2_out$k2), max(r2_out$k2)), 
                     lag2 %in% c(min(r2_out$lag2), max(r2_out$lag2))),
       aes(x = time, y = log10(Density), color = as.factor(log10(a1)))) +
  geom_line() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~lag2*k2)

ggplot(data = filter(r2_out, 
                     Pop == "S1", 
                     k2 %in% c(min(r2_out$k2), max(r2_out$k2)), 
                     lag2 %in% c(min(r2_out$lag2), max(r2_out$lag2))),
       aes(x = time, y = log10(Density), color = as.factor(log10(a1)))) +
  geom_line() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~lag2*k2)

ggplot(data = filter(r2_sum, Pop == "S1"),
       aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(final_dens))) +
  facet_wrap(~round(log10(a1), 2)) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_continuous(name = "log10(Susceptible)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time cost of resistance (%)") +
  #scale_fill_viridis_c(limits = c(0, NA)) +
  NULL

ggplot(data = filter(r2_sum, Pop == "S2"),
       aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(final_dens))) +
  facet_wrap(~round(log10(a1), 3)) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_continuous(name = "log10(Resistant)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time cost of resistance (%)") +
  #scale_fill_viridis_c(limits = c(0, NA)) +
  NULL

ggplot(data = filter(r2_sum, Pop == "P"),
       aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(final_dens))) +
  facet_wrap(~round(log10(a1), 3)) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_continuous(name = "log10(Phage)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time cost of resistance (%)") +
  #scale_fill_viridis_c(limits = c(0, NA)) +
  NULL

ggplot(data = r2_sum_sum,
       aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(S1_S2_ratio))) +
  facet_wrap(~round(log10(a1), 3)) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_viridis_c(name = "log10(S/R)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time cost of resistance (%)") +
  NULL

f5c <- 
  ggplot(data = filter(r2_sum, Pop == "P", log10(a1) > -11.68),
         aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(final_dens))) +
  facet_wrap(~ round(log10(a1), 3), nrow = 4) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_continuous(name = "log10(Phage)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time\ncost of resistance (%)") +
  #scale_fill_viridis_c(limits = c(0, NA)) +
  NULL

f5d <- 
  ggplot(data = filter(r2_sum_sum, log10(a1) > -11.68),
         aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(S1_freq))) +
  facet_wrap(~ round(log10(a1), 3), nrow = 4) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_viridis_c(name = "log10(Sus. freq.)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time\ncost of resistance (%)") +
  NULL

png("./plots/figS_fullmodel.png", width = 20, height = 6, units = "in",
    res = 150)
cowplot::plot_grid(f5c, f5d, labels = c("C", "D"), ncol = 2,
                   align = "v", axis = "lr")
dev.off()

ggplot(data = filter(r2_sum_sum, avgdens_lastpassage_P > 1),
       aes(x = S2_freq, y = log10(avgdens_lastpassage_P))) +
  geom_point(aes(color = round(log10(a1), 3))) +
  geom_smooth() +
  scale_y_continuous(labels = math_format(10^.x),
                     breaks = c(2, 4, 6, 8, 10)) +
  scale_color_continuous(labels = math_format(10^.x),
                         name = "Phage\nadsorption\nrate") +
  labs(subtitle = "Simulation results",
       x = "Resistance frequency",
       y = "Phage population density (pfu/mL)") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16)) +
  NULL

ggplot(data = filter(r2_sum_sum, avgdens_lastpassage_P > 1e2,
                     log10(a1) > -11.52),
       aes(x = 100*(k2/k1-1), y = 100*(lag2/lag1-1))) +
  geom_tile(aes(fill = log10(S1_S2_ratio))) +
  facet_wrap(~round(log10(a1), 2)) +
  scale_x_continuous(transform = "reverse") +
  scale_fill_viridis_c(name = "log10(S/R)") +
  labs(x = "Maximum density cost of resistance (%)",
       y = "Lag time cost of resistance (%)") +
  NULL

ggplot(data = filter(r2_out, 
                     round(log10(a1), 3) == -11.475),
       aes(x = time, y = log10(Density), color = Pop)) +
  geom_line(lwd = 0.75, alpha = 0.8) +
  coord_cartesian(ylim = c(-1, NA)) +
  #scale_y_log10(limits = c(0.1, NA)) +
  theme_bw() +
  scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
                     values = my_cols[c(1, 4, 5, 2, 7)]) +
  facet_grid(lag2 ~ k2)

