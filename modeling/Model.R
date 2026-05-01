#TODO:
#rewrite to use lag and k as params, rather than q0 and c

## Import libraries ----
library(deSolve)
library(ggplot2)
library(dplyr)
library(gcplyr)

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
    
    #Note: in absence of phage, this gives equilibria:
    #  S1 = (k - c2 S2)/c1
    #  S2 = (k - c1 S1)/c2
    
    # -------------------------
    # Implicit resource (with recycling)
    # -------------------------
    N <- k - S1*c1 - I1*c1 - S2*c2 - I2*c2
    
    # -------------------------
    # Strain 1
    # -------------------------
    q1 <- q01 * exp(qr1 * t)
    adj1 <- q1 / (1 + q1)
    
    inf1 <- a1 * S1 * P
    S1grow <- r1 * adj1 * S1 * (N / (k * c1))
    
    dS1 <- S1grow - inf1
    dI1 <- inf1 - I1 / tau1
    
    # -------------------------
    # Strain 2
    # -------------------------
    q2 <- q02 * exp(qr2 * t)
    adj2 <- q2 / (1 + q2)
    
    inf2 <- a2 * S2 * P
    S2grow <- r2 * adj2 * S2 * (N / (k * c2))
    
    dS2 <- S2grow - inf2
    dI2 <- inf2 - I2 / tau2
    
    # -------------------------
    # Phage dynamics
    # -------------------------
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

#Simple example
if(F) {
  parameters <- c(
    # growth
    r1 = 0.75, r2 = 0.75,
    k = 10**10,
    #Note: q0 = 1/(e^(qr*lag)-1)
    #1/(exp(0.75*0.9)-1)
    q01 = .1, q02 = .1,
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
    
    # nutrient conversion (nutrients used per cell)
    c1 = 1,
    c2 = 1.1,
    
    #phage decay rate
    y = 0.01
  )
  
  state <- c(
    S1 = 1e6,
    I1 = 0,
    S2 = 1e6,
    I2 = 0,
    P  = 10**4
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

run_sims_ode <- function(par_list,
                         init_state,
                         times,
                         ode_func = baranyi_phage,
                         combinatorial = TRUE,
                         n_passages = 1,
                         dilution_rate = 1,
                         print_info = TRUE) {
  
  # -----------------------------
  # helper: safe ODE wrapper
  # -----------------------------
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
  
  # -----------------------------
  # build parameter grid
  # -----------------------------
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
  
  # -----------------------------
  # storage
  # -----------------------------
  results <- vector("list", n_sims)
  failures <- vector("list", 0)
  
  # -----------------------------
  # main loop over parameter sets
  # -----------------------------
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
      
      # ---- adjust time so passages are continuous ----
      df$time <- df$time + time_offset
      
      # ---- track metadata ----
      df$run <- i
      df$passage <- p
      
      for (nm in names(pars)) {
        df[[nm]] <- pars[[nm]]
      }
      
      # ---- store ----
      sim_list[[p]] <- df
      
      # ---- prepare next passage ----
      final_state <- as.numeric(df[nrow(df), names(current_state)])
      names(final_state) <- names(current_state)
      
      # apply dilution
      current_state <- final_state * dilution_rate
      
      # update time offset
      time_offset <- max(df$time)
      
      # ---- avoid duplicate timepoints ----
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
  
  # -----------------------------
  # combine all simulations
  # -----------------------------
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


#Note that lag = 1/qr * ln(1 + 1/q0)
#          q0 = 1/(e^(qr*lag)-1)
calc_q0 <- function(qr, lag) {return(1/(exp(qr*lag)-1))}

run1 <- run_sims_ode(
  par_list = list(
    r1 = 0.8,
    r2 = 0.8,
    qr1 = 0.8,
    qr2 = 0.8,
    q01 = calc_q0(0.8, 0.9),
    q02 = signif(calc_q0(0.8, seq(0.9, 1.5, length.out = 2)), 2),
    #Estimated adsorption rate is 8.4e-7/hr
    a1 = signif(10**seq(-10, -12, length.out = 5), 3),
    a2 = 0, #1e-17,
    tau1 = 1,
    tau2 = 1,
    beta1 = 11,
    beta2 = 11,
    c1 = 1,
    c2 = seq(1, 1.4, length.out = 2),
    k = 1e10,
    y = 0.01),
  init_state <- c(
    S1 = 1e6,
    I1 = 0,
    S2 = 1,
    I2 = 0,
    P  = 1e6),
  times = seq(0, 24, by = 0.2),
  n_passages = 20, dilution_rate = 0.01
)

r1_out <- tidyr::pivot_longer(run1[[1]], 
                              cols = c(S1, I1, S2, I2, P), 
                              names_to = "Pop", values_to = "Density")

r1_sum <- summarize(group_by(r1_out, across(-c(time, Density, passage))),
                    final_dens = Density[time == max(time)])
                    

ggplot(data = filter(r1_out, Pop %in% c("S1", "S2", "P"),
                     (q02 == 0.43 & c2 == 1.4) | (q02 == 0.95 & c2 == 1),
                     time %% 24 == 0),
       aes(x = time, y = log10(Density), color = Pop)) +
  geom_line(lwd = 1.5, alpha = 0.8) +
  coord_cartesian(ylim = c(-1, NA)) +
  #scale_y_log10(limits = c(0.1, NA)) +
  theme_bw() +
  scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
                     values = my_cols[c(1, 4, 5, 2, 7)]) +
  facet_grid(q02*c2 ~ a1)

ggplot(data = r1_out,
       aes(x = time, y = Density, color = Pop)) +
  geom_line(lwd = 1.5, alpha = 0.8) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw() +
  scale_color_manual(breaks = c("S1", "S2", "I1", "I2", "P"),
                     values = my_cols[c(1, 4, 5, 2, 7)]) +
  facet_grid(c2 ~ q02)

ggplot(data = mutate(filter(r1_sum, Pop %in% c("S1", "S2", "P")),
                     final_dens = ifelse(final_dens < 0.01, 0, final_dens)),
       aes(x = log10(a1), y = final_dens, color = paste(q02, c2))) +
  geom_line() +
  facet_grid(~ Pop) +
  scale_y_log10(limits = c(0.01, NA))


run2 <- run_sims_ode(
  par_list = list(
    r1 = 0.8,
    r2 = 0.8,
    qr1 = 0.8,
    qr2 = 0.8,
    q01 = calc_q0(0.8, 0.9),
    q02 = signif(calc_q0(0.8, seq(0.9, 1.5, length.out = 7)), 2),
    #Estimated adsorption rate is 8.4e-7/hr
    a1 = signif(10**seq(-10, -12, length.out = 9), 3),
    a2 = 0, #1e-17,
    tau1 = 1,
    tau2 = 1,
    beta1 = 11,
    beta2 = 11,
    c1 = 1,
    c2 = seq(1, 1.4, length.out = 7),
    k = 1e10,
    y = 0.01),
  init_state <- c(
    S1 = 1e6,
    I1 = 0,
    S2 = 1,
    I2 = 0,
    P  = 1e6),
  times = seq(0, 24, by = 0.2),
  n_passages = 10, dilution_rate = 0.01
)

r2_out <- tidyr::pivot_longer(run2[[1]], 
                              cols = c(S1, I1, S2, I2, P), 
                              names_to = "Pop", values_to = "Density")

r2_sum <- summarize(group_by(r2_out, across(-c(time, Density, passage))),
                    final_dens = Density[time == max(time)],
                    avgdens_lastpassage = mean(Density[passage == 10]))

r2_sum_sum <- 
  summarize(group_by(r2_sum, across(-c(Pop, final_dens, avgdens_lastpassage))),
            S1_S2_ratio = 
              avgdens_lastpassage[Pop == "S1"]/avgdens_lastpassage[Pop == "S2"])

ggplot(data = mutate(filter(r2_sum, Pop == "P"),
                     ifelse(avgdens_lastpassage < 1, 1, avgdens_lastpassage)),
       aes(x = c2, y = q02)) +
  geom_tile(aes(fill = log10(avgdens_lastpassage))) +
  facet_wrap(~log10(a1)) +
  scale_fill_viridis_c(limits = c(0, NA))

ggplot(data = r2_sum_sum,
       aes(x = c2, y = q02)) +
  geom_tile(aes(fill = log10(S1_S2_ratio))) +
  facet_wrap(~log10(a1)) +
  scale_fill_viridis_c()

