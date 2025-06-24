
#-------------------------------------------------------------------------------------------
# get confidence intervals on species abundances

# there may be many divergent simulations when calculating the CIs
stats <- abundance_CI(samples, y, b, time, dt, cl)
stats <- cbind(stats, res_long[,2])
names(stats)[6] <- "true"

# visualise simulated data with and without noise
p1 <- ggplot(data = stats, aes(x = time, y = mean, group = species, 
                               colour = species, fill = species)) +
  geom_line(size = 0.7, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_2.5, ymax = upper_2.5), alpha=.4) +
  geom_line(aes(time, true, colour = species), size = 1, linetype = "solid") +
  #scale_color_manual(values = wes_palette("Darjeeling1", n = n_spec)) +
  #scale_fill_manual(values = wes_palette("Darjeeling1", n = n_spec)) +
  labs(title = "gLV Comparison", x = "Time", y = "Abundance") +
  ylim(0, max(stats[6])) +
  theme_classic()
p1
ggsave(paste(pathway, "ts_w_CI.png"), p1, device="png")

#----------------------------------------------------------------------------------------------------------------
# forecasting
# run simulation
extra_time <- 10
res_new <- glv(N = n_spec, full_int_matrix, g_rates[, 1], y, tend = time+extra_time, tstep = dt)

# prepare data
res_new %<>% t() %>% as.data.frame()
names(res_new)      <- paste0("s", seq_len(n_spec))
res_long_new        <- melt(res_new)
res_long_new$t      <- rep(seq(0, time+extra_time, dt), n_spec)
names(res_long_new) <- c("species", "abundance", "time")

# re run original simulation with longer time
res_o <- gLV_sim(n_spec, A, b, y, time+extra_time, dt)

# prepare data
res_long_o        <- melt(res_o)
res_long_o$t      <- rep(seq(0, time+extra_time, dt), n_spec)
names(res_long_o) <- c("species", "abundance", "time")

# append real values from simulation
names(res_long_new)[2] <- "predicted"
res_long_new[4]        <- res_long_o[2]
names(res_long_new)[4] <- "true"

# visualise simulated data with and without noise
p1 <- ggplot(res_long_new) +
  geom_line(aes(time, true, colour = species), size = 1) +
  geom_line(aes(time, predicted, colour = species), size = 1, linetype = "dashed") +
  #scale_color_manual(values = wes_palette("Darjeeling1", n = n_spec)) +
  labs(title = "gLV Comparison", x = "Time", y = "Abundance") +
  geom_vline(xintercept = time) +
  ylim(0, max(stats[6])) +
  theme_classic()
p1
ggsave(paste(pathway, "forecast.png"), p1, device="png")
