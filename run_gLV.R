
# get simulated data ---------------------------------------------------------------------------------------

# run simulation
res <- gLV_sim(n_spec, A, b, y, time, dt)

# prepare data
res_long        <- melt(res)
res_long$t      <- rep(seq(0, time, dt), n_spec)
names(res_long) <- c("species", "abundance", "time")

# visualise simulated data with and without noise
# now it is only without noise
p <- ggplot(res_long) +
     geom_line(aes(time, abundance, colour = species), size = 1) +
     labs(title = "gLV Simulation") +
     theme_classic()
p
ggsave(paste(pathway, "trajectories.png"), p, device="png")

# run Stan model -------------------------------------------------------------------------------------------

# add noise to data
for (i in 1:dim(res)[2])
  res[, i] <- res[, i] + rnorm(dim(res)[1], mean = 0, sd = 0.01)
# need positive values
res <- abs(res)

# compile Stan program
model <- stan_model("R/model/gLVModel.stan")

# feed to Stan model
samples <- sampling(model, 
                    data = list(N = dim(res)[2], 
                                K = dim(res)[1], 
                                x = res, 
                                t = as.numeric(rownames(res)),
                                n_cl = n_clust,
                                cl_gr = cl_assignments,
                                nu_a = nu_a,
                                nu_b = nu_b,
                                nu_w = nu_w),      
                    iter = n_iter, chains = n_chains)

# Get posterior summary of f
summary <- summary(samples)$summary
