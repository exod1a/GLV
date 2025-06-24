
# Posterior predictive checks
# Species posteriors
params <- extract(samples)
pdf(paste(pathway, "posteriors.pdf"))
par(mfrow=c(ceiling(sqrt(n_spec)), ceiling(sqrt(n_spec))))
for (i in 1:n_spec) {
  plot(density(subset(res_long, species == levels(res_long$species)[i])[, 4]),
       main = paste("Species", i, "Abundance"))
  for (j in 1:10) {
    lines(density(params$x_sim[j, , i]), col="lavenderblush2")
  }
}
dev.off()

# re-fit model with the mean posterior draw as the data
pos_sim <- summary[seq(2*n_spec + n_clust^2-n_clust + 4, dim(summary)[1]-n_clust^2-1, n_spec), 1]
for (i in 2:n_spec) {
  pos_sim <- cbind(pos_sim, summary[seq(2*n_spec + n_clust^2-n_clust + 3 + i, dim(summary)[1]-n_clust^2-1, n_spec), 1])
}

samples_sim <- sampling(model, 
                        data = list(N = dim(pos_sim)[2], 
                                    K = dim(pos_sim)[1], 
                                    x = pos_sim, 
                                    t = as.numeric(rownames(res)),
                                    n_cl = n_clust,
                                    cl_gr = melt(cl)[[2]],
                                    nu_a = nu_a, 
                                    nu_b = nu_b,    
                                    nu_w = nu_w),      
                        iter = n_iter, chains = n_chains)

params_sim <- extract(samples_sim)

# plot inferred growth rates from posterior draw data vs inferred from actual data
pdf(paste(pathway, "inf_gr.pdf"))
par(mfrow=c(ceiling(sqrt(n_spec)), ceiling(sqrt(n_spec))))
for (i in 1:n_spec) {
  l <- lims(density(params$a1[,i]), density(params_sim$a1[,i]))
  plot(density(params$a1[,i]), main = paste("Growth Rate Species", i), 
       xlim=c(l[1], l[2]), ylim=c(l[3], l[4]))
  lines(density(params_sim$a1[,i]), col="red")
}
dev.off()

# self-interactions
pdf(paste(pathway, "inf_self_int.pdf"))
par(mfrow=c(ceiling(sqrt(n_spec)), ceiling(sqrt(n_spec))))
for (i in 1:n_spec) {
  l <- lims(density(params$a2[,i]), density(params_sim$a2[,i]))
  plot(density(params$a2[,i]), main = paste("Self-Interaction Species", i), 
       xlim=c(l[1], l[2]), ylim=c(l[3], l[4]))
  lines(density(params_sim$a2[,i]), col="red")
}
dev.off()

# plot inferred matrix elements from posterior draw data vs inferred from actual data
pdf(paste(pathway, "inf_mat.pdf"))
par(mfrow=c(n_clust-1, n_clust))
for (i in 1:(n_clust^2-n_clust)) {
    l <- lims(density(params$b_raw[,i]), density(params_sim$b_raw[,i]))
    plot(density(params$b_raw[,i]), main = paste("b_raw[",i,"]", sep=""), 
         xlim=c(l[1], l[2]), ylim=c(l[3], l[4]))
    lines(density(params_sim$b_raw[,i]), col="red")
}
dev.off()

# hyperparameters
#plot(density(params$var_a), main = "var_a")
#plot(density(params$var_b), main = "var_b")
#plot(density(params$var_b), main = "var_w")
