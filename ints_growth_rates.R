# plot interaction matrices and growth rates and save to files

# check interaction matrix
#for (i in 1:n_spec) {
#  elems <- c(elems, seq(n_spec + i, n_spec + n_spec^2, n_spec))
#}

# get learned matrix from Stan
int_matrix <- get_matrix(summary, n_clust, n_spec)

# 95 % confidence interval
lower_interval <- get_matrix(summary, n_clust, n_spec, rownum = 4)
upper_interval <- get_matrix(summary, n_clust, n_spec, rownum = 8) 

# get self interactions
self_int    <- summary[(n_spec+1):(2*n_spec), 1]
self_int_CI <- data.frame("lower" = summary[(n_spec+1):(2*n_spec), 4], "upper" = summary[(n_spec+1):(2*n_spec), 8])

# full interaction matrix (not just clusters)
full_int_matrix <- build_matrix(cl, int_matrix, self_int)

# compare interaction matrices
p1 <- plotmat(full_int_matrix, "Species j", "Species i", "(Predicted) Interaction Matrix")
p2 <- plotmat(A, "Species j", "Species i", "(True) Interaction Matrix")
p <- plot_grid(p1, p2)
pdf(file = paste(pathway, "full_matrices.pdf"))
print(p)
dev.off()

# create names for the matrix graph
names <- matrix(, nrow = n_clust, ncol = n_clust)
for (i in 1:n_clust)
  for (j in 1:n_clust)
    names[i, j] = paste("[", i, ",", j, "]", sep="")
names <- melt(names)[[3]]

diagonals <- c()
for (i in 1:n_clust) 
  diagonals[i] <- (i - 1) * n_clust + i

`%notin%` <- Negate(`%in%`)

p <- list()
for (i in 1:length(melt(int_matrix)[[3]])) local({
  i <- i
  if (i %notin% diagonals) {
    p[[i]] <<- ggplot() + 
      geom_point(aes(0, melt(int_matrix)[[3]][i]), colour = "skyblue3") +
      geom_errorbar(aes(0), ymin=melt(lower_interval)[[3]][i], 
                    ymax=melt(upper_interval)[[3]][i], colour = "skyblue3") +
      geom_point(aes(0, melt(cl_int)[[3]][i]), colour = "black") +
      ylim(min(melt(lower_interval)[[3]], melt(cl_int)[[3]]), 
           max(melt(upper_interval)[[3]], melt(cl_int)[[3]])) +
      labs(title = names[i], x = "", y = "") +
      theme_classic()
  } else {
    
    p[[i]] <<- ggplot() + 
      geom_point(aes(0, 0), size = 4) +
      theme(axis.line = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = "none",
            panel.background = element_rect(fill = 'white'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }
})

p1 <- plot_grid(plotlist = p, byrow = F, ncol = n_clust)
p1
ggsave(paste(pathway, "matrix.png"), p1, device="png")

# find distance between actual and predicted cluster matrices
cl_dist <- norm(cl_int - int_matrix, type = "2")
cl_dist

# check growth rates
g_rates <- as.data.frame(summary[1:n_spec, c(1, 4, 8)])
g_rates$real <- b

p1 <- ggplot(g_rates) +
  geom_point(aes(rownames(g_rates), mean), colour = "skyblue3") +
  geom_errorbar(aes(rownames(g_rates), ymin=g_rates[, 2], ymax=g_rates[, 3]), colour = "skyblue3") +
  geom_point(aes(rownames(g_rates), real), colour = "black") +
  labs(title = "Growth Rates", x = "Species i", y = "Mean") +
  theme_classic()
p1

p1 <- ggsave(paste(pathway, "growth_rates.png"), p1, device="png")

