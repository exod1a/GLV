
# interaction plot ---------------------------------------------------------------------------------------
cl_int_normal <- normalise_matrix(cl_int)
# names for nodes
names <- c()
for (i in 1:n_clust)
  names[i] <- toString(cl[[i]])

# set non-zero values to one
colnames(cl_int_normal) <- rownames(cl_int_normal) <- names

# build the graph object
network <- graph_from_adjacency_matrix(t(cl_int_normal))

# vertex labels i.e. interaction strengths
labels <- melt(cl_int)[[3]]
labels <- labels[labels != 0]

# plot
pdf(file = paste(pathway, "cluster_ints.pdf"))
set.seed(1)
plot(network, edge.curved=.4, vertex.color = "white", vertex.size = 100, edge.color	= "black", vertex.label.cex = 1.6,
     edge.label	= labels, edge.label.cex	= 1.8, edge.label.color = "blue", vertex.label.color = "black")
dev.off()
