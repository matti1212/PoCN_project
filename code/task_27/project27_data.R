# This is the source code to create data from bond percolation of Erdos-Renyi networks and
# Barabasi-Albert networks. For n > 1e4, it takes some hours to compute the quantities needed 
# for the analysis, hence the running of this script can be avoided because the data collected 
# are saved in files loaded in the main script.

library(igraph)

#function to perform bond percolation with probability of removing an edge = 1-p
#the computed metrics are the diameter and the quantum backbone size, 
#using dynamically updated averages of 15 calculations.
bond_percolation <- function(g, p){
  d <- 0
  QBS <-0
  
  start.time <- Sys.time()
  
  for (n in 1:15){
    keep_edges <- E(g)[runif(ecount(g), 0, 1) < p]
    g_perc <- subgraph.edges(g, keep_edges, delete.vertices = FALSE)
    
    d <- d + (diameter(g_perc, directed = FALSE) - d)/n # diameter: average updated dynamically 
    QBS <- QBS + (max(components(g_perc)$csize) - QBS)/n # Quantum Backbone Size: average updated dynamically
  }
  
  print(paste("n = ", vcount(g), ", p = ", p, ", ", Sys.time() - start.time))
  
  return(list(d, QBS))
}

# E.R. networks: results of bond percolations with avg. degree = 6
#N <- c(1e3, 1e4, 1e5)
N <- round(c(10^2.5, 10^3.5, 10^4.5))

# For all networks, this discretizion of the bond probabilities has been used, in 
# order to have more data around the most interesting points and to speed up the process.
P <- cbind(t(seq(0, 0.4, length.out = 70)), t(seq(0.42, 1, length.out = 30)))
diam <- c()
QBS <- c()
#alpha <- 1
c <- 6

results <- data.frame()

for (n in N){
  for (p in P) {
    
    g <- erdos.renyi.game(n, p.or.m = c/(n-1)) # average degree equal to c = 6
    res <- bond_percolation(g, p)
    results <- rbind(results, data.frame(
      n = n,
      p = p,
      diameter = res[[1]],
      QBS = res[[2]]
    ))
  }
}

saveRDS(results[results$n == N[1],], paste0("project27_N_", N[1], ".rds"))
saveRDS(results[results$n == N[2],], paste0("project27_N_", N[2], ".rds"))
saveRDS(results[results$n == N[3],], paste0("project27_N_", N[3], ".rds"))

# B.A. networks: results of bond percolation
N <- round(c(10^4.5))
diam <- c()
QBS <- c()
#alpha <- 1
c <- 6

results <- data.frame()

for (n in N){
  for (p in P) {
    
    g <- sample_pa(n, m = c/2) # average degree equal to c = 6
    res <- bond_percolation(g, p)
    results <- rbind(results, data.frame(
      n = n,
      p = p,
      diameter = res[[1]],
      QBS = res[[2]]
    ))
  }
}

saveRDS(results[results$n == N[1],], paste0("project27pa_N_", N[1], ".rds"))
#saveRDS(results[results$n == N[2],], paste0("project27pa_N_", N[2], ".rds"))
#saveRDS(results[results$n == N[3],], paste0("project27_N_", N[3], ".rds"))


# E.R. networks: results of bond percolation with avg. degree c = 8
#N <- c(1e3, 1e4, 1e5)
N <- round(c(10^2.5, 10^3.0, 10^3.5, 10^4.0, 10^4.5))
P <- cbind(t(seq(0, 0.4, length.out = 70)), t(seq(0.42, 1, length.out = 30)))
diam <- c()
QBS <- c()
#alpha <- 1
#c <- 6
c <- 8

results <- data.frame()

for (n in N){
  for (p in P) {
    
    g <- erdos.renyi.game(n, p.or.m = c/(n-1)) # average degree equal to c = 6
    res <- bond_percolation(g, p)
    results <- rbind(results, data.frame(
      n = n,
      p = p,
      diameter = res[[1]],
      QBS = res[[2]]
    ))
  }
}

saveRDS(results[results$n == N[1],], paste0("project27c8_N_", N[1], ".rds"))
saveRDS(results[results$n == N[2],], paste0("project27c8_N_", N[2], ".rds"))
saveRDS(results[results$n == N[3],], paste0("project27c8_N_", N[3], ".rds"))
saveRDS(results[results$n == N[4],], paste0("project27c8_N_", N[4], ".rds"))
saveRDS(results[results$n == N[5],], paste0("project27c8_N_", N[5], ".rds"))
