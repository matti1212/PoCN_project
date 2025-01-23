rm(list = ls())
library(igraph)
source("project13_functions.R")

#--- fig 2.a: attraction probability vs degree ----
# Creation of a network following model1, N = 1e5, m = 2, gamma = 2.1 -> beta = 1/(2.1-1) = 1/1.1
N <- 1e5
m <- 2
beta <- 1/1.1

#g <- sample_pop_sim_optimized(N, m, beta)
# Time to grow the network: ~ 1h. It has been saved to speed up the process
#write_graph(g, file = "project13_gForFigure2.graphml", format = "graphml")

# Reading the created graph
g <- read_graph("project13_gForFigure2.graphml", format = "graphml")

# The attraction probability is defined as the probability that an existing node
# of degree k attracts a new link: first of all, let's collect the degrees of 
# the nodes to be connected to all the new nodes.

# Radial coordinate of the new node(s)
t <- N + 1
r_t_log <- log(t)

g_cur <- g
degs <- degree(g_cur)

# Update radial coordinates for popularity fading
V(g_cur)$r[1:(t-1)] <- beta * V(g_cur)$r[1:(t-1)] + (1 - beta) * r_t_log

# Vectors to store results for the current graph
degrees_vec <- c()

# Iterate over angular positions for the new node; here 10000
# new nodes have been selected
for (theta_t in seq(0, 2 * pi, length.out = 10000)) {
  theta_s_t <- pi - abs(pi - abs(V(g_cur)$theta[1:(t-1)] - theta_t))
  x_s_t <- log((1:(t-1)) * t * theta_s_t / 2)
  
  # Find nodes to connect
  to_connect <- order(x_s_t)[1:m]
  
  # Save degrees of connected nodes
  degrees_vec <- c(degrees_vec, degs[to_connect])
  #min_distances_vec <- c(min_distances_vec, min(x_s_t[to_connect]))
}

# To turn these degree counts into probailities, we divide each bin of the histogram
# (centered at each degree value) by the number of nodes who have that degree.
# Histogram for conn_degs[[i]]
h <- hist(degrees_vec,
          breaks = seq(min(degrees_vec) - 0.5, max(degrees_vec) + 0.5),
          plot = FALSE)

hx <- h$mids  # Degree bins
hy <- h$counts  # Counts of selected degrees in conn_degs[[i]]

# Number of nodes in the graph with each degree
node_degree_counts <- table(degs)

# Normalize counts by the number of nodes with each degree
normalized_hy <- sapply(hx, function(d) {
  if (as.character(d) %in% names(node_degree_counts)) {
    hy[which(hx == d)] / as.numeric(node_degree_counts[as.character(d)])
  } else {
    0
  }
})

# Normalization
normalized_hy <- normalized_hy / sum(normalized_hy)

# Filter non-zero values, to display in log coordinates
hx <- hx[normalized_hy > 0]
normalized_hy <- normalized_hy[normalized_hy > 0]

# Plot the result
plot(hx, normalized_hy, 'o', 
     xlab = "Degrees", 
     ylab = "Attraction Probability", 
     log = "xy",
     main = "Attraction probability for model1, N = 1e5, m = 2, beta = 1/1.1")
grid()

#----- fig 2.b: Connection probability vs hyperbolic distance----------
# We discretize the hyperbolic distances and, for each bin, we count how many
# nodes within the that distance from each other are connected. Then, in order to
# normalize, we divide each count by the total (possible) number of edges between 
# nodes within the same bin of the hyperbolic distance

# Warning: this whole calculation is very time consuming, the final plot is 
# shown in the report.

# Function to calculate hyperbolic distance (the exact formula)
hyperbolic_distance <- function(r1, theta1, r2, theta2) {
  theta_st <- pi - abs(pi - abs(theta1 - theta2))
  0.5 * acosh(cosh(2 * r1) * cosh(2 * r2) - sinh(2 * r1) * sinh(2 * r2) * cos(theta_st))
}

num_bins <- 50  # Number of bins for distance
bin_edges <- seq(0, max(V(g)$r) * 2, length.out = num_bins + 1)  # Bin edges for distances

bin_counts <- rep(0, num_bins)  # Total pairs in each bin
bin_connected <- rep(0, num_bins)  # Connected pairs in each bin

# Iterate over all pairs of nodes
for (i in 1:(vcount(g) - 1)) {
  for (j in (i + 1):vcount(g)) {
    # Calculate hyperbolic distance
    r1 <- V(g)$r[i]
    theta1 <- V(g)$theta[i]
    r2 <- V(g)$r[j]
    theta2 <- V(g)$theta[j]
    x <- hyperbolic_distance(r1, theta1, r2, theta2)
    
    # Determine the bin for this distance
    bin_index <- findInterval(x, bin_edges, all.inside = TRUE)
    
    # Increment the total pair count for this bin
    bin_counts[bin_index] <- bin_counts[bin_index] + 1
    
    # Check if the nodes are connected and increment connected count if true
    if (are_adjacent(g, i, j)) {
      bin_connected[bin_index] <- bin_connected[bin_index] + 1
    }
  }
}

# Compute connection probabilities
connection_probs <- bin_connected / bin_counts

# Compute bin centers for plotting
bin_centers <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2

bin_centers <- bin_centers[connection_probs > 0]
connection_probs <- connection_probs[connection_probs > 0]

# Plot the results
plot(bin_centers, connection_probs, type = "o", log = "y", 
     xlab = "Hyperbolic Distance", 
     ylab = "Connection Probability", 
     main = "Connection Probability vs Hyperbolic Distance, N = 1e5, m = 2, beta = 1/1.1",
     col = "blue")
grid()

#----- fig. S3.a distribution of nodes degree -------
#the same graph created above is used.

degs = degree(g)
h <- hist(degs, breaks = seq(min(degs)-0.5, max(degs)+0.5), plot = FALSE)
hx <- h$mids[h$density > 0 & h$mids < 80]
hy <- h$density[h$density > 0 & h$mids < 80]

plot(h$mids[h$density > 0], h$density[h$density > 0],
     log = "xy",
     xlab = "degree distribution",
     ylab = "degree",
     main = "Degree distrib., N = 1e5, m = 2, beta = 1/1.1")
grid()
model <- lm(log(hy) ~ log(hx))
# Extract the coefficients from the linear model
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Define the fitted values for the line in original scale
fitted_line <- exp(intercept) * hx^slope

# Add the fitted line to the existing plot
lines(hx, fitted_line, col = "red", lwd = 2)
legend("topright", legend = "Linear fit, ~ k^(-2.4)", col = "red", lty = 1)
slope

#------- fig S3.b: average clustering vs degree ------

# average clustering of nodes having a certain degree

clust_coeffs <- transitivity(g, type = "local") #clustering coeff of each node
avg_degree_clustering <- vector(, length(hx))
for (i in seq_along(hx)) {
  # for every degree value, I take the nodes having that degree and I average over their clustering coefficients
  avg_degree_clustering[i] <- mean(clust_coeffs[degs == hx[i]])
}

plot(hx, avg_degree_clustering, 
     log = "xy",
     xlab = "degree",
     ylab = "average clustering",
     main = "Average clustering distrib., N = 1e5, m = 2, beta = 1/1.1")
grid()

model <- lm(log(avg_degree_clustering) ~ log(hx))
# Extract the coefficients from the linear model
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Define the fitted values for the line in original scale
fitted_line <- exp(intercept) * hx^slope

# Add the fitted line to the existing plot
lines(hx, fitted_line, col = "red", lwd = 2)

legend("topright", legend = "Fit: ~k^(-1.0)", col = "red", lty = 1)

slope

cat("Average clustering: ", transitivity(g, type="average"))

#-------- Fig. S2.b: degree distribution for model 2 -----------
# Degree distribution of a pop. x sim. network following the model2, with N = 1e4 nodes
# m = 3 and gamma = 3

N2 <- 1e4
m2 <- 3
beta2 <- 1/2
T <- 0.5

g2 <- sample_pop_sim_optimized(1e4, m = 3, beta = 1/2, model = "model2", T = 0.5)

h <- hist(degree(g2), breaks = seq(min(degree(g2))-0.5, max(degree(g2))+0.5), plot = FALSE)
hx <- h$mids[h$density > 0 ]
hy <- h$density[h$density > 0 ]

plot(hx, hy, log = "xy", 
     xlab = "degree",
     ylab = "degree distribution",
     main = "model2, N = 1e4, m = 3, beta = 1/2, T = 0.5")
grid()

hy <- hy[hx < 50]
hx <- hx[hx < 50]

model <- lm(log(hy) ~ log(hx))
# Extract the coefficients from the linear model
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Define the fitted values for the line in original scale
fitted_line <- hx^(slope)*exp(intercept)

# Add the fitted line to the existing plot
lines(hx, fitted_line, col = "red", lwd = 1, lty = 'dashed')

legend("topright", legend = "Fit: ~ k^(-2.7)", col = "red", lty = 'dashed')

slope

#---- Fig S2.a: attraction probability for model 2---

# Radial coordinate of the new node(s)
t <- N2 + 1
r_t_log <- log(t)

g_cur <- g2
degs <- degree(g_cur)

# Update radial coordinates for popularity fading
V(g_cur)$r[1:(t-1)] <- beta2 * V(g_cur)$r[1:(t-1)] + (1 - beta2) * r_t_log

# Vectors to store results for the current graph
degrees_vec <- c()
#min_distances_vec <- c()

# Iterate over angular positions for the new node; here 10000
# new nodes have been selected
for (theta_t in seq(0, 2 * pi, length.out = 10000)) {
  not_conn <- V(g_cur)
  
  theta_s_t <- pi - abs(pi - abs(V(g_cur)$theta[1:(t-1)] - theta_t))
  x_s_t <- log((1:(t-1)) * t * theta_s_t / 2)
  
  Rt <- r_t_log - ifelse(T == 0, ifelse(beta2 == 1, log(2*r_t_log/pi/m2), log(2*(1-exp(-(1-beta2)*r_t_log))/(pi*m*(1-beta2))) )
                         , ifelse(beta2 == 1, log(2*T*r_t_log/(sin(T*pi)*m2)), log(2*T*(1-exp(-(1-beta2)*r_t_log))/(sin(T*pi)*m2*(1-beta2)))))
  
  wiring_prob <- if (T > 0) {
    1 / (1 + exp((x_s_t - Rt) / T))
  } else {
    ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0))
  }
  
  to_connect <- sample(not_conn, size = m, replace = FALSE, prob = wiring_prob)
  
  
  # Find nodes to connect
  #to_connect <- order(x_s_t)[1:m]
  
  # Save degrees of connected nodes
  degrees_vec <- c(degrees_vec, degs[to_connect])
  #min_distances_vec <- c(min_distances_vec, min(x_s_t[to_connect]))
}

# To turn these degree counts into probailities, we divide each bin of the histogram
# (centered at each degree value) by the number of nodes who have that degree.
# Histogram for conn_degs[[i]]
h <- hist(degrees_vec,
          breaks = seq(min(degrees_vec) - 0.5, max(degrees_vec) + 0.5),
          plot = FALSE)

hx <- h$mids  # Degree bins
hy <- h$counts  # Counts of selected degrees in conn_degs[[i]]

# Number of nodes in the graph with each degree
node_degree_counts <- table(degs)

# Normalize counts by the number of nodes with each degree
normalized_hy <- sapply(hx, function(d) {
  if (as.character(d) %in% names(node_degree_counts)) {
    hy[which(hx == d)] / as.numeric(node_degree_counts[as.character(d)])
  } else {
    0
  }
})

# Normalization
normalized_hy <- normalized_hy / sum(normalized_hy)

# Filter non-zero values, to display in log coordinates
hx <- hx[normalized_hy > 0]
normalized_hy <- normalized_hy[normalized_hy > 0]

# Plot the result
plot(hx, normalized_hy, 'o', 
     xlab = "Degrees", 
     ylab = "Attraction Probability", 
     log = "xy",
     main = "Attraction probability for model2, N = 1e4, m = 3, beta = 1/2, T = 0.5")
grid()

#------  fig S4: avg clustering and temperature -------
temperatures <- seq(0.01, 0.99, length.out = 11)
gammas <- c(2.1, 2.5, 3.0)

avg_clustering <- array(dim = c(length(gammas), length(temperatures)))

for (i in seq_along(gammas)) {
  for (j in seq_along(temperatures)) {
    
    g2 <- sample_pop_sim_optimized(N = 1e4, m = 3, beta = 1/(gammas[i] - 1), model = "model2", T = temperatures[j])
    
    avg_clustering[i, j] <- transitivity(g2, type = "average")
    
  }
}

plot(temperatures,avg_clustering[1,], 'o', col = "blue", ylim = c(0, 1.0),
     main = "model2, N = 1e4, m = 3", y_lab = "average clustering")
lines(temperatures,avg_clustering[2,], 'o', col = "red")
lines(temperatures,avg_clustering[3,], 'o', col = "darkred")
legend("topright", legend = c("gamma = 2.1", "gamma = 2.5", "gamma = 3.0"),
       col = c("blue", "red", "darkred"),
       lty = c(1,1,1))
grid()

#----- fig. S5: T=0 (that is, model1) has the maximum possible average clustering coefficient ----
gS5 <- list()
gS5[[1]] <- sample_pop_sim_optimized(1000, m = 3, beta = 1/(2.1-1), model = "model1", T = 0)
gS5[[2]] <- sample_pop_sim_optimized(1000, m = 3, beta = 1/(2.5-1), model = "model1", T = 0)
gS5[[3]] <- sample_pop_sim_optimized(1000, m = 3, beta = 1/(3.0-1), model = "model1", T = 0)

avg_clust <- list(c(transitivity(gS5[[1]], "average")), c(transitivity(gS5[[2]], "average")), c(transitivity(gS5[[3]], "average")))

cat("Average clustering for gamma = 2.1: ", avg_clust[[1]])
cat("Average clustering for gamma = 2.5: ", avg_clust[[2]])
cat("Average clustering for gamma = 3.0: ", avg_clust[[3]])

for (i in 1:3) {
  cur_g <- gS5[[i]]  # current graph in the list
  
  for (rewiring in 1:100) {
    print(rewiring)
    repeat {
      # Randomly select two edges
      edges <- sample(E(cur_g), 2)
      edge1 <- ends(cur_g, edges[1])  # nodes of the first edge
      edge2 <- ends(cur_g, edges[2])  # nodes of the second edge
      
      # Original edge connections (A-B and C-D)
      A <- edge1[1]
      B <- edge1[2]
      C <- edge2[1]
      D <- edge2[2]
      
      # Ensure the new edges are distinct from existing edges
      if (!(are_adjacent(cur_g, A, C) || are_adjacent(cur_g, B, D) || A == C || B == D)) {
        # Create a copy of the graph for testing rewiring
        test_g <- cur_g
        
        # Perform the rewiring (remove old edges and add new edges)
        test_g <- delete_edges(test_g, edges)  # Remove edges A-B and C-D
        test_g <- add_edges(test_g, c(A, C, B, D))  # Add edges A-C and B-D
        
        # Check if the average clustering increases
        old_clust <- transitivity(cur_g, "average")
        new_clust <- transitivity(test_g, "average")
        
        if (new_clust > old_clust) {
          # If clustering improved, accept the rewiring and update the graph
          cur_g <- test_g
          break  # exit the repeat loop
        }
      }
      # If conditions are not met, the loop will continue trying another pair of edges
    }
    
    # Every 100 rewirings, store the average clustering coefficient and print progress
    if (rewiring %% 5 == 0) {
      avg_clust[[i]] <- append(avg_clust[[i]], transitivity(cur_g, "average"))
      cat(i, ": rewiring: ", rewiring, ", avg clustering: ", transitivity(cur_g, "average"), "\n")
    }
  }
}

plot(avg_clust[[1]], col = 'blue', type = 'o', 
     ylim = c(0.7,0.9), ylab = "average clustering",
     xlab = "temperature", main = "model2, N = 1e3, m = 3")
lines(avg_clust[[2]], col = 'red', type = 'o')
lines(avg_clust[[3]], col = 'darkred', type = 'o')
legend("topleft", legend = c("gamma = 2.1", "gamma = 2.5", "gamma = 3.0"),
       lty = c(1,1,1), col = c("blue", "red", "darkred"))
grid()


#--------- Fig S6.a: degree distribution---------
N3 <- 1e4
m3 <- 3
beta3 <- 1/1.1
T <- 0.5
#source("project13_functions_new.R")

g3 <- sample_pop_sim_optimized(10000, 3, beta = 1/1.1, model = "model2'", T = 0.5)

#plot(g3, vertex.size = 5*log(degree(g3)), layout = cbind(V(g3)$r*cos(V(g3)$theta), V(g3)$r*sin(V(g3)$theta)))

h <- hist(degree(g3), breaks = seq(min(degree(g3))-0.5, max(degree(g3))+0.5), plot = FALSE)
hx <- h$mids[h$density > 0 & h$mids > 0]
hy <- h$density[h$density > 0 & h$mids > 0]

plot(hx, hy, log = "xy",
     xlab = "degree",
     ylab = "degree distribution",
     main = "model2', N = 1e4, m = 3, beta = 1/1.1, T = 0.5")
grid()

#---- Fig S6.b: attraction probability ------
# Radial coordinate of the new node(s)
t <- N3 + 1
r_t_log <- log(t)

g_cur <- g3
degs <- degree(g_cur)

# Update radial coordinates for popularity fading
V(g_cur)$r[1:(t-1)] <- beta3 * V(g_cur)$r[1:(t-1)] + (1 - beta3) * r_t_log

# Vectors to store results for the current graph
degrees_vec <- c()
#min_distances_vec <- c()

# Iterate over angular positions for the new node; here 10000
# new nodes have been selected
for (theta_t in seq(0, 2 * pi, length.out = 10000)) {
  theta_s_t <- pi - abs(pi - abs(V(g_cur)$theta[1:(t-1)] - theta_t))
  x_s_t <- log((1:(t-1)) * t * theta_s_t / 2)
  
  Rt <- r_t_log - ifelse(T == 0, ifelse(beta3 == 1, log(2*r_t_log/pi/m3), log(2*(1-exp(-(1-beta3)*r_t_log))/(pi*m3*(1-beta3))) )
                         , ifelse(beta3 == 1, log(2*T*r_t_log/(sin(T*pi)*m3)), log(2*T*(1-exp(-(1-beta3)*r_t_log))/(sin(T*pi)*m3*(1-beta3)))))
  
  not_conn <- V(g3) #set of not-yet-connected nodes
  #theta_s_t <- pi - abs(pi - abs(V(g)$theta[1:(t-1)] - theta_t))
  #x_s_t <- log((1:(t-1)) * t * theta_s_t / 2) # the simplified hyperbolic distance is used
  # again, 3 different cases for the wiring probability, if it is computed at the boundaries of its domain. This formula comes from the Supplementary notes
  #wiring_prob <- ifelse(T > 0, 1 / (1 + exp((x_s_t - Rt) / T)), ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0)))
  wiring_prob <- if (T > 0) {
    1 / (1 + exp((x_s_t - Rt) / T))
  } else {
    ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0))
  }
  
  repeat{
    to_connect <- array(dim = t-1)
    i <- 0
    for (s in 1:(t-1)) {
      if ( runif(1) < wiring_prob[s]) {
        i <- i+1
        to_connect[i] <- s
      }
    }
    to_connect <- to_connect[!is.na(to_connect)]
    
    #if(length(to_connect) > 0) break
    if (TRUE) break
  } 
  
  
  # Find nodes to connect
  #to_connect <- order(x_s_t)[1:m]
  
  # Save degrees of connected nodes
  if (length(to_connect) > 0) degrees_vec <- c(degrees_vec, degs[to_connect])
  #min_distances_vec <- c(min_distances_vec, min(x_s_t[to_connect]))
}

# To turn these degree counts into probailities, we divide each bin of the histogram
# (centered at each degree value) by the number of nodes who have that degree.
# Histogram for conn_degs[[i]]
h <- hist(degrees_vec,
          breaks = seq(min(degrees_vec) - 0.5, max(degrees_vec) + 0.5),
          plot = FALSE)

hx <- h$mids  # Degree bins
hy <- h$counts  # Counts of selected degrees in conn_degs[[i]]

# Number of nodes in the graph with each degree
node_degree_counts <- table(degs)

# Normalize counts by the number of nodes with each degree
normalized_hy <- sapply(hx, function(d) {
  if (as.character(d) %in% names(node_degree_counts)) {
    hy[which(hx == d)] / as.numeric(node_degree_counts[as.character(d)])
  } else {
    0
  }
})

# Normalization
normalized_hy <- normalized_hy / sum(normalized_hy)

# Filter non-zero values, to display in log coordinates
hx <- hx[normalized_hy > 0]
normalized_hy <- normalized_hy[normalized_hy > 0]

# Plot the result
plot(hx, normalized_hy, 'o', 
     xlab = "Degrees", 
     ylab = "Attraction Probability", 
     log = "xy",
     main = "Attraction probability for model2', N = 1e4, m = 3, beta = 1/1.1, T = 0.5")
grid()

#------  fig S7: avg clustering (for model 2') and temperature -------
temperatures <- seq(0.01, 0.99, length.out = 11)
gammas <- c(2.1, 2.5, 3.0)

avg_clustering <- array(dim = c(length(gammas), length(temperatures)))

for (i in seq_along(gammas)) {
  for (j in seq_along(temperatures)) {
    
    gs <- sample_pop_sim_optimized(N = 1e4, m = 3, beta = 1/(gammas[i] - 1), model = "model2'", T = temperatures[j])
    
    avg_clustering[i, j] <- transitivity(gs, type = "average")
    
  }
}

plot(temperatures,avg_clustering[1,], 'o', col = "blue", ylim = c(0, 1.0),
     main = "model2', N = 1e4, m = 3", y_lab = "average clustering")
lines(temperatures,avg_clustering[2,], 'o', col = "red")
lines(temperatures,avg_clustering[3,], 'o', col = "darkred")
legend("topright", legend = c("gamma = 2.1", "gamma = 2.5", "gamma = 3.0"),
       col = c("blue", "red", "darkred"),
       lty = c(1,1,1))
grid()