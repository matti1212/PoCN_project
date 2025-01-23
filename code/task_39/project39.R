rm(list = ls())
library(igraph)
library(deSolve)
library(R6)
main_folder_path <- "C:/Users/matti/OneDrive/Documenti/PoCN_project/code/task_39" # Please, write here the path to the main folder where all the info about cities are gathered
setwd(main_folder_path)
source("project39_functions.R")

nodes_df <- list() # list containing, for each city, info about nodes
lines_df <- list() # list containing, for each city, info about which couple of nodes belongs to which lines
edges_df <- list() # list containing, for each city, info about edges

#cities <- c("Moscow", "New_tokyo", "NewNYC", "NYC", "Osaka", "Paris", "Seoul", "Shanghai")
cities <- c("Moscow", "New_tokyo", "NewNYC", "Osaka", "Paris", "Seoul", "Shanghai")

#---------- Reading nodes -----------

nodes_doc_files <- list(
  Moscow = "moscow-stations-times-ok.txt",
  New_tokyo = "tokyo-stations-positions-years.txt",
  NewNYC = "nyc-stations-times.txt",
  NYC = "nyc-stations-times.txt",
  Osaka = "osaka-stations-positions-years.txt",
  Paris = "paris-stations-positions-years.txt",
  Seoul = "seoul-stations-times-OK.txt",
  Shanghai = "shanghai-stations-positions-years.txt"
)

for (city in cities) {
  
  nodes_df[[city]] <- read_nodes(folder_path = paste0("./", city), nodes_doc_file = nodes_doc_files[[city]])
  setwd(main_folder_path)
}

#---------- Reading lines -----------

lines_doc_files <- list(
  Moscow = "moscow-lines.txt",
  New_tokyo = "tokyo-lignes.txt",
  NewNYC = "nyc-lines.txt",
  NYC = "nyc-lignes.txt",
  Osaka = "osaka-lines.txt",
  Paris = "paris-lignes.txt",
  Seoul = "seoul-lignes.txt",
  Shanghai = "shanghai-lines.txt"
)

for (city in cities) {
  
  lines_df[[city]] <- read_lines(folder_path = paste0("./", city), lines_doc_file = lines_doc_files[[city]])
  setwd(main_folder_path)
}

#---------- Reading topologies (a.k.a. edges) across years -----------

topologies_info <- data.frame( folder = c("Moscow", "Tokyo", "nyc", "osaka", "Paris", "Seoul", "Shanghai"),
                               first_year = c(1936, 1928, 1880, 1934, 1901, 1975, 1996),
                               last_year = c(2009, 2009, 2009, 2010, 2009, 2009, 2010))
topologies_info$folder <- paste0(topologies_info$folder, "-topologies")
rownames(topologies_info) <- cities
topologies_info$edges_file_name <- paste0(cities, "-")
topologies_info["New_tokyo",]$edges_file_name <- "Tokyo-"
topologies_info["NewNYC",]$edges_file_name <- "NewYork-"

load("edges_df.RData") # the computation of adjacency list is very time consuming, better loading the presaved 
                       # object using this line

# for (city in cities) {
#    print(city)
#    edges_df[[city]] <- read_topologies(folder_name = topologies_info[city, ]$folder,
#                                        edges_file_name = topologies_info[city, ]$edges_file_name,
#                                        first_year = topologies_info[city, ]$first_year,
#                                        last_year = topologies_info[city, ]$last_year)
#    setwd(main_folder_path)
# }
# 
# save(edges_df, file = "edges_df.RData")

#----------- Plotting some subways architectures across years --------------
#source("project39_functions.R")
par(mfrow = c(1,1))
g <- graph_from_city_year("New_tokyo", 1948)
plot_subway(g, edge_width = 3, display_legend = TRUE, display_labels = TRUE)

#----------- Some metrics analysis (structural, centrality measures, mixing patterns ...) ----------

# Initialize an empty data frame to store results
metric_df <- data.frame(city = character(),
                        year = integer(),
                        avg_degree = numeric(),
                        avg_betw = numeric(),
                        nodes_count = numeric(),
                        edges_count = numeric(),
                        density = numeric(),
                        diameter = numeric(),
                        avg_path_length = numeric(),
                        clustering_coeff = numeric()) 

# Loop over each city and compute the average degree per year
for (city in cities) {
  years <- unique(edges_df[[city]]$year)  # Get the unique years for the city
  cat("\n", city, " ")
  for (year in years) {
    cat(year, " ")
    g <- graph_from_city_year(city, year)  # Get the graph for the city and year
    m <- length(E(g))
    N <- length(V(g))
    metric_df <- rbind(metric_df, 
                       data.frame(city = city, year = year, 
                                  avg_degree = mean(degree(g)),
                                  avg_betw = mean(betweenness(g)),
                                  nodes_count = N,
                                  edges_count = m,
                                  density = 2*m/(N*(N-1))*100,
                                  diameter = diameter(g, directed = FALSE),
                                  avg_path_length = mean_distance(g, directed = FALSE, unconnected = TRUE),
                                  clustering_coeff = transitivity(g, type = "average")))
  }
}

par(mfrow = c(1,1))
                             # c(division, x_axis, y_axis)
plot_metrics(data = metric_df[,c(1, 2, 3)], xlab = "Year", ylab = "Average degree", main = "Average degree of subways across years", leg_coord = c(1945, 3.8))
plot_metrics(data = metric_df[,c(1, 2, 4)], xlab = "Year", ylab = "Average betweenness", main = "Average betweenness of subways across years", leg_coord = c(1880, 4500))
plot_metrics(data = metric_df[,c(1, 2, 5)], xlab = "Year", ylab = "Number of nodes", main = "Subways' number of nodes across years", leg_coord = c(1870, 450))
plot_metrics(data = metric_df[,c(1, 2, 6)], xlab = "Year", ylab = "Number of edges", main = "Subways' number of edges across years", leg_coord = c(1945, 800))
plot_metrics(data = metric_df[,c(1, 2, 7)], xlab = "Year", ylab = "Density (%)", main = "Subways' density percentage across years", leg_coord = c(1945, 50))
plot_metrics(data = metric_df[,c(1, 2, 8)], xlab = "Year", ylab = "Maximum path length ", main = "Diameter of subways across years", leg_coord = c(1960, 10))
plot_metrics(data = metric_df[,c(1, 2, 9)], xlab = "Year", ylab = "Average path length", main = "Average path length of subways across years", leg_coord = c(1960, 4.3))
plot_metrics(data = metric_df[,c(1, 2, 10)], xlab = "Year", ylab = "Average clustering coefficient", main = "Average clustering coefficient of subways across years", leg_coord = c(1880, 0.039))

mix_pattern_df <- data.frame(city = character(),
                        year = integer(),
                        deg_assort = numeric())

for (city in cities) {
  years <- unique(edges_df[[city]]$year)  # Get the unique years for the city
  cat("\n", city, " ")
  for (year in years) {
    cat(year, " ")
    g <- graph_from_city_year(city, year)  # Get the graph for the city and year
    mix_pattern_df <- rbind(mix_pattern_df, 
                       data.frame(city = city, year = year, 
                                  deg_assort = assortativity_degree(g)))
  }
}

plot_metrics(data = mix_pattern_df[,c(1, 2, 3)], xlab = "Year", ylab = "Degree assortativity", main = "Degree assortativity od subways across years", leg_coord = c(1940, -0.29))

#------------- Modularity and community detection ---------

g <- graph_from_city_year(city = "Shanghai", 2009)
par(mfrow = c(2, 3))

# Walktrap Algorithm
communities <- cluster_walktrap(g)
# Plot the graph with community colors
plot(communities, g, vertex.label = NA, vertex.size = 3, vertex.label.cex = 0.8, layout = cbind(V(g)$latitude, V(g)$longitude), main = paste0(g$name, ", Walktrap algorithm"))
modularity(communities)

# Louvein Method
communities <- cluster_louvain(g)
plot(communities, g, vertex.label = NA, vertex.size = 3, vertex.label.cex = 0.8, layout = cbind(V(g)$latitude, V(g)$longitude), main = paste0(g$name, ", Louvein Method"))
modularity(communities)

# Spin-glass-like techinque
communities <- cluster_spinglass(g)
plot(communities, g, vertex.label = NA, vertex.size = 3, vertex.label.cex = 0.8, layout = cbind(V(g)$latitude, V(g)$longitude), main = paste0(g$name, ", Spin-glass"))
modularity(communities)

# Infomap method
communities <- cluster_infomap(g)
plot(communities, g, vertex.label = NA, vertex.size = 3, vertex.label.cex = 0.8, layout = cbind(V(g)$latitude, V(g)$longitude), main = paste0(g$name, ", Infomap"))
modularity(communities)

# Leading eigenvalue
communities <- cluster_leading_eigen(g)
plot(communities, g, vertex.label = NA, vertex.size = 3, vertex.label.cex = 0.8, layout = cbind(V(g)$latitude, V(g)$longitude), main = paste0(g$name, ", Leading Eigenvalue"))
modularity(communities)

plot_subway(g, edge_width = 2, display_legend = FALSE, display_labels = FALSE)

#--------------- Robustness ----------------
source("project39_functions.R")
g <- graph_from_city_year("NewNYC", 2009)
#plot_subway(g, 3, TRUE, FALSE)
df <- node_percolation(g)
par(mfrow = c(1, 1))
#---------- What if...? ----------
# In a post-apocalyptic world, humanity is forced to live underground.
# Little human communities are established in each station of the subways, and the only way 
# people have for moving is using the subway lines. What if some of them are infected by a 
# very contagious disease? The following is a simulation of such an epidemic.
# A possible reinfection scenario is also implemented. 

g <- graph_from_city_year("Osaka", 1980)
M <- length(g)
plot_subway(g, 3, FALSE, TRUE)

N <- round(rnorm(M, 500, 200)) # initial inhabitants of each spatial patch

# Let's start with some infected people in one random spatial patch
I_init <- rep(0, length.out = M)
I_init[sample(1:M, 1)] <- runif(1, min = 1, max = 10)

# No recovered (yet)
R_init <- rep(0, length.out = M)

# Susceptible are the total population per patch minus the infected
S_init <- N - R_init
S_init <- abs(S_init) # to prevent negative random values

# Mobility flow: in this scenario there are no planes, cars or other ways to reach
# different patches rather than the subway tunnels. Hence, the matrix Phi is positive only
# where an edge exists. 
# Get the adjacency matrix
adj_mat <- as_adjacency_matrix(g, sparse = FALSE)

# Generate a random symmetric matrix based on the adjacency matrix
random_mat <- matrix(runif(M * M), ncol = M)
Phi <- adj_mat * (random_mat + t(random_mat)) / 2  # Symmetrize

# Normalize the matrix to make it doubly stochastic
row_sums <- rowSums(Phi)
Phi <- sweep(Phi, 1, row_sums, FUN = "/")  # Normalize rows to sum to 1
col_sums <- colSums(Phi)
Phi <- sweep(Phi, 2, col_sums, FUN = "/")  # Normalize columns to sum to 1

# Define parameters
beta <- 0.5        # Transmission rate
gamma <- 0.3         # Recovery rate
epsilon <- 0.01      # Mobility scaling factor
rho <- 0.1        # Re-susceptibility rate
  
# Combine initial states
init_states <- c(S = S_init, I = I_init, R = R_init)

# ODE function
sir_model <- function(time, state, parameters) {
  S <- state[1:M]
  I <- state[(M+1):(2*M)]
  R <- state[(2*M+1):(3*M)]
  
  dS <- numeric(M)
  dI <- numeric(M)
  dR <- numeric(M)
  
  for (i in 1:M) {
    dS[i] <- -beta * S[i] * I[i] / N[i] + rho * R[i] + (1 / epsilon) * (sum(Phi[, i] * S) - sum(Phi[i, ] * S[i]))
    dI[i] <- beta * S[i] * I[i] / N[i] - gamma * I[i] + (1 / epsilon) * (sum(Phi[, i] * I) - sum(Phi[i, ] * I[i]))
    dR[i] <- gamma * I[i] - rho * R[i] + (1 / epsilon) * (sum(Phi[, i] * R) - sum(Phi[i, ] * R[i]))
  }
  
  return(list(c(dS, dI, dR)))
}

# Time sequence for simulation
times <- seq(0, 100, by = 1)

# Solve the system
result <- ode(y = init_states, times = times, func = sir_model, parms = NULL)

# Convert result to a data frame
result_df <- as.data.frame(result)

# Extract time and SIR columns for a sample of nodes
sample_nodes <- sample(1:M, min(5, M))  # Select up to 10 random nodes for clarity
time <- result_df$time
S_columns <- result_df[, sample_nodes + 1]  # Columns for sampled S variables
I_columns <- result_df[, sample_nodes + M + 1]  # Columns for sampled I variables
R_columns <- result_df[, sample_nodes + 2 * M + 1]  # Columns for sampled R variables

# Combine sampled states for plotting
all_states <- cbind(S_columns, I_columns, R_columns)

# Assign colors to S, I, R populations
colors <- c(rep("yellow", length(sample_nodes)), 
            rep("green", length(sample_nodes)), 
            rep("blue", length(sample_nodes)))

# Create labels for the legend (S, I, R only)
legend_labels <- c("S (Susceptible)", "I (Infected)", "R (Recovered)")

# Plot all sampled states using matplot
matplot(time, all_states, type = "l", lty = 1, lwd = 2, col = colors,
        xlab = "Time", ylab = "Population (5 cities)", main = paste("Mobility-Driven SIR Model with reinfection,",g$name ), xlim = c(-2, max(time) + 15))

# Add a legend with colors specific to S, I, R
legend("topright", legend = legend_labels, col = c("yellow", "green", "blue"), 
       lty = 1, cex = 0.7, ncol = 1)

# Add city names as labels next to the respective lines for all S, I, R populations
for (i in seq_along(sample_nodes)) {
  # Susceptible (S)
  text(x = time[nrow(result_df)], 
       y = S_columns[nrow(result_df), i], 
       labels = V(g)$city_name[sample_nodes[i]], 
       pos = 4, cex = 0.8, col = "yellow")
  
  # Infected (I)
  text(x = time[nrow(result_df)], 
       y = I_columns[nrow(result_df), i], 
       labels = V(g)$city_name[sample_nodes[i]], 
       pos = 4, cex = 0.8, col = "green")
  
  # Recovered (R)
  text(x = time[nrow(result_df)], 
       y = R_columns[nrow(result_df), i], 
       labels = V(g)$city_name[sample_nodes[i]], 
       pos = 4, cex = 0.8, col = "blue")
}
grid()





library(gifski)

time <- result_df$time
S_columns <- result_df[, 2:(M + 1)]       # Columns for S variables
I_columns <- result_df[, (M + 2):(2 * M + 1)]  # Columns for I variables
R_columns <- result_df[, (2 * M + 2):(3 * M + 1)]  # Columns for R variables


# Define the function to generate a frame for each time stamp
generate_frame <- function(time_idx) {
  # Get populations for the current time step
  S <- S_columns[time_idx, ]
  I <- I_columns[time_idx, ]
  R <- R_columns[time_idx, ]
  N <- S + I + R
  max_N <- max(N)
  
  #print(S)
  
  # Define node colors based on population values
  S_colors <- rgb(1, 1, 0, alpha = S/max_N)#log(10*S+1)/log(10*max_N +1 ))  # Shades of yellow
  I_colors <- rgb(0, 1, 0, alpha = I/max_N)#log(10*I+1)/log(10*max_N +1 )) # Shades of green
  R_colors <- rgb(0, 0, 1, alpha = R/max_N)#log(10*R+1)/log(10*max_N +1 )) # Shades of blue
  
  # Open a graphics device for combined plot
  png(sprintf("frame_%03d.png", time_idx), width = 1200, height = 400)
  par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))  # Arrange plots horizontally
  
  layout <- cbind(V(g)$latitude, V(g)$longitude)
  
  # Plot the S graph
  plot(g, layout = layout, 
       vertex.color = S_colors, 
       vertex.label = round(S), 
       main = "Susceptible (S)", 
       vertex.size = (5 * N / max_N) + 5,
       edge.width = 4)
  
  # Plot the I graph
  plot(g, layout = layout, 
       vertex.color = I_colors, 
       vertex.label = round(I), 
       main = "Infected (I)", 
       vertex.size = (5 * N / max_N) + 5,
       edge.width = 4)
  
  # Plot the R graph
  plot(g, layout = layout, 
       vertex.color = R_colors, 
       vertex.label = round(R), 
       main = "Recovered (R)", 
       vertex.size = (5 * N / max_N) + 5,
       edge.width = 4)
  
  # Close the graphics device
  dev.off()
}

# Generate frames for all time steps
for (time_idx in 1:nrow(result_df)) {
  generate_frame(time_idx)
}

# Create GIF from the frames
png_files <- sprintf("frame_%03d.png", 1:nrow(result_df))
gifski(png_files, gif_file = paste(g$name, "sir_simulation.gif"), width = 1200, height = 400, delay = 0.2)

# Cleanup temporary frame files
file.remove(png_files)



#-------- files creation -------
# Define the path to the "data" folder two levels back
output_path <- file.path("..", "..", "data")

# Ensure the "data" folder exists
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Combine all dataframes in the nodes_df list into a single dataframe
combined_df <- do.call(rbind, nodes_df)

# Write the combined dataframe to a CSV file in the "data" folder
write.csv(combined_df, file.path(output_path, "nodes.csv"), row.names = FALSE)

# Combine all dataframes in the edges_df list into a single dataframe
combined_edges_df <- do.call(rbind, edges_df)

# Write the combined dataframe to a CSV file in the "data" folder
write.csv(combined_edges_df, file.path(output_path, "edges.csv"), row.names = FALSE)

