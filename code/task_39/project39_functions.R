read_nodes <- function(folder_path, nodes_doc_file) {
  setwd(folder_path)
  
  # Read the file line by line
  lines <- readLines(nodes_doc_file)
  first_four_columns <- list()
  
  for (line in lines) {
    # Split the line using a regular expression that handles different delimiters
    values <- strsplit(line, "\\s+|,|;|:|\\t")[[1]]
    
    # Extract nodeLabel, geographical coordinates and year
    first_value <- values[1]  # The first value is a string
    numeric_values <- as.numeric(values[2:4])  # The next three are numeric
        combined_values <- c(first_value, numeric_values)
    
    # Store these values in the list
    first_four_columns <- append(first_four_columns, list(combined_values))
  }
  nodes_df <- do.call(rbind, first_four_columns)
  nodes_df <- as.data.frame(nodes_df, stringsAsFactors = FALSE)
  nodes_df <- cbind(seq(from = 1, to = nrow(nodes_df)), nodes_df)
  colnames(nodes_df) <- c("nodeID", "nodeLabel", "latitude", "longitude", "year")
  
  return(nodes_df[,1:5])
}

read_lines <- function(folder_path, lines_doc_file) {
  # moving to folder with lines informations and reading of file names having them
  setwd(folder_path)
  lines_table <- read.table(lines_doc_file)[-1,]
  lines_table
  
  # creating the dataframe linking each couple of nodes to the line it belong to, using the single files in the folder
  lines_df <- data.frame(node_from = character(), node_to = character(), line = character())
  
  total_lines <- nrow(lines_table)
  for (l in 1:total_lines) {
    line_name <- lines_table[l, 1]
    line_file <- lines_table[l, 2]
    #cat(line_name, "\t", line_file, "\n")
    current_line <- read.table(line_file)
    current_line <- cbind(current_line, rep(line_name, nrow(current_line)))
    lines_df <- rbind(lines_df, current_line)
  }
  colnames(lines_df) <- c("node_from", "node_to", "line")
  
  return(lines_df)
}

read_topologies <- function(folder_name, edges_file_name, first_year, last_year) {
  #print(paste0("./", city, "/", folder_name))
  setwd(paste0("./", city, "/", folder_name))
  edges_df <- data.frame(nodeID_from = integer(), nodeID_to = integer(), line = character(), year = integer())
  
  nodes <- nodes_df[[city]]
  lines <- lines_df[[city]]
  
  for (year in first_year:last_year) {
    print(year)
    current_topology <- read.table(paste0(edges_file_name, year, "-adjacency.txt"))
    current_topology <- data.frame(node_from = pmin(current_topology[,1], current_topology[,2]),
                                   node_to = pmax(current_topology[,1], current_topology[,2]))
    current_topology <- current_topology[!duplicated(current_topology), ]
    for (i in 1:nrow(current_topology)) {
      nodeID_from <- nodes[nodes$nodeLabel == current_topology[i, 1], 1][1]  #qui potrebbero esserci dei problemi, controllare
      nodeID_to <- nodes[nodes$nodeLabel == current_topology[i, 2], 1][1]
      line <- lines[(lines$node_from == current_topology[i, 1] & lines$node_to == current_topology[i, 2]) | (lines$node_from == current_topology[i, 2] & lines$node_to == current_topology[i, 1]), "line"]
      if (length(line) == 0) line <- NA
      edges_df <- rbind(edges_df, data.frame(nodeID_from = nodeID_from, nodeID_to = nodeID_to, line = line, year = year))
    }
  }
  
  return(edges_df)
}

graph_from_city_year <- function(city, year) {
  my_df <- edges_df[[city]]
  my_df <- my_df[my_df$year == year,]
  g <- graph_from_data_frame(my_df[,1:2], directed = FALSE)
  E(g)$line <- my_df$line
  nod_df <- nodes_df[[city]]
  V(g)$latitude <- as.numeric(nod_df[V(g)$name, "latitude"])
  V(g)$longitude <- as.numeric(nod_df[V(g)$name, "longitude"])
  V(g)$city_name <- nod_df[V(g)$name, "nodeLabel"]
  g$name <- paste(city, "-", year)
  return(g)
}

plot_subway <- function(g, edge_width, display_legend, display_labels) {
  lines <- unique(E(g)$line)
  color_palette <- rainbow(length(lines))
  lay <- as.matrix(cbind(V(g)$latitude, V(g)$longitude))
  edge_colors <- color_palette[match(E(g)$line, lines)]
  # Plot the graph
  plot(g,
       vertex.size = log(degree(g) + 1),
       vertex.label = if (display_labels) V(g)$city_name else NA,
       vertex.color = "#FFFFFFBF",
       vertex.label.cex = 0.8,
       edge.curved = 0,
       edge.color = edge_colors,  # Apply edge colors
       edge.width = edge_width,            # Apply edge width
       layout = lay,
       asp = 1)
  title(main = g$name)
  
  if (display_legend) {
    legend("bottomright",                 # Position of the legend
           legend = lines,             # Labels (line names)
           col = color_palette,        # Colors corresponding to lines
           lty = 1,                    # Line type (solid)
           lwd = edge_width,           # Line width (matches edge width)
           title = "Lines",
           cex = 0.8)            # Title of the legend
  }
}

plot_metrics <- function(data, xlab, ylab, main, leg_coord) {
  
  plot(NULL, xlim = range(data[,2]), ylim = range(data[,3]),
       xlab = xlab, ylab = ylab, main = main)
  
  colors <- rainbow(length(cities))
  
  for (i in seq_along(cities)) {
    city_data <- subset(data, city == cities[i])
    lines(city_data[,2], city_data[,3], type = "o", col = colors[i], 
          pch = 19, cex = 0.8, lwd = 1)
  }
  
  legend(x = leg_coord[1], y = leg_coord[2], legend = cities, col = colors, 
         pch = 19, lty = 1, lwd = 2, cex = 0.7, ncol = 3,title = "City")
  grid()
  
}

# Function to perform node percolation and compute metrics
node_percolation <- function(g) {
  set.seed(2131065)
  
  # Number of nodes in the graph
  num_nodes <- vcount(g)
  
  # Initialize results
  random_removal_sizes <- numeric(num_nodes)
  targeted_removal_sizes <- numeric(num_nodes)
  random_diameters <- numeric(num_nodes)
  targeted_diameters <- numeric(num_nodes)
  
  # Copy the graph for random and targeted removal
  g_random <- g
  g_targeted <- g
  
  # Random removal simulation
  random_order <- sample(V(g)) # Randomly shuffle nodes
  for (i in seq_along(random_order)) {
    # Ensure the node to remove is valid for the current graph
    node_to_remove <- V(g_random)[name == random_order[i]$name]
    g_random <- delete_vertices(g_random, node_to_remove)
    
    # Check if the graph is empty
    if (vcount(g_random) == 0) {
      largest_cc_size <- 0
      network_diameter <- 0
    } else {
      largest_cc_size <- max(components(g_random)$csize, na.rm = TRUE)
      network_diameter <- diameter(g_random, directed = FALSE, unconnected = TRUE)
    }
    random_removal_sizes[i] <- largest_cc_size
    random_diameters[i] <- network_diameter
  }
  
  # Targeted removal simulation
  for (i in 1:num_nodes) {
    # Check if the graph is empty
    if (vcount(g_targeted) == 0) {
      largest_cc_size <- 0
      network_diameter <- 0
    } else {
      # Find the node with the maximum degree in the current graph
      max_degree_node <- V(g_targeted)[which.max(degree(g_targeted))]
      g_targeted <- delete_vertices(g_targeted, max_degree_node)
      
      # Compute the largest connected component size and the network diameter
      largest_cc_size <- max(components(g_targeted)$csize, na.rm = TRUE)
      network_diameter <- diameter(g_targeted, directed = FALSE, unconnected = TRUE)
    }
    targeted_removal_sizes[i] <- largest_cc_size
    targeted_diameters[i] <- network_diameter
  }
  
  # Create a dataframe for results
  percolation_results <- data.frame(
    Removed_Nodes = 1:num_nodes,
    Random_Removal_Size = random_removal_sizes,
    Targeted_Removal_Size = targeted_removal_sizes,
    Random_Diameter = random_diameters,
    Targeted_Diameter = targeted_diameters
  )
  
  # Plot the results
  par(mfrow = c(2, 1))
  
  # Largest Connected Component Size
  plot(percolation_results$Removed_Nodes, percolation_results$Random_Removal_Size, 
       type = "l", col = "blue", lwd = 2,
       xlab = "Number of Removed Nodes", ylab = "LCC Size",
       main = paste("LCC size:", g$name))
  lines(percolation_results$Removed_Nodes, percolation_results$Targeted_Removal_Size, 
        col = "red", lwd = 2)
  legend("topright", legend = c("Random Removal", "Targeted Removal"), 
         col = c("blue", "red"), lty = 1, lwd = 2)
  grid()
  
  # Diameter
  plot(percolation_results$Removed_Nodes, percolation_results$Random_Diameter, 
       type = "l", col = "blue", lwd = 2,
       xlab = "Number of Removed Nodes", ylab = "Network Diameter",
       main = paste("Network Diameter:", g$name))
  lines(percolation_results$Removed_Nodes, percolation_results$Targeted_Diameter, 
        col = "red", lwd = 2)
  legend("topright", legend = c("Random Removal", "Targeted Removal"), 
         col = c("blue", "red"), lty = 1, lwd = 2)
  grid()
  
  return(percolation_results)
}
