library(igraph)

# This is the main code where the calculations are performed
#N <- c(1e3, 1e4, 1e5)
P <- cbind(t(seq(0, 0.4, length.out = 70)), t(seq(0.42, 1, length.out = 30)))

# operational distance as a function of p
l <- function(p, N, alpha, n_mean = 15){
  n <- (n_mean * log(N))^alpha
  #print("In l function: ")
  #print(paste(p, n, alpha, -n*log(p), abs(-n * log(p))^(1/alpha)))
  if (is.na((-n * log(p))^(1/alpha))){
    #return((abs(-n * log(p)))^(1/alpha))
    return(0)
  }
  return(((-n * log(p)))^(1/alpha))
}

# function to solve with bisection method
func <- function(p_index, p_ext, N, alpha, n_mean = 15){
  p <- P[p_index]
  #print(paste(p, p_ext))
  #print((diameter_p[p_index] - l(p/p_ext, N, alpha))^2)
  return (diameter_p[p_index] - l(p/p_ext, N, alpha, n_mean))
}

# bisection method implementation
Bisection_method <- function(P, n, alpha, n_mean = 15, thresh = 1e-9, v = FALSE){
  p_ext_of_p <- c()
  
  for (p_index in 1:length(P)) {
    if (p_index < 7) {
      p_ext_of_p[p_index] <- NA
      next
    }
    
    a <- 0.05
    b <- 1.1
    c <- (a+b)/2
    
    counter <- 0 
    while( abs(func(p_index, p_ext = c, n, alpha, n_mean)) > thresh) {
      if (func(p_index, p_ext = c, n, alpha, n_mean)*func(p_index, p_ext = b, n, alpha, n_mean) < 0){
        a <- c
      } else
        b <- c
      
      c <- (a+b)/2
      counter <- counter + 1
      if (counter > 3000) break
    }
    
    p_ext <- c
    p_ext_of_p[p_index] <- p_ext
    if (v) print(paste(p_index, P[p_index], p_ext, func(p_index, p_ext, n, alpha, n_mean)))
  }
  
  return (p_ext_of_p)
}

# function to find the two probabilities setting the range in which the 
# critical one occurs
stable_solutions_search <- function(QBS, p_ext_of_p) {
  max_x <- 0
  pc1 <- p_ext_of_p[which.max(diameter_p)]
  pc2 <- p_ext_of_p[which.max(diameter_p)]
  for (i in order(seq_along(p_ext_of_p), decreasing = TRUE)){
    if (i == 1) next
    x <- p_ext_of_p[i-1] - p_ext_of_p[i]
    if (is.na(x)) next
    if (x > max_x){
      max_x <- x
      if (pc1 > p_ext_of_p[i]) pc1 <- p_ext_of_p[i]
      if (pc2 < p_ext_of_p[i-1]) pc2 <- p_ext_of_p[i-1]
    }
    
  }
  
  return (c(pc1, pc2))
}

plot_p_ext_of_p <- function(P, p_ext_of_p, n, alpha, ens, c = 6) {
  par(mar = c(5, 5, 4, 2) + 0.1)  # Bottom, Left, Top, Right margins
  
  # Generate the plot
  plot(P, p_ext_of_p, type = "o", 
       xlab = expression(p), 
       ylab = expression(p^ext~(p)), 
       main = bquote(p^ext ~ ": D(p) = l(p/p^ext), " * .(ens) * ", n =" ~ .(n) ~ ", alpha =" ~ .(alpha) ~ ", c =" ~ .(c)))
  # Add a reference line (y = p)
  lines(P, P, col = "red", lty = 2) 
  
  # Add grid
  grid()
  
  # Add legend
  legend("topleft", 
         legend = c(expression(p^ext(p)), "y = p"), 
         col = c("black", "red"), 
         lty = c(1, 2), 
         pch = c(1, NA), 
         bg = "white")
}

plot_QBS <- function(p_ext_of_p, results, n, alpha, ens, c = 6) {
  QBS_values <- results[results$n == n, ]$QBS
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(p_ext_of_p, QBS_values, type = "o", 
       xlab = expression(p^ext), 
       ylab = "QBS", 
       main = bquote("QBS for " * .(ens) * ", n =" ~ .(n) ~ ", alpha =" ~ .(alpha) ~ ", c =" ~ .(c)))
  grid()
}

plot_Diameter1 <- function(P, diameter_p, n, alpha, ens, c = 6, P_ext = seq(0.2, 1, by = 0.2)) {
  # Adjust margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Plot Diameter(p) normalized by log(n)
  plot(P, diameter_p / log(n), type = "o", 
       xlab = expression(p), 
       ylab = expression(D(p) / log(n)), 
       main = bquote("Diameter(p)/log(n) for " * .(ens) * ", n =" ~ .(n) ~ ", alpha =" ~ .(alpha) ~ ", c =" ~ .(c)))
  
  # Overlay l(p/p_ext) for each p_ext
  for (p_ext in P_ext) {
    l_values <- sapply(P, function(p) l(p / p_ext, n, alpha))  # Compute l values
    lines(P, l_values / log(n))  # Add normalized l values as solid lines
  }
  
  # Add grid
  grid()
  
  # Add legend
  legend("topright", 
         legend = c(expression(D(p) / log(n)), 
                    expression(l^op~(p / p^ext) ~ "for some values of" ~ p^ext)), 
         col = c("black", "black"), 
         lty = c(1, 1), 
         pch = c(1, NA), 
         bg = "white")
}

plot_QBS_with_criticals <- function(p_ext_of_p, results, n, alpha, ens = "E.R.", c = 6) {
  # Extract QBS values for the given network size
  QBS <- results[results$n == n, ]$QBS
  
  # Find critical points pc1 and pc2
  pcs <- stable_solutions_search(QBS, p_ext_of_p)
  pc1 <- pcs[1]
  pc2 <- pcs[2]
  
  # Adjust margins to ensure labels fit well
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Generate the QBS plot
  plot(p_ext_of_p, QBS, type = "o", 
       xlab = expression(p^ext), 
       ylab = "QBS", 
       main = bquote("QBS with Critical Points, " * .(ens) * ", n =" ~ .(n) ~ ", alpha =" ~ .(alpha) ~ ", c =" ~ .(c)))
  
  # Add vertical lines for pc1 and pc2
  lines(c(pc1, pc1), c(0, max(QBS)), col = "red", lty = 1)
  lines(c(pc2, pc2), c(0, max(QBS)), col = "darkgreen", lty = 1)
  
  # Add grid
  grid()
  
  # Add legend
  legend("bottomright", 
         legend = c("QBS", expression(p[c]^1), expression(p[c]^2)), 
         col = c("black", "red", "darkgreen"), 
         lty = c(1, 1, 1), 
         pch = c(1, NA, NA), 
         bg = "white")
}

plot_Diameter2 <- function(P, diameter_p, n, alpha, ens = "E.R.", c = 6, pc1, pc2) {
  # Adjust margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Plot Diameter(p) normalized by log(n)
  plot(P, diameter_p / log(n), type = "o", 
       xlab = expression(p), 
       ylab = expression(D(p) / log(n)), 
       main = bquote("Diameter(p)/log(n) for " * .(ens) * ", n =" ~ .(n) ~ ", alpha =" ~ .(alpha) ~ ", c =" ~ .(c)))
  
  # Overlay l(p/pc1) as a red line
  l_values <- sapply(P, function(p) l(p / pc1, n, alpha))
  lines(P, l_values / log(n), col = "red", lwd = 2)  # Thicker red line
  
  # Overlay l(p/pc2) as a dark green line
  l_values <- sapply(P, function(p) l(p / pc2, n, alpha))
  lines(P, l_values / log(n), col = "darkgreen", lwd = 2)  # Thicker green line
  
  # Add grid
  grid()
  
  # Add legend
  legend("topright", 
         legend = c(expression(D(p) / log(n)), 
                    expression(l^op~(p / p[c]^1)), 
                    expression(l^op~(p / p[c]^2))), 
         col = c("black", "red", "darkgreen"), 
         lty = c(1, 1, 1), 
         lwd = c(1, 2, 2), 
         pch = c(1, NA, NA), 
         bg = "white")
}

# For the following E.R. networks, we are loading the data from bond percolations, 
# mapping the p_ext to the bond percolation probabilities, and, on the new domain for the
# QBS, we are computing the critical probabilities. This is carried out by varying the number 
# of nodes in the networks and the alpha parameter in {1, 2}. 
# --------- 1e3, alpha = 1 ----------
n <- 1e3
results <- readRDS("project27_N_1000.rds")
alpha <- 1
diameter_p <- results[results$n == n,]$diameter

# Now I want for each p, the value of p_ext such that: D(p) = l(p/p_ext)
p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)

plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "E.R.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS(p_ext_of_p, results, n, alpha, ens = "E.R.")

plot_Diameter1(P, diameter_p, n, alpha, ens = "E.R.")

# --------- 1e3, alpha = 2 --------------
alpha <- 2
p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)

plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "E.R.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS(p_ext_of_p, results, n, alpha, ens = "E.R.")

plot_Diameter1(P, diameter_p, n, alpha, ens = "E.R.")
# --------- 1e4, alpha = 1 ----------
# case alpha = 1
n <- 1e4
results <- readRDS("project27_N_10000.rds")
alpha <- 1
diameter_p <- results[results$n == n,]$diameter

p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "E.R.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "E.R.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter2(P, diameter_p, n, alpha, ens = "E.R.", c = 6, pc1, pc2)
# --------- 1e4, alpha = 2 --------------
alpha <- 2
p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "E.R.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "E.R.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter2(P, diameter_p, n, alpha, ens = "E.R.", c = 6, pc1, pc2)
#----- 5*1e4 alpha = 1---
n <- 5*1e4
results <- readRDS("project27_N_50000.rds")
alpha <- 1
diameter_p <- results[results$n == n,]$diameter

p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "E.R.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "E.R.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter2(P, diameter_p, n, alpha, ens = "E.R.", c = 6, pc1, pc2)

# --------- 5*1e4, alpha = 2 --------------
alpha <- 2

p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "E.R.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "E.R.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter2(P, diameter_p, n, alpha, ens = "E.R.", c = 6, pc1, pc2)

#----- phase diagrams E.R. c = 6 ------
alpha <- 1
n_means <- seq(20, 300, length.out = 150)
n <- 5*1e4
results <- readRDS("project27_N_50000.rds")
diameter_p <- results[results$n == n,]$diameter
QBS <- results[results$n == n,]$QBS

pc1s <- c()
pc2s <- c()

for (n_mean in n_means) {
  
  print(n_mean)
  p_ext_of_p <- Bisection_method(P, n, alpha, n_mean, v = FALSE)
  pcs <- stable_solutions_search(QBS, p_ext_of_p)

  pc1s <- append(pc1s, pcs[1])
  pc2s <- append(pc2s, pcs[2])
  print(paste("n_mean = ", n_mean, " pc1 = ", pcs[1], " pc2 = ", pcs[2]))
}

plot(pc1s, n_means / log(n), type = 'l', col = 'red', xlim = c(0.1, 0.5), 
     ylab = "Mean number of entangled pairs / log(n)", 
     xlab = "Critical probabilities", 
     main = "Phase Diagram: E.R., alpha = 1, n = 5*10^4")
lines(pc2s, n_means / log(n), col = 'darkgreen')
grid()
legend("topright", 
       legend = c(expression(p[c]^1), expression(p[c]^2)), 
       col = c("red", "darkgreen"), 
       lty = 1, 
       bg = "white")

n <- 1e4
results <- readRDS("project27_N_10000.rds")
pc1s <- c()
pc2s <- c()

for (n_mean in n_means) {
  
  print(n_mean)
  p_ext_of_p <- Bisection_method(P, n, alpha, n_mean, v = FALSE)
  pcs <- stable_solutions_search(QBS, p_ext_of_p)
  
  pc1s <- append(pc1s, pcs[1])
  pc2s <- append(pc2s, pcs[2])
  print(paste("n_mean = ", n_mean, " pc1 = ", pcs[1], " pc2 = ", pcs[2]))
}

plot(pc1s, n_means / log(n), type = 'l', col = 'red', xlim = c(0.1, 0.5), 
     ylab = "Mean number of entangled pairs / log(n)", 
     xlab = "Critical probabilities", 
     main = "Phase Diagram: E.R., alpha = 1, n = 10^4")
lines(pc2s, n_means / log(n), col = 'darkgreen')
grid()
legend("topright", 
       legend = c(expression(p[c]^1), expression(p[c]^2)), 
       col = c("red", "darkgreen"), 
       lty = 1, 
       bg = "white")


# n_mean <- 50
# p_ext_of_p <- Bisection_method(P, n, alpha, n_mean, v = TRUE)
# plot(P, p_ext_of_p)
# plot(P, diameter_p)
# for (p_ext in P_ext){
#   l_values <- c()
#   for (p in P){
#     l_values <- append(l_values, l(p/p_ext, n, alpha))
#   }
#   lines(P, l_values/log(n))
# }
# grid()
# pcs <- stable_solutions_search(QBS, p_ext_of_p)
# pc1 <- pcs[1]
# pc2 <- pcs[2]
# plot(p_ext_of_p, QBS, 'o')
# lines(c(pc1, pc1), c(0, max(QBS)), col = 'red')
# lines(c(pc2, pc2), c(0, max(QBS)), col = 'darkgreen')

# Same results, but for the Barabasi-Albert networks.
#------ B.A. 1e3 alpha = 1 --------
# case alpha = 1
n <- 1e3
results <- readRDS("project27pa_N_1000.rds")
alpha <- 1
diameter_p <- results[results$n == n,]$diameter

p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "B.A.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "B.A.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter1(P, diameter_p, n, alpha, ens = "B.A.", c = 6)

#-- ba 1e3 alpha = 2----
alpha <- 2
p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "B.A.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "B.A.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter1(P, diameter_p, n, alpha, ens = "B.A.", c = 6)

#------ ba 1e4 alpha = 1 --------
n <- 1e4
results <- readRDS("project27pa_N_10000.rds")
alpha <- 1
diameter_p <- results[results$n == n,]$diameter
p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "B.A.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "B.A.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter1(P, diameter_p, n, alpha, ens = "B.A.", c = 6)
#-- ba 1e4 alpha = 2---
alpha <- 2
p_ext_of_p <- Bisection_method(P, n, alpha, v = TRUE)
plot_p_ext_of_p(P, p_ext_of_p, n, alpha, ens = "B.A.")

plot(results[results$n == n,]$p, results[results$n == n,]$QBS, 'o')
grid()

plot_QBS_with_criticals(p_ext_of_p, results, n, alpha, ens = "B.A.")

pcs <- stable_solutions_search(QBS, p_ext_of_p)
pc1 <- pcs[1]
pc2 <- pcs[2]

plot_Diameter1(P, diameter_p, n, alpha, ens = "B.A.", c = 6)

#-- critical probability and max diameter -----
res <- data.frame()

n_values <- list(
  "ER_6" = c(316, 1e3, 3162, 1e4, 5*1e4),
  "ER_8" = c(316, 1e3, 3162, 1e4, 31623),
  "BA"   = c(316, 1e3, 3162, 1e4, 31623)
)
c_values <- list(
  "ER_6" = 6,
  "ER_8" = 8,
  "BA"   = 6
)
file_patterns <- list(
  "ER_6" = "project27_N_",
  "ER_8" = "project27c8_N_",
  "BA"   = "project27pa_N_"
)

for (pattern_key in names(file_patterns)) {
  ensemble <- ifelse(grepl("ER", pattern_key), "ER", "BA")
  c <- c_values[[pattern_key]]
  file_pattern <- file_patterns[[pattern_key]]
  
  for (n in n_values[[pattern_key]]) {
    file_name <- paste0(file_pattern, n, ".rds")
    if (!file.exists(file_name)) next
    
    # Load data and compute values
    results <- readRDS(file_name)
    diameter_p <- results[results$n == n, ]$diameter
    if (length(diameter_p) == 0) next # Skip if no data for this `n`
    
    res <- rbind(res, data.frame(
      Ensemble = ensemble, N = n, c = c,
      max_diam = max(diameter_p),
      pc = P[which.max(diameter_p)]
    ))
  }
}

print(res)

# Critical Probability (Pc) vs Network Size
plot(
  log10(res$N), res$pc, type = "n",
  xlab = "log10(Number of Nodes)", ylab = "Critical Probability (Pc)",
  main = "Critical Probability vs Network Size"
)
grid()

ensemble_colors <- c("ER_6" = "blue", "ER_8" = "blue", "BA" = "red")
ensemble_symbols <- c("ER_6" = 15, "ER_8" = 19, "BA" = 17)

for (ensemble in unique(res$Ensemble)) {
  for (avg_degree in unique(res$c)) {
    subset_data <- res[res$Ensemble == ensemble & res$c == avg_degree, ]
    
    key <- if (ensemble == "ER") paste0(ensemble, "_", avg_degree) else ensemble
    
    if (nrow(subset_data) > 0) {
      lines(log10(subset_data$N), subset_data$pc, col = ensemble_colors[key], lwd = 2)
      points(log10(subset_data$N), subset_data$pc, col = ensemble_colors[key], pch = ensemble_symbols[key])
    }
  }
}

legend(
  "topright", legend = c("ER (c=6)", "ER (c=8)", "BA (c = 6)"),
  col = c("blue", "blue", "red"), pch = c(15, 19, 17), lwd = 2, bty = "n"
)

# Maximum Diameter vs Network Size
plot(
  log10(res$N), res$max_diam, type = "n",
  xlab = "log10(Number of Nodes)", ylab = "Maximum Diameter",
  main = "Maximum Diameter vs Network Size"
)
grid()

for (ensemble in unique(res$Ensemble)) {
  for (avg_degree in unique(res$c)) {
    subset_data <- res[res$Ensemble == ensemble & res$c == avg_degree, ]
    
    key <- if (ensemble == "ER") paste0(ensemble, "_", avg_degree) else ensemble
    
    if (nrow(subset_data) > 0) {
      lines(log10(subset_data$N), subset_data$max_diam, col = ensemble_colors[key], lwd = 2)
      points(log10(subset_data$N), subset_data$max_diam, col = ensemble_colors[key], pch = ensemble_symbols[key])
    }
  }
}

legend(
  "topleft", legend = c("ER (c=6)", "ER (c=8)", "BA (c=6)"),
  col = c("blue", "blue", "red"), pch = c(15, 19, 17), lwd = 2, bty = "n"
)

