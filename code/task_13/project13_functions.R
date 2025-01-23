sample_pop_sim_optimized <- function(N, m, beta = 1, model = "model1", T = 0) {
  
  g <- make_empty_graph(n = 1, directed = FALSE)  # Start with 1 node
  V(g)$r <- 0
  V(g)$theta <- runif(1, 0, 2 * pi)
  
  for (t in 2:N) {
    # Create one node with its coordinates r and theta
    r_t <- log(t)
    theta_t <- runif(1, 0, 2 * pi)
    g <- add_vertices(g, 1, r = r_t, theta = theta_t)
    
    r_t_log <- V(g)$r[t]
    
    #popularity fading
    if (beta < 1) {
      V(g)$r[1:(t-1)] <- beta * V(g)$r[1:(t-1)] + (1-beta) * r_t_log
    }
    
    if (model == "model2") {
      
      # to compute Rt we need to take into account the four cases: T=0 and (beta = 1 or beta != 1) or T != 0 and (beta = 1 or beta != 1)
      Rt <- r_t_log - ifelse(T == 0, ifelse(beta == 1, log(2*r_t_log/pi/m), log(2*(1-exp(-(1-beta)*r_t_log))/(pi*m*(1-beta))) )
                             , ifelse(beta == 1, log(2*T*r_t_log/(sin(T*pi)*m)), log(2*T*(1-exp(-(1-beta)*r_t_log))/(sin(T*pi)*m*(1-beta)))))
      
      # Wiring phase for model 2
      if (t > m + 1) {
        
        not_conn <- V(g)[-t] #set of not-yet-connected nodes
        theta_s_t <- pi - abs(pi - abs(V(g)$theta[1:(t-1)] - theta_t))
        x_s_t <- log((1:(t-1)) * t * theta_s_t / 2) # the simplified hyperbolic distance is used
        # 3 different cases for the wiring probability, if it is computed at the boundaries of its domain. This formula comes from the Supplementary notes
        #wiring_prob <- ifelse(T > 0, 1 / (1 + exp((x_s_t - Rt) / T)), ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0)))
        wiring_prob <- if (T > 0) {
          1 / (1 + exp((x_s_t - Rt) / T))
        } else {
          ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0))
        }
        
        # As aforementioned, we use the sample function to extract m nodes following the wiring probability of each node to be sampled. Note: there's no need for probabilities to be normalized to 1
        # Moreover, as an observation, when T = 0 (boundary case), it occurs that all of the wiring probabilities are 0 but just 1 or, however, a number smaller than m. To solve this problem,
        # when T=0, sampling with replacement is allowed, that is, creating more than one edge with the same node. The alternative, could be allowing links with less than m nodes per iteration.
        #to_connect <- ifelse(T==0 && sum(wiring_prob == 1) < m, sample(not_conn, size = m, replace = TRUE, prob = wiring_prob), sample(not_conn, size = m, replace = FALSE, prob = wiring_prob))
        to_connect <- if( T == 0) {
          sample(not_conn, size = m, replace = TRUE, prob = wiring_prob)
        } else {
          sample(not_conn, size = m, replace = FALSE, prob = wiring_prob)
        }
        
        new_edges <- rep(t, 2 * length(to_connect))
        new_edges[seq(2, 2 * length(to_connect), 2)] <- to_connect
        g <- add_edges(g, new_edges, mode = "undirected")
        
      } else {
        
        # If t <= m, connect t to all existing nodes
        to_connect <- 1:(t-1)
        new_edges <- rep(t, 2 * (t-1))
        new_edges[seq(2, 2 * (t-1), 2)] <- to_connect
        g <- add_edges(g, new_edges, mode = "undirected")
      }
      
    } else if (model == "model2'") {
      # to compute Rt we need to take into account the four cases: T=0 and (beta = 1 or beta != 1) or T != 0 and (beta = 1 or beta != 1)
      Rt <- r_t_log - ifelse(T == 0, ifelse(beta == 1, log(2*r_t_log/pi/m), log(2*(1-exp(-(1-beta)*r_t_log))/(pi*m*(1-beta))) )
                             , ifelse(beta == 1, log(2*T*r_t_log/(sin(T*pi)*m)), log(2*T*(1-exp(-(1-beta)*r_t_log))/(sin(T*pi)*m*(1-beta)))))
      
      not_conn <- V(g)[-t] #set of not-yet-connected nodes
      theta_s_t <- pi - abs(pi - abs(V(g)$theta[1:(t-1)] - theta_t))
      x_s_t <- log((1:(t-1)) * t * theta_s_t / 2) # the simplified hyperbolic distance is used
      # again, 3 different cases for the wiring probability, if it is computed at the boundaries of its domain. This formula comes from the Supplementary notes
      #wiring_prob <- ifelse(T > 0, 1 / (1 + exp((x_s_t - Rt) / T)), ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0)))
      wiring_prob <- if (T > 0) {
        1 / (1 + exp((x_s_t - Rt) / T))
      } else {
        ifelse(x_s_t - Rt < 0, 1, ifelse(x_s_t - Rt == 0, 0.5, 0))
      }
      
      # for model 2' we loop once along every node s and, with probability wiring_prob[s], we
      # connect that node to the new one. m new connections are allowed on average: there could
      # be more than or less than m connections at a given loop.
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
        break
      } 
      
      if(length(to_connect) > 0){
      new_edges <- rep(t, 2 * length(to_connect))
      new_edges[seq(2, 2 * length(to_connect), 2)] <- to_connect
      g <- add_edges(g, new_edges, mode = "undirected")
      }
    } else {
      # Model 1: search for m hyperbolically closest nodes
      if (t > m + 1) {
        theta_s_t <- pi - abs(pi - abs(V(g)$theta[1:(t-1)] - theta_t))
        x_s_t <- log((1:(t-1)) * t * theta_s_t / 2)
        
        to_connect <- order(x_s_t)[1:m]
        new_edges <- rep(t, 2 * m)
        new_edges[seq(2, 2 * m, 2)] <- to_connect
        g <- add_edges(g, new_edges, mode = "undirected")
        
      } else {
        to_connect <- 1:(t-1)
        new_edges <- rep(t, 2 * (t-1))
        new_edges[seq(2, 2 * (t-1), 2)] <- to_connect
        g <- add_edges(g, new_edges, mode = "undirected")
      }
    }
  }
  
  return(g)
}

# 
# sample_pop_sim <- function(N, m, beta = 1, model = "model1", T = 0) {
#   
#   g <- make_empty_graph(directed = FALSE)
#   g <- add_vertices(g, 1, r = 0, theta = runif(1, 0, 2*pi))
#   for (t in 2:N) {
#     #cat(t, "\n")
#     # create one node and its coordinate
#     g <- add_vertices(g, 1, r = log(t), theta = runif(1, 0, 2*pi)) 
#     
#     # popularity fading
#     if (beta < 1) {
#       for (s in 1:(t-1)){
#         V(g)$r[s] <- beta * V(g)$r[s] + (1-beta) * V(g)$r[t]
#       }
#     }
#     
#     # wiring method selection
#     if (model == "model2") {
#       
#       # if model 2 is selected, there exists a wiring probability (Fermi-Dirac distributed) of node t.
#       if (beta != 1) {
#         Rt <- V(g)$r[t] - log(2*(1-exp(-(1-beta)*V(g)$r[t]))/(pi*m*(1-beta))) #radius of hyperbolic disc
#       } else Rt <- V(g)$r[t]
#       #Rt not defined in beta = 1, so I take the limit of Rt for beta -> 1
#       #cat("t: ", t, ", r[t]: ", V(g)$r[t], ", Rt: ", Rt, "\n")
#       if (t>m+1) {
#         
#         not_conn <- V(g)[-t]
#         to_connect <- c()
#         while (length(to_connect) < m & length(not_conn) > 0) { 
#           
#           # sample a node s among all nodes not connected to t
#           s <- sample(not_conn, 1)
#           # compute wiring probability
#           theta_s_t <- pi - abs(pi - abs(V(g)$theta[s] - V(g)$theta[t]))
#           #x_s_t <- 0.5 * acosh(cosh(2*V(g)$r[s]) * cosh(2*V(g)$r[t]) - sinh(2*V(g)$r[s]) * sinh(2*V(g)$r[t]) * cos(theta_s_t))
#           x_s_t <- log(s*t*theta_s_t/2)
#           # wiring probability function is not defined for T = 0; instead, limits for T->0 by varying x_s_t - Rt are used
#           if (T > 0) {
#             wiring_prob <- 1/(1+exp((x_s_t - Rt)/T))
#           } else if (x_s_t - Rt < 0) {
#             wiring_prob <- 1
#           } else if (x_s_t - Rt == 0) {
#             wiring_prob <- 0.5
#           } else {
#             wiring_prob <- 0
#           }
#           
#           if (wiring_prob == 0 ) { # se per quel nodo la prob di legarsi a t è nulla, lo tolgo dai possibili nodi da connettere
#             not_conn <- setdiff(not_conn, s)
#             #cat("wiring prob of ", s, " is 0 \n")
#             next
#           }
#           
#           sample_unif <- runif(1, min = 0, max = 1)
#           #cat("s: ", s, ", x_s_t: ", x_s_t, ", x_s_t-Rt: ", x_s_t - Rt,", wiring prob: ", wiring_prob, ", sample_unif: ", sample_unif,"\n")
#           
#           if(sample_unif <= wiring_prob) {
#             not_conn <- setdiff(not_conn, s)
#             to_connect <- c(to_connect, s)
#           }
#           
#           #cat("not_conn: ", not_conn, ", to_connect: ", to_connect, "\n \n")
#           
#         }
#         
#         if (length(to_connect) > 0) {
#           new_edges <- rep(t, 2*length(to_connect))
#           new_edges[seq(2, 2*length(to_connect), 2)] <- to_connect
#         } 
#         
#       } else {
#         
#         #if t < m, connect t to all present nodes
#         to_connect <- 1:(t-1)
#         new_edges <- rep(t, 2*(t-1))
#         new_edges[seq(2, 2*(t-1), 2)] <- to_connect
#       }
#     } else { 
#       # search for the m hyperbolically closest nodes to connect to node t. If t<=m, t is connected to every node alredy in the graph
#       #to_connect <- c()
#       if (t>m+1) {
#         
#         x <- c()
#         for (s in 1:(t-1)){
#           theta_s_t <- pi - abs(pi - abs(V(g)$theta[s] - V(g)$theta[t]))
#           #x[s] <- 0.5 * acosh(cosh(2*V(g)$r[s]) * cosh(2*V(g)$r[t]) - sinh(2*V(g)$r[s]) * sinh(2*V(g)$r[t]) * cos(theta_s_t))
#           x[s] <- log(s*t*theta_s_t/2)
#         }
#         #to_connect <- sort(x, index.return = TRUE)$ix[1:m]
#         to_connect <- order(x)[1:m]
#         new_edges <- rep(t, 2*m)
#         new_edges[seq(2, 2*m, 2)] <- to_connect
#       } else {
#         
#         #if t < m, connect t to all present nodes
#         to_connect <- 1:(t-1)
#         new_edges <- rep(t, 2*(t-1))
#         new_edges[seq(2, 2*(t-1), 2)] <- to_connect
#       }
#     }
#     
#     g <- add_edges(g, new_edges, mode = "undirected")
#   }
#   
#   return (g)
#   
# }
