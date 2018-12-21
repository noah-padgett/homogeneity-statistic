

# The following function is used to find the squared
#   euclidean distance between a persons score vector and the 
#   total_centroids vector found above.
dist_func <- function(
  score, centroid, C, 
  dist.measure = "euclidean", mat.opt = T)
{
  centroid <- ifelse(mat.opt == T, centroid[C,], centroid)
  
  x <- rbind(score, centroid)
  y <- (dist(x, method = dist.measure))**2
  return(y)
}


# ========================================= #
#   EV_Calc()
# ========================================= #
# Inputs (arguments)
# data = raw data frame
# cluster_var = vector of variables in data used for clustering
# cluster_id = variable with class assignments 
# class_centroids = (optional)
# 
# Outputs
# a single matrix(#class + 1, 1)
#   - EV (explained error sum of squares %)
#   - Ec for each cluster 
EV_Calc <- function(
  data, cluster_var, cluster_id, 
  class_centroids = NULL,...)
{
  # data <- mydata
  # cluster_var <- c(paste0("var",1:6))
  # cluster_id <- "k"
  # class_centroids = NULL
  
  dat <- data[,c(cluster_var, cluster_id)]
  nV <- length(cluster_var) # number of cluster indicators
  classSize <- table(dat[,cluster_id])
  nClass <- length(classSize)
  
  
  # Calculate the centroid: mean of each variable
  total_centroid <- apply(dat[,1:nV], 2, mean, na.rm = T)
  
  
  # Calculate the distance for each individual
  total_dists <- apply(dat[,1:nV], 1, dist_func, centroid = total_centroid, mat.opt = F)
  Et <- sum(total_dists, na.rm = T) # SST = sum of squared deviations from the centroid
  
  
  # Check is class centroids were supplied, if there were, they don't need to be computed.
  if(is.null(class_centroids) == T)
  {
    # Compute centroids
    class_centroids <- matrix(0, ncol = nV, nrow = nClass)
    class_centroids <- aggregate(dat[,1:nV], 
                                 by = list(class = dat[,cluster_id]), 
                                 mean, na.rm = T)
    class_centroids <- class_centroids[,-1]
  } else {
    cat("Class specific means were supplied and were not computed.\n 
        Double check input statement if this is correct.\n")
  }
  
  # Calculate the within class error sum of squares
  class_dists <- list()
  subdat <- numeric()
  for(c in 1:nClass)
  {
    subdat <- dat[dat[,cluster_id] == c,1:nV]
    class_dists[[c]] <- apply(subdat, 1, dist_func, 
                              centroid = class_centroids, C = c)
  }
  
  # Next, sum across individuals within a class, then unpack the lists
  #   in order to get the within class errors into 1 vector
  class_error <- unlist(lapply(class_dists, sum, na.rm = T))
  
  # compute class specific Ec's as per Bergman's chapter:
  # This weight's class Ec by class sample size and number of variables
  Adj_Ec <- numeric(nClass)
  Adj_Ec <- class_error / (classSize * nV)
  
  # Then, sum across classes for total Ec
  Ec <- sum(class_error)/nV
  
  # Now finally, compute EV
  EV <- 100*(Et - Ec)/Et
  
  output <- matrix(0, nrow=(nClass+1), ncol=1)
  output[,1] <- c(EV, Adj_Ec)
  rownames(output) <- c("EV", paste0("Adj_Ec_Class_",1:nClass))
  colnames(output) <- "Est"
  
  return(output)
}

