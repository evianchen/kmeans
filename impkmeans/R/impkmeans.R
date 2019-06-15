#' Perform improvement of k-means clustering on a data matrix.
#' 
#' Use K-means++ to find the beginning center point to improve K-means.
#' @usage imp_kmeans(datasets, k, convergence = 1e-6)
#' @param dataset numeric matrix of data, or an object that can be coerced to such a matrix(such as a numeric vector or a data frame with all numeric columns).
#' @param k the number of clusters.
#' @param convergence the convergence tolerance allowed.
#' @param begining_centroids The beginning center point.
#' @param final_centroids The final center point.
#' @param label The cluster to which each point is assigned.
#' @param loss  Total sum of squared.
#' @param iter The number of iterations.
#' @export 
#' @seealso \code{kmeans}
#' @references D. Arthur, S. Vassilvitskii, "K-Means++: The advantages of careful seeding", Proc. Symp. Discrete Algorithms, pp. 1027-1035, 2007
#' @examples
#' ## a 2-dimensional example
#' data(iris)
#' iris = iris[, -5] #remove label
#' (k_means = imp_kmeans(iris, 3))
#' plot(iris$Petal.Length, iris$Petal.Width, col = k_means$label)
#'        
#' ## b 2-dimensional example
#' x = rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#' colnames(x) = c("x", "y")
#' (cl = imp_kmeans(x, 2))
#' plot(x, col = cl$label)


imp_kmeans = function(datasets, k, convergence = 1e-6){
  datasets = data.matrix(datasets) 
  obs = nrow(datasets) 
  features = ncol(datasets) 
  prob_points = function(){
    init_pts = matrix(0, k, features) 
    index = sample(obs, 1)
    init_pts[1, ] = datasets[index, ]
    for(i in 1:k){
      null = rep(NA, obs)
      for(item in 1:obs){
        pts_dist = c()
        for(pts in 1:i){
          pts_dist = rbind(pts_dist, c(dist(rbind(init_pts[pts, ], datasets[item, ])))^2)
        }
        null[item] = min(pts_dist)
      }
      cdf = cumsum(null/sum(null))
      uniform = runif(1)
      for(interval in 1:obs){
        if(uniform >= cdf[interval]){
          next
        }
        else{
          interval = interval
          break
        }
      }
      init_pts[i, ] = datasets[interval, ]
    }
    return(init_pts)
  }
  cluster_centroids = function(){
    cluster_centroids = aggregate(datasets[, c(1:features)], list(null_clusterAssign), mean)[,-1]
    return(cluster_centroids)
  }
  begining_centroids = init_cluster_centroids = prob_points()
  begining_centroids = data.frame(begining_centroids)
  colnames(begining_centroids) = colnames(datasets)
  datasets = cbind(datasets, clusterAssign = rep(NA, obs))
  iter = 1
  while(TRUE){
    null_clusterAssign = rep(NA, obs)
    SSE = rep(NA, obs)
    for(value in 1:obs){
      null_dist = matrix(0, k, 1)
      for(j in 1:k){
        null_dist[j, ] = dist(rbind(datasets[, -(features+1)][value, ], init_cluster_centroids[j, ]))^2
      }
      SSE[value] = min(null_dist)
      null_clusterAssign[value] = which.min(null_dist)
    }
    if(length(unique(null_clusterAssign)) != k){
      init_cluster_centroids = prob_points(datasets[, -(features+1)], k)
      next
    }
    else{
      datasets[, features + 1] = null_clusterAssign
    }
    if(sqrt(sum((cluster_centroids() - init_cluster_centroids)^2)) < convergence){
      break
    }
    else{
      init_cluster_centroids = cluster_centroids()
      iter = iter + 1
    }
  }
  result = list("begining_centroids" = begining_centroids,
                "final_centroids" = cluster_centroids(),
                "label" = null_clusterAssign,
                "loss" = sum(SSE),
                "iter" = iter)
  return(result)
}

