merge_design_points <- function(res, Kv, threshold = 0.1, digits = 4) {
  # res is the output object of the locally() function
  des <- res$design[-1]
  x <- as.vector(unlist(des[1:Kv]))
  w <- as.vector(unlist(des[(Kv + 1):(2 * Kv)]))
  
  # Keep points with positive weight (to guard against numerical zeros)
  keep <- w > 1e-4
  x <- x[keep]
  w <- w[keep]
  
  # Sort by x
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  
  # Greedily merge points whose distance is smaller than the threshold
  merged_x <- c()
  merged_w <- c()
  i <- 1
  n <- length(x_sorted)
  while (i <= n) {
    group <- i
    j <- i + 1
    while (j <= n && (x_sorted[j] - x_sorted[i]) < threshold) {
      group <- c(group, j)
      j <- j + 1
    }
    new_x <- mean(x_sorted[group])
    new_w <- sum(w_sorted[group])
    merged_x <- c(merged_x, new_x)
    merged_w <- c(merged_w, new_w)
    i <- j
  }
  
  # Build a data frame
  DW <- data.frame(x = round(merged_x, digits = digits), w = merged_w)
  # Aggregate (to ensure that identical x values are merged; kept for robustness,
  # although the previous step has already done so)
  agg <- aggregate(w ~ x, data = DW, FUN = sum)
  # Generate an xtable
  library(xtable)
  xtable(t(agg), digits = digits)
}


merge_2d_design_points <- function(res, Kv, threshold = 0.1, digits = 4) {
  # Extract design points and weights
  des <- res$design[-1]
  x1 <- as.vector(unlist(des[1:Kv]))
  x2 <- as.vector(unlist(des[(Kv + 1):(2 * Kv)]))
  w  <- as.vector(unlist(des[(2 * Kv + 1):(3 * Kv)]))
  
  # Keep points with positive weight (to guard against numerical zeros)
  keep <- w > 1e-4
  x1 <- x1[keep]
  x2 <- x2[keep]
  w  <- w[keep]
  
  if (length(x1) == 0) {
    stop("No design points with positive weight.")
  }
  
  # Build the point matrix
  points <- cbind(x1, x2)
  
  # Compute the distance matrix (Euclidean distance)
  d <- dist(points)
  
  # Perform hierarchical clustering with complete linkage and cut at height = threshold
  hc <- hclust(d, method = "complete")
  clusters <- cutree(hc, h = threshold)
  
  # Merge points within the same cluster
  merged_x1 <- tapply(x1, clusters, mean)
  merged_x2 <- tapply(x2, clusters, mean)
  merged_w  <- tapply(w,  clusters, sum)
  
  # Build a data frame and sort it (by the first coordinate, optional)
  DW <- data.frame(x1 = round(merged_x1, digits),
                   x2 = round(merged_x2, digits),
                   w  = merged_w)
  DW <- DW[order(DW$x1, DW$x2), ]
  
  # Generate a summary table (combine x1 and x2 into a label, or keep two columns)
  # For display, use (x1, x2) as row names and the weights as values
  agg <- data.frame(
    label = paste0("(", DW$x1, ",", DW$x2, ")"),
    weight = DW$w
  )
  # Transpose into one row, suitable for xtable output
  result <- t(agg$weight)
  colnames(result) <- agg$label
  
  library(xtable)
  xtable(result, digits = digits)
}