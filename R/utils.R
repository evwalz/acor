# Shared utility functions used by both AKC and AGC implementations

# Helper function to check if Y is binary
is_binary <- function(Y) {
  unique_vals <- unique(Y)
  length(unique_vals) == 2
}

# Fenwick Tree (Binary Indexed Tree) implementation
# Supports point updates and prefix sum queries in O(log n)

fenwick_create <- function(n) {
  numeric(n)
}

fenwick_update <- function(tree, i, delta = 1) {
  n <- length(tree)
  while (i <= n) {
    tree[i] <- tree[i] + delta
    i <- i + bitwAnd(i, -i)
  }
  tree
}

fenwick_query <- function(tree, i) {
  if (i <= 0) return(0)
  s <- 0
  while (i > 0) {
    s <- s + tree[i]
    i <- i - bitwAnd(i, -i)
  }
  s
}