sanitise <- function(data) {
  bool <- apply(data, c(1, 2), function(x) {x > 0})
  active_genes <- apply(bool, 2, sum)
  return(data[, active_genes > 2])
}