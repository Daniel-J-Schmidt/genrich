#' Population-specific FST with sample rarefaction.
#'
#' Randomly subsample each population in the dataset to a value (rar.n) set by
#' user, then calculate population-specific FST for each population sample using
#' hierfstat::betas. Expects data in hierfstat format.
#' @param rar.n Rarefied sample size for each population. Must be equal or less
#'   than the minimum population sample size. If not specified, rar.n defaults
#'   to a value equal to the smallest sample size in the dataset.
#' @param reps Number of iterations to repeat the rarefied FST calculation. If
#'   not specified, defaults to 100.
#' @return Dataframe of FST estimates, where number of rows equal to iterations,
#'   and columns represent populations.
#' @examples
#' test <- rarefy_fst(dat.sim, rar.n = 5, reps = 100)
#' @export
#'
rarefy_fst <- function(data, rar.n = NULL, reps = 100) {
  if (is.null(rar.n)) {
    rar.n <- min(table(data[, 1]))
  }
  if (min(table(data[, 1])) < rar.n |
      rar.n > max(table(data[, 1]))) {
    stop("rar.n must be equal to, or smaller than, the lowest population size")
  }
  loop <- lapply(1:reps, function(y) {
    samp.rows <- lapply(unique(data[, 1]), function(x) {
      data[as.numeric(sample(as.list(which(data[, 1] == x)), rar.n)), ]
    })
    dat <- do.call(rbind, samp.rows)
    beta <- hierfstat::betas(dat, nboot = 0, diploid = TRUE)
    return(beta$betaiovl)
  })
  out <- data.frame(do.call(rbind, loop))
  names(out) <- unique(dat.sim[, 1])
  return(out)
}
