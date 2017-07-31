#' Plot difference in sample size for each locus.
#'
#' Provides a quick visual check to compare missing data between two populations
#' and to verify that function bal_loci() has equalised sample sizes correctly.
#' Expects data in hierfstat format. #' If more than two population samples are
#' present in the input file, then only the first two listed pop samples are
#' processed. Plot can be modified using any additional arguments to plot() e.g.
#' add a title with the argument: main = "title text".
#' @return A lineplot showing difference in sample size for each locus.
#' @examples
#' test <- plot_bal_loci(dat.sim, main = "With missing data")
#' @export
#'
plot_bal_loci <- function(data, ...) {
  pop1 <- unique(data[, 1])[1]
  pop2 <- unique(data[, 1])[2]
  a <- apply(data[data[, 1] == pop1, ], 2, function(x) {
    sum(!is.na(x))
  })
  b <- apply(data[data[, 1] == pop2, ], 2, function(x) {
    sum(!is.na(x))
  })
  score.diff <- a - b
  plot(score.diff, type = "l", ylab = "diff. n pop1 & pop2", xlab = "Locus index", las=1, ...)
}
