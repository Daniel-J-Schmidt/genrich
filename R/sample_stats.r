#' Balance per-locus sample size between two populations.
#'
#' Takes a hierfstat input file and equalizes sample size independently for each
#' locus by randomly removing genotypes from the sample with higher n. Returns
#' the modified hierfstat input file. If more than two population samples are
#' present in the input file, then only the first two listed pop samples are
#' processed.
#' @return A modified version of the input dataframe in hierfstat format.
#' @examples
#' test <- bal_loci(dat.sim)
#' @export
#'
bal_loci <- function(data) {
  pop1       <- unique(data[, 1])[1]
  pop2       <- unique(data[, 1])[2]
  pop1.index <- which(data[, 1] == pop1)
  pop2.index <- which(data[, 1] == pop2)
  data.bal   <- lapply(data[-1], function(x) {
    n.ind.samp1 <- sum(!is.na(x[pop1.index]))
    n.ind.samp2 <- sum(!is.na(x[pop2.index]))
    if (n.ind.samp1 > n.ind.samp2) {
      lose.n <- n.ind.samp1 - n.ind.samp2
      x[pop1.index][as.numeric(sample(as.list(which(!is.na(
        x[pop1.index]
      ))), lose.n))] <- NA
      return(x)
    } else if (n.ind.samp1 < n.ind.samp2) {
      lose.n <- n.ind.samp2 - n.ind.samp1
      x[pop2.index][as.numeric(sample(as.list(which(!is.na(
        x[pop2.index]
      ))), lose.n))] <- NA
      return(x)
    } else {
      return(x)
    }
  })
  data.bal <- data.frame(Pop = data[, 1], do.call(cbind, data.bal))
}

#' Gene diversity (Hs).
#'
#' Calculates Nei's (1987) gene diversity (= Hs = expected heterozygosity
#' corrected for small sample size) for each locus and each population sample.
#' Code is modified from function hierfstat::basic.stats and relies on several
#' functions from the hierfstat library: pop.freq(), ind.count(), getal.b().
#' @return A matrix of calculated Hs values, with one row per locus and one
#'   column per population sample.
#' @examples
#' test <- Hs_calc(dat.sim)
#' @export
#'
Hs_calc <- function(data){
  p <- pop.freq(data)
  n <- t(ind.count(data))
  dum <- getal.b(data[, -1])
  Ho <- dum[, , 1] == dum[, , 2]
  sHo <- (1 - t(apply(Ho, 2, fun <- function(x) tapply(x, data[, 1], mean, na.rm = TRUE))))
  mHo <- apply(sHo, 1, mean, na.rm = TRUE)
  sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
  sp2 <- matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
  Hs <- (1 - sp2 - sHo/2/n)
  Hs <- n/(n - 1) * Hs
  return(Hs)
}

#' Calculate null distribution of difference in gene diversity (Hs) between two
#' populations.
#'
#' Randomly permutes individuals between sample groups and calculates difference
#' in gene diversity (Hs) between the two randomised groups. Used to generate
#' null distribution of difference in gene diversity between two samples for
#' permutation testing
#' @param bal.loci Option to incorporate balancing of per-locus sample size into
#'   the permutation (see function bal_loci). When bal.loci = TRUE, each
#'   replicate of the permutation first balances sample sizes before gene
#'   diversity calculation is performed. When bal.loci = FALSE, sample balancing
#'   is skipped.
#' @param reps Number of iterations to repeat the permutation. If not specified,
#'   defaults to 100.
#' @return A vector of values of length "reps". Each value represents the
#'   difference in gene diversity between two populations after individuals have
#'   been randomly permuted between the populations. This difference is
#'   calculated as Hs(pop2) - Hs(pop1).
#' @examples
#' test <- shuff_data(data = dat.sim, reps = 100, bal.loci = FALSE)
#' @export
#'
shuff_data <- function(data, reps = 100, bal.loci = TRUE) {
  if (bal.loci) {
    bal <- lapply(1:reps, function(x) {
      data[, 1] <- sample(data[, 1])
      data.bal <- bal_loci(data)
      Hs <- Hs_calc(data.bal)
      mHs <- colMeans(Hs, na.rm = TRUE)
      Hs_mdiff <- mHs[2] - mHs[1]
    })
    return(as.vector(do.call(rbind, bal)))
  } else {
    not.bal <- lapply(1:reps, function(x) {
      data[, 1] <- sample(data[, 1])
      Hs <- Hs_calc(data)
      mHs <- colMeans(Hs, na.rm = TRUE)
      Hs_mdiff <- mHs[2] - mHs[1]
    })
    return(as.vector(do.call(rbind, not.bal)))
  }
}
