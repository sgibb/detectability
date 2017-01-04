#' @param x \code{character}, amino acid sequence(s).
#' @param arg \code{double}, base ChemScore for Arginine (default: 100).
#' @param cys \code{double}, base ChemScore for Cysteine (default: 0).
#' @param lys \code{double}, base ChemScore for Lysine (default: 10).
#' @param detrimentalCys \code{double}, score for detrimental cysteine.
#' @param metOxF \code{double}, methionine oxidation factor (default: 0.2).
#' @param nMetOx \code{integer}, number of oxidized methionine (default: 1).
#' @param \ldots further arguments passed to internal functions.
#' @export
chemScore <- function(x, arg=100, cys=0, lys=10,
                      detrimentalCys=10,
                      metOxF=0.2, nMetOx=1L,
                      verbose=interactive()) {
  score <- .baseChemScore(x, arg=arg, cys=cys, lys=lys,
                          verbose=verbose)
  score <- score / .cysteine(x, detrimentalCys=detrimentalCys,
                             verbose=verbose)
  score <- score * .methionine(x, metOxF=metOxF, nMetOx=nMetOx,
                               verbose=verbose)
  score / .chemScorePartialFactor(x)
}

#' rule 1-5
#' @param x \code{character}, amino acid sequence(s)
#' @param arg \code{double}, base ChemScore for Arginine (default: 100).
#' @param cys \code{double}, base ChemScore for Cysteine (default: 0).
#' @param lys \code{double}, base ChemScore for Lysine (default: 10).
#' @noRd
.baseChemScore <- function(x, arg=100, cys=0, lys=10, verbose=interactive()) {
  n <- length(x)

  if (any(!is.double(arg)) || any(arg < 0L)) {
    stop(sQuote("arg"), " has to be a ", sQuote("double"),
         " greater than or equal 0.")
  }

  if (any(!is.double(cys)) || any(cys < 0L)) {
    stop(sQuote("cys"), " has to be a ", sQuote("double"),
         " greater than or equal 0.")
  }

  if (any(!is.double(lys)) || any(lys < 0L)) {
    stop(sQuote("lys"), " has to be a ", sQuote("double"),
         " greater than or equal 0.")
  }

  arg <- rep_len(arg, n)
  cys <- rep_len(cys, n)
  lys <- rep_len(lys, n)

  scores <- cbind(arg, cys, lys, 1)

  r <- lapply(c("R", "C", "K"), grepl, x=x)
  r <- do.call(cbind, c(r, TRUE))
  mode(r) <- "integer"
  cs <- r * scores
  cs <- cs[cbind(1L:nrow(cs), max.col(cs, ties.method="first"))]

  .msg(verbose, paste0("rule 1-5: ", x,
                       ", RCK=", apply(r[, 1L:3L, drop=FALSE], 1, paste0, collapse=""),
                       ", arg=", arg, ", cys=", cys, ", lys=", lys,
                       ", score=", cs, collapse="\n"))
  cs
}

#' rule 6-7
#' @param x \code{character}, amino acid sequence(s).
#' @param detrimentalCys \code{double}, score for detrimental cysteine.
#' @param verbose \code{logical}, verbose ouput?
#' @noRd
.cysteine <- function(x, detrimentalCys=10, verbose=interactive()) {
  if (any(!is.double(detrimentalCys)) || any(detrimentalCys < 0L)) {
    stop(sQuote("detrimentalCys"), " has to be a ", sQuote("double"),
         " greater than or equal 0.")
  }

  if (is.character(x)) {
    x <- AAStringSet(x)
  }

  nC <- as.integer(letterFrequency(x, letters="C"))
  detrimentalCys <- (nC > 0L) * detrimentalCys
  detrimentalCys[detrimentalCys == 0L] <- 1L # avoid division by zero

  .msg(verbose, paste0("rule 6 and 7: ", x, ", nC=", nC,
                       ", detrimentalCys=", detrimentalCys, collapse="\n"))

  detrimentalCys
}

#' rule 8-13
#' @param x \code{character}, amino acid sequence(s).
#' @param metOxF \code{double}, methionine oxidation factor (default: 0.2).
#' @param nMetOx \code{integer}, number of oxidized methionine (default: 1).
#' @param verbose \code{logical}, verbose ouput?
#' @importFrom Biostrings AAStringSet letterFrequency
#' @noRd
.methionine <- function(x, metOxF=0.2, nMetOx=1L, verbose=interactive()) {
  if (any(!is.double(metOxF)) || any(metOxF < 0.2) || any(metOxF > 5.0)) {
    stop(sQuote("metOxF"), " has to be a ", sQuote("double"),
         " value between 0.2 and 5.0.")
  }

  if (any(!is.numeric(nMetOx)) || any(nMetOx < 0L)) {
    stop(sQuote("nMetOx"), " has to be an ", sQuote("integer"),
         " greater than or equal 0.")
  }

  if (is.character(x)) {
    x <- AAStringSet(x)
  }

  nM <- as.integer(letterFrequency(x, letters="M"))
  nR <- pmax(nM - nMetOx, 0L)
  nX <- pmin(nM, nMetOx)

  metOxF <- rep_len(metOxF, length(x))
  gt1 <- which(metOxF > 1L)
  lt1 <- which(metOxF < 1L)
  metOxF[metOxF == 1L] <- 1 / 2
  metOxF[gt1] <- 1 / metOxF[gt1]^(nR[gt1])
  metOxF[lt1] <- metOxF[lt1]^(nX[lt1])

  .msg(verbose, paste0(sprintf("rule 8-13: %s, nM=%i, nR=%i, nX=%i, metOxF=%0.3f",
                               x, nM, nR, nX, metOxF), collapse="\n"))

  metOxF
}

#' @param x \code{character}, amino acid sequence(s)
#' @param bmcf \code{double}, Basal Missed Cleavage Factor, typical 100
#' @importFrom cleaver cleavageSites
#' @noRd
.chemScorePartialFactor <- function(x, bmcf=100,
                                    rules=data.frame(pattern=c("[KR]P",       # 16
                                                               "^[KR].",      # 17
                                                               "[DE][KR].",   # 18
                                                               "[KR][DE]",    # 19
                                                               "[KR][ILV]",   # 20
                                                               "[KR].$",      # 21
                                                               "[DE].[KR].",  # 22
                                                               "[KR].[DE]",   # 23
                                                               ".[KR].",      # 24
                                                               "[KR]..$"),    # 25
                                                     score=c(100,   # 16
                                                             30,    # 17
                                                             20,    # 18
                                                             20,    # 19
                                                             5,     # 20
                                                             3,     # 21
                                                             2,     # 22
                                                             2,     # 23
                                                             2,     # 24
                                                             1.5),  # 25
                                                     stringsAsFactors=FALSE)) {
  r <- lapply(rules$pattern, function(p)lengths(cleavageSites(x, custom=p)))
  r <- do.call(rbind, r) # rows == pattern, columns == x
  r <- r * rules$score
  r[r == 0L] <- 1L # avoid multiplication by zero
  r <- apply(r, 2L, prod)
  r <- (bmcf + r)/r
  r[r == 101L] <- 1L # no cleavage rule matched
  r
}