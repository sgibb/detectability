#' ChemScore
#'
#' This function calculates the ChemScore proposed in Parker 2001.
#'
#' @param x \code{character}, amino acid sequence(s).
#' @param arg \code{double}, base ChemScore for Arginine (default: 100).
#' @param cys \code{double}, base ChemScore for Cysteine (default: 0).
#' @param lys \code{double}, base ChemScore for Lysine (default: 10).
#' @param detrimentalCys \code{double}, score for detrimental cysteine.
#' @param metOxF \code{double}, methionine oxidation factor (default: 0.2).
#' @param nMetOx \code{integer}, number of oxidized methionine (default: 1).
#' @param verbose \code{logical}, verbose output?
#' @references
#' Parker KC. Scoring methods in MALDI peptide mass fingerprinting: ChemScore,
#' and the ChemApplex program. \cr
#' Journal of the American Society for Mass Spectrometry. 2002 Jan 31;13(1):22-39.
#' @export
#' @examples
#' library("detectability")
#'
#' # Example sequences taken from Table 4 in Parker 2001.
#' sequences <- c("HGLDNYR",
#'                "RHGLDNYR",
#'                "WWCNDGR",
#'                "GTDVQAWIR",
#'                "GYSLGNWVCAAK",
#'                "FESNFNTQATNR",
#'                "IVSDGNGMNAWVAWR",
#'                "NTDGSTDYGILQINSR",
#'                "KIVSDGNGMNAWVAWR")
#' chemScore(sequences)
chemScore <- function(x, arg=100, cys=0, lys=10,
                      detrimentalCys=10,
                      metOxF=0.2, nMetOx=1L,
                      verbose=interactive()) {
  score <- .baseChemScore(x, arg=arg, cys=cys, lys=lys,
                          verbose=verbose)
  score <- score / .cysteine(x, detrimentalCys=detrimentalCys,
                             verbose=verbose)
  score <- score / .methionine(x, metOxF=metOxF, nMetOx=nMetOx,
                               verbose=verbose)
  score <- score / .proline(x, verbose)
  setNames(score / .chemScorePartialFactor(x, verbose=verbose), x)
}

#' rule 1-5
#' @param x \code{character}, amino acid sequence(s)
#' @param arg \code{double}, base ChemScore for Arginine (default: 100).
#' @param cys \code{double}, base ChemScore for Cysteine (default: 0).
#' @param lys \code{double}, base ChemScore for Lysine (default: 10).
#' @param verbose \code{logical}, verbose output?
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
  metOxF[metOxF == 1L] <- 2
  metOxF[gt1] <- metOxF[gt1]^(nR[gt1])
  metOxF[lt1] <- 1 / metOxF[lt1]^(nX[lt1])
  metOxF[nM == 0L] <- 1L

  .msg(verbose, paste0(sprintf("rule 8-13: %s, nM=%i, nR=%i, nX=%i, metOxF=%0.1f",
                               x, nM, nR, nX, metOxF), collapse="\n"))

  metOxF
}

#' rule 14
#' @param x \code{character}, amino acid sequence(s).
#' @param verbose \code{logical}, verbose output?
#' @noRd
.proline <- function(x, verbose=interactive()) {
  p <- rep.int(1, length(x))
  p[substr(x, 1, 1) == "P"] <- 100
  .msg(verbose, paste0("rule 14: ", x, ", score=", p, collapse="\n"))
  p
}

#' rule 16-25
#' @param x \code{character}, amino acid sequence(s).
#' @param bmcf \code{double}, Basal Missed Cleavage Factor, typical 100.
#' @param cleavageRule \code{character}, regular trypsin cleavage.
#' @param missedCleavageRules \code{data.frame}, two-columns (pattern, score).
#' @param verbose \code{logical}, verbose output?
#' @importFrom cleaver cleavageSites
#' @noRd
.chemScorePartialFactor <- function(x, bmcf=100, cleavageRule="[KR].",
                                    missedCleavageRules=data.frame(
                                      pattern=c("[KR]P",              # 16
                                                "^[KR].",             # 17
                                                "(?<=[DE])[KR].",     # 18
                                                "[KR][DE]",           # 19
                                                "[KR][ILV]",          # 20
                                                "[KR][KR]|[KR].$",    # 21
                                                "(?<=[DE].)[KR].",    # 22
                                                "[KR].[DE]",          # 23
                                                "(?<=^.)[KR].",       # 24
                                                "[KR].[KR]|[KR]..$"), # 25
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
                                      stringsAsFactors=FALSE),
                                    verbose=interactive()) {
  sites <- cleavageSites(x, custom=cleavageRule)

  isMissingCleavage <- which(lengths(sites) != 0L)

  chpf <- rep.int(1, length(x))

  r <- lapply(missedCleavageRules$pattern,
              function(p)cleavageSites(x[isMissingCleavage], custom=p))

  for (i in seq(along=isMissingCleavage)) {
    curSites <- sites[isMissingCleavage[i]]
    mcf <- rep.int(1, length(curSites))

    for (j in seq(along=r)) {
      m <- match(curSites, r[[j]][i])
      if (verbose && !is.na(m)) {
        .msg(verbose, paste0("rule ", j + 15, ": ", x[isMissingCleavage[i]],
                             ", pattern=", sQuote(missedCleavageRules$pattern[j]),
                             ", score=", missedCleavageRules$score[j]))
      }
      mcf[m] <- mcf[m] * missedCleavageRules$score[j]
    }
    chpf[isMissingCleavage[i]] <- prod(chpf[isMissingCleavage[i]], mcf/(bmcf + mcf))
  }
  chpf = 1 / chpf
  .msg(verbose, paste0("ChPF: ", x, ", score=", chpf))
  chpf
}
