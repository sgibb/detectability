fetchAAindex <- function(url="ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1",
                         verbose=interactive()) {
  tmp <- tempfile()

  if (download.file(url, tmp, quiet=!verbose) != 0L) {
    stop("Could not download ", sQuote(url))
  }
  tmp
}

parseAAindex1 <- function(file, verbose=interactive()) {
  lines <- readLines(file)

  key <- c(AccessionNumber="H",
           DataDescription="D",
           LitDbEntryNumber="R",
           Authors="A",
           Title="T",
           JournalReference="J",
           Correlation="C",
           Index="I",
           End="\\/\\/")
  nkey <- length(key)

  # find lines
  pos <- matrix(grep(paste0("^", key, collapse="|"), lines),
                nrow=nkey, dimnames=list(names(key), NULL))
  pos["Index", ] <- pos["Index", ] + 1L
  pos["End", ] <- pos["End", ] - 1L

  # remove keys
  lines <- substring(lines, 3L, nchar(lines))

  # fetch metadata
  metadata <- lapply(1:6, function(i) {
    mapply(function(from, to) paste0(lines[from:to], collapse=""),
           from=pos[i,], to=pos[i + 1L,] - 1L)
  })
  names(metadata) <- names(key)[seq_along(metadata)]

  # fetch index
  index <- matrix(scan(text=lines[pos[c("Index", "End"),]], quiet=TRUE),
                  nrow=ncol(pos), ncol=20L, byrow=TRUE,
                  dimnames=list(accession, .aaNames1))

}
