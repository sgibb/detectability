context("chemScore")

sequences <- c("HGLDNYR",
               "RHGLDNYR",
               "WWCNDGR",
               "GTDVQAWIR",
               "GYSLGNWVCAAK",
               "FESNFNTQATNR",
               "IVSDGNGMNAWVAWR",
               "NTDGSTDYGILQINSR",
               "KIVSDGNGMNAWVAWR",
               "DKLDAAAK")

test_that("chemScore generates same scores as in Parker 2002", {
  scores <- setNames(c(100,
                       23.3,
                       5,
                       100,
                       0.5,
                       100,
                       100,
                       100,
                       30.1,
                       8), sequences)

  expect_equal(chemScore(sequences[-c(3, 5)], metOxF=1.0), scores[-c(3, 5)])
  expect_equal(chemScore(sequences[c(3, 5)], cys=100, detrimentalCys=10), scores[c(3, 5)])
})

test_that(".baseChemScore", {
  scoreNone <- c(100, 100, 100, 100, 10, 100, 100, 100, 100, 10)
  scoreVp <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 10)

  expect_error(detectability:::.baseChemScore(sequences, arg=c(100, -1)),
               ".*arg.* has to be a .*double.* greater than or equal 0")
  expect_error(detectability:::.baseChemScore(sequences, cys=-1),
               ".*cys.* has to be a .*double.* greater than or equal 0")
  expect_error(detectability:::.baseChemScore(sequences, lys=-1),
               ".*lys.* has to be a .*double.* greater than or equal 0")

  expect_message(detectability:::.baseChemScore("ACE", verbose=TRUE),
                 "rule 1-5: ACE, RCK=010, arg=100, cys=0, lys=10, score=1")
  expect_silent(detectability:::.baseChemScore(sequences, verbose=FALSE))

  expect_equal(detectability:::.baseChemScore(sequences), scoreNone)
  expect_equal(detectability:::.baseChemScore(sequences, cys=100), scoreVp)
  expect_equal(detectability:::.baseChemScore("ACE"), 1)
  expect_equal(detectability:::.baseChemScore("ACE", cys=100), 100)
  expect_equal(detectability:::.baseChemScore("ACE", cys=10), 10)
  expect_equal(detectability:::.baseChemScore(c("RAC", "RACR"),
                                              arg=c(200, 100), cys=c(0, 150)),
               c(200, 150))
})

test_that(".cysteine", {
  s <- c("ACE", "AKE")
  expect_error(detectability:::.cysteine(s, detrimentalCys=-1),
               paste0(".*detrimentalCys.* has to be a .*double.* ",
                      "greater than or equal 0"))

  expect_message(detectability:::.cysteine(s, verbose=TRUE),
                 paste0("rule 6 and 7: ", s, ", nC=", 1:0,
                        ", detrimentalCys=", c(10, 1), collapse="\n"))
  expect_silent(detectability:::.cysteine(s, verbose=FALSE))

  expect_equal(detectability:::.cysteine(s), c(10, 1))
})

test_that(".methionine", {
  s <- c("RMR", "RMM", "RMMM")
  m <- c(1, 0.5, 2)

  expect_error(detectability:::.methionine(s, metOxF=0.1),
               ".*metOxF.* has to be a .*double.* between 0.2 and 5.0")
  expect_error(detectability:::.methionine(s, metOxF=5.1),
               ".*metOxF.* has to be a .*double.* between 0.2 and 5.0")
  expect_error(detectability:::.methionine(s, metOxF=c(1, 0.1, 2)),
               ".*metOxF.* has to be a .*double.* between 0.2 and 5.0")
  expect_error(detectability:::.methionine(s, metOxF=0.2, nMetOx=-1),
               ".*nMetOx.* has to be an .*integer.* greater than or equal 0")

  expect_message(detectability:::.methionine(s, verbose=TRUE),
                 paste0("rule 8-13: ", s,
                        ", nM=", 1:3, ", nR=", 0:2, ", nX=", 1, ", metOxF=5.0",
                        collapse="\n"))
  expect_silent(detectability:::.methionine(s, verbose=FALSE))

  expect_equal(detectability:::.methionine(s), rep(5, 3))
  expect_equal(detectability:::.methionine(s, nMetOx=1:3), 5^(1:3))
  expect_equal(detectability:::.methionine(s, metOxF=m), c(2, 2, 4))

  expect_equal(detectability:::.methionine(rep("ACE", 3),
                                           metOxF=c(0.2, 1, 5)), rep(1, 3))
})

test_that(".proline", {
  s <- c("PACE", "AKE")
  expect_message(detectability:::.proline(s, verbose=TRUE),
                 paste0("rule 14: ", s, ", score=", c(100, 1), collapse="\n"))
  expect_silent(detectability:::.proline(s, verbose=FALSE))

  expect_equal(detectability:::.proline(s), c(100, 1))
})

test_that(".chemScorePartialFactor", {
  scores <- c(1,
              (100 + 30)/30,
              1,
              1,
              1,
              1,
              1,
              1,
              (100 + 30 * 5)/(30 * 5),
              (100 + 20 * 5 * 2 * 2)/(20 * 5 * 2 * 2))

  expect_equal(detectability:::.chemScorePartialFactor(sequences), scores)
})
