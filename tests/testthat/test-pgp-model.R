# Tests for the P-glycoprotein (P-gp/ABCB1) substrate prediction model.
# These exercise the pure-R Random Forest stored in R/pgp_model.R without
# requiring rcdk or any external ML package.

## ---- ..predict_pgp basic behavior ----------------------------------------

test_that("..predict_pgp returns numeric(0) for empty input", {
  expect_equal(..predict_pgp(numeric(0), numeric(0), numeric(0),
                              numeric(0), numeric(0), numeric(0),
                              numeric(0), numeric(0), numeric(0)),
               numeric(0))
})

test_that("..predict_pgp returns NA when any descriptor is NA", {
  p <- ..predict_pgp(300, NA, 80, 2, 4, 3, 20, 6, 70)
  expect_true(is.na(p))
})

test_that("..predict_pgp returns a probability in [0, 1] for valid input", {
  p <- ..predict_pgp(300, 2, 80, 2, 4, 3, 20, 6, 70)
  expect_false(is.na(p))
  expect_true(p >= 0 && p <= 1)
})

## ---- Pharmacological sanity checks ---------------------------------------
## The model was trained on Metrabase P-gp substrate data. We verify that
## known P-gp substrates get higher probabilities than known non-substrates.

test_that("small simple molecules are predicted as non-substrates", {
  ## Ethanol: CCO — very small, very polar, not a P-gp substrate
  p_ethanol <- ..predict_pgp(46, -0.1, 20, 1, 1, 0, 3, 0, 11)
  expect_true(p_ethanol < 0.5)
  expect_equal(..predict_pgp(46, -0.1, 20, 1, 1, 0, 3, 0, 11) >= 0.5, FALSE)
})

test_that("large lipophilic molecules are predicted as substrates", {
  ## A large, lipophilic molecule with many aromatic atoms and high MW
  ## (vinblastine-like profile: MW ~800, LogP ~4, TPSA ~180, many rings)
  p_large <- ..predict_pgp(811, 4.0, 180, 4, 8, 10, 60, 12, 200)
  expect_true(p_large >= 0.5)
})

test_that("caffeine and aspirin are predicted as non-substrates", {
  ## Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C — MW 194, known non-substrate
  p_caff <- ..predict_pgp(194, -0.1, 58, 0, 6, 0, 14, 9, 50)
  expect_true(p_caff < 0.5)

  ## Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O — MW 180, known non-substrate
  p_asp <- ..predict_pgp(180, 1.2, 64, 1, 4, 3, 13, 6, 45)
  expect_true(p_asp < 0.5)
})

## ---- Vectorized input -----------------------------------------------------

test_that("..predict_pgp handles vectorized input correctly", {
  mw   <- c(46, 180, 194, 811)
  logp <- c(-0.1, 1.2, -0.1, 4.0)
  tpsa <- c(20, 64, 58, 180)
  hbd  <- c(1, 1, 0, 4)
  hba  <- c(1, 4, 6, 8)
  rb   <- c(0, 3, 0, 10)
  ha   <- c(3, 13, 14, 60)
  arom <- c(0, 6, 9, 12)
  mr   <- c(11, 45, 50, 200)

  probs <- ..predict_pgp(mw, logp, tpsa, hbd, hba, rb, ha, arom, mr)
  expect_length(probs, 4)
  expect_true(all(probs >= 0 & probs <= 1))
  ## The large molecule (index 4) should have the highest probability
  expect_true(probs[4] > probs[1])
  expect_true(probs[4] >= 0.5)
  ## The small molecules (indices 1-3) should be non-substrates
  expect_true(probs[1] < 0.5)
  expect_true(probs[2] < 0.5)
  expect_true(probs[3] < 0.5)
})

test_that("..predict_pgp handles mixed NA and valid input", {
  mw   <- c(300, 400, NA)
  logp <- c(2, 3, 2)
  tpsa <- c(80, 60, 50)
  hbd  <- c(2, 1, 2)
  hba  <- c(4, 3, 4)
  rb   <- c(3, 5, 3)
  ha   <- c(20, 25, 18)
  arom <- c(6, 0, 6)
  mr   <- c(70, 90, 60)

  probs <- ..predict_pgp(mw, logp, tpsa, hbd, hba, rb, ha, arom, mr)
  expect_length(probs, 3)
  expect_false(is.na(probs[1]))
  expect_false(is.na(probs[2]))
  expect_true(is.na(probs[3]))  ## third molecule has NA MW
})

## ---- computeADMETProperties integration ----------------------------------

test_that("computeADMETProperties adds Pgp substrate and probability columns", {
  d <- data.frame(
    LogP = c(2.5, 6),
    TPSA = c(70, 220),
    MW = c(300, 650),
    "#H-bond donors" = c(2, 7),
    "#H-bond acceptors" = c(4, 12),
    "#Rotatable bonds" = c(3, 14),
    "#Heavy atoms" = c(20, 80),
    "#Aromatic heavy atoms" = c(6, 12),
    MR = c(70, 150),
    check.names = FALSE
  )
  out <- computeADMETProperties(d)
  expect_true("Pgp substrate" %in% names(out))
  expect_true("Pgp probability" %in% names(out))
  expect_true(all(out[["Pgp probability"]] >= 0 & out[["Pgp probability"]] <= 1))
  ## Pgp substrate should be "Yes" or "No"
  expect_true(all(out[["Pgp substrate"]] %in% c("Yes", "No")))
})

test_that("computeADMETProperties returns NA Pgp when descriptors are missing", {
  ## Only LogP and TPSA — no MW, HBD, HBA, etc. — P-gp prediction should be NA
  d <- data.frame(LogP = c(2, 5), TPSA = c(60, 90), check.names = FALSE)
  out <- computeADMETProperties(d)
  expect_true(all(is.na(out[["Pgp substrate"]])))
  expect_true(all(is.na(out[["Pgp probability"]])))
})
