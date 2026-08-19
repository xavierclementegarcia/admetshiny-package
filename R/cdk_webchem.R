# ---------------------------------------------------------------------------
# cdk_webchem.R
# Retrieval of canonical SMILES from PubChem (via webchem) and local
# calculation of molecular descriptors with CDK (via rcdk), plus the
# drug-likeness violation and BOILED-Egg ADMET property calculators shared by
# all modules that use CDK.
# ---------------------------------------------------------------------------

#' Retrieve canonical SMILES from PubChem
#'
#' Queries PubChem (via the suggested package \pkg{webchem}) to obtain CIDs and
#' canonical SMILES for a vector of chemical identifiers (names, CAS numbers,
#' InChIKeys or PubChem CIDs).
#'
#' @param ids Character vector of identifiers.
#' @param from Character. Type of identifier: one of \code{"name"},
#'   \code{"cas"}, \code{"inchikey"}, \code{"cid"}. Default \code{"name"}.
#' @return A data.frame with columns \code{query}, \code{cid},
#'   \code{CanonicalSMILES}, \code{IsomericSMILES}, \code{MolecularFormula},
#'   \code{IUPACName}.
#' @examples
#' \dontrun{
#' ids <- c("aspirin", "ibuprofen")
#' smiles <- getSmilesFromIdentifiers(ids, from = "name")
#' }
#' @export
#' @seealso \code{\link{calcCDKDescriptors}}.
getSmilesFromIdentifiers <- function(ids, from = "name") {

  if (!requireNamespace("webchem", quietly = TRUE)) {
    stop("This function requires the 'webchem' package. Install it with: ",
         "install.packages('webchem')", call. = FALSE)
  }

  ids <- trimws(as.character(ids))
  ids <- ids[ids != "" & !is.na(ids)]
  if (length(ids) == 0) {
    stop("No identifiers were provided.", call. = FALSE)
  }

  cid_df <- webchem::get_cid(ids, from = from, match = "first")

  ok <- !is.na(cid_df$cid)
  if (!any(ok)) {
    stop("No CID was found in PubChem for the given identifiers.",
         call. = FALSE)
  }

  props <- webchem::pc_prop(
    cid_df$cid[ok],
    properties = c("CanonicalSMILES", "IsomericSMILES",
                   "MolecularFormula", "IUPACName")
  )

  out <- merge(
    cid_df[ok, c("query", "cid")],
    props,
    by.x = "cid", by.y = "CID",
    all.x = TRUE
  )

  out[order(match(out$cid, cid_df$cid[ok])), ]
}

## ---------------------------------------------------------------------------

#' Calculate molecular descriptors with CDK
#'
#' Parses SMILES strings and calculates the requested molecular descriptors
#' locally using the Chemistry Development Kit (CDK) through the suggested
#' package \pkg{rcdk}. Requires a working Java JDK (a JRE alone is not
#' sufficient).
#'
#' @param smiles Character vector of SMILES strings.
#' @param which Character vector of descriptor short names. Any subset of
#'   \code{c("mw", "alogp", "tpsa", "hbd", "hba", "rotb", "heavy", "aroma",
#'   "mr")}.
#'   By default all are calculated.
#' @return A data.frame with one row per valid molecule and a \code{SMILES}
#'   column.
#' @examples
#' \dontrun{
#' smiles <- c("CCO", "CC(=O)Oc1ccccc1C(=O)O")
#' desc <- calcCDKDescriptors(smiles)
#' }
#' @export
#' @seealso \code{\link{mapCDKDescriptors}}, \code{\link{getSmilesFromIdentifiers}}.
calcCDKDescriptors <- function(smiles, which = c(
  "mw", "alogp", "tpsa", "hbd", "hba", "rotb", "heavy", "aroma", "mr")) {

  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("This function requires the 'rcdk' package (it depends on rJava + a ",
         "JDK, not just a JRE). Install it with: install.packages('rcdk').",
         call. = FALSE)
  }

  which <- match.arg(which, c(
    "mw", "alogp", "tpsa", "hbd", "hba", "rotb", "heavy", "aroma", "mr"
  ), several.ok = TRUE)

  ## Core descriptors: use CDK eval.desc with the standard Java classes.
  ## NOTE: "heavy" is NOT in core_classes because there is no portable CDK
  ## descriptor class for heavy-atom count. AtomCountDescriptor counts ALL
  ## atoms (including H), and HeavyAtomCountDescriptor is not present in all
  ## CDK versions. Instead we compute the heavy-atom count directly from the
  ## parsed IAtomContainer via getAtomCount() BEFORE convert.implicit.to.explicit
  ## is called (when the container holds only the heavy atoms from the SMILES).
  core_classes <- c(
    mw    = "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
    alogp = "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
    mr    = "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
    tpsa  = "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
    hbd   = "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
    hba   = "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
    rotb  = "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
    aroma = "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor"
  )

  core_which <- intersect(which, names(core_classes))

  if (length(core_which) == 0 && !"heavy" %in% which) {
    stop("No valid descriptor was selected.", call. = FALSE)
  }

  smiles_chr <- trimws(as.character(smiles))
  valid_idx <- !is.na(smiles_chr) & nzchar(smiles_chr)
  smiles_valid <- smiles_chr[valid_idx]
  if (length(smiles_valid) == 0) {
    stop("No valid SMILES strings provided.", call. = FALSE)
  }

  mols <- rcdk::parse.smiles(smiles_valid)
  ok <- !vapply(mols, is.null, logical(1))
  if (!any(ok)) {
    stop("No SMILES could be interpreted by CDK.", call. = FALSE)
  }

  ## Heavy-atom count: computed BEFORE convert.implicit.to.explicit(m), so the
  ## IAtomContainer only contains the atoms written in the SMILES. We use
  ## rcdk::get.atoms(m) and filter out any "H" atoms explicitly, which correctly
  ## handles both the common case (SMILES with implicit H, where getAtomCount()
  ## would also work) and the edge case of SMILES with explicit [H] atoms.
  ## After convert.implicit.to.explicit(m), all implicit H would become explicit
  ## atoms and a naive count would over-count by the number of H.
  heavy_counts <- if ("heavy" %in% which) {
    vapply(mols[ok], function(m) {
      tryCatch({
        atoms <- rcdk::get.atoms(m)
        if (length(atoms) == 0) return(NA_integer_)
        syms <- vapply(atoms, function(a) a$getSymbol(), character(1))
        sum(syms != "H")
      }, error = function(e) NA_integer_)
    }, integer(1))
  } else {
    NULL
  }

  ## Mandatory molecule configuration for the other descriptors
  invisible(lapply(mols[ok], function(m) {
    tryCatch({
      rcdk::convert.implicit.to.explicit(m)
      rcdk::do.aromaticity(m)
      rcdk::set.atom.types(m)
    }, error = function(e) {
      warning("Error configuring molecule: ", e$message)
    })
  }))

  ## Core descriptors via eval.desc (only if at least one non-heavy descriptor
  ## was requested)
  if (length(core_which) > 0) {
    classes <- unique(core_classes[core_which])
    desc_df <- rcdk::eval.desc(mols[ok], classes)
  } else {
    ## Only "heavy" was requested; start from an empty data.frame
    desc_df <- data.frame(row.names = seq_len(sum(ok)))
  }

  ## Fix ALOGPDescriptor column names. The ALOGPDescriptor class returns three
  ## columns (ALogP, AMR, ALogP2) as a single atomic unit; we cannot request
  ## one without getting the others. The fallback rename below ensures the
  ## canonical names exist regardless of the casing/variant returned by CDK.
  if ("alogp" %in% which && !"ALogP" %in% names(desc_df)) {
    alogp_col <- grep("alogp", names(desc_df), ignore.case = TRUE, value = TRUE)
    if (length(alogp_col) > 0) names(desc_df)[names(desc_df) == alogp_col[1]] <- "ALogP"
  }
  if ("mr" %in% which && !"AMR" %in% names(desc_df)) {
    amr_col <- grep("amr", names(desc_df), ignore.case = TRUE, value = TRUE)
    if (length(amr_col) > 0) names(desc_df)[names(desc_df) == amr_col[1]] <- "AMR"
  }

  ## Respect the user's descriptor selection: drop side-effect columns that
  ## were not requested. Because ALOGPDescriptor returns ALogP + AMR + ALogP2
  ## together, the user could end up with MR even if "mr" was deselected, or
  ## with LogP even if "alogp" was deselected. Drop these to make the
  ## selection in the UI actually meaningful.
  if (!"alogp" %in% which) {
    drop_cols <- grep("alogp", names(desc_df), ignore.case = TRUE)
    if (length(drop_cols) > 0) desc_df <- desc_df[, -drop_cols, drop = FALSE]
  }
  if (!"mr" %in% which) {
    drop_cols <- grep("^amr$", names(desc_df), ignore.case = TRUE)
    if (length(drop_cols) > 0) desc_df <- desc_df[, -drop_cols, drop = FALSE]
  }
  ## ALogP2 is never used by any plot/filter; always drop it
  if ("ALogP2" %in% names(desc_df)) {
    desc_df$ALogP2 <- NULL
  }

  ## Attach the heavy-atom count (computed manually above). This is the
  ## canonical "nAtom" column expected by mapCDKDescriptors / mapADMETColumns.
  if ("heavy" %in% which) {
    desc_df$nAtom <- heavy_counts
  }

  desc_df$SMILES <- smiles_valid[ok]
  rownames(desc_df) <- NULL
  desc_df
}

## ---------------------------------------------------------------------------

#' Map CDK descriptors to the application's standard schema
#'
#' Translates the raw CDK descriptor names into the application's canonical
#' column names and computes the drug-likeness violation columns and the
#' BOILED-Egg ADMET properties (GI absorption, BBB permeability, P-gp
#' substrate).
#'
#' The CDK \code{ALogP} descriptor is mapped to the generic \code{LogP}
#' column. No artificial \code{MLOGP}/\code{WLOGP}/\code{XLOGP3} columns are
#' created; the application uses a single \code{LogP} column for all
#' drug-likeness filters, with thresholds from the original publications.
#'
#' @param cdk_df A data.frame as returned by \code{\link{calcCDKDescriptors}}.
#' @return A data.frame with renamed columns, violation columns and ADMET
#'   properties.
#' @examples
#' \dontrun{
#' smiles <- c("CCO", "CC(=O)Oc1ccccc1C(=O)O")
#' desc <- calcCDKDescriptors(smiles)
#' mapped <- mapCDKDescriptors(desc)
#' }
#' @export
#' @seealso \code{\link{calcCDKDescriptors}},
#'   \code{\link{computeViolationColumns}}, \code{\link{computeADMETProperties}}.
mapCDKDescriptors <- function(cdk_df) {

  rename_map <- c(
    MW         = "MW",
    ALogP      = "LogP",
    AMR        = "MR",
    TopoPSA    = "TPSA",
    nHBDon     = "#H-bond donors",
    nHBAcc     = "#H-bond acceptors",
    nRotB      = "#Rotatable bonds",
    nAtom      = "#Heavy atoms",
    naAromAtom = "#Aromatic heavy atoms"
  )

  for (old in names(rename_map)) {
    if (old %in% names(cdk_df)) {
      names(cdk_df)[names(cdk_df) == old] <- rename_map[[old]]
    }
  }

  ## Drop ALogP2 (returned by ALOGPDescriptor but not used by any plot/filter)
  if ("ALogP2" %in% names(cdk_df)) {
    cdk_df$ALogP2 <- NULL
  }

  cdk_df <- computeViolationColumns(cdk_df)
  cdk_df <- computeADMETProperties(cdk_df)

  cdk_df
}

## ---------------------------------------------------------------------------

#' Compute drug-likeness violation columns
#'
#' Computes the five drug-likeness \code{"#violations"} columns (Lipinski,
#' Ghose, Veber, Egan, Muegge) from the physicochemical properties, using the
#' thresholds from the original publications. Missing columns are safely
#' filled with \code{NA}.
#'
#' The generic \code{LogP} column is used for all LogP-dependent rules. The
#' original publication thresholds are used:
#' \itemize{
#'   \item Lipinski: LogP > 5 (original Rule-of-Five)
#'   \item Ghose: LogP < -0.4 or LogP > 5.6
#'   \item Egan: LogP > 5.88
#'   \item Muegge: LogP < -2 or LogP > 5
#' }
#'
#' The Veber violation is binary: a compound violates if it has more than 10
#' rotatable bonds OR if both polarity conditions fail (TPSA > 140 AND
#' HBA + HBD > 12).
#'
#' @param data A data.frame with physicochemical property columns.
#' @return A data.frame with the added violation columns.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = 300, LogP = 2, TPSA = 80, MR = 90,
#'   "#H-bond acceptors" = 4, "#H-bond donors" = 2,
#'   "#Rotatable bonds" = 3, "#Heavy atoms" = 22, check.names = FALSE)
#' d <- computeViolationColumns(d)
#' }
#' @export
#' @seealso \code{\link{lipinskiFilter}}, \code{\link{ghoseFilter}},
#'   \code{\link{veberFilter}}, \code{\link{eganFilter}},
#'   \code{\link{mueggeFilter}}.
computeViolationColumns <- function(data) {

  safe_get <- function(col) {
    if (col %in% names(data)) {
      as.numeric(data[[col]])
    } else {
      rep(NA_real_, nrow(data))
    }
  }

  mw    <- safe_get("MW")
  logp  <- safe_get("LogP")
  tpsa  <- safe_get("TPSA")
  mr    <- safe_get("MR")
  hba   <- safe_get("#H-bond acceptors")
  hbd   <- safe_get("#H-bond donors")
  rb    <- safe_get("#Rotatable bonds")
  ha    <- safe_get("#Heavy atoms")

  ## Lipinski (1997): MW > 500, LogP > 5, HBA > 10, HBD > 5
  data[["Lipinski #violations"]] <-
    (mw > 500) + (logp > 5) + (hba > 10) + (hbd > 5)

  ## Ghose (1999): MW 160-480, MR 40-130, LogP -0.4 to 5.6, heavy atoms 20-70
  data[["Ghose #violations"]] <-
    (mw < 160 | mw > 480) + (mr < 40 | mr > 130) +
    (logp < -0.4 | logp > 5.6) + (ha < 20 | ha > 70)

  ## Veber (2002): RB <= 10 AND (TPSA <= 140 OR HBA+HBD <= 12).
  ## Violation = RB > 10 OR (TPSA > 140 AND HBA+HBD > 12). Binary (0 or 1).
  data[["Veber #violations"]] <-
    as.integer((rb > 10) | (tpsa > 140 & (hba + hbd) > 12))

  ## Egan (2000): TPSA > 131.6, LogP > 5.88
  data[["Egan #violations"]] <-
    (tpsa > 131.6) + (logp > 5.88)

  ## Muegge (2001): MW 200-600, LogP -2 to 5, HBA > 10, HBD > 5,
  ##   RB > 15, TPSA > 150, pharmacophore points < 4.
  ## Pharmacophore points are approximated as HBA + HBD.
  pharma_points <- hba + hbd
  pharma_violation <- as.integer(pharma_points < 4)

  data[["Muegge #violations"]] <-
    (mw < 200 | mw > 600) + (logp < -2 | logp > 5) + (hba > 10) +
    (hbd > 5) + (rb > 15) + (tpsa > 150) + pharma_violation

  data
}

## ---------------------------------------------------------------------------

#' Compute BOILED-Egg ADMET properties
#'
#' Computes the gastrointestinal absorption (\code{GI absorption}), blood-brain
#' barrier permeability (\code{BBB permeant}) and P-glycoprotein substrate
#' (\code{Pgp substrate}) categorical columns from the \code{LogP} and
#' \code{TPSA} values, using the official BOILED-Egg polygon coordinates
#' (Daina & Zoete, 2016, Data S3) for GI/BBB classification via
#' point-in-polygon testing, and a Random Forest classifier for P-gp
#' substrate prediction.
#'
#' The BOILED-Egg model was originally calibrated with WLOGP; here the
#' application's generic \code{LogP} column is used (which may be WLOGP,
#' Consensus Log P, ALogP, or the source platform's own LogP). This is an
#' acceptable approximation for exploratory visualization.
#'
#' @param data A data.frame with \code{LogP}, \code{TPSA}, \code{MW},
#'   \code{#H-bond donors} and \code{#H-bond acceptors} columns.
#' @return A data.frame with the added ADMET property columns.
#' @examples
#' \dontrun{
#' d <- data.frame(LogP = 2, TPSA = 80, MW = 300,
#'   "#H-bond donors" = 2, "#H-bond acceptors" = 4, check.names = FALSE)
#' d <- computeADMETProperties(d)
#' }
#' @export
#' @references Daina, A., & Zoete, V. (2016). A boiled egg to predict
#'   gastrointestinal absorption and brain penetration of small molecules.
#'   \emph{ChemMedChem}, 11(11), 1117-1121.
#' @references Sedykh, A., Fourches, D., Duan, J., et al. (2013). Human
#'   intestinal transporter database: QSAR modeling and virtual excretion
#'   experiments. \emph{J. Cheminformatics} 2015, 7:21 (Metrabase P-gp data).
computeADMETProperties <- function(data) {

  safe_get <- function(col) {
    if (col %in% names(data)) {
      as.numeric(data[[col]])
    } else {
      rep(NA_real_, nrow(data))
    }
  }

  logp <- safe_get("LogP")
  tpsa <- safe_get("TPSA")
  mw   <- safe_get("MW")
  hbd  <- safe_get("#H-bond donors")
  hba  <- safe_get("#H-bond acceptors")
  rb   <- safe_get("#Rotatable bonds")
  ha   <- safe_get("#Heavy atoms")
  arom <- safe_get("#Aromatic heavy atoms")
  mr   <- safe_get("MR")

  ## Select the appropriate BOILED-Egg polygon based on the LogP source.
  ## If the data has a WLOGP column, use the official WLOGP
  ## polygons from Daina & Zoete (2016). Otherwise (CDK, ADMET Master,
  ## which use ALogP or their own logP), use the ALogP-trained polygons.
  use_alogp <- !"WLOGP" %in% names(data)
  hia_poly <- if (use_alogp) .hia_polygon_alogp else .hia_polygon
  bbb_poly <- if (use_alogp) .bbb_polygon_alogp else .bbb_polygon

  ## Classify points using point-in-polygon
  ## Polygons are stored as (TPSA, LogP) = (x, y) in R/boiled_egg_data.R
  hia_inside <- ..point_in_polygon(tpsa, logp, hia_poly)
  data[["GI absorption"]] <- ifelse(
    is.na(hia_inside), NA,
    ifelse(hia_inside, "High", "Low")
  )

  bbb_inside <- ..point_in_polygon(tpsa, logp, bbb_poly)
  data[["BBB permeant"]] <- ifelse(
    is.na(bbb_inside), NA,
    ifelse(bbb_inside, "Yes", "No")
  )

  ## 3. P-gp substrate prediction
  ##    Uses a Random Forest classifier (100 trees, max depth = 10) trained
  ##    from 882 experimental P-gp (ABCB1) substrate/non-substrate
  ##    classifications curated from Metrabase (J. Cheminformatics 2015, 7:21).
  ##    The model uses the 9 CDK descriptors and achieves 69.4% cross-validated
  ##    accuracy with MCC = 0.383 (5-fold stratified CV). Stored as pure R
  ##    code in R/pgp_model.R; no external ML package required at runtime.
  pgp_prob <- ..predict_pgp(mw, logp, tpsa, hbd, hba, rb,
                            ha, arom, mr)
  data[["Pgp substrate"]] <- ifelse(
    is.na(pgp_prob), NA,
    ifelse(pgp_prob >= 0.5, "Yes", "No")
  )
  data[["Pgp probability"]] <- pgp_prob

  data
}
