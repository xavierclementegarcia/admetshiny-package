# admetshiny

<!-- badges: start -->
<!-- badges: end -->

An open-source R package and Shiny application for the management,
calculation, filtering, visualization and exploratory analysis of molecular
descriptors and ADMET properties of small molecules.

## Overview

**admetshiny** integrates cheminformatics and bioinformatics workflows with an
intuitive dashboard to support the prioritization of compounds in early-stage
drug discovery. The application is organised in two complementary modules:

- **CDK & webchem** — retrieve canonical SMILES from PubChem (via the
  `webchem` package), enter them manually, or upload them as a CSV, then
  compute nine physicochemical descriptors locally with the Chemistry
  Development Kit (CDK).
- **ADMET Master Manager** — upload any CSV or Excel (`.xlsx`) ADMET dataset
  and manually map its columns to the application's 20-field standard schema.
  Missing descriptors are back-filled from SMILES via CDK when available.

Both modules share the same drug-likeness filters, BOILED-Egg model, P-gp
substrate classifier and 14-chart catalogue.

### Features

- Drug-likeness filters with the original publication thresholds: **Lipinski**
  (1997), **Veber** (2002), **Ghose** (1999), **Egan** (2000), **Muegge**
  (2001).
- **BOILED-Egg** model (Daina & Zoete, 2016) for gastrointestinal absorption
  and blood-brain barrier permeability, with dual WLOGP / ALogP polygon
  systems.
- **P-glycoprotein (ABCB1) substrate** prediction via an embedded Random
  Forest classifier (100 trees, 5-fold CV: 69.4% accuracy, MCC = 0.383).
- Chemical-space visualizations: **PCA**, **t-SNE**, **UMAP**, **parallel
  coordinates**, **violin**, **radar**, **correlation heatmap**, **cluster
  heatmap** and **Tanimoto / AGNES** structural-similarity clustering.
- Dark mode, custom colour palettes and high-resolution figure export
  (PNG, PDF, SVG, TIFF up to 1200 DPI).
- Comprehensive HTML / PDF / Word report generation.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("xavierclementegarcia/admetshiny")
```

The optional CDK-based descriptor calculation requires the `rcdk` package,
which in turn requires **Java JDK** (a JRE alone is not sufficient).

## Usage

Launch the interactive application:

```r
admetshiny::run_app()
```

Use the computational functions programmatically:

```r
library(admetshiny)

# Compute CDK descriptors for a few SMILES and apply the Lipinski filter
smiles <- c("CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
desc  <- calcCDKDescriptors(smiles)
desc  <- mapCDKDescriptors(desc)               # adds #violations and ADMET cols
filtered <- applyFilters(desc, filters = c("Lipinski", "Veber", "Ghose"))

# Plot the BOILED-Egg
plotBoiledEgg(filtered)

# Or normalize any external ADMET dataset via manual column mapping
# (returns the standard schema with #violations + ADMET properties):
# d <- read.csv("my_admet.csv", check.names = FALSE)
# mapping <- setNames(c("SMILES", "MW", "LogP", "TPSA"),
#                     c("CanonicalSMILES", "MW", "iLOGP", "TPSA"))
# d <- mapADMETColumns(d, mapping, calculate_cdk = TRUE)
```

## References

- Lipinski, C. A., et al. (1997). *Advanced Drug Delivery Reviews*, 23(1-3),
  3-25.
- Ghose, A. K., et al. (1999). *J. Combinatorial Chemistry*, 1(1), 55-68.
- Veber, D. F., et al. (2002). *J. Medicinal Chemistry*, 45(12), 2615-2623.
- Egan, W. J., et al. (2000). *J. Medicinal Chemistry*, 43(21), 3867-3877.
- Muegge, I., et al. (2001). *J. Medicinal Chemistry*, 44(12), 1841-1846.
- Daina, A., & Zoete, V. (2016). *ChemMedChem*, 11(11), 1117-1121.
  (BOILED-Egg model)

## License

CC BY 4.0
