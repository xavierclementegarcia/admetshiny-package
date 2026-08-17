# admetshiny

<!-- badges: start -->
<!-- badges: end -->

An open-source R package and Shiny application for the management,
calculation, filtering, visualization and exploratory analysis of molecular
descriptors and ADMET properties of small molecules.

## Overview

**admetshiny** integrates cheminformatics and bioinformatics workflows with an
intuitive dashboard to support the prioritization of compounds in early-stage
drug discovery. It integrates data exported from:

### Features

- Drug-likeness filters: **Lipinski**, **Veber**, **Ghose**, **Egan**,
  **Muegge**.
- **BOILED-Egg** model (Daina & Zoete, 2016) for GI absorption and BBB
  permeability.
- Chemical-space visualizations: **PCA**, **t-SNE**, **parallel coordinates**,
  **violin**, **radar** and **correlation heatmap**.
- **Tanimoto / AGNES** structural-similarity clustering.
- Dark mode, custom colour palettes and high-resolution figure export
  (PNG, PDF, SVG, TIFF up to 1200 DPI).

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("xavierclementegarcia/admetshiny-package")
```

The optional CDK-based descriptor calculation requires the `rcdk` package,
which in turn requires **Java JDK** (a JRE alone is not sufficient).

## Usage

Launch the interactive application:

```r
admetshiny::run_app()
```

## References

- Daina, A., Michielin, O., & Zoete, V. (2017). *Scientific Reports*, 7(1), 42717.
- Daina, A., & Zoete, V. (2016). *ChemMedChem*, 11(11), 1117-1121.
- Fu, L., et al. (2024). *Nucleic Acids Research*, 52(W1), W422-W431.
- Myung, Y., de Sá, A. G., & Ascher, D. B. (2024). *Nucleic Acids Research*, 52(W1), W469-W475.

## License

MIT License
