# ExpoRiskR

ExpoRiskR is a Bioconductor package for **exposure‑aware multi‑omics integration**.  
It provides a unified framework to quantify exposure perturbation effects, build exposure‑adjusted networks, and derive interpretable risk scores from integrated omics data.

---

## Installation

Install the released version from **Bioconductor**:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ExpoRiskR")
```

You can install the development version from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("ppchaudhary/ExpoRiskR")
```

---

## Overview

ExpoRiskR enables:

- Exposure‑aware multi‑omics integration  
- Exposure perturbation scoring  
- Network construction adjusted for environmental effects  
- Feature importance ranking  
- Individual risk profiling  
- Visualization tools for interpretable results  

The package is designed to work seamlessly with **Bioconductor data structures** such as `SummarizedExperiment`.

---

## Quick Example

```r
library(ExpoRiskR)

# Generate example data
data <- generate_dummy_exporisk()

# Align omics data
aligned <- align_omics(
  data$omics,
  data$metadata,
  exposure_col = "exposure"
)

# Compute exposure perturbation score
scores <- compute_exposure_scores(aligned)

# Build exposure‑adjusted network
net <- build_exposure_network(scores)

plot_network(net)
```

---

## Documentation

Full tutorials and workflow examples are available in the package vignette:

```r
browseVignettes("ExpoRiskR")
```

---

## Citation

If you use ExpoRiskR in your research, please cite:

> Chaudhary PP et al. *ExpoRiskR: Exposure‑aware multi‑omics risk modeling framework.*

---

## Support

For questions, issues, or feature requests:

- Bioconductor Support Site: https://support.bioconductor.org  
- GitHub Issues: https://github.com/ppchaudhary/ExpoRiskR/issues
