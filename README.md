# ImpactEffectsize

An R package for calculating and visualizing the Impact effect size measure, which combines central tendency differences with morphological differences in data distribution shapes between two groups.

## Overview

The Impact effect size measure provides a comprehensive assessment of differences between two groups by considering both:
- **Central tendency component (CT)**: Differences in group medians normalized by the Gini Mean Difference
- **Morphological component (Morph)**: Differences in the shapes of probability density functions using Pareto Density Estimation

This dual approach makes Impact particularly effective for detecting group differences that traditional effect size measures might miss, especially when groups have similar central tendencies but different distribution shapes.

## Features

- **Dual-component analysis**: Combines location and shape-based differences
- **Robust median-based calculations**: Less sensitive to outliers than mean-based measures
- **Pareto Density Estimation**: Non-parametric approach for capturing distribution shape differences
- **Built-in visualization**: Optional density plot generation with customizable aesthetics
- **Automatic data validation**: Handles missing values and validates input requirements
- **Kolmogorov-Smirnov integration**: Uses KS test to assess distributional differences


## Installation

You can install ImpactEffectsize from CRAN:

```r
install.packages("ImpactEffectsize")
```
Or install the development version directly from GitHub:
``` r
devtools::install_github("JornLotsch/ImpactEffectsize")
```
For more information, visit the CRAN package page: [https://cran.r-project.org/package=ImpactEffectsize](https://cran.r-project.org/package=ImpactEffectsize)

## Dependencies

The package requires:
- `caTools` (for numerical integration)
- `methods` (for argument checking)
- `matrixStats` (for matrix operations)
- `stats` (for statistical functions)

## Quick start

### Basic usage
```r 
library(ImpactEffectsize)
# Load example data
data("FeatureselectionData")
# Calculate Impact effect size
result <- Impact(Data = FeatureselectionDataVar0011, Cls = FeatureselectionDataClasses)
# View results
print(resultImpact) # Main effect size measure print(resultCTDiff) # Central tendency component print(result$MorphDiff) # Morphological component
``` 

### With visualization
```r
# Calculate Impact with density plot
result <- Impact(Data = FeatureselectionDataVar0011, Cls = FeatureselectionDataClasses, PlotIt = TRUE, col = c("red", "blue"), meanLines = TRUE, medianLines = TRUE)
``` 

## Main function: `Impact()`

### Call
```r 
result <- Impact(Data, Cls, PlotIt = FALSE, pde = TRUE, col = c("red", "blue"), meanLines = FALSE, medianLines = FALSE, ...)
``` 

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `Data` | numeric vector | - | The data of both groups as a vector |
| `Cls` | vector | - | Class information vector of similar length as Data |
| `PlotIt` | logical | FALSE | Whether to plot the probability density functions |
| `pde` | logical | TRUE | Use Pareto density estimation (TRUE) vs standard density (FALSE) |
| `col` | character vector | c("red", "blue") | Colors for the two density lines in plots |
| `meanLines` | logical | FALSE | Add perpendicular lines at group means |
| `medianLines` | logical | FALSE | Add perpendicular lines at group medians |
| `...` | - | - | Additional graphical parameters for plotting |

### Output

| Output | Type | Description |
|--------|------|-------------|
| `Impact` | numeric | The main Impact effect size measure |
| `CTDiff` | numeric | Central tendency difference component |
| `MorphDiff` | numeric | Morphological difference component |

## Examples

### Example 1: Basic Impact calculation
```r 
data("FeatureselectionData") result <- Impact(Data = FeatureselectionDataVar0011, Cls = FeatureselectionDataClasses, PlotIt = TRUE)
cat("Impact Effect Size:", resultImpact, "\n") cat("CT Component:", resultCTDiff, "\n") cat("Morphological Component:", result$MorphDiff, "\n")
``` 

### Example 2: Custom visualization
```r 
data("BcellLymphomaCD79") result <- Impact(Data = BcellLymphomaCD79SomeVariable, Cls = BcellLymphomaCD79Classes, PlotIt = TRUE, col = c("darkred", "darkblue"), meanLines = TRUE, medianLines = TRUE, main = "Impact Analysis with Reference Lines")
```

## Interpretation

The Impact effect size measure ranges from 0 to positive values:
- **0**: No meaningful difference between groups
- **> 0**: Increasing difference between groups
- **Positive/Negative values**: Direction indicates which group has higher central tendency

The measure combines:
- **CTDiff**: Captures location differences (similar to standardized mean difference)
- **MorphDiff**: Captures shape differences not detected by location-based measures

## License

GPL-3

## Citation

Lötsch J, Ultsch A. A non-parametric effect-size measure capturing changes in central tendency and data distribution shape. PLoS One. 2020 Sep 24;15(9):e0239623. doi: 10.1371/journal.pone.0239623. PMID: 32970758; PMCID: PMC7514071.

