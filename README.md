# Coastal Erosion Prediction - Machine Learning & Statistical Modeling

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Predictive modeling of coastal rockfall events using meteorological and oceanographic indicators. This project develops statistical models to estimate annual erosion rates and assess risk levels on coastal cliffs.

## Project Overview

This study analyzes 15 years of coastal erosion data (2001-2022) across 5 temporal periods with varying durations. The objective is to identify which meteorological, oceanic, and atmospheric indicators best explain rockfall frequency and build predictive models capable of estimating future erosion rates.

**Key Challenge:** With only 5 observation periods and 70+ potential predictors, the primary methodological challenge is avoiding overfitting while extracting meaningful relationships.

## Objectives

1. **Feature Selection:** Identify the most predictive indicators among 70+ meteorological and marine variables
2. **Predictive Modeling:** Develop robust models to estimate annual rockfall rates
3. **Risk Classification:** Classify periods into risk categories (low/medium/high)
4. **Statistical Rigor:** Handle small sample size, temporal normalization, and multicollinearity issues

## Methodology

### Part 1: Feature Selection & Correlation Analysis

**Problem:** Variables representing cumulative counts (days of frost, number of storms) mechanically increase with period duration, creating spurious correlations.

**Solution:** Temporal normalization by dividing count variables by period duration to obtain annual rates.

**Metrics Computed:**
- **Pearson correlation:** Linear relationship strength
- **Spearman correlation:** Monotonic relationship (rank-based)
- **Poisson regression (univariate):** Direct effect estimation
- **Pseudo-R² (McFadden):** Explanatory power of each indicator

**Results:** Identified top 15 predictive variables combining correlation and pseudo-R² scores.

### Part 2: Predictive Modeling with LOOCV

**Three approaches tested to avoid overfitting:**

1. **LASSO Regularization (L1):** Automatic feature selection through coefficient penalization
2. **Simple Model (2 variables max):** Following the empirical rule n/3 for small samples
3. **Composite Score:** Single synthetic indicator combining multiple standardized variables

**Validation:** Leave-One-Out Cross-Validation (LOOCV) to assess generalization capability with limited data.

## Key Results

### Top 15 Predictive Variables

1. nb_jours_pmermin - Low pressure days
2. rafale_max_kmh - Maximum wind gust
3. nb_seq_depression_3j - 3-day depression sequences
4. jours_vent_fort_60 - Strong wind days (>60 km/h)
5. plus_longue_serie_tempete_consecutive - Longest consecutive storm sequence
6. t02_max - Maximum wave period
7. marnage_std_m - Tidal range standard deviation
8. energie_vent_cumulee - Cumulative wind energy
9. nb_seq_seche_10j - 10-day dry sequences
10. pression_min_hpa - Minimum atmospheric pressure
11. nb_jours_vent_dir_ouest - West wind days
12. t02_moy - Average wave period
13. ifm - High tide index
14. jours_pluie - Rainy days
15. nb_combinaisons_critiques - Critical event combinations

### Model Performance Comparison

| Model | MAE (rockfalls/year) | Relative Error | R² | Classification Accuracy |
|-------|---------------------|----------------|-----|------------------------|
| LASSO | 6.4 | 23.5% | 0.72 | 40% |
| Simple (2 var) | 12.0 | 45.7% | 0.90 | 20% |
| Composite | 6.7 | 23.8% | 0.79 | 0% |

**Recommended Model:** Composite score achieves the best balance with 6.7 rockfalls/year MAE and 24% relative error.

### Key Findings

- Meteorological and marine indicators explain approximately 75% of rockfall variability despite limited sample size
- Storm-related variables are dominant predictors: wind gusts, storm sequences, and atmospheric depressions show the strongest relationships
- Prediction accuracy is reasonable (24% relative error) for estimating annual erosion rates
- Fine-grained classification into 3 risk categories fails with only 5 periods (accuracy 0-40%)
- Binary classification (low vs high risk) would be more realistic given the data constraints

## Technical Stack

**Languages & Tools:**
- R (ggplot2, glmnet, readxl, dplyr, corrplot)
- Statistical modeling (Poisson regression, LASSO)
- Machine Learning (Cross-validation, regularization)
- Data visualization

**Key Techniques:**
- Feature engineering and temporal normalization
- Correlation analysis (Pearson, Spearman)
- Regularization methods (L1 penalty)
- Leave-One-Out Cross-Validation (LOOCV)
- Statistical bias diagnosis (multicollinearity, spurious correlations)

## Getting Started

### Prerequisites

```r
install.packages(c("readxl", "ggplot2", "dplyr", "tidyr", 
                   "corrplot", "glmnet"))
```

### Running the Analysis

```r
# Step 1: Feature selection and correlation analysis
source("scripts/01_selection_variables.R")

# Step 2: LOOCV predictive modeling
source("scripts/02_loocv_regularise.R")
```

### Outputs

The analysis generates:
- 11 figures from feature selection (correlation distributions, top variables, relationships)
- 8 figures from LOOCV models (comparative performance, classification, error analysis)
- CSV files with detailed metrics and predictions

All outputs are saved in the `outputs/` directory.

## Limitations & Future Work

**Current Limitations:**
- Small sample size (n=5 periods) limits statistical power
- Unable to capture non-linear interactions or threshold effects
- Geomorphological and geotechnical variables not included
- Temporal dependencies between successive periods not modeled

**Potential Improvements:**
- Increase observations by subdividing periods (with autocorrelation management)
- Integrate spatial and geotechnical features
- Explore non-linear models (Random Forests, Gradient Boosting)
- Investigate critical thresholds and variable interactions
- Implement binary classification (low vs high risk) instead of 3-class system

## License

This project is licensed under the MIT License - see the LICENSE file for details.


*This project demonstrates advanced statistical modeling, feature engineering, and machine learning techniques applied to a real-world environmental prediction problem with inherent data limitations.*
