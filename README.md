# acor

R package for **Asymmetric Correlation Measures**: AKC, AGC, CID, and CMA.

## Installation

```r
devtools::install_github("evwalz/acor")
```

## Usage

```r
# Compute correlation
x <- rnorm(100)
y <- rnorm(100)
acor(x, y, method = "akc")   # Asymmetric Kendall Correlation
acor(x, y, method = "cma")   # Coefficient of Monotone Association

# Statistical test
acor.test(x, y, method = "akc", alternative = "two.sided")

# Multiple predictors
X <- cbind(rnorm(100), rnorm(100))
acor.test(X, y, method = "agc")
```

## Methods

| Method | Scale | Independence value |
|--------|-------|---------------------|
| AKC, AGC | [-1, 1] | 0 |
| CID, CMA | [0, 1] | 0.5 |

- **AKC** = Asymmetric Kendall Correlation (Kendall framework)
- **AGC** = Asymmetric Grade Correlation (Spearman framework)
- **CID** = (AKC + 1) / 2
- **CMA** = (AGC + 1) / 2

## Structure

- `R/acor_functions.R`: Main API (`acor`, `acor.test`)
- `R/acor_internal.R`: Internal helpers
- `R/akc_functions.R`: AKC/Kendall implementation
- `R/agc_functions.R`: AGC/Spearman implementation

