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

### Asymmetric measures (directional, Y is the outcome)

| Method | Description | Scale | Independence value |
|--------|-------------|-------|-------------------|
| AKC | Asymmetric Kendall Correlation | [-1, 1] | 0 |
| AGC | Asymmetric Grade Correlation | [-1, 1] | 0 |
| CID | Concordance-Discordance Index = (AKC + 1) / 2 | [0, 1] | 0.5 |
| CMA | Coefficient of Monotone Association = (AGC + 1) / 2 | [0, 1] | 0.5 |

### Symmetric measures

| Method | Description | Scale | Independence value |
|--------|-------------|-------|-------------------|
| tau_a | Kendall's tau-a (no tie correction) | [-1, 1] | 0 |
| tau_b | Kendall's tau-b (pair-based tie correction) | [-1, 1] | 0 |
| tau_b_mod | Modified tau-b (triple-based tie correction) | [-1, 1] | 0 |
| gamma | Goodman-Kruskal gamma | [-1, 1] | 0 |
| rho_a | Spearman's rho (no tie correction) | [-1, 1] | 0 |
| rho_b | Spearman's rho (with tie correction) | [-1, 1] | 0 |
| pearson | Pearson correlation | [-1, 1] | 0 |

