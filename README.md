# mancova_hocking

R code for the paper

> Yoon, S.-H., Busababodhin, P., and Park, J.-S.
> *Variable selection theorem for the multivariate analysis of covariance model.*

The paper extends the variable selection theorem of Hocking (1976) to the
multivariate analysis of covariance (MANCOVA) model with a
less-than-full-rank treatment design matrix. Omitting covariates reduces the
covariance matrices of the estimators of the retained covariate coefficients,
of every estimable function of the treatment parameters, and of the predicted
response vector, and the mean squared error (MSE) matrices also decrease when
tr(Delta) <= 1, where Delta is the noncentrality matrix of the multivariate
test of the omitted coefficients.

## Files

| File | Content |
|---|---|
| `R/rohwer_example.R` | Reproduces Section 5 and Table 1 (Rohwer data): model checks, the test of H0: B_r = 0, the unbiased, positive-part, and maximum likelihood estimates of tr(Delta), the traces of the MSE matrices, the crossing values, and the estimators of Sigma. |
| `R/formula_check.R` | Additional check, not reported in the paper: exact traces of the MSE matrices for a small simulated one-way MANCOVA design, compared with a Monte Carlo estimate. |

## Requirements

- R (>= 4.0)
- Packages: `heplots` (Rohwer data), `car` (MANOVA tests), `MASS` (generalized inverse)

```r
install.packages(c("heplots", "car"))   # MASS ships with R
```

## Usage

From the top directory of the repository:

```r
source("R/rohwer_example.R")   # Section 5 and Table 1
source("R/formula_check.R")    # optional check of the formulas
```

## Data

The Rohwer data (Timm, 1975, Examples 4.3, 4.7, and 4.13) are available in the
R package `heplots` (Friendly, 2007).

## References

- Friendly, M. (2007). HE plots for multivariate linear models. *Journal of Computational and Graphical Statistics*, 16, 421-444.
- Hocking, R. R. (1976). The analysis and selection of variables in linear regression. *Biometrics*, 32, 1-49.
- Timm, N. H. (1975). *Multivariate Analysis with Applications in Education and Psychology*. Brooks/Cole, Monterey, CA.

## License

GPL-3. See `LICENSE`.
