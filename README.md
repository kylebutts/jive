
<!-- README.md is generated from README.Rmd. Please edit that file -->

## \[🚧 WIP 🚧\] `jive`

<!-- badges: start -->
<!-- badges: end -->

The goal of jive is to implement jackknife instrumental-variable
estimators (JIVE) and various alternatives.

## Installation

You can install the development version of jive like so:

``` r
devtools::install_github("kylebutts/jive") 
```

## Example Usage

We are going to use the data from Stevenson (2018). Stevenson leverages
the quasi-random assignment of 8 judges (magistrates) in Philadelphia to
study the effects pretrial detention on several outcomes, including
whether or not a defendant subsequently pleads guilty.

``` r
library(jive)
data(stevenson)
```

### Juke n’ JIVE

``` r
jive(
  guilt ~ i(black) + i(white) | bailDate | jail3 ~ 0 | judge_pre,
  data = stevenson
)
#> Coefficients: 
#>         Estimate  Robust SE Z value   Pr(>z)   
#> jail3 -0.0218451 -0.0075172   2.906 0.003661 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 331,971 observations, 7 instruments, 2,352 covariates
#> First-stage F: stat = 32.626
#>        Sargan: stat = 3.342, p = 0.765
#>            CD: stat = 3.319, p = 0.768
```

``` r
ujive(
  guilt ~ i(black) + i(white) | bailDate | jail3 ~ 0 | judge_pre,
  data = stevenson
)
#> Coefficients: 
#>       Estimate Robust SE Z value  Pr(>z)  
#> jail3 0.159091  0.070564  2.2546 0.02416 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 331,971 observations, 7 instruments, 2,352 covariates
#> First-stage F: stat = 32.626
#>        Sargan: stat = 3.342, p = 0.765
#>            CD: stat = 3.319, p = 0.768
```

``` r
ijive(
  guilt ~ i(black) + i(white) | bailDate | jail3 ~ 0 | judge_pre,
  data = stevenson
)
#> Coefficients: 
#>       Estimate Robust SE Z value  Pr(>z)  
#> jail3 0.159529  0.070533  2.2618 0.02371 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 331,971 observations, 7 instruments, 2,352 covariates
#> First-stage F: stat = 32.626
#>        Sargan: stat = 3.342, p = 0.765
#>            CD: stat = 3.319, p = 0.768
```

### (Leave-out) Leniency Measures

The package will allow you to estimate (leave-out) leniency measures:

``` r
# TODO :-) 
# I have the function, just thinking of the API
```

## Econometric Details on JIVE, UJIVE, IJIVE, and CJIVE

Consider the following instrumental variables setup

$$
  y_i = T_i \beta + W_i' \psi + \varepsilon_i,
$$

where $T_i$ is a scalar endogenous variable, $W_i$ is a vector of
exogenous covariates, and $\varepsilon_i$ is an error term.
Additionally, assume there is a set of valid instruments, $Z_i$. The
intuition of two-stage least squares is to use the instruments to
“predict” $T_i$:

$$
  T_i = Z_i' \pi + W_i \gamma + \eta_i.
$$

Then, the prediction, $T_i$, is used in place of $T_i$ in the original
regression.

When the dimension of $Z_i$ grows with the number of observations,
two-stage least squares is biased (Kolesar, 2013). Without getting into
the details, the problem comes because you’re predicting $T_i$ using
$i$’s own observation in the first-stage. For this reason, Angrist,
Imbens, and Krueger (1999) developed the jackknife instrumental-variable
estimator (JIVE). In short, for each $i$, $T_i$ is predicted using the
first-stage equation leaving out $i$’s own observation.

In general, the JIVE estimator (and variants) are given by

$$
  \frac{P' Y}{P' T}
$$

where $P$ is a function of $W$ and $T_i$. The particulars differ across
the JIVE, the unbiased JIVE (UJIVE), the improved JIVE (IJIVE), and the
cluster JIVE (CJIVE).

### JIVE definition

**Source:** Kolesar (2013) and Angrist, Imbens, and Kreuger (1999)

The original JIVE estimate produces $T$ given by

$$
  T_{JIVE} = (I - D_{(Z,W)})^{-1} (H_{(Z,W)} - D_{(Z,W)}) T,
$$

where $H_{(Z,W)}$ is the hat/projection matrix for $(Z,W)$ and
$D_{(Z,W)}$ is the diagonal matrix with the diagonal elements from
$H_{Z, W}$.

$$
  P_{JIVE} = M_W T_{JIVE}
$$

### UJIVE definition

**Source:** Kolesar (2013)

$$
  T_{UJIVE} = (I - D_{(Z,W)})^{-1} (H_{(Z,W)} - D_{(Z,W)}) T = T_{JIVE}
$$

$$
  P_{UJIVE} = T_{UJIVE} - (I - D_{W})^{-1} (H_{W} - D_{W}) T
$$

### IJIVE definition

**Source:** Ackerberg and Devereux (2009)

Start residualizing $T$, $Y$, and $Z$ by the covariates $W$. The
definition they give (assuming things have already been residualized) is

$$
  T_{IJIVE} = P_{IJIVE} = (I - D_{\tilde{Z}})^{-1} (H_{\tilde{Z}} - D_{\tilde{Z}}) \tilde{T}
$$

Note that $P = T$ because you have already residualized by $W$.

### CJIVE definition

**Source:** Frandsen, Leslie, and McIntyre (2023)

This is a modified version of IJIVE as proposed by Frandsen, Leslie, and
McIntyre (2023). This is necessary if the errors are correlated within
clusters (e.g. court cases assigned on the same day to the same judge).
The modified version is given by:

$$
  T_{CJIVE} = P_{CJIVE} = (I - \mathbb{D}(P_Z, \{ n_1, \dots, n_G \}))^{-1} (H_{\tilde{Z}} - D_{\tilde{Z}}) \tilde{T},
$$

where $\mathbb{D}(P_Z, \{ n_1, \dots, n_G \})$ is a block-diagonal
matrix equal to the projection matrix of $Z$ zerod out except for each
cluster.

In this package, the same adjustment (replacing the diagonal $D$ with
the cluster block-diagonal $\mathbb{D}$ version) can be done for all
three estimators, though the details have only been worked out for IJIVE
in particular (free research idea).

### Standard errors

Heteroskedastic-robust standard errors are given by

$$ 
  \frac{\sqrt{\sum_i P_i^2 \hat{\epsilon}_i^2}}{\sum_i P_i T_i},
$$

where $\hat{\epsilon}_i = M_W Y_i - M_W T_i * \hat{\beta}$. See
<https://github.com/kolesarm/ivreg/blob/master/ivreg.pdf> for
derivations.

### Advantages of alternative estimators:

Quoting from the papers that propose each:

- UJIVE: “UJIVE is consistent for a convex combination of local average
  treatment effects under many instrument asymptotics that also allow
  for many covariates and heteroscedasticity”
- IJIVE: “We introduce two simple new variants of the jackknife
  instrumental variables (JIVE) estimator for overidentified linear
  models and show that they are superior to the existing JIVE estimator,
  significantly improving on its small-sample-bias properties”
- CJIVE: “In settings where inference must be clustered, however, the
  \[IJIVE\] fails to eliminate the many-instruments bias. We propose a
  cluster-jackknife approach in which first-stage predicted values for
  each observation are constructed from a regression that leaves out the
  observation’s entire cluster, not just the observation itself. The
  cluster-jackknife instrumental variables estimator (CJIVE) eliminates
  many-instruments bias, and consistently estimates causal effects in
  the traditional linear model and local average treatment effects in
  the heterogeneous treatment effects framework.”

## Computation Efficiency

This package uses a ton of algebra tricks to speed up the computation of
the JIVE (and variants). First, the package is written using `fixest`
which estiamtes high-dimensional fixed effects very quickly. All the
credit to @lrberge for this.

All uses of the projection matrix, $H$, can be efficiently calculated
using the fitted values of regressiones. The diagonal of the hat matrix
needs to be calculated using $X_i (X'X)^{-1} X_i'$ which could be quite
slow when there are high-dimensional fixed effects. However, we use
sparse matrices and the following trick. Partition the covariates
$W = [A, B]$ where $B$ are the fixed effects. Then the projection matrix
can be decomposed as $P_W = P_B + P_{M_B A}$. $P_B$ is quick to
calculate since it’s a bunch of dummy variables. $P_{M_B A}$ is quick to
calculate since you $M_B$ “sweeps out” the fixed effects and is a
byproduct of `fixest` estimates.

For IJIVE, another trick is used. Instead of actually residualizing
everything, it’s useful to use the original $Z$ for computational
effiicency (e.g. if $Z$ is high-dimensional FEs). If we actually
residualized the $Z$ fixed-effects, then we couldn’t use `fixest` to
estimate them. For this reason, we can use [block projection matrix
decomposition](https://en.wikipedia.org/wiki/Projection_matrix#Blockwise_formula)
to transform the hat matrix into two components that are easier to
calculate with `fixest` regressions:

$$
  H_{\tilde{Z}} = H_{M_W Z} = H_{[W Z]} - H_{W}.
$$

Then, through some algebra $H_{\tilde{Z}} M_W T$ can be written as
$H_{[W Z]} T - H_W T$, i.e. the fitted values of the unresidualized $T$
from a regression on $W$ and $Z$ and from a regression on $W$. Likewise,
$D_{\tilde{Z}}$ can be written as $D_{[W Z]} - D_W$.

## Sources

[Ackerberg and Deverux (2009)](literature/Ackerberg_Deverux_2009.pdf)

[Angrist, Imbens, and Krueger
(1999)](literature/Angrist_Imbens_Krueger_1999.pdf)

[Frandsen, Leslie, and McIntyre
(2023)](literature/Frandsen_Leslie_McIntyre_2023.pdf)

[Kolesar (2013)](literature/Kolesar_2013.pdf)

[kolesarm/ivreg](https://github.com/kolesarm/ivreg)

[ghpk-metrics/stata-manyiv](https://github.com/gphk-metrics/stata-manyiv)
