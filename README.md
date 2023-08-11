
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
  data = stevenson, y = ~ guilt, exogenous = ~ i(black) + i(white) | bailDate, 
  endogenous = ~ jail3, instruments = ~ 0 | judge_pre
)
#> $beta
#> [1] -0.02184514
#> 
#> $se
#> [1] -0.007517231
#> 
#> $F
#> [1] 32.62642
#> 
#> $Omega
#>              [,1]         [,2]
#> [1,]  0.243162932 -0.001992143
#> [2,] -0.001992143  0.233459364
#> 
#> $Xi
#>              [,1]         [,2]
#> [1,] 1.099707e-06 2.453454e-05
#> [2,] 2.453454e-05 1.556895e-04
#> 
#> $Sargan
#> $Sargan$statistic
#> [1] 3.342377
#> 
#> $Sargan$pvalue
#> [1] 0.7648087
#> 
#> 
#> $CD
#> $CD$statistic
#> [1] 3.31866
#> 
#> $CD$pvalue
#> [1] 0.7679163
#> 
#> 
#> attr(,"class")
#> [1] "UJIVE"    "jive_est"
```

``` r
ujive(
  data = stevenson, y = ~ guilt, exogenous = ~ i(black) + i(white) | bailDate, 
  endogenous = ~ jail3, instruments = ~ 0 | judge_pre
)
#> $beta
#> [1] 0.1590913
#> 
#> $se
#> [1] 0.07056371
#> 
#> $F
#> [1] 32.62642
#> 
#> $Omega
#>              [,1]         [,2]
#> [1,]  0.243162932 -0.001992143
#> [2,] -0.001992143  0.233459364
#> 
#> $Xi
#>              [,1]         [,2]
#> [1,] 1.099707e-06 2.453454e-05
#> [2,] 2.453454e-05 1.556895e-04
#> 
#> $Sargan
#> $Sargan$statistic
#> [1] 3.342377
#> 
#> $Sargan$pvalue
#> [1] 0.7648087
#> 
#> 
#> $CD
#> $CD$statistic
#> [1] 3.31866
#> 
#> $CD$pvalue
#> [1] 0.7679164
#> 
#> 
#> attr(,"class")
#> [1] "UJIVE"    "jive_est"
```

TBD. I’m not sure if this syntax is better. The advantage is the syntax
is a lot less verbose. I’m not sure if it could be confusing though.

``` r
# This doesn't run at the moment
jive(
  guilt ~ i(black) + i(white) | bailDate | jail3 ~ 0 | judge_pre,
  data = stevenson, 
)
```

### (Leave-out) Leniency Measures

he package will allow you to estimate (leave-out) leniency measures:

``` r
# TODO :-) 
# I have the function, just thinking of the API
```

## Econometric Details on JIVE, UJIVE, IJIVE, and CJIVE

Consider the following instrumental variables setup

$$
y_i = T_i \beta + W_i' \psi + \varepsilon_i,
$$ where $T_i$ is a scalar endogenous variable, $W_i$ is a vector of
exogenous covariates, and $\varepsilon_i$ is an error term.
Additionally, assume there is a set of valid instruments, $Z_i$. The
intuition of two-stage least squares is to use the instruments to
“predict” $T_i$: $$
T_i = Z_i' \pi + W_i \gamma + \eta_i.
$$ Then, the prediction, $\hat{T}_i$, is used in place of $T_i$ in the
original regression.

When the dimension of $Z_i$ grows with the number of observations,
two-stage least squares is biased (Kolesar, 2013). Without getting into
the details, the problem comes because you’re predicting $T_i$ using
$i$’s own observation in the first-stage. For this reason, Angrist,
Imbens, and Krueger (1999) developed the jackknife instrumental-variable
estimator (JIVE). In short, for each $i$, $T_i$ is predicted using the
first-stage equation leaving out $i$’s own observation.

In general, the JIVE estimator (and variants) are given by $$
\frac{\hat{P}' \tilde{Y}}{\hat{P}' \tilde{T}}
$$ where $\hat{P}$ is a function of $W$ and $\hat{T}_i$. The particulars
differ across the JIVE, the unbiased JIVE (UJIVE), the improved JIVE
(IJIVE), and the cluster JIVE (CJIVE).

### JIVE definition

**Source:** Kolesar (2013) and Angrist, Imbens, and Kreuger (1999)

The original JIVE estimate produces $\hat{T}$ given by $$
\hat{T}_{\text{JIVE}} = (I - D_{(Z,W)})^{-1} (H_{(Z,W)} - D_{(Z,W)}) T,
$$ where $H_{(Z,W)}$ is the hat/projection matrix for $(Z,W)$ and
$D_{(Z,W)}$ is the diagonal matrix with the diagonal elements from
$H_{Z, W}$.

$$
\begin{align*}
\hat{P}_{\text{JIVE}} &= M_W \hat{T}_{\text{JIVE}}
\end{align*}
$$

### UJIVE definition

**Source:** Kolesar (2013)

$$
\hat{T}_{\text{UJIVE}} = (I - D_{(Z,W)})^{-1} (H_{(Z,W)} - D_{(Z,W)}) T = \hat{T}_{\text{JIVE}}
$$ $$
\hat{P}_{\text{UJIVE}} = \hat{T}_{\text{UJIVE}} - (I - D_{W})^{-1} (H_{W} - D_{W}) T
$$

### IJIVE definition

**Source:** Ackerberg and Devereux (2009)

Start residualizing $T$, $Y$, and $Z$ by the covariates $W$. The
definition they give (assuming things have already been residualized) is
$$
\hat{T}_{IJIVE} = \hat{P}_{IJIVE} = (I - D_{\tilde{Z}})^{-1} (H_{\tilde{Z}} - D_{\tilde{Z}}) \tilde{T}
$$ Note that $P = T$ because you have already residualized by $W$.

### CJIVE definition

**Source:** Frandsen, Leslie, and McIntyre (2023)

This is a modified version of IJIVE as proposed by Frandsen, Leslie, and
McIntyre (2023). This is necessary if the errors are correlated within
clusters (e.g. court cases assigned on the same day to the same judge).
The modified version is given by:

$$
\hat{T}_{CJIVE} = \hat{P}_{IJIVE} = (I - \mathbb{D}(P_Z, \{ n_1, \dots, n_G \}))^{-1} (H_{\tilde{Z}} - D_{\tilde{Z}}) \tilde{T},
$$ where $\mathbb{D}(P_Z, \{ n_1, \dots, n_G \})$ is a block-diagonal
matrix equal to the projection matrix of $Z$ zerod out except for each
cluster.

In this package, the same adjustment (replacing the diagonal $D$ with
the cluster block-diagonal $\mathbb{D}$ version) can be done for all
three estimators, though the details have only been worked out for IJIVE
in particular (free research idea).

### Standard errors

Heteroskedastic-robust standard errors are given by $$ 
  \frac{\sqrt{\sum_i \hat{P}_i^2 \hat{\epsilon}_i^2}}{\sum_i \hat{P}_i T_i},
$$ where $\hat{\epsilon}_i = M_W Y_i - M_W T_i * \hat{\beta}$. See
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

## Sources

[Ackerberg and Deverux (2009)](literature/Ackerberg_Deverux_2009.pdf)

[Angrist, Imbens, and Krueger
(1999)](literature/Angrist_Imbens_Krueger_1999.pdf)

[Frandsen, Leslie, and McIntyre
(2023)](literature/Frandsen_Leslie_McIntyre_2023.pdf)

[Kolesar (2013)](literature/Kolesar_2013.pdf)

[kolesarm/ivreg](https://github.com/kolesarm/ivreg)

[ghpk-metrics/stata-manyiv](https://github.com/gphk-metrics/stata-manyiv)
