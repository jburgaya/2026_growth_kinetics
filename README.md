# 2026 Growth Kinetics

Code for the analysis and visualisation of growth curves and extraction of growth kinetics parameters.

---

## Input and output

**Input:** Excel file with optical density (OD) measurements across time (hours) for multiple samples and replicates. Negative OD values are treated as 0.

**Output:**

| File | Description |
|---|---|
| `Growth_features_estimates_{batch}.tsv` | Fitted parameters + 95% CI for every replicate |
| `Growth_features_estimates_failed_{batch}.tsv` | Replicates that could not be fitted (OD < 0.01 or optimisation failure) |
| `figures/<strain>_growth_curves.pdf/.png` | Fitted curves with CI bands per strain |

---

## Scripts

| Script | Purpose |
|---|---|
| `1_growth_curve_fitting-{batch}` | Fit models, extract parameters and CIs |
| `2_growth_curves_plot-{batch}` | Plot fitted curves with shaded CI bands |

---

## Growth kinetics pipeline

The pipeline is adapted from [Chaguza et al.](https://github.com/ChrispinChaguza/SpnGrowthKinetics/blob/main/growth_curves.R), extended with a multi-model fitting framework, AICc-based model selection, and Jacobian-derived confidence intervals.

### Step 1 — preprocessing

For each replicate, data are masked to retain only the rising phase (up to and including the observed OD maximum) plus a 5% plateau tail. This focuses the optimiser on the biologically meaningful growth window and reduces the influence of post-stationary decline.

### Step 2 — model fitting

The pipeline always attempts the **4-parameter logistic first**, using three biologically-informed starting points refined by Trust Region Reflective (TRF) nonlinear least-squares. The logistic fit is accepted if it passes six sanity checks:

- growth rate `r > 0.05`
- no runaway extrapolation: `OD(t=40h) < 3.0`
- fit explains > 50% of total variance
- inflection time `t₀` falls within ±5 h of the observation window
- amplitude `L > 0.005`
- baseline `L₀ ≥ 0`

If the smoothed first derivative of the growth curve shows two distinct peaks, a **biphasic model** is also fitted and compared to the logistic by AICc. The biphasic model is selected only if it achieves a lower AICc — meaning the improvement in fit justifies its three additional parameters.

If the logistic fails the sanity checks, a **fallback cascade** is run: Gompertz → Richards → Exponential saturation, all scored by AICc. The model with the lowest AICc is selected.

### Step 3 — confidence intervals

95% confidence intervals on all parameters are derived analytically from the Jacobian covariance matrix at the TRF solution:

```
C = σ² (JᵀJ)⁻¹
```

where σ² is the residual variance and J is the Jacobian returned by `scipy.optimize.least_squares` with `jac="3-point"`. This requires no resampling and adds negligible runtime.

---

## Models

### 4-parameter logistic (primary model)

```
OD(t) = L₀ + L / (1 + exp(−r · (t − t₀)))
```

Symmetric S-shaped curve. Used for the large majority of well-behaved growth curves.

### Gompertz (fallback)

```
OD(t) = L₀ + A · exp(−exp(−r · (t − t₀)))
```

Asymmetric sigmoid with a longer lag phase and sharper exponential rise.

### Richards (fallback)

```
OD(t) = L₀ + A / (1 + v · exp(−r · (t − t₀)))^(1/v)
```

Generalised logistic with a shape parameter `v` that moves the inflection freely. Recovers the logistic when `v = 1` and approaches Gompertz as `v → 0`.

### Exponential saturation (fallback)

```
OD(t) = L₀ + A · (1 − exp(−r · t))
```

No inflection point. Applied when growth never plateaus within the 22-hour observation window.

### Biphasic — double logistic (AICc-selected)

```
OD(t) = L₀ + A₁ / (1 + exp(−r₁ · (t − t₁))) + A₂ / (1 + exp(−r₂ · (t − t₂)))
```

Two sequential logistic phases sharing a common baseline. Models diauxic growth — two distinct substrate utilisations, or a secondary growth phase following partial inhibition.

---

## Output parameters

All parameters are extracted into a unified output regardless of which model was selected.

| Column | Symbol | Description |
|---|---|---|
| `L` | *L* or *A* | Total OD gain (amplitude); difference between plateau and baseline |
| `r` | *r* | Maximum growth rate parameter (h⁻¹); related to µ_max by model-specific formula |
| `t0` | *t₀* | Inflection time (h); time at which growth rate is highest |
| `L0` | *L₀* | Baseline OD before growth begins |
| `delta_H` | *ΔH* | Maximum OD reached: *L₀ + L* |
| `lag` | *λ* | Lag time (h): first time fitted curve exceeds *L₀ + 5% × L* |
| `L_lo / L_hi` | | 95% CI on amplitude |
| `r_lo / r_hi` | | 95% CI on growth rate |
| `t0_lo / t0_hi` | | 95% CI on inflection time |
| `L0_lo / L0_hi` | | 95% CI on baseline |
| `model` | | Winning model: `logistic4`, `gompertz4`, `richards`, `exponential`, `biphasic` |
| `AICc` | | Small-sample corrected Akaike Information Criterion of the winning model |
| `converged` | | 1 = converged, 0 = fallback result |
| `raw_params` | | JSON with exact optimised parameter values (used for curve plotting) |

### Lag time definition

The original Chaguza script used `1.5 × OD_min` as the lag threshold. When initial OD values are near zero, this collapses to lag ≈ 0. We instead define lag as:

```
threshold = L₀ + 0.05 × L
```

evaluated on a 5000-point dense grid of the fitted curve. This is model-agnostic, robust to near-zero baselines, and consistent across all five model families.

---

## Model selection rationale

AICc was chosen over AIC because sample sizes per replicate are small relative to the number of parameters in the biphasic model (7), making the finite-sample correction non-negligible. AICc is defined as:

```
AICc = −2·log(L̂) + 2k + 2k(k+1)/(n−k−1)
```

where *k* is the number of free parameters plus one (for σ²) and *n* is the number of fitted data points. A lower AICc indicates a better balance of fit and parsimony. The biphasic model must achieve a substantially lower AICc (typically ΔAIC < −6) to be preferred over the 4-parameter logistic, preventing overfitting of noisy curves.

---

## Dependencies

```
python >= 3.9
numpy
pandas
scipy
matplotlib
seaborn
openpyxl
```

Install with:

```bash
pip install numpy pandas scipy matplotlib seaborn openpyxl
```

---

## References

- Chaguza C. *SpnGrowthKinetics*. GitHub. https://github.com/ChrispinChaguza/SpnGrowthKinetics
- Burnham KP, Anderson DR. *Model Selection and Multimodel Inference*. 2nd ed. Springer; 2002.
- Richards FJ. A flexible growth function for empirical use. *J Exp Bot*. 1959;10(2):290–301.
