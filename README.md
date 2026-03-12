# 2026 Growth Kinetics
Code for the analysis and plots of the growth curves and growth kinetics params obtention.

# Input and Output data
Data is OD measurements of differnet *S. pneumoniae* strains across different conditions.

**Input** : excel file with optical density (OD) measurements across time (h) for multiple samples and replicates.
**Output** : `Growth_features_estimates.tsv`, estimated parameters for each replicate. `Growth_features_estimates_failed_samples.tsv`, list of replicates that failed.

Values with negaitve OD were treated as 0.

# Growth Kinetics analysis

The script from [Chaguza](https://github.com/ChrispinChaguza/SpnGrowthKinetics/blob/main/growth_curves.R) has been adapted.
The pipeline processes microbial growth curve data, fits logistic growth models, and extracts growth parameters for each sample replicate.

## Model

4-parameters logistic growth functions:

`OD(t) = L0 + L / (1 + exp(-r * (t - t0)))`

Where:

```
L0 : initial OD
L : max OD increase
r : growth rate
t0 : inflection point 
```

For the fitting we masked to focus on main growth and plateau phases, and use multiple optimization methods: `BFGS`, `Nelder-Mead`, `CG`, `Powell`, `L-BFGS-B`, finally refined the best fit using `least_quares`.

Feature extraction includes: `delta.H`, max OD increae and `lag`, lag phase (OD first exceeds 1.5 * initial OD).

Modifications from the original script:
* **lag calculation**: The Chaguza script was using the observed OD data to extract lag, but because we have values close to 0 or 0 at the begining this leads to lag == 0. We instead used the fitted logisitic curve params and ocmpute threshold as 5% above the baseline avoiding values too close to zero.
* **allow maxOD to 1.5**



