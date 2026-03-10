# 2026 Growht Kinetics
Code for the analysis and plots of the growht curves and growht kinetics params obtention.

# Input and Output data
Data is OD measurements of differnet *S. pneumoniae* strains across different conditions.

**Input** : excel file with optical density (OD) measurements across time (h) for multiple samples and replicates.
**Output** : `Growth_features_estimates.tsv`, estimated parameters for each replicate. `Growth_features_estimates_failed_samples.tsv`, list of replicates that failed.

# Growth Kinetics analysis

The script from (Chaguza)[https://github.com/ChrispinChaguza/SpnGrowthKinetics/blob/main/growth_curves.R] has been adapted.
The pipeline processes microbial growth curve data, fits logistic growhth models, and extracts groeth parameters for each sample replicate.

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



