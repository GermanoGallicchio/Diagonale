# _PhysioExplorer_ for MATLAB (soon to be rebranded and changed substantially to include new features)

Author: Germano Gallicchio, [Bangor University](https://www.bangor.ac.uk/)

Code developed on MATLAB R2025b on a Linux OS (Kubuntu).



## Overview

Not yet available.




## Documentation

Documentation for PhysioExplorer is available at this (link)[https://germanogallicchio.github.io/PhysioExplorer_documentation/index.html]



## Tutorials

Not yet available.



## Installation

Not yet available.




## How to contribute

Not yet available.



## How to cite

Gallicchio, G. (2025). PhysioExplorer. Zenodo. [https://doi.org/10.5281/zenodo.16808782](https://doi.org/10.5281/zenodo.16808782)

---

## Overview
**PhysioExplorer (PE)** is a set of functions to extract patterns from multivariate physiological data. 

Testing effects on datasets with many and correlated variables (i.e., multivariate)--or even much larger than the number of observations (i.e., megavariate, Eriksson et al., 2013)--can be challenging due to the multiple comparison problem. (See the [green jelly bean comic](https://xkcd.com/882).) This challenge can be overcome through various approaches.
1. One solution is to run mass (i.e., a lot of) univariate tests and then correct for False Discovery Rate (e.g., Benjamini & Hochberg, 1995).
2. Another solution is to still run mass univariate tests and then cluster their results where there is contiguity in some physical dimension (e.g., time, frequency, sensor space). Then compute cluster metrics (e.g., their statistical mass) and perform inference on them (Groppe et al., 2011; Maris & Oostenveld, 2007).
3. Another solution is to find the combination of the whole set of features that best describes behavioral data or experimental design group/condition (Add reference).

If hypothesis testing is not the goal, but rather stability of the statistical metric across sampling variability, the bootstrap framework provides such metrics.

## What PE's current version can do
[X] = PE can do it
<br>
[ ] = PE cannot _yet_ do it
<br>
| analysis &<br>objective | symmetric association between variables | compare groups | compare levels of one repeated-measure factor
| ---: | :---: | :---: | :---: | 
| empiricalL1_FDR<br>permutation             | [X]<br>(2 variables)     | [X]<br>(2 groups)   | [X]<br>(2 levels) |
| empiricalL1_FDR<br>bootstrap               | [X]<br>(2 variables)     | [X]<br>(2 groups)   | [X]<br>(2 levels) |
| theoreticalL1_clusterMaxT<br>permutation   | [X]<br>(2 variables)     | [X]<br>(2 groups)   | [X]<br>(2 levels) |
| theoreticalL1_clusterMaxT<br>bootstrap     | [ ]<br>(2 variables)     | [ ]<br>(2 groups)   | [ ]<br>(2 levels) |
| PLS_SVD<br>permutation                     | [ ]<br>(2 variable sets) | [ ]<br>(2+ groups)  | [ ]<br>(2+ levels) |
| PLS_SVD<br>bootstrap                       | [ ]<br>(2 variable sets) | [ ]<br>(2+ groups)  | [ ]<br>(2+ levels) |

---


## Referencing: <br>
If you use PhysioExplorer, please cite the Zenodo DOI

Example:

Gallicchio, G. (2025). PhysioExplorer. Zenodo. [https://doi.org/10.5281/zenodo.16808782](https://doi.org/10.5281/zenodo.16808782)

---

## Documentation
_(this documentation is an unpolished draft)_
PhysioExplorer can perform any combination of _analysis_ and _objective_ described below in both multivariate and megavariate contexts (with no distinction). 



```mermaid
  graph LR;
    A(pe_cfg.analysis)
    B(empiricalL1_FDR)
    C(theoreticalL1_clusterMaxT)
    D(PLS_SVD)
    E(pe_cfg.objective)
    F(permutationH0testing)
    G(bootstrapStability)

    A-->B;
    A-->C;
    A-->D;
    E-->F;
    E-->G;

```

## Analysis (pe_cfg.analysis)
### 'empiricalL1_FDR'
...upcoming description...

### 'theoreticalL1_clusterMaxT'
**Cluster-level analysis** (Groppe et al., 2011; Maris & Oostenveld, 2007) is a two-step procedure: (1) compute univariate test statistics, evaluate them against the associated theoretical distribution to get p-values, and threshold them, (2) form spatial/temporal/spectral clusters of suprathreshold points. Clusters can be defined in a 3-dimensional space (e.g, time-frequency-channel, frequency-frequency-channel) or a lower-dimensional subset (e.g., time-channel, time-frequency, frequency-channel, time). At the heart of the code is a_cluster forming algorithm that combines adjacency criteria (e.g., spatial-temporal-spectral) with the results of univariate statistical testing (e.g., p-values). The code forms clusters on the observed data and, depending on the _objective_ many sets of surrogate data artificially created under the null hypothesis of exchangeability of group/condition labels (permutataion) or many replicates, each with sampling variability, of the original data (bootstrap). The surrogate data are sampled through the Monte-Carlo approach. 

### 'PLD_SVD'
**SVD-based Partial Least Squares** is a form of symmetric covariance mapping (Note: SVD stands for singular value decomposition.) It handles multi/megavariate data structures natively (in one step) to find combinations of features that best describe the linear associations between two sets of variables. The number of combinations found depends on how much linear independence is in the combined data (the rank). Each combination is characterized by the singular value, informing on how much this combination explains of the covariance, and two singular vectors (one for each variable set), telling how the original variables should be weighted to form that specific combination. Resampling statistics are then used to evaluate whether a certain mapping has a magnitude larger than noise (permutation testing on the singular value based metrics) or whether a certain combination's weights are stable under sampling variability (bootstrap evaluation on the weights).

## Objective (pe_cfg.objective)
### 'permutation' (pe_cfg.objective = permutationH0testing)
**Permutation** is for null-hypothesis testing. In each Monte-Carlo iteration, group/condition labels are shuffled, and the statistics are recomputed. The code compares the observed cluster metrics (e.g., cluster mass, singular value) with the distribution of the same metrics under the null hypothesis to evaluate their statistical significance. (Note: for cluster analysis, inference is done at the cluster level and not at the point level.)
### 'bootstrap' (pe_cfg.objective = bootstrapStability)
**Bootstrap** is for stability estimation.





## Wish list (maybe future updates)
- cluster descriptives: effect size for group/condition comparison 
- visualize 3d results: pe_x1y2z3View
- complete code PLS_SVD
- improve own version of topoplot to allow spherical interpolation
- write tutorials on how to use PhysioExplorer



## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal statistical society: series B (Methodological), 57(1), 289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Eriksson, L., Byrne, T., Johansson, E., Trygg, J., & Vikström, C. (2013). Multi-and megavariate data analysis basic principles and applications. Umetrics Academy.

Groppe, D. M., Urbach, T. P., & Kutas, M. (2011). Mass univariate analysis of event‐related brain potentials/fields I: A critical tutorial review. Psychophysiology, 48(12), 1711-1725.

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of neuroscience methods, 164(1), 177-190.


