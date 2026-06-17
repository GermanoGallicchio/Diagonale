# _Diagonale_ for MATLAB

Author: Germano Gallicchio, [Bangor University](https://www.bangor.ac.uk/)

Code developed on MATLAB R2025b on a Linux OS (Kubuntu).

Note: Diagonale is an expansive project. It will expand to improve and develop analytic features. 

## How to cite / Referencing
If you use Diagonale, please cite this software and mention the specific version you have used.

Gallicchio, G. (2025). Diagonale. Zenodo. https://doi.org/10.5281/zenodo.16808782

---

## Overview
**Diagonale** is a set of functions to extract patterns from multivariate physiological data.

Testing effects on datasets with many and correlated variables (i.e., multivariate), or even much larger than the number of observations (i.e., megavariate, Eriksson et al., 2013), can be challenging due to the multiple comparison problem. (See the [green jelly bean comic](https://xkcd.com/882).) This challenge can be overcome through various approaches.
1. One solution is to run mass (i.e., a lot of) univariate tests (non-parametric via permutation or parametric) and then correct for False Discovery Rate (e.g., Benjamini and Hochberg, 1995) or most extreme statistic (maxT algorithm).
2. Another solution is to still run mass univariate tests (parametric for convenience) and then cluster significant results that show contiguity in some physical dimension (e.g., time, frequency, sensor space). Then compute cluster metrics (e.g., cluster mass) and perform inference on them (Groppe et al., 2011; Maris and Oostenveld, 2007) and correct for multiple comparisonsvia extreme statistics (maxT algorithm).
3. Another solution is to find the linear combination of the whole set of features that best describes behavioral data or experimental design group/condition.

If hypothesis testing is not the goal, but rather stability of the statistical metric across sampling variability, the bootstrap framework provides such metrics.

Diagonale offers a command line interface (i.e., no graphic user interface) to perform these analyses. While scripting is required, the script would consist of preparing the data in a standardized format and making choices about the type of analysis and its parameters, so no actual programming is required.

## What Diagonale's current version can do

| objective         | inferenceLevel | supported | 
| --- | --- | --- |
| `permutationH0testing` | `feature` | yes | 
| `permutationH0testing` | `cluster` | yes | 
| `permutationH0testing` | `latent`  | yes | 
| `bootstrapStability`   | `feature` | yes |
| `bootstrapStability`   | `cluster` | no  |
| `bootstrapStability`   | `latent`  | yes |

## Documentation
Documentation for Diagonale is in preparation and will be available at this [link](https://germanogallicchio.github.io/Diagonale_documentation/index.html)

## Installation
An how to install guide will be hosted in the Documentation page. But in short, for now, installation is manual and simple:
1. Download (or clone) this repository.
2. Add the Diagonale folder to the MATLAB path.

## Tutorials
Tutorials are in preparation and will be hosted in the Documentation page. 

---

## Wish list (maybe future updates)
- cross validation for generalizability evaluation (high priority, highest effort). 
- cluster descriptives: effect sizes at descriptive level (low priority, very low effort)
- enable 3d di_view functions (e.g., channels, time, frequency) using imagesc() instead of plot() waveforms (medium priority, low effort)
- write tutorials on how to use Diagonale (medium priority, low effort)
- improve own version of topoplot to allow spherical interpolation (very low priority, high effort)




## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B (Methodological), 57(1), 289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Eriksson, L., Byrne, T., Johansson, E., Trygg, J., & Vikström, C. (2013). Multi- and megavariate data analysis: basic principles and applications. Umetrics Academy.

Groppe, D. M., Urbach, T. P., & Kutas, M. (2011). Mass univariate analysis of event‐related brain potentials/fields I: A critical tutorial review. Psychophysiology, 48(12), 1711-1725.

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.


