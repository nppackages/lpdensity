# LPDENSITY

The `lpdensity` package provides tools for local polynomial distribution and density estimation, bandwidth selection, robust bias-corrected inference, and graphical procedures.

- `lpdensity`: implements local polynomial distribution and density estimation with pointwise and uniform inference.
- `lpbwdensity`: implements data-driven bandwidth selectors for local polynomial distribution and density estimation.

See Cattaneo, Jansson and Ma (2020, 2022, 2024) for references.

Website: [https://nppackages.github.io/](https://nppackages.github.io/).

Source code: [https://github.com/nppackages/lpdensity](https://github.com/nppackages/lpdensity).

## Authors

Matias D. Cattaneo, maintainer (<matias.d.cattaneo@gmail.com>)

Rajita Chandak (<rajita.chandak@gmail.com>)

Michael Jansson (<michael.jansson.berkeley@gmail.com>)

Xinwei Ma (<xinweima.pku@gmail.com>)

## Installation

To install/update use pip
```
pip install lpdensity
```

# Usage
```
from lpdensity import lpdensity, lpbwdensity
```

- Replication: [lpdensity illustration](https://github.com/nppackages/lpdensity/blob/main/Python/lpdensity_illustration.py).

## Dependencies

- numpy
- pandas
- scipy
- sympy
- plotnine

## References

For overviews and introductions, see [NP Packages website](https://nppackages.github.io/).

### Software and Implementation

- Cattaneo, Jansson and Ma (2022): [lpdensity: Local Polynomial Density Estimation and Inference](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf).<br>
_Journal of Statistical Software_ 101(2): 1-25.

### Technical and Methodological

- Cattaneo, Jansson and Ma (2020): [Simple Local Polynomial Density Estimators](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf).<br>
_Journal of the American Statistical Association_ 115(531): 1449-1455.<br>
[Supplemental appendix](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA--Supplement.pdf).

- Cattaneo, Jansson and Ma (2024): [Local Regression Distribution Estimators](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2024_JoE.pdf).<br>
_Journal of Econometrics_ 240(2): 105074.<br>
[Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2024_JoE--Supplement.pdf).
