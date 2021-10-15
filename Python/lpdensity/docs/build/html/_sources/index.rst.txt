.. lpdensity documentation master file, created by
   sphinx-quickstart on Wed Oct 13 09:51:29 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

lpdensity Homepage 
=====================================
``lpdensity`` implements the local polynomial regression based density (and derivatives) estimator and bandwidth selection tools proposed in Cattaneo, Jansson and Ma (2020).  Robust bias-corrected inference methods,
both pointwise (confidence intervals) and uniform (confidence bands), are also implemented
following the results in Cattaneo, Jansson and Ma (2020, 2021a).
See Cattaneo, Jansson and Ma (2021b) for more implementation details and illustrations.

Install the ``lpdensity`` package by running
``pip install lpdensity``.

Import density estimation and bandwidth selection functions by running the following lines

>>> from lpdensity import lpdensity

>>> from lpdensity import lpbwdensity


Related ``Python``, ``R`` and ``Stata`` packages useful for nonparametric estimation and inference are
available at `https://nppackages.github.io/ <https://nppackages.github.io/>`_.

For source code and replication files, visit the `lpdensity repository <https://github.com/nppackages/lpdensity/>`_.

References
----------
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. 
`On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference <https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf>`_
*Journal of the American Statistical Association*, 113(522): 767-779.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020. 
`Coverage Error Optimal Confidence Intervals for Local Polynomial Regression <https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2020_CEopt.pdf>`_
. Working paper.

Cattaneo, M. D., M. Jansson, and X. Ma. 2020.
`Simple Local Polynomial Density Estimators <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf>`_.
*Journal of the American Statistical Association*, 115(531): 1449-1455.

Cattaneo, M. D., M. Jansson, and X. Ma. 2021a.
`Local Regression Distribution Estimators <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2021_JoE.pdf>`_
*Journal of Econometrics*, forthcoming.

Cattaneo, M. D., M. Jansson, and X. Ma. 2021b.
`lpdensity: Local Polynomial Density Estimation and Inference <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2021_JSS.pdf>`_
*Journal of Statistical Software*, forthcoming.

Authors
-------
Matias D. Cattaneo, Princeton University. (cattaneo@princeton.edu).

Rajita Chandak (maintainer), Princeton University. (rchandak@princeton.edu). 

Michael Jansson, University of California Berkeley. (mjansson@econ.berkeley.edu).

Xinwei Ma (maintainer), University of California San Diego. (x1ma@ucsd.edu).


.. toctree::
   :maxdepth: 2
   :caption: Contents:
  
   modules


