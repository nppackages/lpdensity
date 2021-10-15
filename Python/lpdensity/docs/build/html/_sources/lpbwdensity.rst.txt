lpbwdensity
===========
Data-driven Bandwidth Selection for Local Polynomial Density Estimators

Description
-----------
``lpbwdensity`` implements the bandwidth selection methods for local
polynomial based density (and derivatives) estimation proposed and studied
in Cattaneo, Jansson and Ma (2020, 2021a).
See Cattaneo, Jansson and Ma (2021b) for more implementation details and illustrations.

Companion command: ``lpdensity`` for estimation and robust bias-corrected inference.


References
----------
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


.. automodule:: lpdensity.lpbwdensity
   :members:
   :undoc-members:

Example
-------
>>> import numpy as np
>>> from lpdensity import lpbwdensity
>>> data = np.random.normal(0,1,500)
>>> grid = np.linspace(min(data), max(data), 10)
>>> est = lpbwdensity(data=data, grid=grid, bwselect="mse-dpi")
>>> print(repr(est))
