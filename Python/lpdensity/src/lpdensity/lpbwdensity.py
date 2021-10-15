#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import math
from . import lpbwdensity_fn

pd.options.mode.chained_assignment = None

def lpbwdensity(data, grid=None, p=None, v=None, kernel="triangular", bwselect="mse-dpi", massPoints=True, stdVar=True, regularize=True, nLocalMin=None, nUniqueMin=None, Cweights=None, Pweights=None):
	"""
	Parameters
	----------
		 data : vector
			Numeric vector or one dimensional matrix/data frame, the raw data.
		 grid : vector
			Numeric, specifies the grid of evaluation points. When set to default, grid points will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05.
		 p : int
			Nonnegative integer, specifies the order of the local polynomial used to construct point estimates. (Default is *2*.)
		 v : int
			Nonnegative integer, specifies the derivative of the distribution function to be estimated. *0* for the distribution function, *1* (default) for the density funtion, etc.
		 kernel : string
			Specifies the kernel function, should be one of *"triangular"*, *"uniform"* or
		   *"epanechnikov"*.
		 bwselect : string
			Specifies the method for data-driven bandwidth selection. This option will be
			ignored if *bw* is provided. Can be (1) *"mse-dpi"* (default, mean squared error-optimal
			bandwidth selected for each grid point); or (2) *"imse-dpi"* (integrated MSE-optimal bandwidth,
			common for all grid points); (3) *"mse-rot"* (rule-of-thumb bandwidth with Gaussian
			reference model); and (4) *"imse-rot"* (integrated rule-of-thumb bandwidth with Gaussian
			reference model).
		 massPoints : boolean
			*True* (default) or *False*, specifies whether point estimates and standard errors
		   should be adjusted if there are mass points in the data.
		 stdVar : boolean
			*True* (default) or *False*, specifies whether the data should be standardized for
		   bandwidth selection.
		 regularize : boolean
			*True* (default) or *False*, specifies whether the bandwidth should be
		   regularized. When set to *True*, the bandwidth is chosen such that the local region includes
		   at least *nLocalMin* observations and at least *nUniqueMin* unique observations.
		 nLocalMin : int
			Nonnegative integer, specifies the minimum number of observations in each local neighborhood. This option
		   will be ignored if *regularize=False*. Default is *20+p+1*.
		 nUniqueMin : int
			Nonnegative integer, specifies the minimum number of unique observations in each local neighborhood. This option
		   will be ignored if *regularize=False*. Default is *20+p+1*.
		 Cweights : vector
			Numeric vector, specifies the weights used
		   for counterfactual distribution construction. Should have the same length as the data.
		   This option will be ignored if *bwselect* is *"mse-rot"* or *"imse-rot"*.
		 Pweights : vector
			Numeric vector, specifies the weights used
		   in sampling. Should have the same length as the data.
		   This option will be ignored if *bwselect* is *"mse-rot"* or *"imse-rot"*.
		   
	Returns
	-------
	BW : object
		A :py:meth:`~bw_output` class object containing a matrix of (1) *grid* (grid point), (2) *bw* (bandwidth),
		(3) *nh* (number of observations in each local neighborhood),
		(4) *nhu* (number of unique observations in each local neighborhood), and
		(5) *opt* additional parameters.
	"""

	# data
	data = pd.DataFrame(data)
	if data.isnull().values.any():
		raise Exception(data.isnull().sum() + 'missing observation(s) are ignored.\n')
		data = data.dropna(axis=0)
	
	n = len(data)
	data = np.array(data).reshape(n,)
	if n==0 or not np.isscalar(n):
		raise Exception('Data should be numeric, and cannot be empty.')
	
	# grid
	if grid is None:
		grid = np.quantile(data, np.arange(0.05, 0.95, 0.05))
		ng = len(grid)
	else:
		ng = len(grid)
		if not np.isreal(grid).any():
			raise Exception('Grid points should be numeric.')
	
	# bwselect
	if len(bwselect) == 0:
		bwselect = "mse-dpi"
	else:
		bwselect = bwselect.lower()
	
	# p
	if p is None:
		p = 2
	elif not np.isscalar(p) or p not in range(21):
		raise Exception('Polynomial order p incorrectly specified.')
		
	# v
	if v is None:
		v = min(1, p)
	elif not np.isscalar(v) or v not in range(21) or v>p:
		raise Exception('Derivative order v incorrectly specified')
		
	#kernel
	if kernel is None:
		flag_no_kernel = True
		kernel = "triangular"
	else:
		kernel = kernel.lower()
		if kernel not in ["triangular", "uniform", "epanechnikov"]:
			raise Exception('Kernel function incorrectly specified.')
	
	# Cweights
	if Cweights is None:
		Cweights = np.repeat(1,n)
	elif not np.isreal(Cweights).any():
		raise Exception('Counterfactual weights incorrectly specified.')
	elif len(Cweights) !=n:
		raise Exception('Counterfactual weights should have the same length as sample.')
	
	if isinstance(Cweights, pd.Series):
		Cweights = Cweights.to_numpy()
	
	# Pweights
	if Pweights is None:
		Pweights = np.repeat(1,n)
	elif not np.isreal(Pweights).any():
		raise Exception('Probability weights incorrectly specified.')
	elif len(Pweights) !=n:
		raise Exception('Probability weights should have the same length as sample.')
	elif any(x<0 for x in Pweights)==True:
		raise Exception('Probability weights should be nonnegative.')
		
	if isinstance(Pweights, pd.Series):
		Pweights = Pweights.to_numpy()
	
	# massPoints
	if massPoints is None:
		massPoints = True
	elif type(massPoints) != bool:
		raise Exception('Option massPoints incorrectly specified.')
	
	
	# stdVar
	if stdVar is None:
		stdVar = True
	elif type(stdVar) != bool:
		raise Exception('Option stdVar incorrectly specified.')
	
	
	# regularize
	if regularize is None:
		regularize = True
	elif type(regularize) != bool:
		raise Exception('Regularization parameter incorrectly specified')
	
	
	# nLocalMin
	if nLocalMin is None:
		nLocalMin = 20 + p + 1
	if not np.isreal(nLocalMin) or np.isnan(nLocalMin):
		raise Exception('Option nLocalMin incorrectly specified.')
	elif math.ceil(nLocalMin)<0:
		raise Exception('Option nLocalMin incorrectly specified.')
	else:
		nLocalMin = math.ceil(nLocalMin)

	# nUniqueMin
	if nUniqueMin is None:
		nUniqueMin = 20 + p + 1
	if not np.isreal(nUniqueMin) or np.isnan(nUniqueMin):
		raise Exception('Option nUniqueMin incorrectly specified.')
	elif math.ceil(nUniqueMin)<0:
		raise Exception('Option nUniqueMin incorrectly specified.')
	else:
		nUniqueMin = math.ceil(nUniqueMin)
		
	#########################
	##   Sample Trimming   ##
	#########################
#	trim_index = (Pweights==0)
#	if np.all(trim_index==True):
#		raise Exception('All weights are zero.')
#	else:
#		data = data[Pweights!=0]
#		Cweights = Cweights[Pweights!=0]
#		Pweights = Pweights[Pweights!=0]
	
	if abs(sum(Cweights*Pweights)) <= 10*np.finfo(float).eps:
		raise Exception('Composited weights (Cweights * Pweights) are numerically zero.')
		
	#############################
	##   Bandwidth Selection   ##
	#############################
	
	if bwselect == "mse-dpi":
		bw = lpbwdensity_fn.__bw_MSE( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
	elif bwselect == "imse-dpi":
		bw =  lpbwdensity_fn.__bw_IMSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
		bw = np.repeat(bw, len(grid))
	elif bwselect == "mse-rot":
		bw = lpbwdensity_fn.__bw_ROT( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
	else:
		bw = lpbwdensity_fn.__bw_IROT(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
		bw = np.repeat(bw, len(grid))
		
	BW = pd.DataFrame(columns = ["grid", "bw", "nh", "nhu"])
	BW['grid'] = grid
	BW['bw'] = bw
	
	dataUnique =  lpbwdensity_fn.__lpdensityUnique(np.sort(data))['unique']
	
	for i in range(ng):
		BW['nh'][i] = sum(abs(data-BW['grid'][i])<= BW['bw'][i])
		BW['nhu'][i] = sum(abs(dataUnique-BW['grid'][i]) <= BW['bw'][i])
	
	return(bw_output(BW, bw, p, v, kernel, n, ng, bwselect, massPoints, stdVar, regularize, nLocalMin, nUniqueMin, data_min=min(data), data_max=max(data), grid_min=min(grid), grid_max=max(grid)))


class bw_output:
	"""
	Class of lpbwdensity function outputs.
	
	Object type returned by :py:meth:`~lpbwdensity`.
	"""
	def __init__(self, BW, bws, p, v, kernel, n, ng, bwselect, massPoints, stdVar, regularize, nLocalMin, nUniqueMin, data_min, data_max, grid_min, grid_max):
		self.table = BW
		self.bws = bws
		self.p = p
		self.v = v
		self.kernel = kernel
		self.n = n
		self.ng = ng
		self.bwselect = bwselect
		self.massPoints = massPoints
		self.stdVar = stdVar
		self.regularize = regularize
		self.nLocalMin = nLocalMin
		self.nUniqueMin = nUniqueMin
		self.data_min = data_min
		self.data_max = data_max
		self.grid_min = grid_min
		self.grid_max = grid_max
	
	def __repr__(self):
		print('Call: lpbwdensity')
		print('')
		fw = 30
		fw_r = 14
		print('Sample Size:'.ljust(fw), str(self.n).rjust(25))
		print('Polynomial order for point estimation (p=):'.ljust(fw), str(self.p).rjust(12))
		print('Order of derivative estimated         (v=):'.ljust(fw), str(self.v).rjust(12))
		print('Kernel function:'.ljust(fw), str(self.kernel).rjust(25))
		print('Bandwidth Method:'.ljust(fw), str(self.bwselect).rjust(25))
		print('')
		
		#Table
		fw_l = 15
		fw_c = 8
		fw_ci = 18
		n_dec = 3
		print('='*65)
		print('Index'.ljust(fw_c), 'Grid'.rjust(fw_c),
			'Bandwidths'.rjust(fw_l),
			'Eff.n'.rjust(fw_c), 'Uniq.n'.rjust(fw_c))
		print('-'*65)
		for j in range(self.ng):
			print(str(j+1).ljust(fw_c), str(round(self.table['grid'][j], n_dec)).rjust(fw_c),
			str(round(self.table['bw'][j], n_dec)).rjust(fw_l), str(int(self.table['nh'][j])).rjust(fw_c), str(int(self.table['nhu'][j])).rjust(fw_c))
		print('='*65)
		return('')
		
	def __str__(self):
		print('Call: lpbwdensity')
		print('')
		fw = 30
		fw_r = 14
		print('Sample Size:'.ljust(fw), str(self.n).rjust(25))
		print('Polynomial order for point estimation (p=):'.ljust(fw), str(self.p).rjust(12))
		print('Order of derivative estimated         (v=):'.ljust(fw), str(self.v).rjust(12))
		print('Kernel function:'.ljust(fw), str(self.kernel).rjust(25))
		print('Bandwidth Method:'.ljust(fw), str(self.bwselect).rjust(25))
		print('')
		print('Use .repr() to show bandwidths.')
		return ('')
	
	def coef(self):
		"""
		Returns estimate of bandwidths.
		"""
		coef_data = self.table[['grid', 'bw']]
		return(coef_data)
