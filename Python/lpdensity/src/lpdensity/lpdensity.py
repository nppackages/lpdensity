#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import warnings
import math
import pandas as pd
import scipy.stats as spstat
from . import lpdensity_fn
from . import lpbwdensity_fn
import plotnine as pn

pd.options.mode.chained_assignment = None

warnings.filterwarnings("ignore")


def lpdensity(data, grid=None, bw=None, p=None, q=None, v=None,
	kernel='triangular', scale=None,
	massPoints=True, bwselect='mse-dpi',
	stdVar=True, regularize=True,nLocalMin=None, nUniqueMin=None,
	Cweights=None, Pweights=None):
	
	"""
	Parameters
	----------
	data : vector
		Numeric vector or one dimensional matrix/data frame, the raw data.
	grid : number or vector
		Numeric, specifies the grid of evaluation points. When set to default, grid points will be chosen as 0.05-0.95 percentiles of the data, with a step size of 0.05.
	bw : number or vector
		Numeric, specifies the bandwidth used for estimation. Can be (1) a positive scalar (common bandwidth for all grid points); or (2) a positive numeric vector specifying bandwidths for each grid point (should be the same length as *grid*).
	p : int
		Nonnegative integer, specifies the order of the local polynomial used to construct point estimates. (Default is *2*.)
	q : int
		Nonnegative integer, specifies the order of the local polynomial used to construct confidence intervals/bands (a.k.a. the bias correction order). Default is *p+1*. When set to be the same as *p*, no bias correction will be performed. Otherwise it should be strictly larger than *p*.
	v : int
		Nonnegative integer, specifies the derivative of the distribution function to be estimated. *0* for the distribution function, *1* (default) for the density funtion, etc.
	kernel : string
		Specifies the kernel function, should be one of *"triangular"*, *"uniform"*, and *"epanechnikov"*.
	scale : number
		Numeric, specifies how estimates are scaled. For example, setting this parameter to 0.5 will scale down both the point estimates and standard errors by half. Default is *1*. This parameter is useful if only part of the sample is employed for estimation, and should not be confused with *Cweights* or *Pweights*.
	massPoints : boolean
		*True* (default) or *False*, specifies whether point estimates and standard errors should be adjusted if there are mass points in the data.
	bwselect : string
		String, specifies the method for data-driven bandwidth selection. This option will be ignored if *bw* is provided. Options are (1) *"mse-dpi"* (default, mean squared error-optimal bandwidth selected for each grid point); (2) *"imse-dpi"* (integrated MSE-optimal bandwidth, common for all grid points); (3) *"mse-rot"* (rule-of-thumb bandwidth with Gaussian reference model); and (4) *"imse-rot"* (integrated rule-of-thumb bandwidth with Gaussian reference model).
	stdVar : boolean
		*True* (default) or *False*, specifies whether the data should be standardized for bandwidth selection.
	regularize : boolean
		*True* (default) or *False*, specifies whether the bandwidth should be regularized. When set to *True*, the bandwidth is chosen such that the local region includes at least *nLocalMin* observations and at least *nUniqueMin* unique observations.
	nLocalMin : int
		Nonnegative integer, specifies the minimum number of observations in each local neighborhood. This option will be ignored if *regularize=False*. Default is *20+p+1*.
	nUniqueMin : int
		Nonnegative integer, specifies the minimum number of unique observations in each local neighborhood. This option will be ignored if *regularize=False*. Default is *20+p+1*.
	Cweights : numeric vector
		Numeric, specifies the weights used for counterfactual distribution construction. Should have the same length as the data.
	Pweights : numeric vector
		Numeric, specifies the weights used in sampling. Should have the same length as the data.
	
	Returns
	-------
	Estimate:
		A matrix containing (1) *grid* (grid points), (2) *bw* (bandwidths), (3) *nh* (number of observations in each local neighborhood), (4) *nhu* (number of unique observations in each local neighborhood), (5) *f_p* (point estimates with p-th order local polynomial), (6) *f_q* (point estimates with q-th order local polynomial, only if option *q* is nonzero), (7) *se_p* (standard error corresponding to *f_p*), and (8) *se_q* (standard error corresponding to *f_q*).
	CovMat_p :
		The variance-covariance matrix corresponding to *f_p*.
	CovMat_q :
		The variance-covariance matrix corresponding to *f_q*.
	"""
	
	##############################
	##   INPUT ERROR HANDLING   ##
	##############################
	
	#data
	data = pd.DataFrame(data)
	if data.isnull().values.any():
		raise Exception(data.isnull().sum() + 'missing observation(s) are ignored.\n')
		data = data.dropna(axis=0)
	
	n = len(data)
	data = np.array(data).reshape(n,)
	if n==0 or not np.isscalar(n):
		raise Exception("Data should be numeric, and cannot be empty.")
	
	#grid
	if grid is None:
		grid = np.quantile(data, np.arange(0.05, 0.95, 0.05))
		ng = len(grid)
	else:
		ng = len(grid)
		if not np.isreal(grid).any():
			raise Exception("Grid points should be numeric.")

	if p is None:
		p = 2
	elif not np.isscalar(p) or p not in range(21):
		raise Exception("Polynomial order p incorrectly specified.")

	if q is None:
		q = p + 1
	elif not np.isscalar(q) or q not in range(21) or q<p:
		raise Exception("Polynomial order (for bias correction) q incorrectly specified")

	if v is None:
		v = min(1, p)
	elif not np.isscalar(q) or v not in range(21) or v>p:
		raise Exception("Derivative order v incorrectly specified")

	bw_list = ['mse-dpi','imse-dpi', 'mse-rot', 'imse-rot']
	if bw is None:
		if len(bwselect) == 0:
			bwselect = 'mse-dpi'
		else:
			bwselect = bwselect.lower()
			if bwselect not in bw_list:
				raise Exception("bwselect incorrectly specified")
	elif len(bw) == 1:
		if not np.isreal(bw) or any(bw) < 0:
			raise Exception("Bandwidth incorrectly specified.")
		else:
			bw = np.repeat(bw, ng)
			bwselect = 'user provided'
	else:
		bw = np.asarray(bw)
		if not np.isreal(bw):
			raise Exception("Bandwidth incorrectly specified.")
		elif len(bw) != ng:
			raise Exception("Bandwidth has to be the same length as grid.")
		else:
			bwselect = 'user provided'
			
	
	#kernel
	if kernel is None:
		flag_no_kernel = True
		kernel = 'epanechnikov'
	else:
		kernel = kernel.lower()
		if kernel not in ['triangular', 'uniform', 'epanechnikov']:
			raise Exception("Kernel function incorrectly specified.")

	# Cweights
	if Cweights is None:
		Cweights = np.repeat(1,n)
	elif not np.isreal(Cweights).any():
		raise Exception("Counterfactual weights incorrectly specified.")
	elif len(Cweights) !=n:
		raise Exception("Counterfactual weights should have the same length as sample.")
	
	if isinstance(Cweights, pd.Series):
		Cweights = Cweights.to_numpy()
	
	# Pweights
	if Pweights is None:
		Pweights = np.repeat(1,n)
	elif not np.isreal(Pweights).any():
		raise Exception("Probability weights incorrectly specified.")
	elif len(Pweights) !=n:
		raise Exception("Probability weights should have the same length as sample.")
	elif any(x<0 for x in Pweights)==True:
		raise Exception("Probability weights should be nonnegative.")
	
	if isinstance(Pweights, pd.Series):
		Pweights = Pweights.to_numpy()
	
	# scale
	if scale is None:
		scale = 1
	elif len(scale) != 1 or not np.isreal(scale) or scale <= 0:
		raise Exception("Scale incorrectly specified.")

	# regularize
	if regularize is None:
		regularize = True
	elif type(regularize) != bool:
		raise Exception("Regularization parameter incorrectly specified")

	# massPoints
	if massPoints is None:
		massPoints = True
	elif type(massPoints) != bool:
		raise Exception("Option massPoints incorrectly specified.")


	# stdVar
	if stdVar is None:
		stdVar = True
	elif type(stdVar) != bool:
		raise Exception("Option stdVar incorrectly specified.")

	# nLocalMin
	if nLocalMin is None:
		nLocalMin = 20 + p + 1
	if not np.isreal(nLocalMin) or np.isnan(nLocalMin):
		raise Exception("Option nLocalMin incorrectly specified.")
	elif math.ceil(nLocalMin)<0:
		raise Exception("Option nLocalMin incorrectly specified.")
	else:
		nLocalMin = math.ceil(nLocalMin)


	# nUniqueMin
	if nUniqueMin is None:
		nUniqueMin = 20 + p + 1
	if not np.isreal(nUniqueMin) or np.isnan(nUniqueMin):
		raise Exception("Option nUniqueMin incorrectly specified.")
	elif math.ceil(nUniqueMin)<0:
		raise Exception("Option nUniqueMin incorrectly specified.")
	else:
		nUniqueMin = math.ceil(nUniqueMin)

	#########################
	##   Sample Trimming   ##
	#########################
	
	#	trim_index = (Pweights == 0)
#	if np.all(trim_index==True):
#		raise Exception("All weights are zero.")
#	else:
#		data = data[trim_index == False]
#		Cweights = Cweights[trim_index == False]
#		Pweights = Pweights[trim_index == False]

		if abs(sum(Cweights * Pweights)) <=np.finfo(float).eps * 10:
			raise Exception("Composited weights (Cweights * Pweights) are numerically zero.")

	#############################
	##   Bandwidth Selection   ##
	#############################
	if bwselect == 'mse-dpi':
		bw = lpbwdensity_fn.__bw_MSE( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
	elif bwselect == 'imse-dpi':
		bw =  lpbwdensity_fn.__bw_IMSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
		bw = np.repeat(bw, len(grid))
	elif bwselect == 'mse-rot':
		bw = lpbwdensity_fn.__bw_ROT( data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
	else:
		bw = lpbwdensity_fn.__bw_IROT(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, stdVar=stdVar, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin)
		bw = np.repeat(bw, len(grid))
		
	##########################
	##   Point Estimation   ##
	##########################
	Temp_Result = lpdensity_fn.__lpdensity_fn(data=data, grid=grid, bw=bw, p=p, q=q, v=v, kernel=kernel,Cweights=Cweights, Pweights=Pweights, massPoints=massPoints, showSE=True)
	
	#scaling output
	Temp_Result['Estimate'][['f_p', 'f_q', 'se_p', 'se_q']] = Temp_Result['Estimate'][['f_p', 'f_q', 'se_p', 'se_q']] * scale
	Temp_Result['CovMat_p'] = Temp_Result['CovMat_p']* np.power(scale, 2)
	Temp_Result['CovMat_q'] = Temp_Result['CovMat_q']* np.power(scale, 2)

	
	#################################
	##   Return lpdensity Object   ##
	#################################
	
	lpoutput = lpdensity_output(Temp_Result['Estimate'], Temp_Result['CovMat_p'], Temp_Result['CovMat_q'], p, q, v, kernel, scale, massPoints, n, ng, bwselect, stdVar, regularize, nLocalMin, nUniqueMin, data_min = min(data), data_max = max(data), grid_min = min(grid), grid_max = max(grid))
	
	return(lpoutput)



#class of lpdensity output
class lpdensity_output:
	"""
	Class of lpdensity function outputs.
	
	Object type returned by :py:meth:`~lpdensity`.
	"""
	def __init__(self, Estimate, CovMat_p, CovMat_q, p, q, v, kernel, scale, massPoints, n, ng, bwselect, stdVar, regularize, nLocalMin, nUniqueMin, data_min, data_max, grid_min, grid_max):
		self.Estimate = Estimate
		self.CovMat_p = CovMat_p
		self.CovMat_q = CovMat_q
		self.p = p
		self.q = q
		self.v = v
		self.kernel = kernel
		self.scale = scale
		self.massPoints = massPoints
		self.n = n
		self.ng = ng
		self.bwselect = bwselect
		self.stdVar = stdVar
		self.regularize = regularize
		self.nLocalMin = nLocalMin
		self.nUniqueMin = nUniqueMin
		self.data_min = data_min
		self.data_max = data_max
		self.grid_min = grid_min
		self.grid_max = grid_max
	
	#summary output
	def __repr__(self, CIuniform=False, alpha=0.05):
		print('Call: lpdensity')
		print('')
		fw = 30
		fw_r = 14
		print('Sample Size:'.ljust(fw), str(self.n).rjust(25))
		print('Polynomial order for point estimation    (p=):'.ljust(fw), str(self.p).rjust(9))
		print('Order of derivative estimated            (v=):'.ljust(fw), str(self.v).rjust(9))
		print('Polynomial order for confidence interval (q=):'.ljust(fw), str(self.q).rjust(9))
		print('Kernel function:'.ljust(fw), str(self.kernel).rjust(25))
		print('Scaling factor:'.ljust(fw), str(self.scale).rjust(25))
		print('Bandwidth Method:'.ljust(fw), str(self.bwselect).rjust(25))
		print('')
		
		#critical value
		ci_df = pd.DataFrame(columns=['CI_l_q', 'CI_r_q'])

		if CIuniform==True:
			if np.isnan(CIsimul):
				raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
				z_val = spstat.norm.ppf(1-alpha/2)
			elif math.ceil(CIsimul)<2:
				raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
				z_val = spstat.norm.ppf(1-alpha/2)
			else:
				corrMat = self.CovMat_q.mul(1/self.Estimate['se_q'], axis=0).mul(1/self.Estimate['se_q'], axis=1)
				#cehck PSD matrix
				if np.all(np.linalg.eigvals(corrMat)>=0):
					#simulate multivariate normal
					normalSimul = np.random.multivariate_normal(mean = np.repeat(0, len(self.CovMat_q)), cov = corrMat, size=2000)
					#comput estimated quantile
					z_val = np.quantile(list(map(lambda x: max(abs(x)), normalSimul)), 1-alpha)
					
				else:
					#exception for non-PSD covariance matrix
					raise warnings.warn('Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.')
					z_val = spstat.norm.ppf(1-alpha/2)
		else:
			z_val = spstat.norm.ppf(1-alpha/2)
		
		ci_df['CI_l_q'] = self.Estimate['f_q'] - z_val * self.Estimate['se_q']
		ci_df['CI_r_q'] = self.Estimate['f_q'] + z_val * self.Estimate['se_q']
		
		#Table
		fw_l = 15
		fw_c = 8
		fw_ci = 18
		n_dec = 3
		print('='*81)
		print('Index'.ljust(fw_c), 'Grid'.ljust(fw_c),
			 'B.W.'.ljust(fw_c), 'Eff.n'.ljust(fw_c), 'Uniq.n'.ljust(fw_c), 'Point.Est.'.rjust(fw_c),  'Std.Err'.rjust(fw_c), 'C.I.'.rjust(fw_c))
		print('-'*81)
		for j in range(self.ng):
			print(str(j+1).ljust(fw_c), str(round(self.Estimate['grid'][j], n_dec)).ljust(fw_c),
			str(round(self.Estimate['bw'][j], n_dec)).ljust(fw_c), str(int(self.Estimate['nh'][j])).ljust(fw_c), str(int(self.Estimate['nhu'][j])).ljust(fw_c), str(round(self.Estimate['f_p'][j], n_dec)).rjust(fw_c), str(round(self.Estimate['se_p'][j], n_dec)).rjust(fw_c), str('['+str(round(ci_df['CI_l_q'][j], n_dec)) + ',' + str(round(ci_df['CI_r_q'][j], n_dec))+']').rjust(fw_l))
		print('='*81)
		return('')
		
	#print call
	def __str__(self):
		print('Call: lpdensity')
		print('')
		fw = 30
		fw_r = 14
		print('Sample Size:'.ljust(fw), str(self.n).rjust(25))
		print('Polynomial order for point estimation    (p=):'.ljust(fw), str(self.p).rjust(9))
		print('Order of derivative estimated            (v=):'.ljust(fw), str(self.v).rjust(9))
		print('Polynomial order for confidence interval (q=):'.ljust(fw), str(self.q).rjust(9))
		print('Kernel function:'.ljust(fw), str(self.kernel).rjust(25))
		print('Scaling factor:'.ljust(fw), str(self.scale).rjust(25))
		print('Bandwidth Method:'.ljust(fw), str(self.bwselect).rjust(25))
		print('')
		print('Use .repr() to show estimates.')
		return ('')
	
	#coefficients of estimation
	def coef(self):
		"""
		Returns estimate coefficients.
		"""
		coef_data = self.Estimate[['grid', 'f_p', 'f_q']]
		return(coef_data)
	
	#variance-covariance
	def vcov(self):
		"""
		Returns estimate standard error and covariance matrices.
		"""
		vcov_data = self.Estimate[['grid', 'se_p', 'se_q']]
		return({'stdErr': vcov_data, 'CovMat_p': self.CovMat_p, 'CovMat_q': self.CovMat_q})
	
	#confidence interval function
	def confint(self, alpha=0.05, CIuniform=False, CIsimul=2000):
		"""
		Returns confindence intervals/bands for prespecified confidence level.
		
		alpha : number
			Confindence level, must be between 0 and 1.
		CIuniform : boolean
			Boolean on wehther to construct uniform confidence bands, *True* or *False (default)*.
		CIsimul : int
			Number of simulations used to generate confidence intervals.
		"""
		if alpha<0 or alpha>1:
			raise Exception("Significance level incorrectly specified.")
		
		#initialized output DF
		estimate_df = pd.DataFrame(columns=['grid', 'f_p', 'CI_l_p', 'CI_r_p', 'f_q', 'CI_l_q', 'CI_r_q'])
		
		estimate_df['grid'] = self.Estimate['grid']
		estimate_df['f_p'] = self.Estimate['f_p']
		estimate_df['f_q'] = self.Estimate['f_q']
		
		#critical value
		if CIuniform==True:
			if np.isnan(CIsimul):
				raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
				z_val = spstat.norm.ppf(1-alpha/2)
			elif math.ceil(CIsimul)<2:
				raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
				z_val = spstat.norm.ppf(1-alpha/2)
			else:
				CIsimul = math.ceil(CIsimul)
				corrMat = self.CovMat_q.mul(1/self.Estimate['se_q'], axis=0).mul(1/self.Estimate['se_q'], axis=1)
				#cehck PSD matrix
				if np.all(np.linalg.eigvals(corrMat)>=0):
					#simulate multivariate normal
					normalSimul = np.random.multivariate_normal(mean = np.repeat(0, len(self.CovMat_q)), cov = corrMat, size=CIsimul)
					#comput estimated quantile
					z_val = np.quantile(list(map(lambda x: max(abs(x)), normalSimul)), 1-alpha)
					
				else:
					#exception for non-PSD covariance matrix
					raise warnings.warn('Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.')
					z_val = spstat.norm.ppf(1-alpha/2)
		else:
			z_val = spstat.norm.ppf(1-alpha/2)
		
		estimate_df['CI_l_p'] = self.Estimate['f_p'] - z_val * self.Estimate['se_p']
		estimate_df['CI_r_p'] = self.Estimate['f_p'] + z_val * self.Estimate['se_p']
		
		estimate_df['CI_l_q'] = self.Estimate['f_q'] - z_val * self.Estimate['se_q']
		estimate_df['CI_r_q'] = self.Estimate['f_q'] + z_val * self.Estimate['se_q']
		
		return(estimate_df)
	
	#plotting function
	def plot(self, alpha=0.05, type='line', CItype='region', CIuniform=False, CIsimul=2000, hist=False, histData=None, histBins=None, histFillCol=3, histFillShade=0.2, histLineCol="white", title=None, xlabel=None, ylabel=None, CIshade=0.2):
	
		"""
		Method to plot estimate and confidence bands.
		Requires ggplot.
		
		alpha : number
			Confindence level, must be between 0 and 1.
		CIuniform : boolean
			Boolean on wehther to construct uniform confidence bands, *True* or *False (default)*.
		CIsimul : int
			Number of simulations used to generate confidence intervals.
		type : string
			type of estimate plot, *line (default)*, or *points*, or *all*.
		CItype : string
			type of confidence interval plot, *region (default)*, *lines*, *ebar*, or *all*.
		"""
		if self.ng <2:
			raise Exception("At least two grid points are needed to plot input.")
			
		if alpha<0 or alpha>1:
			raise Exception("Significance level incorrectly specified.")
		#TODO: error handling
		if hist==True and histData is not None:
			histData = pd.DataFrame(histData, columns='d1')
			if histBins is None:
				histBreaks = 25
				temp_plot = (pn.ggplot(histData, pn.aes(x='d1')) + pn.geom_histogram(bins=histBins) + pn.theme_bw())
		else:
			temp_plot = pn.ggplot() + pn.theme_bw()
		
		
		data_x = self.Estimate[['grid', 'f_p', 'f_q', 'se_p', 'se_q']]
		
		#critical value
		if CIuniform==True:
			if np.isnan(CIsimul):
				raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
				z_val = spstat.norm.ppf(1-alpha/2)
			elif math.ceil(CIsimul)<2:
				raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
				z_val = spstat.norm.ppf(1-alpha/2)
			else:
				CIsimul = math.ceil(CIsimul)
				corrMat = self.CovMat_q.mul(1/self.Estimate['se_q'], axis=0).mul(1/self.Estimate['se_q'], axis=1)
				#cehck PSD matrix
				if np.all(np.linalg.eigvals(corrMat)>=0):
					#simulate multivariate normal
					normalSimul = np.random.multivariate_normal(mean = np.repeat(0, len(self.CovMat_q)), cov = corrMat, size=CIsimul)
					#comput estimated quantile
					z_val = np.quantile(list(map(lambda x: max(abs(x)), normalSimul)), 1-alpha)
					
				else:
					#exception for non-PSD covariance matrix
					raise warnings.warn('Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.')
					z_val = spstat.norm.ppf(1-alpha/2)
		else:
			z_val = spstat.norm.ppf(1-alpha/2)
				
		data_x['CI_l'] = data_x['f_q'] - z_val * data_x['se_q']
		data_x['CI_r'] = data_x['f_q'] + z_val * data_x['se_q']
		
		if CItype in ['region', 'all']:
			temp_plot = temp_plot + pn.geom_ribbon(data=data_x, mapping=pn.aes(x='grid', ymin='CI_l', ymax='CI_r'), alpha=CIshade)
		
		if CItype in ['lines', 'all']:
			temp_plot = temp_plot + pn.geom_line(data=data_x, mapping=pn.aes(x='grid', y='CI_l')) + pn.geom_line(data=data_x, mapping=pn.aes(x='grid', y='CI_r'))
			
		if CItype in ['ebar', 'all']:
			temp_plot = temp_plot + pn.geom_errorbar(data=data_x, mapping=pn.aes(x='grid', ymin='CI_l', ymax='CI_r'))
		
		if type in ['line', 'both']:
			temp_plot = temp_plot + pn.geom_line(data=data_x, mapping=pn.aes(x='grid', y='f_p'))
		
		if type in ['points', 'both']:
			temp_plot = temp_plot + pn.geom_point(data=data_x, mapping=pn.aes(x='grid', y='f_p'))
		
		if xlabel is None:
			xlabel = ''
		
		if ylabel is None:
			ylabel = ''
		
		if title is None:
			title = ''
			
		temp_plot = temp_plot + pn.xlab(xlabel) + pn.ylab(ylabel)
		
		return(temp_plot)
