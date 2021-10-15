#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ignore future warnings from pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy.stats import norm
import sympy
from sympy.abc import x, y
import math
import scipy.integrate as integrate
from numpy.linalg import inv
import scipy.optimize as optimize
import statistics as stat
import pandas as pd


def __kernel_fn(a, h, kernel):
	if kernel == "triangular":
		kh = (1-abs(a))/h
	elif kernel == "uniform":
		kh = np.repeat(0.5/h, len(a))
	else:
		kh = 0.75*(1-(a)**2)/h
	return(kh)

def __lpdensityUnique(x):
	n = len(x)
	x = pd.DataFrame(x)
	# if x has less than 2 elements
	if n==0:
		return({'unique':None, 'freq':[None], 'index':[None]})
	if n==1:
		return({'unique':x, 'freq':[1], 'index':[0]})
		
	#if x has more than 1 element
	unique, uniqueIndex, nUniquelist = np.unique(x, return_index=True, return_counts=True, axis=0)
	nUnique = sum(unique)
#	uniqueIndex = x != x.shift()
#	unique = np.array(x)[uniqueIndex]
#	nUnique = len(unique)
	
	#if all are distinct
	if nUnique==n:
		return({'unique':unique, 'freq':np.repeat(1, len(x)), 'index':list(range(n))})
		
	#if all are the same
	if nUnique ==1:
		return({'unique':unique, 'freq':[n], 'index':[n-1]})
	
#	freq = np.cumsum(uniqueIndex==False)[uniqueIndex].iloc[:,0]
#	freq = freq - pd.concat([pd.Series([0]), freq[:-1]]) + 1
#	freq = pd.DataFrame(freq).reset_index(drop=True).dropna()[0]
	
	return({'unique':unique, 'freq':nUniquelist, 'index':uniqueIndex})

# Internal function
# Gaussian derivatives
# Output matches R
def __normal_dgps(x_val, v, mean, sd):
	if v==0:
		return(norm.cdf(x_val, loc=mean, scale=sd))
	else:
		
		temp = sympy.exp(-(x-mean)**2/(2*sd**2))/sympy.sqrt(2*math.pi*sd**2)
		while v>1:
			temp = sympy.diff(temp, x)
			v = v-1
		return(temp.evalf(subs={x:x_val}))

# Generate S matrix
# Output matches R
def __Sgenerate(p, low=-1, up=1, kernel="triangular"):
	S = np.zeros((p+1, p+1))
	for i in range(1,p+2):
		for j in range(1,p+2):
			if kernel=="uniform":
				integrand = lambda l: 0.5*np.power(l, i+j-2)
			elif kernel=="epanechnikov":
				integrand = lambda l: 0.75*np.power(l, i+j-2)*(1-l**2)
			else:
				integrand = lambda l: (1-abs(l))*np.power(l, i+j-2)
			
			S[i-1, j-1] = integrate.quad(integrand, low, up)[0]
	
	return(S)
	

# Generate T matrix
# Output matches R
def __Tgenerate(p, low=-1, up=1, kernel="triangular"):
	S = np.zeros((p+1, p+1))
	for i in range(1,p+2):
		for j in range(1,p+2):
			if kernel=="uniform":
				integrand = lambda l: np.power(0.5, 2)*np.power(l, i+j-2)
			elif kernel=="epanechnikov":
				integrand = lambda l: np.power(0.75*(1-l**2), 2)*np.power(l, i+j-2)
			else:
				integrand = lambda l: np.power(1-abs(l), 2)*np.power(l, i+j-2)
			
			S[i-1, j-1] = integrate.quad(integrand, low, up)[0]
	
	return(S)

# Generate C matrix
# Output matches R
def __Cgenerate(k, p, low=-1, up=1, kernel="triangular"):
	C = np.zeros((p+1, 1))
	for i in range(1, p+2):
		if kernel=="uniform":
			integrand = lambda l: 0.5*np.power(l, i+k-1)
		elif kernel=="epanechnikov":
			integrand = lambda l: 0.75*np.power(l, i+k-1)*(1-l**2)
		else:
			integrand = lambda l: (1-abs(l))*np.power(l, i+k-1)
			
		C[i-1, 0] = integrate.quad(integrand, low, up)[0]
	
	return(C)

# Generate G matrix
#TODO: Fix triangular kernel
def __Ggenerate(p, low=-1, up=1, kernel="uniform"):
	G = np.zeros((p+1, p+1))
	for i in range(1, p+2):
		for j in range(1, p+2):
			def integrand_1(x,y):
				return(np.power(x,i) * np.power(y, (j-1))*(1-abs(x))*(1-abs(y)))
				
			def integrand_2(x,y):
				return(np.power(x,i-1) * np.power(y, j)*(1-abs(x))*(1-abs(y)))
				
			def x_integral_1(y):
				return (integrate.quad(integrand_1, low, y, args=(y))[0])
			
			def x_integral_2(y):
				return (integrate.quad(integrand_2, y, up, args=(y))[0])
			
			G[i-1, j-1] = integrate.quad(x_integral_1, low, up)[0] + integrate.quad(x_integral_2, low, up)[0]
				
	return(G)



# MSE-ROT Bandwidth function
# matches R output
def __bw_ROT(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin):
	n = len(data)
	ng = len(grid)
	
	dataUnique = np.unique(data)
	nUnique = len(dataUnique)
	
	if stdVar==True:
		center_temp = stat.mean(data)
		scale_temp = stat.stdev(data)
		data = (data-center_temp)/scale_temp
		dataUnique = (dataUnique-center_temp)/scale_temp
		grid = (grid-center_temp)/scale_temp
		
	#estimate normal reference model
	mean_hat = sum(Cweights*Pweights*data)/sum(Cweights*Pweights)
	sd_hat = np.sqrt(sum(Cweights*Pweights*(data-mean_hat)**2)/sum(Cweights*Pweights))
	
	#normal quantities
	temp_1 = sympy.exp(-(x-mean_hat)**2/(2*sd_hat**2))/sympy.sqrt(2*math.pi*sd_hat**2)
	temp_2 = sympy.diff(temp_1, x)
	temp_3 = sympy.exp(-(x-mean_hat)**2/(2*sd_hat**2))/sympy.sqrt(2*math.pi*sd_hat**2)
	temp_4 = sympy.exp(-(x-mean_hat)**2/(2*sd_hat**2))/sympy.sqrt(2*math.pi*sd_hat**2)
	
	j = p+1
	while j>1:
		temp_3 = sympy.diff(temp_3, x)
		j = j-1
	
	j=p+2
	while j >1:
		temp_4 = sympy.diff(temp_4, x)
		j = j-1
		
	#bias estimate, no rate added, DGP constant
	bias_dgp = np.zeros((ng, 2))
	for j in range(ng):
		bias_dgp[j, 0] = temp_3.evalf(subs={x:grid[j]})/math.factorial(p+1) * math.factorial(v)
		bias_dgp[j, 1] = temp_4.evalf(subs={x:grid[j]})/math.factorial(p+2) * math.factorial(v) + bias_dgp[j,0] * temp_2.evalf(subs={x:grid[j]})/temp_1.evalf(subs={x:grid[j]})
	
	#bias estimate, no rate added, kernel constant
	S = __Sgenerate(p=p, low=-1, up=1, kernel=kernel)
	C1 = __Cgenerate(k=p+1, p=p, low=-1, up=1, kernel=kernel)
	C2 = __Cgenerate(k=p+2, p=p, low=-1, up=1, kernel=kernel)
	G = __Ggenerate(p=p, low=-1, up=1, kernel=kernel)
	S2 = __Tgenerate(p=p, low=-1, up=1, kernel=kernel)
	bias_dgp[:,0] = bias_dgp[:,0] * (np.matmul(inv(S), C1)[v,:])
	bias_dgp[:,1] = bias_dgp[:,1] * (np.matmul(inv(S), C2)[v,:])
	
	#variance estimate, sample size added
	sd_dgp = np.zeros((ng, 1))
	if v>0:
		for j in range(ng):
			temp_1eval = np.array(temp_1.evalf(subs={x:grid[j]}), dtype=np.float64)
			sd_dgp[j, 0] = math.factorial(v) * np.sqrt(temp_1eval/n)
			
		sd_dgp = sd_dgp * np.sqrt(abs((inv(S) * G * inv(S))[v, v]))
		
	else:
		for j in range(ng):
			#this comes from a higher-order variance expansion. See Lemma 4 in the Appendix of Cattaneo, Jansson and Ma (2019a)
			sd_dgp[j, 0] = np.sqrt(norm.cdf(grid[j], loc=mean_hat, scale=sd_hat) * (1-norm.cdf(grid[j], loc=mean_hat, scale=sd_hat))/norm.pdf(grid[j], loc=mean_hat, scale=sd_hat)/(0.5*(n**2)))
			
		sd_dgp = sd_dgp * np.sqrt(abs((inv(S) * G * inv(S))[v, v]))
	
	# bandwidth
	h = np.zeros((ng, 1))
	for j in range(ng):
		if v>0:
			opt_f = lambda a: np.power(a, 2*p+2-2*v) * (bias_dgp[j, 0] + a * bias_dgp[j, 1])**2 + sd_dgp[j,0]**2 / np.power(a, 2*v - 1)
		else:
			opt_f = lambda a: np.power(a, 2*p+2-2*v) * (bias_dgp[j, 0] + a * bias_dgp[j, 1])**2 + sd_dgp[j,0]**2 / a
		h[j] = optimize.minimize_scalar(opt_f, bounds=(0, max(data)-min(data)), method='bounded').x
		
	
	for j in range(ng):
		if np.isnan(h[j])==True:
			h[j] = np.sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))]
			
		if regularize==True:
			if nLocalMin>0:
				h[j] = max(h[j], np.sort(abs(data-grid[j]))[min(n, nLocalMin)-1])
			if nUniqueMin>0:
				h[j] = max(h[j], np.sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)-1])
			
			h[j] = min(h[j], max(abs(dataUnique-grid[j])))
			
	if stdVar==True:
		h = h*scale_temp
		
	return(h)


# MSE-IROT Bandwidth function
# matches R output
#slightly off only on optimize function
def __bw_IROT(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin):
	n = len(data)
	ng = len(grid)
	
	dataUnique = np.unique(data)
	nUnique = len(dataUnique)
	
	if stdVar==True:
		center_temp = stat.mean(data)
		scale_temp = stat.stdev(data)
		data = (data-center_temp)/scale_temp
		dataUnique = (dataUnique-center_temp)/scale_temp
		grid = (grid-center_temp)/scale_temp
		
	#estimate normal reference model
	mean_hat = sum(Cweights*Pweights*data)/sum(Cweights*Pweights)
	sd_hat = np.sqrt(sum(Cweights*Pweights*(data-mean_hat)**2)/sum(Cweights*Pweights))
	
	#normal quantities
	temp_1 = sympy.exp(-(x-mean_hat)**2/(2*sd_hat**2))/sympy.sqrt(2*math.pi*sd_hat**2)
	temp_2 = sympy.diff(temp_1, x)
	temp_3 = sympy.exp(-(x-mean_hat)**2/(2*sd_hat**2))/sympy.sqrt(2*math.pi*sd_hat**2)
	temp_4 = sympy.exp(-(x-mean_hat)**2/(2*sd_hat**2))/sympy.sqrt(2*math.pi*sd_hat**2)
	
	j = p+1
	while j>1:
		temp_3 = sympy.diff(temp_3, x)
		j = j-1
	
	j=p+2
	while j >1:
		temp_4 = sympy.diff(temp_4, x)
		j = j-1
		
	#bias estimate, no rate added, DGP constant
	bias_dgp = np.zeros((ng, 2))
	for j in range(ng):
		bias_dgp[j, 0] = temp_3.evalf(subs={x:grid[j]})/math.factorial(p+1) * math.factorial(v)
		bias_dgp[j, 1] = temp_4.evalf(subs={x:grid[j]})/math.factorial(p+2) * math.factorial(v) + bias_dgp[j,0] * temp_2.evalf(subs={x:grid[j]})/temp_1.evalf(subs={x:grid[j]})
	
	#bias estimate, no rate added, kernel constant
	S = __Sgenerate(p=p, low=-1, up=1, kernel=kernel)
	C1 = __Cgenerate(k=p+1, p=p, low=-1, up=1, kernel=kernel)
	C2 = __Cgenerate(k=p+2, p=p, low=-1, up=1, kernel=kernel)
	G = __Ggenerate(p=p, low=-1, up=1, kernel=kernel)
	S2 = __Tgenerate(p=p, low=-1, up=1, kernel=kernel)
	bias_dgp[:,0] = bias_dgp[:,0] * (np.matmul(inv(S), C1)[v,:])
	bias_dgp[:,1] = bias_dgp[:,1] * (np.matmul(inv(S), C2)[v,:])
	
	#variance estimate, sample size added
	sd_dgp = np.zeros((ng, 1))
	if v>0:
		for j in range(ng):
			temp_1eval = np.array(temp_1.evalf(subs={x:grid[j]}), dtype=np.float64)
			sd_dgp[j, 0] = math.factorial(v) * np.sqrt(temp_1eval/n)
			
		sd_dgp = sd_dgp * np.sqrt(abs((np.matmul(inv(S), np.matmul(G, inv(S))))[v, v]))
		
	else:
		for j in range(ng):
			#this comes from a higher-order variance expansion. See Lemma 4 in the Appendix of Cattaneo, Jansson and Ma (2019a)
			sd_dgp[j, 0] = np.sqrt(norm.cdf(grid[j], loc=mean_hat, scale=sd_hat) * (1-norm.cdf(grid[j], loc=mean_hat, scale=sd_hat))/norm.pdf(grid[j], loc=mean_hat, scale=sd_hat)/(0.5*(n**2)))
			
		sd_dgp = sd_dgp * np.sqrt(abs((np.matmul(inv(S), np.matmul(S2, inv(S))))[v, v]))
		
	# bandwidth
	if v>0:
		opt_f = lambda a: np.power(a, 2*p+2-2*v) * sum((bias_dgp[:, 0] + a * bias_dgp[:, 1])**2) + sum(sd_dgp[:,0]**2) / np.power(a, 2*v - 1)
	else:
		opt_f = lambda a: np.power(a, 2*p+2-2*v) * sum((bias_dgp[:, 0] + a * bias_dgp[:, 1])**2) + sum(sd_dgp[:,0]**2) / a
			
	h = optimize.minimize(opt_f, x0=0.5).x
	
	if np.isnan(h)==True:
		h = np.sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))]
			
	if regularize==True:
		for j in range(ng):
			if nLocalMin>0:
				h = max(h, np.sort(abs(data-grid[j]))[min(n, nLocalMin)-1])
			if nUniqueMin>0:
				h = max(h, np.sort(abs(dataUnique-grid[j]))[min(nUnique, nUniqueMin)-1])
			
	h = min(h, max(abs(max(dataUnique)-min(grid)), abs(min(dataUnique)-max(grid))))
			
	if stdVar==True:
		h = h*scale_temp
		
	return(h)
	
	
# MSE-optimal bandwidth

def __bw_MSE(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin):
	ii = np.argsort(data)
	data = data[ii]
	Cweights = Cweights[ii]
	Pweights = Pweights[ii]
	n = len(data)
	ng = len(grid)
	
	dataUnique = __lpdensityUnique(data)
	freqUnique = dataUnique["freq"]
	indexUnique = dataUnique["index"]
	dataUnique = dataUnique["unique"]
	nUnique = len(dataUnique)
	
	if stdVar==True:
		center_temp = stat.mean(data)
		scale_temp = stat.stdev(data)
		data = (data-center_temp)/scale_temp
		dataUnique = (dataUnique-center_temp)/scale_temp
		grid = (grid-center_temp)/scale_temp
	
	# whether considering mass points when constructing the empirical distribution function
	if massPoints == True:
		Fn = np.repeat((np.cumsum(Cweights * Pweights) / sum(Cweights * Pweights))[indexUnique], repeats = freqUnique)
	else:
		Fn = np.cumsum(Cweights * Pweights) / sum(Cweights * Pweights)
	
	weights_normal = Cweights * Pweights / sum(Cweights * Pweights) * n
	Cweights = Cweights / sum(Cweights) * n
	Pweights = Pweights / sum(Pweights) * n
	
	weights_normalUnique = np.cumsum(weights_normal)[indexUnique]
	if nUnique > 1:
		weights_normalUnique = weights_normalUnique - np.insert(weights_normalUnique[0:(nUnique-1)],0,0)
	
	CweightsUnique = np.cumsum(Cweights)[indexUnique]
	if nUnique > 1:
		CweightsUnique = CweightsUnique - np.insert(CweightsUnique[0:(nUnique-1)],0,0)
	
	PweightsUnique = np.cumsum(Pweights)[indexUnique]
	if nUnique > 1:
		PweightsUnique = PweightsUnique - np.insert(PweightsUnique[0:(nUnique-1)],0,0)
		
	# Preliminary Bandwidth estimate
	# bandwidth for preasymptotic matrics
	h1 = __bw_IROT(data=data, grid=grid, p=2, v=1, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=True, stdVar=True, regularize=True, nLocalMin=20+2+1,  nUniqueMin=20+2+1)
	#bandwidth for F_p+1
	hp1 = __bw_IROT(data=data, grid=grid, p=p+2, v=p+1, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=True, stdVar=True, regularize=True, nLocalMin=20+p+2+1,  nUniqueMin=20+p+2+1)
	
	#bandwidth fro F_p+2
	hp2 = __bw_IROT(data=data, grid=grid, p=p+3, v=p+2, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=True, stdVar=True, regularize=True, nLocalMin=20+p+3+1,  nUniqueMin=20+p+3+1)
	
	#density estimates with normalization constants
	dgp_hat = np.zeros((ng, 2))
	const_hat = np.zeros((ng, 3))
	h = np.repeat(np.nan, ng)
	
	for j in range(ng):
		#estimate F_p+2
		index_temp = abs(data-grid[j]) <= hp2
		Xh_temp = pd.DataFrame(data[index_temp]-grid[j])/hp2
		Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+3+1)]
		Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+3+1)), axis=1)
		
		
		#kernel function calculations
		Kh_temp = pd.DataFrame(__kernel_fn(Xh_temp, hp2, kernel))
		Kh_temp = Kh_temp.mul(Pweights[index_temp], axis=0)
		Y_temp = pd.DataFrame(Fn[index_temp])

		#try solving matrix inverse
		try:
			temp = (np.matmul(np.matmul(inv(np.matmul(Xh_p_temp.transpose(), Xh_p_temp.mul(Kh_temp, axis=0))), Xh_p_temp.transpose()), Y_temp.mul(Kh_temp, axis=0))).loc[p+2, 0]/np.power(hp2, p+2)
		except np.linalg.LinAlgError:
				continue
		dgp_hat[j,1] = temp
		
		

		#estimate F_p+1
		index_temp = abs(data-grid[j]) <= hp1
		Xh_temp = pd.DataFrame(data[index_temp]-grid[j])/hp1
		Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+2+1)]
		Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+2+1)), axis=1)
		
		#kernel function calculations
		Kh_temp = pd.DataFrame(__kernel_fn(Xh_temp, hp1, kernel))
#		Kh_temp = Kh_temp.mul(Pweights[index_temp], axis=0)
		Y_temp = pd.DataFrame(Fn[index_temp])
		
		#try solving matrix inverse
		try:
			temp = (np.matmul(np.matmul(inv(np.matmul(Xh_p_temp.transpose(), Xh_p_temp.mul(Kh_temp, axis=0))), Xh_p_temp.transpose()), Y_temp.mul(Kh_temp, axis=0))).loc[p+1, 0] / (hp1**(p+1))
		except np.linalg.LinAlgError:
				continue
		
		dgp_hat[j,0] = temp
		
		#estimating matrices
		index_temp = abs(data-grid[j]) <= h1
		Xh_temp = pd.DataFrame(data[index_temp]-grid[j])/h1
		Kh_temp = pd.DataFrame(__kernel_fn(Xh_temp, h1, kernel))
		
		#Cp matrix
		if p>0:
			C_p_hat = pd.DataFrame(Xh_temp[np.repeat(Xh_temp.columns, p+1)].apply(lambda x: np.power(x, range(p+1, 2*p+1+1)), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)
		else:
			C_p_hat = pd.DataFrame(Xh_temp.apply(lambda x: np.power(x, 1), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)
		
		#Cp+1 matrix
		if p>0:
			C_p1_hat = pd.DataFrame(Xh_temp[np.repeat(Xh_temp.columns, p+1)].apply(lambda x: np.power(x, range(p+2, 2*p+2+1)), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)
		else:
			C_p1_hat = pd.DataFrame(Xh_temp.apply(lambda x: np.power(x, 2), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)

		#S matrix
		Xh_p_temp = pd.DataFrame(Xh_temp[np.repeat(Xh_temp.columns, p+1)].apply(lambda x: np.power(x, range(p+1)), axis=1))
		
		
		S_hat = np.matmul(Xh_p_temp.transpose(),Xh_p_temp.mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0))/n
		
		try:
			S_hat_inv = inv(S_hat)
		except np.linalg.LinAlgError:
			continue
			
		#G matrix
		if v==0:
			G_hat = np.matmul(Xh_p_temp, Xh_p_temp.mul(Kh_temp.mul(Pweights[index_temp], axis=0)**2)/n)
		else:
			if massPoints==True:
				Y_temp = pd.DataFrame(Fn[indexUnique])
				Xh_temp = pd.DataFrame(dataUnique-grid[j])/h1
				Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+1)]
				Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+1)), axis=1)
				if kernel == "triangular":
					Kh_temp = (1-abs(Xh_temp))/h1
					Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
				elif kernel == "uniform":
					Kh_temp = index_temp[indexUnique] * 0.5/h1
				else:
					Kh_temp = 0.75*(1-np.power(Xh_temp,2))/h1
					Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
				
				Xh_p_Kh_temp = Xh_p_temp.mul(Kh_temp, axis=0)
				Xh_p_Kh_Pweights_temp = Xh_p_Kh_temp.mul(PweightsUnique, axis=0)
				
				F_Xh_p_Kh_temp = np.matmul(Y_temp.transpose(), Xh_p_Kh_Pweights_temp)/n

				G = np.zeros((n,  len(Xh_p_Kh_Pweights_temp.columns)))
				for jj in range(len(Xh_p_Kh_Pweights_temp.columns)):
					G[:,jj] = (np.repeat((np.cumsum(Xh_p_Kh_Pweights_temp.iloc[::-1, jj])/n)[::-1], freqUnique) - F_Xh_p_Kh_temp.loc[0, jj]) * weights_normal
				
				G = pd.DataFrame(G)
				G_hat = np.matmul(G.transpose(), G)/n
				
			else:
				Y_temp = pd.DataFrame(Fn)
				Xh_temp = pd.DataFrame(data-grid[j])/h1
				Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+1)]
				Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+1)), axis=1)
				Kh_temp = __kernel_fn(Xh_temp, h1, kernel)
				
				Xh_p_Kh_temp = Xh_p_temp.mul(Kh_temp, axis=0)
				Xh_p_Kh_Pweights_temp = Xh_p_Kh_temp.mul(Pweights, axis=0)
				
				F_Xh_p_Kh_temp = np.matmul(Y_temp.transpose(), Xh_p_Kh_Pweights_temp)/n
				
				
				G = np.zeros((n,  len(Xh_p_Kh_Pweights_temp.columns)))
				for jj in range(len(Xh_p_Kh_Pweights_temp.columns)):
					G[:,jj] = (np.repeat(np.cumsum(Xh_p_Kh_Pweights_temp.iloc[:, jj]/n), freqUnique) - F_Xh_p_Kh_temp.loc[0, jj]) * weights_normal
				
				G = pd.DataFrame(G)
				G_hat = (np.matmul(G.transpose(), G))/n
				
		#Constants
		const_hat[j, 0] = math.factorial(v) * (np.matmul(S_hat_inv, C_p_hat)).loc[v, 0]
		const_hat[j, 1] = math.factorial(v) * (np.matmul(S_hat_inv, C_p1_hat)).loc[v, 0]
		
		if v>0:
			const_hat[j, 2] = math.factorial(v) * np.sqrt(abs(np.matmul(S_hat_inv, np.matmul(G_hat, S_hat_inv)).loc[v,v]) / (n*h1))
		else:
			temp_ii = min(max(mean(data <= grid[j]), 1/n), 1 - 1/n )
			const_hat[j, 2] = math.factorial(v) * np.sqrt(abs(np.matmul(S_hat_inv, np.matmul(G_hat, S_hat_inv)).loc[v,v] / (0.5*n**2)*h1 * temp_ii*(1-temp_ii)))

		#optimal bandwidth
		if v>0:
			opt_f = lambda a: np.power(a, 2*p+2-2*v) * (dgp_hat[j, 0]*const_hat[j,0] + a * dgp_hat[j, 1]*const_hat[j,1])**2 + const_hat[j,2]**2 / np.power(a, 2*v - 1)
		else:
			opt_f = lambda a: np.power(a, 2*p+2-2*v) * (dgp_hat[j, 0]*const_hat[j,0] + a * dgp_hat[j, 1]*const_hat[j,1])**2 + const_hat[j,2]**2 / a
		
		h[j] = optimize.minimize_scalar(opt_f, bounds=(0, max(data)-min(data)), method='bounded').x
	

	for j in range(ng):
		if np.isnan(h[j])==True:
			h[j] = np.sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))]
			
		if regularize==True:
			if nLocalMin>0:
				h[j] = max(h[j], np.sort(abs(data-grid[j]))[min(n-1, nLocalMin-1)])
			if nUniqueMin>0:
				h[j] = max(h[j], sorted(abs(dataUnique-grid[j]))[min(nUnique-1, nUniqueMin-1)])
			
			h[j] = min(h[j], max(abs(dataUnique-grid[j])))
			
	if stdVar==True:
		h = h*scale_temp
		
	return(h)

#IMSE optimal bandwidth

def __bw_IMSE(data, grid, p, v, kernel, Cweights, Pweights, massPoints, stdVar, regularize, nLocalMin, nUniqueMin):
	ii = np.argsort(data)
	data = data[ii]
#	Cweights = Cweights.reset_index(drop=True)
	Cweights = Cweights[ii]
	Pweights = Pweights[ii]
	n = len(data)
	ng = len(grid)
	
	dataUnique = __lpdensityUnique(data)
	freqUnique = dataUnique["freq"]
	indexUnique = dataUnique["index"]
	dataUnique = dataUnique["unique"]
	nUnique = len(dataUnique)
	
	if stdVar==True:
		center_temp = stat.mean(data)
		scale_temp = stat.stdev(data)
		data = (data-center_temp)/scale_temp
		dataUnique = (dataUnique-center_temp)/scale_temp
		grid = (grid-center_temp)/scale_temp
	
	# whether considering mass points when constructing the empirical distribution function
	if massPoints == True:
		Fn = np.repeat((np.cumsum(Cweights * Pweights) / sum(Cweights * Pweights))[indexUnique], repeats = freqUnique)
	else:
		Fn = np.cumsum(Cweights * Pweights) / sum(Cweights * Pweights)
		
	weights_normal = Cweights * Pweights / sum(Cweights * Pweights) * n
	Cweights = Cweights / sum(Cweights) * n
	Pweights = Pweights / sum(Pweights) * n
	
	weights_normalUnique = np.cumsum(weights_normal)[indexUnique]
	if nUnique > 1:
		weights_normalUnique = weights_normalUnique - np.insert(weights_normalUnique[:-1],0,0)
	
	CweightsUnique = np.cumsum(Cweights)[indexUnique]
	if nUnique > 1:
		CweightsUnique = CweightsUnique - np.insert(CweightsUnique[:-1],0,0)
	
	PweightsUnique = np.cumsum(Pweights)[indexUnique]
	if nUnique > 1:
		PweightsUnique = PweightsUnique - np.insert(PweightsUnique[:-1],0,0)
		
	# Preliminary Bandwidth estimate
	# bandwidth for preasymptotic matrics
	h1 = __bw_IROT(data=data, grid=grid, p=2, v=1, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=True, stdVar=True, regularize=True, nLocalMin=20+2+1,  nUniqueMin=20+2+1)
	
	#bandwidth for F_p+1
	hp1 = __bw_IROT(data=data, grid=grid, p=p+2, v=p+1, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=True, stdVar=True, regularize=True, nLocalMin=20+p+2+1,  nUniqueMin=20+p+2+1)

	#bandwidth fro F_p+2
	hp2 = __bw_IROT(data=data, grid=grid, p=p+3, v=p+2, kernel=kernel, Cweights=Cweights, Pweights=Pweights, massPoints=True, stdVar=True, regularize=True, nLocalMin=20+p+3+1,  nUniqueMin=20+p+3+1)

	#density estimates with normalization constants
	dgp_hat = np.zeros((ng, 2))
	const_hat = np.zeros((ng, 3))
	
	for j in range(ng):
		#estimate F_p+2
		index_temp = abs(data-grid[j]) <= hp2
		Xh_temp = pd.DataFrame(data[index_temp]-grid[j])/hp2
		Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+3+1)]
		Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+3+1)), axis=1)
		
		
		#kernel function calculations
		Kh_temp = pd.DataFrame(__kernel_fn(Xh_temp, hp2, kernel))
		Kh_temp = Kh_temp.mul(Pweights[index_temp], axis=0)
		Y_temp = pd.DataFrame(Fn[index_temp])

		#try solving matrix inverse
		try:
			temp = (np.matmul(np.matmul(inv(np.matmul(Xh_p_temp.transpose(), Xh_p_temp.mul(Kh_temp, axis=0))), Xh_p_temp.transpose()), Y_temp.mul(Kh_temp, axis=0))).loc[p+2, 0]/np.power(hp2, p+2)
		except np.linalg.LinAlgError:
				continue
		dgp_hat[j,1] = temp
		
		

		#estimate F_p+1
		index_temp = abs(data-grid[j]) <= hp1
		Xh_temp = pd.DataFrame(data[index_temp]-grid[j])/hp1
		Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+2+1)]
		Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+2+1)), axis=1)
		
		#kernel function calculations
		Kh_temp = pd.DataFrame(__kernel_fn(Xh_temp, hp1, kernel))
#		Kh_temp = Kh_temp.mul(Pweights[index_temp], axis=0)
		Y_temp = pd.DataFrame(Fn[index_temp])
		
		#try solving matrix inverse
		try:
			temp = (np.matmul(np.matmul(inv(np.matmul(Xh_p_temp.transpose(), Xh_p_temp.mul(Kh_temp, axis=0))), Xh_p_temp.transpose()), Y_temp.mul(Kh_temp, axis=0))).loc[p+1, 0] / (hp1**(p+1))
		except np.linalg.LinAlgError:
				continue
		
		dgp_hat[j,0] = temp
		
		#estimating matrices
		index_temp = abs(data-grid[j]) <= h1
		Xh_temp = pd.DataFrame(data[index_temp]-grid[j])/h1
		Kh_temp = pd.DataFrame(__kernel_fn(Xh_temp, h1, kernel))
		
		#Cp matrix
		if p>0:
			C_p_hat = pd.DataFrame(Xh_temp[np.repeat(Xh_temp.columns, p+1)].apply(lambda x: np.power(x, range(p+1, 2*p+1+1)), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)
		else:
			C_p_hat = pd.DataFrame(Xh_temp.apply(lambda x: np.power(x, 1), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)
		
		#Cp+1 matrix
		if p>0:
			C_p1_hat = pd.DataFrame(Xh_temp[np.repeat(Xh_temp.columns, p+1)].apply(lambda x: np.power(x, range(p+2, 2*p+2+1)), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)
		else:
			C_p1_hat = pd.DataFrame(Xh_temp.apply(lambda x: np.power(x, 2), axis=1).mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0).sum(axis=0)/n)

		#S matrix
		Xh_p_temp = pd.DataFrame(Xh_temp[np.repeat(Xh_temp.columns, p+1)].apply(lambda x: np.power(x, range(p+1)), axis=1))
		
		
		S_hat = np.matmul(Xh_p_temp.transpose(),Xh_p_temp.mul(Kh_temp.mul(Pweights[index_temp], axis=0), axis=0))/n
		
		try:
			S_hat_inv = inv(S_hat)
		except np.linalg.LinAlgError:
			continue
			
		#G matrix
		if v==0:
			G_hat = np.matmul(Xh_p_temp, Xh_p_temp.mul(Kh_temp.mul(Pweights[index_temp], axis=0)**2)/n)
		else:
			if massPoints==True:
				Y_temp = pd.DataFrame(Fn[indexUnique])
				Xh_temp = pd.DataFrame(dataUnique-grid[j])/h1
				Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+1)]
				Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+1)), axis=1)

				if kernel == "triangular":
					Kh_temp = (1-abs(Xh_temp))/h1
					Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
				elif kernel == "uniform":
					Kh_temp = index_temp[indexUnique] * 0.5/h1
				else:
					Kh_temp = 0.75*(1-np.power(Xh_temp,2))/h1
					Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
				
				Xh_p_Kh_temp = Xh_p_temp.mul(Kh_temp, axis=0)
				Xh_p_Kh_Pweights_temp = Xh_p_Kh_temp.mul(PweightsUnique, axis=0)
				
				F_Xh_p_Kh_temp = np.matmul(Y_temp.transpose(), Xh_p_Kh_Pweights_temp)/n

				G = np.zeros((n,  len(Xh_p_Kh_Pweights_temp.columns)))
				for jj in range(len(Xh_p_Kh_Pweights_temp.columns)):
					G[:,jj] = (np.repeat((np.cumsum(Xh_p_Kh_Pweights_temp.iloc[::-1, jj])/n)[::-1], freqUnique) - F_Xh_p_Kh_temp.loc[0, jj]) * weights_normal
				
				G = pd.DataFrame(G)
				G_hat = np.matmul(G.transpose(), G)/n
				
			else:
				Y_temp = pd.DataFrame(Fn)
				Xh_temp = pd.DataFrame(data-grid[j])/h1
				Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+1)]
				Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+1)), axis=1)
				Kh_temp = __kernel_fn(Xh_temp, h1, kernel)
				
				Xh_p_Kh_temp = Xh_p_temp.mul(Kh_temp, axis=0)
				Xh_p_Kh_Pweights_temp = Xh_p_Kh_temp.mul(Pweights, axis=0)
				
				F_Xh_p_Kh_temp = np.matmul(Y_temp.transpose(), Xh_p_Kh_Pweights_temp)/n
				
				
				G = np.zeros((n,  len(Xh_p_Kh_Pweights_temp.columns)))
				for jj in range(len(Xh_p_Kh_Pweights_temp.columns)):
					G[:,jj] = (np.cumsum(Xh_p_Kh_Pweights_temp.iloc[:, jj]/n) - F_Xh_p_Kh_temp.loc[0, jj]) * weights_normal
				
				G = pd.DataFrame(G)
				G_hat = (np.matmul(G.transpose(), G))/n
				
		#Constants
		const_hat[j, 0] = math.factorial(v) * (np.matmul(S_hat_inv, C_p_hat)).loc[v, 0]
		const_hat[j, 1] = math.factorial(v) * (np.matmul(S_hat_inv, C_p1_hat)).loc[v, 0]
		
		if v>0:
			const_hat[j, 2] = math.factorial(v) * np.sqrt(abs(np.matmul(S_hat_inv, np.matmul(G_hat, S_hat_inv)).loc[v,v]) / (n*h1))
		else:
			temp_ii = min(max(mean(data <= grid[j]), 1/n), 1 - 1/n )
			const_hat[j, 2] = math.factorial(v) * np.sqrt(abs(np.matmul(S_hat_inv, np.matmul(G_hat, S_hat_inv)).loc[v,v] / (0.5*n**2)*h1 * temp_ii*(1-temp_ii)))
	

	#optimal bandwidth
	if v>0:
		opt_f = lambda a: np.power(a, 2*p+2-2*v) * sum((dgp_hat[:, 0]* const_hat[:, 0] + a * dgp_hat[:, 1]* const_hat[:, 1])**2) + sum(const_hat[:,2]**2) / np.power(a, 2*v - 1)
	else:
		opt_f = lambda a: np.power(a, 2*p+2-2*v) * sum((dgp_hat[:, 0]* const_hat[:, 0] + a * dgp_hat[:, 1]* const_hat[:, 1])**2) + sum(const_hat[:,2]**2) / a
			
	h = optimize.minimize_scalar(opt_f, bounds=(0, max(data)-min(data)), method='bounded').x
	
	
	if np.isnan(h)==True:
		h = max(h, np.sort(abs(data-grid[j]))[min(n, max(nLocalMin, 20+p+1))])
	
	if regularize==True:
		for j in range(ng):
			if nLocalMin>0:
				h = max(h, np.sort(abs(data-grid[j]))[min(n-1, nLocalMin-1)])
			if nUniqueMin>0:
				h = max(h, sorted(abs(dataUnique-grid[j]))[min(nUnique-1, nUniqueMin-1)])

	h = min(h, max(abs(max(dataUnique)-min(grid)), abs(min(dataUnique)-max(grid))))
			
	if stdVar==True:
		h = h*scale_temp
		
	return(h)
