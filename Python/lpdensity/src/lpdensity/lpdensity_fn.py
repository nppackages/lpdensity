#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ignore future warnings from pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np
import pandas as pd
import math
from numpy.linalg import inv
from . import lpbwdensity_fn

def __lpdensity_fn(data, grid, bw, p, q, v, kernel, Cweights, Pweights, massPoints, showSE=True):
	# data prep
	ii = np.argsort(data)
	data = data[ii]
	Cweights = Cweights[ii]
	Pweights = Pweights[ii]
	n = len(data)
	ng = len(grid)
	
	dataUnique = lpbwdensity_fn.__lpdensityUnique(data)
	freqUnique = dataUnique["freq"]
	indexUnique = dataUnique["index"]
	dataUnique = dataUnique["unique"]
	nUnique = len(dataUnique)

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
	
	# initialization
	hat_p = np.repeat(float("NaN"), ng)
	hat_q = np.repeat(float("NaN"), ng)
	se_p = np.repeat(float("NaN"), ng)
	se_q = np.repeat(float("NaN"), ng)
	iff_p = np.empty((n, ng))
	iff_p[:] = np.NaN
	iff_q = np.empty((n, ng))
	iff_q[:] = np.NaN
	nh = np.repeat(float("NaN"), ng)
	nhu = np.repeat(float("NaN"), ng)
	
	# looping over grid points
	for j in range(ng) :
		index_temp = abs(data-grid[j]) <= bw[j]
		nh[j] = sum(index_temp)
		nhu[j] = sum(index_temp[indexUnique])
			
		# LHS and RHS variables
		if massPoints == True:
			Y_temp = pd.DataFrame(Fn[indexUnique])
			Xh_temp = pd.DataFrame(data[indexUnique]-grid[j])/bw[j]
			Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+1)]
			Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+1)), axis=1)
			
			# weights
			if kernel=="triangular":
				Kh_temp = (1-abs(Xh_temp))/bw[j]
				Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
			elif kernel=="uniform":
				Kh_temp = (0.5/bw[j]) * index_temp[indexUnique]
			else :
				Kh_temp = 0.75*(1-np.power(Xh_temp,2))/bw[j]
				Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
			
			Xh_p_Kh_temp = Xh_p_temp.mul(Kh_temp, axis=0)
			Xh_p_Kh_Pweights_temp = Xh_p_Kh_temp.mul(PweightsUnique, axis=0)
			
			try:
				XhKhXh_inv = inv(np.matmul(Xh_p_temp[index_temp[indexUnique]].transpose(), Xh_p_Kh_Pweights_temp[index_temp[indexUnique]])/ n)
			except np.linalg.LinAlgError:
				continue
			
			#point estimate
			hat_p[j] = math.factorial(v) * np.matmul(Y_temp.transpose(), np.matmul(Xh_p_Kh_Pweights_temp, XhKhXh_inv))[v]/(np.power(bw[j],v)) /n
		else:
			Y_temp = pd.DataFrame(Fn)
			Xh_temp = pd.DataFrame(data-grid[j])/bw[j]
			Xh_p_temp = Xh_temp[np.repeat(Xh_temp.columns, p+1)]
			Xh_p_temp = Xh_p_temp.apply(lambda x: np.power(x, range(p+1)), axis=1)
			
			# weights
			if kernel=="triangular":
				Kh_temp = ((1-abs(Xh_temp))/bw[j]) * index_temp
			elif kernel=="uniform":
				Kh_temp = (0.5/bw[j]) * index_temp
			else :
				Kh_temp = (0.75*(1-np.power(Xh_temp, 2))/bw[j]) * index_temp
				
			Xh_p_Kh_temp = Xh_p_temp.mul(Kh_temp, axis=0)
			Xh_p_Kh_Pweights_temp = Xh_p_temp.mul(PweightsUnique, axis=0)
			
			try:
				XhKhXh_inv = inv(np.matmul(Xh_p_temp[index_temp].transpose(), Xh_p_Kh_Pweights_temp[index_temp])/ n)
			except np.linalg.LinAlgError:
				continue
			
			#point estimate
			hat_p[j] = math.factorial(v) * np.matmul(XhKhXh_inv, np.matmul(Xh_p_Kh_Pweights_temp[index_temp,:].transpose(), Y_temp[index_temp]))[v]/(np.power(bw[j], v)) /n
			
		if showSE==True:
			if massPoints==True:
				F_Xh_p_Kh_temp = np.matmul(Y_temp.transpose(), Xh_p_Kh_Pweights_temp)/n
				
				G = np.zeros((n,  len(Xh_p_Kh_Pweights_temp.columns)))
				for jj in range(len(Xh_p_Kh_Pweights_temp.columns)):
					G[:,jj] = (np.repeat((np.cumsum(Xh_p_Kh_Pweights_temp.iloc[::-1, jj])/n)[::-1], freqUnique) - F_Xh_p_Kh_temp.loc[0, jj]) * weights_normal
				
				G = pd.DataFrame(G)
				
				iff_p[:,j] = np.matmul(XhKhXh_inv, G.transpose()).iloc[v, :] * math.factorial(v)/np.sqrt(n*np.power(bw[j], 2*v))
			else:
				F_Xh_p_Kh_temp = np.matmul(Y_temp.transpose(), Xh_p_Kh_Pweights_temp)/n
				
				G = np.zeros((n,  len(Xh_p_Kh_Pweights_temp.columns)))
				for jj in range(len(Xh_p_Kh_Pweights_temp.columns)):
					G[:,jj] = (np.repeat(np.cumsum(Xh_p_Kh_Pweights_temp.iloc[:, jj]/n), freqUnique) - F_Xh_p_Kh_temp.loc[0, jj]) * weights_normal
				
				G = pd.DataFrame(G)
				
				iff_p[:,j] = np.matmul(XhKhXh_inv, G.transpose()).iloc[v, :] * math.factorial(v)/np.sqrt(n*np.power(bw[j], 2*v))
				
		if q>p:
			if massPoints == True:
				Y_temp = pd.DataFrame(Fn[indexUnique])
				Xh_temp = pd.DataFrame(data[indexUnique]-grid[j])/bw[j]
				Xh_q_temp = Xh_temp[np.repeat(Xh_temp.columns, q+1)]
				Xh_q_temp = Xh_q_temp.apply(lambda x: np.power(x, range(q+1)), axis=1)
			
				# weights
				if kernel=="triangular":
					Kh_temp = (1-abs(Xh_temp))/bw[j]
					Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
				elif kernel=="uniform":
					Kh_temp = (0.5/bw[j]) * index_temp[indexUnique]
				else :
					Kh_temp = 0.75*(1-np.power(Xh_temp,2))/bw[j]
					Kh_temp = Kh_temp.mul(index_temp[indexUnique], axis=0)
				
				Xh_q_Kh_temp = Xh_q_temp.mul(Kh_temp, axis=0)
				Xh_q_Kh_Pweights_temp = Xh_q_Kh_temp.mul(PweightsUnique, axis=0)
				
				try:
					XhKhXh_inv = inv(np.matmul(Xh_q_temp[index_temp[indexUnique]].transpose(), Xh_q_Kh_Pweights_temp[index_temp[indexUnique]])/ n)
				except np.linalg.LinAlgError:
					continue
				
				#point estimate
				hat_q[j] = math.factorial(v) * np.matmul(Y_temp.transpose(), np.matmul(Xh_q_Kh_Pweights_temp, XhKhXh_inv))[v]/(np.power(bw[j],v)) /n
			else:
				Y_temp = pd.DataFrame(Fn)
				Xh_temp = pd.DataFrame(data-grid[j])/bw[j]
				Xh_q_temp = Xh_temp[np.repeat(Xh_temp.columns, q+1)]
				Xh_q_temp = Xh_q_temp.apply(lambda x: np.power(x, range(q+1)), axis=1)
				
				# weights
				if kernel=="triangular":
					Kh_temp = ((1-abs(Xh_temp))/bw[j]) * index_temp
				elif kernel=="uniform":
					Kh_temp = (0.5/bw[j]) * index_temp
				else :
					Kh_temp = (0.75*(1-np.power(Xh_temp, 2))/bw[j]) * index_temp
					
				Xh_q_Kh_temp = Xh_q_temp.mul(Kh_temp, axis=0)
				Xh_q_Kh_Pweights_temp = Xh_q_temp.mul(PweightsUnique, axis=0)
				
				try:
					XhKhXh_inv = inv(np.matmul(Xh_q_temp[index_temp].transpose(), Xh_q_Kh_Pweights_temp[index_temp])/ n)
				except np.linalg.LinAlgError:
					continue
				
				#point estimate
				hat_q[j] = math.factorial(v) * np.matmul(XhKhXh_inv, np.matmul(Xh_q_Kh_Pweights_temp[index_temp,:].transpose(), Y_temp[index_temp]))[v]/(np.power(bw[j], v)) /n
				
			if showSE==True:
				if massPoints==True:
					F_Xh_q_Kh_temp = np.matmul(Y_temp.transpose(), Xh_q_Kh_Pweights_temp)/n
					
					G = np.zeros((n,  len(Xh_q_Kh_Pweights_temp.columns)))
					for jj in range(len(Xh_q_Kh_Pweights_temp.columns)):
						G[:,jj] = (np.repeat((np.cumsum(Xh_q_Kh_Pweights_temp.iloc[::-1, jj])/n)[::-1], freqUnique) - F_Xh_q_Kh_temp.loc[0, jj]) * weights_normal
					
					G = pd.DataFrame(G)
					
					iff_q[:,j] = np.matmul(XhKhXh_inv, G.transpose()).iloc[v, :] * math.factorial(v)/np.sqrt(n*np.power(bw[j], 2*v))
				else:
					F_Xh_q_Kh_temp = np.matmul(Y_temp.transpose(), Xh_q_Kh_Pweights_temp)/n
					
					G = np.zeros((n,  len(Xh_q_Kh_Pweights_temp.columns)))
					for jj in range(len(Xh_q_Kh_Pweights_temp.columns)):
						G[:,jj] = (np.repeat(np.cumsum(Xh_q_Kh_Pweights_temp.iloc[:, jj]/n), freqUnique) - F_Xh_q_Kh_temp.loc[0, jj]) * weights_normal
					
					G = pd.DataFrame(G)
					
					iff_q[:,j] = np.matmul(XhKhXh_inv, G.transpose()).iloc[v, :] * math.factorial(v)/np.sqrt(n*np.power(bw[j], 2*v))
					

	# putting together all the data in the right dataframes
	
	CovMat_p = pd.DataFrame(data=np.matmul(iff_p.transpose(), iff_p)/n)
	CovMat_q = pd.DataFrame(data=np.matmul(iff_q.transpose(), iff_q)/n)

	se_p = np.sqrt(abs(np.diag(CovMat_p)))
	se_q = np.sqrt(abs(np.diag(CovMat_q)))
	
	Estimate = np.empty((ng, 8))
	Estimate[:,0] = grid
	Estimate[:,1] = bw
	Estimate[:,2] = nh
	Estimate[:,3] = nhu
	Estimate[:,4] = hat_p
	Estimate[:,5] = hat_q
	Estimate[:,6] = se_p
	Estimate[:,7] = se_q
	
	est_data = pd.DataFrame(data=Estimate, columns=["grid", "bw", "nh", "nhu", "f_p", "f_q", "se_p", "se_q"])
	
	return({'Estimate': est_data, 'CovMat_p': CovMat_p, 'CovMat_q': CovMat_q})
