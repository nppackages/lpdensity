#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import math

import numpy as np
import pandas as pd
from numpy.linalg import inv

from . import lpbwdensity_fn


def _poly_matrix(x, degree):
    x = np.asarray(x, dtype=float)
    return np.vander(x, N=degree + 1, increasing=True)


def _kernel_values(xh, h, kernel):
    xh = np.asarray(xh, dtype=float)
    if kernel == "triangular":
        return (1 - np.abs(xh)) / h
    if kernel == "uniform":
        return np.full(len(xh), 0.5 / h)
    return 0.75 * (1 - np.power(xh, 2)) / h


def _group_increments(values, index):
    cumulative = np.cumsum(values)[index]
    if len(cumulative) > 1:
        cumulative = cumulative - np.insert(cumulative[:-1], 0, 0)
    return cumulative


def _local_fit(
    data,
    grid_point,
    bandwidth,
    degree,
    derivative,
    kernel,
    n,
    mass_points,
    index_temp,
    index_unique,
    freq_unique,
    fn,
    weights_normal,
    pweights,
    pweights_unique,
):
    if mass_points:
        y = np.asarray(fn[index_unique], dtype=float)
        xh = (data[index_unique] - grid_point) / bandwidth
        weights = pweights_unique
        mask = index_temp[index_unique]
    else:
        y = np.asarray(fn, dtype=float)
        xh = (data - grid_point) / bandwidth
        weights = pweights
        mask = index_temp

    design = _poly_matrix(xh, degree)
    kernel_weights = _kernel_values(xh, bandwidth, kernel) * mask
    design_kh_weights = design * kernel_weights[:, None] * weights[:, None]

    try:
        gram_inv = inv((design[mask].T @ design_kh_weights[mask]) / n)
    except np.linalg.LinAlgError:
        return np.nan, None

    estimate = (
        math.factorial(derivative)
        * ((y @ design_kh_weights) @ gram_inv)[derivative]
        / np.power(bandwidth, derivative)
        / n
    )

    moment_center = (y @ design_kh_weights) / n
    influence_raw = np.zeros((n, degree + 1))
    for column in range(degree + 1):
        if mass_points:
            cumulative = (np.cumsum(design_kh_weights[::-1, column]) / n)[::-1]
            expanded = np.repeat(cumulative, freq_unique)
        else:
            expanded = np.cumsum(design_kh_weights[:, column] / n)
        influence_raw[:, column] = (expanded - moment_center[column]) * weights_normal

    influence = (
        (gram_inv @ influence_raw.T)[derivative, :]
        * math.factorial(derivative)
        / np.sqrt(n * np.power(bandwidth, 2 * derivative))
    )
    return estimate, influence


def __lpdensity_fn(
    data, grid, bw, p, q, v, kernel, Cweights, Pweights, massPoints, showSE=True
):
    data = np.asarray(data, dtype=float).reshape(-1)
    grid = np.asarray(grid, dtype=float).reshape(-1)
    bw = np.asarray(bw, dtype=float).reshape(-1)
    Cweights = np.asarray(Cweights, dtype=float).reshape(-1)
    Pweights = np.asarray(Pweights, dtype=float).reshape(-1)

    ii = np.argsort(data)
    data = data[ii]
    Cweights = Cweights[ii]
    Pweights = Pweights[ii]
    n = len(data)
    ng = len(grid)

    dataUnique = lpbwdensity_fn.__lpdensityUnique(data)
    freqUnique = np.asarray(dataUnique["freq"], dtype=int)
    indexUnique = np.asarray(dataUnique["index"], dtype=int)
    dataUnique = np.asarray(dataUnique["unique"], dtype=float).reshape(-1)

    if massPoints == True:
        fn_unique = (np.cumsum(Cweights * Pweights) / sum(Cweights * Pweights))[
            indexUnique
        ]
        Fn = np.repeat(fn_unique, repeats=freqUnique)
    else:
        Fn = np.cumsum(Cweights * Pweights) / sum(Cweights * Pweights)

    weights_normal = Cweights * Pweights / sum(Cweights * Pweights) * n
    Cweights = Cweights / sum(Cweights) * n
    Pweights = Pweights / sum(Pweights) * n

    weights_normalUnique = _group_increments(weights_normal, indexUnique)
    CweightsUnique = _group_increments(Cweights, indexUnique)
    PweightsUnique = _group_increments(Pweights, indexUnique)

    hat_p = np.repeat(float("NaN"), ng)
    hat_q = np.repeat(float("NaN"), ng)
    iff_p = np.empty((n, ng))
    iff_p[:] = np.nan
    iff_q = np.empty((n, ng))
    iff_q[:] = np.nan
    nh = np.repeat(float("NaN"), ng)
    nhu = np.repeat(float("NaN"), ng)

    for j in range(ng):
        index_temp = abs(data - grid[j]) <= bw[j]
        nh[j] = sum(index_temp)
        nhu[j] = sum(index_temp[indexUnique])

        hat_p[j], influence_p = _local_fit(
            data=data,
            grid_point=grid[j],
            bandwidth=bw[j],
            degree=p,
            derivative=v,
            kernel=kernel,
            n=n,
            mass_points=massPoints,
            index_temp=index_temp,
            index_unique=indexUnique,
            freq_unique=freqUnique,
            fn=Fn,
            weights_normal=weights_normal,
            pweights=Pweights,
            pweights_unique=PweightsUnique,
        )
        if showSE == True and influence_p is not None:
            iff_p[:, j] = influence_p

        if q > p:
            hat_q[j], influence_q = _local_fit(
                data=data,
                grid_point=grid[j],
                bandwidth=bw[j],
                degree=q,
                derivative=v,
                kernel=kernel,
                n=n,
                mass_points=massPoints,
                index_temp=index_temp,
                index_unique=indexUnique,
                freq_unique=freqUnique,
                fn=Fn,
                weights_normal=weights_normal,
                pweights=Pweights,
                pweights_unique=PweightsUnique,
            )
            if showSE == True and influence_q is not None:
                iff_q[:, j] = influence_q

    CovMat_p = pd.DataFrame(data=np.matmul(iff_p.transpose(), iff_p) / n)
    CovMat_q = pd.DataFrame(data=np.matmul(iff_q.transpose(), iff_q) / n)

    se_p = np.sqrt(abs(np.diag(CovMat_p)))
    se_q = np.sqrt(abs(np.diag(CovMat_q)))

    Estimate = np.empty((ng, 8))
    Estimate[:, 0] = grid
    Estimate[:, 1] = bw
    Estimate[:, 2] = nh
    Estimate[:, 3] = nhu
    Estimate[:, 4] = hat_p
    Estimate[:, 5] = hat_q
    Estimate[:, 6] = se_p
    Estimate[:, 7] = se_q

    est_data = pd.DataFrame(
        data=Estimate, columns=["grid", "bw", "nh", "nhu", "f_p", "f_q", "se_p", "se_q"]
    )

    return {"Estimate": est_data, "CovMat_p": CovMat_p, "CovMat_q": CovMat_q}
