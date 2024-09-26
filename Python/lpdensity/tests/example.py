import numpy as np
import random
import pandas as pd
from lpdensity import lpbwdensity
from lpdensity import lpdensity
import statistics as stat
import plotnine as pn

jtpa = pd.read_csv("jtpa.csv")
X = jtpa[
    [
        "income",
        "hsorged",
        "male",
        "nonwhite",
        "married",
        "wkless13",
        "afdc",
        "age2225",
        "age2629",
        "age3035",
        "age3644",
        "age4554",
    ]
]

n = len(jtpa["income"])
ones_df = pd.DataFrame(np.repeat(1, n))
X = pd.concat([X, ones_df], axis=1)

jtpa_notreat = jtpa.loc[jtpa["treatment"] == 0]
grid = np.linspace(2, 5, 10)

# uniquedata = lpdensityUnique(jtpa_notreat['logincome'])
# print(len(uniquedata['unique']))
# print(len(uniquedata['index']))
# print(len(uniquedata['freq']))

# n = len(data)
# grid = np.linspace(min(data), max(data), 19)
p = 2
q = p + 1
v = 1
kernel = "uniform"
Cweights = np.repeat(1, n)
Pweights = np.repeat(1, n)
massPoints = True
stdVar = True
regularize = True
nLocalMin = 20 + p + 1
nUniqueMin = 20 + p + 1
alpha = 0.05

est = lpdensity(
    data=jtpa_notreat["logincome"], bwselect="imse-dpi", grid=grid, kernel="triangular"
)
est2 = lpdensity(
    data=jtpa_notreat["logincome"],
    bwselect="imse-dpi",
    grid=grid,
    Cweights=jtpa_notreat["hsorged"],
)
est3 = lpdensity(
    data=jtpa_notreat["logincome"],
    bwselect="imse-dpi",
    grid=grid,
    Cweights=jtpa_notreat["hsorged"],
)
print(type(est3))
print(est3.plot(CIuniform=True))
