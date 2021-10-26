#--------------------------------------------------------------------------------
# lpdensity: Local Polynomial Density Estimation and Inference
# Matias D. Cattaneo, Rajita Chandak, Michael Jansson, and Xinwei Ma
# Replication Code
#--------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from lpdensity import lpdensity
from lpdensity import lpbwdensity
from scipy.stats import norm
import plotnine as plt

#----------------------------------------
# Generate data
#----------------------------------------
data = pd.DataFrame(data=np.random.normal(-1, 1, 4000))
data = data[data<0].dropna()
data = -1*data[0:2000].reset_index(drop=True)
data.columns = ['v1']

#----------------------------------------
# Figure 1
#----------------------------------------
data['pdf'] = norm.pdf(data, 1, 1)/norm.sf(0, 1, 1)
hist_plot = plt.ggplot(data, plt.aes(x='v1', y=plt.after_stat('density')))+ plt.geom_histogram(alpha=0.6, fill='red') + plt.theme_bw() + plt.geom_line(plt.aes(x='v1', y='pdf'))

model2 = lpdensity(data['v1'], bw=[0.5], grid = np.linspace(0, 4, num=int(4/0.5)))
model2.plot()

#-------------------------------------------------------------------
# lpdensity(): Estimation with bandwidth 0.5 on provided grid points
#-------------------------------------------------------------------
model1 = lpdensity(data['v1'], bw=[0.5], grid = np.linspace(0, 4, num=int(4/0.5)))
model1.plot()
print(repr(model1))

#----------------------------------------
# lpdensity(): extracting estimation
#   results
#----------------------------------------
print(model1.Estimate)

#----------------------------------------
# lpdensity(): conventional inference
#----------------------------------------
print(repr(lpdensity(data=data['v1'], bw=[0.5], p=2, q=2)))

model1.confint(alpha=0.01, CIuniform=True)

#----------------------------------------
# lpdensity(): inconsistent density
#   estimation using partial sample
#----------------------------------------
lpdensity(data['v1'][data['v1'] < 1.5].reset_index(drop=True), bw=[0.5], grid=[1.5]).Estimate['f_p']
lpdensity(data['v1'][data['v1'] > 1.5].reset_index(drop=True), bw=[0.5], grid=[1.5]).Estimate['f_p']
norm.pdf(1.5, 1, 1)/norm.sf(0, 1, 1)

#----------------------------------------
# lpdensity(): consistent density
#   estimation using partial sample and
#   option "scale"
#----------------------------------------
lpdensity(data['v1'][data['v1'] < 1.5].reset_index(drop=True), bw=[0.5], grid=[1.5], scale = [sum(data['v1'] < 1.5)/2000]).Estimate['f_p']
lpdensity(data['v1'][data['v1'] > 1.5].reset_index(drop=True), bw=[0.5], grid=[1.5], scale = [sum(data['v1'] > 1.5)/2000]).Estimate['f_p']

#----------------------------------------
# plot(): customization
#----------------------------------------
model2.plot(CItype='line')
model2.plot(type='points', CItype='ebar')
model2.plot(alpha=0.1, CIuniform=True)

#----------------------------------------
# lpbwdensity(): illustration
#----------------------------------------
model1bw = lpbwdensity(data['v1'], grid = np.linspace(0, 4, num=int(4/0.5)))
print(repr(model1bw))

#----------------------------------------
# lpdensity(): automatic bandwidth
#   selection
#----------------------------------------
model5 = lpdensity(data['v1'], grid = np.linspace(0, 4, num=int(4/0.5)), bwselect='imse-dpi')
print(repr(model5))

#----------------------------------------
# lpdensity(): undersmoothing
#----------------------------------------

# Estimation and plot using IMSE bandwidth
model6bwIMSE = lpbwdensity(data['v1'], grid = np.linspace(0, 4, num=int(4/0.05)), bwselect='imse-dpi')
model6 = lpdensity(data['v1'], grid = np.linspace(0, 4, num=int(4/0.05)), bwselect='imse-dpi')
model6.plot()

# Estimation and plot using half IMSE bandwidth
model7 = lpdensity(data['v1'], grid = np.linspace(0, 4, num=int(4/0.5)), bw=[model6bwIMSE.bws[0]/2])
model7.plot()
