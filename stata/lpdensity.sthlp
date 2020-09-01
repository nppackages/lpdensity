{smcl}
{* *! version 2.2 2020-09-01}{...}

{title:Title}

{p 4 8}{cmd:lpdensity} {hline 2} Local Polynomial Density Estimation and Inference.{p_end}


{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:lpdensity} {it:Var} {ifin} 
[{cmd:,} {p_end}
{p 14 18}{cmd:grid(}{it:Var}{cmd:)} 
{cmd:bw(}{it:Var} or {it: #}{cmd:)} 
{cmd:p(}{it:#}{cmd:)}
{cmd:q(}{it:#}{cmd:)}
{cmd:v(}{it:#}{cmd:)}
{cmd:kernel(}{it:KernelFn}{cmd:)}
{cmd:scale(}{it:#}{cmd:)}
{cmd:nomasspoints}{p_end}
{p 14 18}{cmd:bwselect(}{it:BwMethod}{cmd:)}
{cmd:nlocalmin(}{it:#}{cmd:)}
{cmd: nuniquemin(}{it:#}{cmd:)}
{cmd:noregularize}
{cmd:nostdvar}
{p_end}
{p 14 18}{cmd:cweights(}{it:Var}{cmd:)}
{cmd:pweights(}{it:Var}{cmd:)}{p_end}
{p 14 18}{cmd:genvars(}{it:NewVarName}{cmd:)}{p_end}
{p 14 18}{cmd:rgrid(}{it:Var}{cmd:)} 
{cmd:rindex(}{it:Var}{cmd:)} 
{cmd:level(}{it:#}{cmd:)}
{cmd:ciuniform}
{cmd:cisimul(}{it:#}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}{p_end}
{p 14 18}{cmd:plot}{p_end}
{p 14 18}{cmd:estype(}{it:ESOpts}{cmd:)}
{cmd:esline_options(}{it:ESLineOpts}{cmd:)}
{cmd:espoint_options(}{it:ESPointOpts}{cmd:)}{p_end}
{p 14 18}{cmd:citype(}{it:CIOpts}{cmd:)}
{cmd:ciregion_options(}{it:CIRegionOpts}{cmd:)}
{cmd:ciline_options(}{it:CILineOpts}{cmd:)}
{cmd:ciebar_options(}{it:CIEbarOpts}{cmd:)}{p_end}
{p 14 18}{cmd:histogram}
{cmd:hiplot_options(}{it:HistOpts}{cmd:)}{p_end}
{p 14 18}{cmd:graph_options(}{it:GraphOpts}{cmd:)}
{p_end}
{p 14 18}]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:lpdensity} implements the local polynomial regression based density (and derivatives) estimator proposed in 
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf":Cattaneo, Jansson and Ma (2020a)}.
Robust bias-corrected inference, both pointwise (confidence intervals) and uniform (confidence bands) are also implemented following the results in
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf":Cattaneo, Jansson and Ma (2020a)}
and
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JoE.pdf":Cattaneo, Jansson and Ma (2020b)}. 
See {browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JSS.pdf":Cattaneo, Jansson and Ma (2020c)} for more 
implementation details and illustrations.{p_end}

{p 8 8} Companion command: {help lpbwdensity:lpbwdensity} for bandwidth selection.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} functions are also available {browse "https://nppackages.github.io/lpdensity":here}.{p_end}

{p 4 8} Related Stata and R packages are available in the following website:{p_end}

{p 8 8} {browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimation}

{p 4 8}{opt gri:d}({it:var}) specifies the grid on which density is estimated. When set to default, grid points will be chosen as 0.05-0.95
percentiles of the data, with 0.05 step size.{p_end}

{p 4 8}{opt bw}({it:var} or {it:#}) specifies the bandwidth (either a variable containing bandwidth for each grid point or a single number) used for estimation. When omitted, bandwidth will be computed by method specified 
in {cmd:bwselect(}{it:BwMethod}{cmd:)}.{p_end}

{p 4 8}{opt p}({it:#}) specifies the local polynomial order for constructing point estimates.
Default is {cmd:p(2)} (local quadratic regression).{p_end}

{p 4 8}{opt q}({it:#}) specifies the local polynomial order for constructing
confidence intervals/bands (a.k.a. the bias correction order). Default is {cmd:p(}{it:#}{cmd:)+1}. When specified
the same as {cmd:p(}{it:#}{cmd:)}, no bias correction will be performed. Otherwise it should be
strictly larger than {cmd:p(}{it:#}{cmd:)}. {p_end}

{p 4 8}{opt v}({it:#}) specifies the derivative of distribution function to be estimated. {cmd:v(0)} for
the distribution function, {cmd:v(1)} (default) for the density funtion, etc.{p_end}

{p 4 8}{opt ker:nel}({it:KernelFn}) specifies the kernel function used to construct the local-polynomial estimator(s). {p_end}
{p 8 12}{opt triangular}{bind:  } {it:K(u) = (1 - |u|) * (|u|<=1)}. This is the default option. {p_end}
{p 8 12}{opt epanechnikov}{bind:}  {it:K(u) = 0.75 * (1 - u^2) * (|u|<=1)}. {p_end}
{p 8 12}{opt uniform}{bind:     }  {it:K(u) = 0.5 * (|u|<=1)}. {p_end}

{p 4 8}{opt sca:le}({it:#}) controls how estimates are scaled. For example, setting this parameter to 0.5 will scale down both the
 point estimates and standard errors by half. Default is {cmd:scale(1)}. This parameter is useful when only
 a subsample is employed for estimation.{p_end}

{p 4 8}{opt nomass:points} will not adjust point estimates or standard errors even if there are mass points in the data.{p_end}

{dlgtab:Bandwidth Selection}

{p 4 8}{opt bws:elect}({it:BwMethod}) specifies method for data-driven bandwidth selection. This option will be
ignored if {cmd:bw(}{it:Var}{cmd:)} is provided.{p_end}
{p 8 12}{opt mse-dpi}{bind: } mean squared error optimal bandwidth for each grid point. This is the default option.{p_end}
{p 8 12}{opt imse-dpi}{bind:} integrated mean squared error optimal bandwidth which is common for all grid points.{p_end}
{p 8 12}{opt mse-rot}{bind: } rule-of-thumb bandwidth based on a Gaussian reference model.{p_end}
{p 8 12}{opt imse-rot}{bind:} integrated rule-of-thumb bandwidth based on a Gaussian reference model which is common for all grid points.{p_end}

{p 4 8}{opt nloc:almin}({it:#}) specifies the minimum number of observations in each local neighborhood. This option will be ignored if set to 0, or if {cmd:noregularize} is used. The default value is {cmd:20+p(}{it:#}{cmd:)+1}. {p_end}

{p 4 8}{opt nuni:quemin}({it:#}) specifies the minimum number of unique observations in each local neighborhood. This option will be ignored if set to 0, or if {cmd:noregularize} is used. The default value is {cmd:20+p(}{it:#}{cmd:)+1}. {p_end}

{p 4 8}{opt noreg:ularize} suppresses local sample size checking.{p_end}

{p 4 8}{opt nostdv:ar} will not standardize the data for bandwidth selection. Note that this may lead to unstable performance of the numerical optimization procedure. {p_end}

{dlgtab:Weights}

{p 4 8}{opt cw:eights}({it:Var}) specifies weights used for counterfactual distribution construction.{p_end}
 
{p 4 8}{opt pw:eights}({it:Var}) specifies weights used in sampling. Should be nonnegative.{p_end}

{dlgtab:Storing and displaying results}

{p 4 8}{opt gen:vars}({it:NewVarName}) specifies if new varaibles should be generated to store estimation results. If {it:NewVarName} is provided, the following new varaibles will be
generated: {p_end}
{p 8 12}{it:NewVarName_grid}{bind:} grid points, {p_end}
{p 8 12}{it:NewVarName_bw}{bind:  } bandwidth, {p_end}
{p 8 12}{it:NewVarName_nh}{bind:  } local/effective sample sizes, {p_end}
{p 8 12}{it:NewVarName_f_p} and {it:NewVarName_se_p}{bind:}
point estimates with polynomial order {cmd:p(}{it:#}{cmd:)} and the corresponding standard errors, {p_end}
{p 8 12}{it:NewVarName_f_q} and {it:NewVarName_se_q}{bind:}
point estimates with polynomial order {cmd:q(}{it:#}{cmd:)} and the corresponding standard errors, only available if different from {cmd:p(}{it:#}{cmd:)},{p_end}
{p 8 12}{it:NewVarName_CI_l} and {it:NewVarName_CI_r}{bind:} confidence intervals/bands.{p_end}

{p 4 8}{opt rgri:d}({it:var}) specifies a set of grid points to display the results. When omitted, this will be the same as {cmd:grid(}{it:Var}{cmd:)}. 

{p 4 8}{opt rind:ex}({it:var}) specifies a set of indices to display the results. This option will be ignored if {cmd:rgrid(}{it:Var}{cmd:)} is provided. 

{p 4 8}{opt l:evel}({it:#}) controls the level of the confidence interval, and should be between 0 and 100. Default is {cmd:level(95)}.{p_end}

{p 4 8}{opt ciu:iform} computes a uniform confidence band instead of pointwise confidence intervals. {p_end}

{p 4 8}{opt cis:imul}({it:#}) specifies the number of simulations used to construct critical values. Default is {cmd:cisimul(2000)}. This option will be ignored unless {cmd:ciuniform} is provided.{p_end}

{p 4 8}{opt sep:arator}({it:#}) draw a seperation line after every {it:#} variables; default is {cmd:separator(5)}.{p_end}

{dlgtab:Plotting}

{p 4 8}{opt pl:ot} if specified, point estimates and confidence intervals will be plotted.{p_end}

{p 4 8}{opt esty:pe}({it:ESOpts}) specifies the plotting style of point estimates.{p_end}
{p 8 12}{opt line}{bind:  } a curve. This is the default option.{p_end}
{p 8 12}{opt points}{bind:} individual points.{p_end}
{p 8 12}{opt both}{bind:  } both of the above.{p_end}
{p 8 12}{opt none}{bind:  } will not plot point estimates.{p_end}

{p 4 8}{opt esli:ne_options}({it:ESlineOpts}){bind:  } specifies additional {cmd:twoway line}{bind:   } options for plotting point estimates. {p_end}

{p 4 8}{opt espo:int_options}({it:ESPointOpts}){bind:} specifies additional {cmd:twoway scatter}{bind:} options for plotting point estimates. {p_end}

{p 4 8}{opt city:pe}({it:CIOpts}) specifies the plotting style of confidence intervals/bands.{p_end}
{p 8 12}{opt region}{bind:} shaded region. This is the default option.{p_end}
{p 8 12}{opt line}{bind:  } upper and lower bounds.{p_end}
{p 8 12}{opt ebar}{bind:  } error bars.{p_end}
{p 8 12}{opt all}{bind:   } all of the above.{p_end}
{p 8 12}{opt none}{bind:  } will not plot confidence intervals/bands.{p_end}

{p 4 8}{opt cire:gion_options}({it:CIRegionOpts}){bind:} specifies additional {cmd:twoway rarea}{bind:}  options for plotting confidence intervals/regions. {p_end}

{p 4 8}{opt cili:ne_options}({it:CILineOpts}){bind:    } specifies additional {cmd:twoway rline}{bind:} options for plotting confidence intervals/regions. {p_end}

{p 4 8}{opt cieb:r_options}({it:CIEbarOpts}){bind:     } specifies additional {cmd:twoway rcap}{bind: } options for plotting confidence intervals/regions. {p_end}

{p 4 8}{opt hist:gram} if specified, a histogram will be included in the background.{p_end}

{p 4 8}{opt hipl:ot_options}({it:HistOpts}){bind:} specifies additional {cmd:twoway histogram}{bind:}  options for the histogram. {p_end}

{p 4 8}{opt gra:ph_options}({it:GraphOpts}) specifies additional options for plotting, such as legends and labels. {p_end}

{marker remarks}{...}
{title:Remarks}

{p 4 8}
Bias correction is only used for the construction of confidence intervals/bands, but not for point estimation. The point estimates, denoted by f_p, are constructed using local polynomial estimates of order
{cmd:p(}{it:#}{cmd:)},
while the centering of the confidence intervals/bands, denoted by f_q, are constructed using local polynomial estimates of order
{cmd:q(}{it:#}{cmd:)}.
The confidence intervals/bands take the form:
[f_q - cv * SE(f_q) , f_q + cv * SE(f_q)],
where cv denotes the appropriate critical value and SE(f_q) denotes an standard error estimate for the centering of the confidence interval/band.
As a result, the confidence intervals/bands may not be centered at the point estimates because they have been bias-corrected. Setting
{cmd:q(}{it:#}{cmd:)}
and
{cmd:p(}{it:#}{cmd:)}
to be equal results on centered at the point estimate confidence intervals/bands, but requires undersmoothing for valid inference (i.e., (I)MSE-optimal bandwdith for the density point estimator cannot be used). Hence the bandwidth would need to be specified manually when
{cmd:q(}{it:#}{cmd:)} = {cmd:p(}{it:#}{cmd:)},
and the point estimates will not be (I)MSE optimal. See Cattaneo, Jansson and Ma
({browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf":2020a},
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JoE.pdf":2020b})
for details, and also Calonico, Cattaneo, and Farrell
({browse "https://nppackages.github.io/Calonico-Cattaneo-Farrell_2018_JASA.pdf":2018},
{browse "https://nppackages.github.io/Calonico-Cattaneo-Farrell_2020_CEopt.pdf":2020}) 
for robust bias correction methods. {p_end}

{p 4 8}
Sometimes the density point estimates may lie outside of the confidence intervals/bands, which can happen if the underlying distribution exhibits high curvature at some evaluation point(s). One possible solution in this case is to increase the polynomial order {cmd:p(}{it:#}{cmd:)} or to employ a smaller bandwidth. {p_end}	
		
{marker examples}{...}
{title:Examples}

{p 4 8}Generate artifitial data:{p_end}
{p 8 8}{cmd:. set obs 2000}{p_end}
{p 8 8}{cmd:. set seed 42}{p_end}
{p 8 8}{cmd:. gen lpd_data = rnormal()}{p_end}

{p 4 8}Density estimation at empirical quantiles: {p_end}
{p 8 8}{cmd:. lpdensity lpd_data}{p_end}

{p 4 8}Density estimation at empirical quantiles with the IMSE-optimal bandwidth: {p_end}
{p 8 8}{cmd:. lpdensity lpd_data, bwselect(imse-dpi)}{p_end}

{p 4 8}Density estimation on a fixed grid (0.1, 0.2, ..., 1):{p_end}
{p 8 8}{cmd:. gen lpd_grid = _n / 10 if _n <= 10}{p_end}
{p 8 8}{cmd:. lpdensity lpd_data, grid(lpd_grid)}{p_end}

{p 4 8}Report uniform confidence bands (instead of pointwise confidence intervals): {p_end}
{p 8 8}{cmd:. lpdensity lpd_data, ciuniform}{p_end}
{p 8 8}{cmd:. lpdensity lpd_data, ciuniform level(99)}{p_end}

{p 4 8}Save estimation results to new variables:{p_end}
{p 8 8}{cmd:. capture drop temp_*}{p_end}
{p 8 8}{cmd:. lpdensity lpd_data, genvars(temp)}{p_end}

{p 4 8}Density plot:{p_end}
{p 8 8}{cmd:. lpdensity lpd_data, plot}{p_end}
{p 8 8}{cmd:. lpdensity lpd_data, plot histogram}{p_end}
{p 8 8}{cmd:. lpdensity lpd_data, plot histogram ciuniform level(90)}{p_end}

{marker saved_results}{...}
{title:Saved results}

{p 4 8}{cmd:lpdensity} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(p)}}option {cmd:p(}{it:#}{cmd:)}{p_end}
{synopt:{cmd:e(q)}}option {cmd:q(}{it:#}{cmd:)}{p_end}
{synopt:{cmd:e(v)}}option {cmd:v(}{it:#}{cmd:)}{p_end}
{synopt:{cmd:e(scale)}}option {cmd:scale(}{it:#}{cmd:)}{p_end}
{synopt:{cmd:e(level)}}option {cmd:level(}{it:#}{cmd:)}{p_end}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(bwselect)}}option {cmd:bwselect(}{it:BwMethod}{cmd:)}{p_end}
{synopt:{cmd:e(kernel)}}option {cmd:kernel(}{it:KernelFn}{cmd:)}{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(result)}}estimation result{p_end}

{marker references}{...}
{title:References}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018.
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf":On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}.{p_end}
{p 8 8}{it:Journal of the American Statistical Association} 113(522): 767-779.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020.
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2020_CEopt.pdf":Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}.{p_end}
{p 8 8}Working paper.{p_end}

{p 4 8}Cattaneo, M. D., Michael Jansson, and Xinwei Ma. 2020a.
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf":Simple Local Polynomial Density Estimators}.{p_end}
{p 8 8}{it:Journal of the American Statistical Association} 115(531): 1449-1455.{p_end}

{p 4 8}Cattaneo, M. D., Michael Jansson, and Xinwei Ma. 2020b.
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JoE.pdf":Local Regression Distribution Estimators}.{p_end}
{p 8 8}Working paper.{p_end}

{p 4 8}Cattaneo, M. D., Michael Jansson, and Xinwei Ma. 2020c.
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JSS.pdf":lpdensity: Local Polynomial Density Estimation and Inference}.{p_end}
{p 8 8}Working paper.{p_end}

{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Michael Jansson, University of California Berkeley, Berkeley, CA.
{browse "mailto:mjansson@econ.berkeley.edu":mjansson@econ.berkeley.edu}.{p_end}

{p 4 8}Xinwei Ma, University of California San Diego, La Jolla, CA.
{browse "mailto:x1ma@ucsd.edu":x1ma@ucsd.edu}.{p_end}


