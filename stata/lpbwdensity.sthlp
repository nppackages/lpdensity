{smcl}
{* *!version 2.3.1 2022-05-03}{...}

{title:Title}

{p 4 8}{cmd:lpbwdensity} {hline 2} Bandwidth Selection for Local Polynomial Density Estimation and Inference.{p_end}


{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:lpbwdensity} {it:Var} {ifin} 
[{cmd:,} {p_end}
{p 16 20}{cmd:grid(}{it:Var}{cmd:)} 
{cmd:p(}{it:#}{cmd:)}
{cmd:v(}{it:#}{cmd:)}
{cmd:kernel(}{it:KernelFn}
{cmd:)}{cmd:bwselect(}{it:BwMethod}{cmd:)}
{cmd:nomasspoints}
{cmd:nostdvar}{p_end}
{p 16 20}
{cmd:nlocalmin(}{it:#}{cmd:)}
{cmd:nuniquemin(}{it:#}{cmd:)}
{cmd:noregularize}
{p_end}
{p 16 20}{cmd:cweights(}{it:Var}{cmd:)}
{cmd:pweights(}{it:Var}{cmd:)}{p_end}
{p 16 20}{cmd:genvars(}{it:NewVarName}{cmd:)}{p_end}
{p 16 20}{cmd:rgrid(}{it:Var}{cmd:)} 
{cmd:rindex(}{it:Var}{cmd:)} 
{cmd:separator(}{it:#}{cmd:)}{p_end}
{p 16 20}]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:lpbwdensity} implements the bandwidth selection methods for local polynomial based density (and derivatives) estimation proposed and studied in 
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf":Cattaneo, Jansson and Ma (2020)}
and
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE.pdf":Cattaneo, Jansson and Ma (2022a)}. 
See {browse "https://doi.org/10.18637/jss.v101.i02":Cattaneo, Jansson and Ma (2022b)} for more 
implementation details and illustrations.{p_end}

{p 8 8} Companion command: {help lpdensity:lpdensity} for estimation and robust bias-corrected inference.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} functions are also available {browse "https://nppackages.github.io/lpdensity":here}.{p_end}

{p 4 8} Related Stata and R packages are available in the following website:{p_end}

{p 8 8} {browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Bandwidth Selection}

{p 4 8}{opt gri:d}({it:Var}) specifies the grid on which density is estimated. When set to default, grid points will be chosen as 0.05-0.95
percentiles of the data, with 0.05 step size.{p_end}

{p 4 8}{opt p}({it:#}) specifies the the order of the local-polynomial used to construct point estimates.
Default is {cmd:p(2)} (local quadratic regression).{p_end}

{p 4 8}{opt v}({it:#}) specifies the derivative of distribution function to be estimated. {cmd:v(0)} for
the distribution function, {cmd:v(1)} (default) for the density funtion, etc.{p_end}

{p 4 8}{opt ker:nel}({it:KernelFn}) specifies the kernel function used to construct the local-polynomial estimator(s). {p_end}
{p 8 12}{opt triangular}{bind:  } {it:K(u) = (1 - |u|) * (|u|<=1)}. This is the default option. {p_end}
{p 8 12}{opt epanechnikov}{bind:}  {it:K(u) = 0.75 * (1 - u^2) * (|u|<=1)}. {p_end}
{p 8 12}{opt uniform}{bind:     }  {it:K(u) = 0.5 * (|u|<=1)}. {p_end}

{p 4 8}{opt bws:elect}({it:BwMethod}) specifies method for data-driven bandwidth selection. This option will be
ignored if {cmd:bw(}{it:Var}{cmd:)} is provided.{p_end}
{p 8 12}{opt mse-dpi}{bind: } mean squared error optimal bandwidth for each grid point. This is the default option.{p_end}
{p 8 12}{opt imse-dpi}{bind:} integrated mean squared error optimal bandwidth which is common for all grid points.{p_end}
{p 8 12}{opt mse-rot}{bind: } rule-of-thumb bandwidth based on a Gaussian reference model.{p_end}
{p 8 12}{opt imse-rot}{bind:} integrated rule-of-thumb bandwidth based on a Gaussian reference model which is common for all grid points.{p_end}

{p 4 8}{opt nomass:points} will not adjust for mass points in the data.{p_end}

{p 4 8}{opt nostdv:ar} will not standardize the data for bandwidth selection. Note that this may lead to unstable performance of the numerical optimization procedure. {p_end}

{dlgtab:Local Sample Size Checking}

{p 4 8}{opt nloc:almin}({it:#}) specifies the minimum number of observations in each local neighborhood. This option will be ignored if set to 0, or if {cmd:noregularize} is used. The default value is {cmd:20+p(}{it:#}{cmd:)+1}. {p_end}

{p 4 8}{opt nuni:quemin}({it:#}) specifies the minimum number of unique observations in each local neighborhood. This option will be ignored if set to 0, or if {cmd:noregularize} is used. The default value is {cmd:20+p(}{it:#}{cmd:)+1}. {p_end}

{p 4 8}{opt noreg:ularize} suppresses local sample size checking.{p_end}

{dlgtab:Weights}

{p 4 8}{opt cw:eights}({it:Var}) specifies weights used for counterfactual distribution construction.{p_end}
 
{p 4 8}{opt pw:eights}({it:Var}) specifies weights used in sampling. Should be nonnegative.{p_end}

{dlgtab:Storing and displaying results}

{p 4 8}{opt gen:vars}({it:NewVarName}) specifies if new varaibles should be generated to store estimation results. If {it:newVarName} is provided, the following new varaibles will be
generated: {p_end}
{p 8 12}{it:NewVarName_grid}{bind:} grid points, {p_end}
{p 8 12}{it:NewVarName_bw}{bind:  } bandwidths, {p_end}
{p 8 12}{it:newVarName_nh}{bind:  } local/effective sample sizes.{p_end}

{p 4 8}{opt rgri:d}({it:var}) specifies a set of grid points to display the results. When omitted, this will be the same as {cmd:grid(}{it:Var}{cmd:)}.{p_end}

{p 4 8}{opt rind:ex}({it:var}) specifies a set of indices to display the results. This option will be ignored if {cmd:rgrid(}{it:Var}{cmd:)} is provided.{p_end}

{p 4 8}{opt sep:arator}({it:#}) draw a seperation line after every {it:#} variables; default is {cmd:separator(5)}.{p_end}
	
		
{marker examples}{...}
{title:Examples}

{p 4 8}Generate artifitial data:{p_end}
{p 8 8}{cmd:. set obs 1000}{p_end}
{p 8 8}{cmd:. set seed 42}{p_end}
{p 8 8}{cmd:. gen lpd_data = rnormal()}{p_end}

{p 4 8}MSE-optimal bandwidths for empirical quantiles: {p_end}
{p 8 8}{cmd:. lpbwdensity lpd_data}{p_end}

{p 4 8}Save estimation results to variables:{p_end}
{p 8 8}{cmd:. capture drop temp_*}{p_end}
{p 8 8}{cmd:. lpbwdensity lpd_data, genvars(temp)}{p_end}


{marker saved_results}{...}
{title:Saved results}

{p 4 8}{cmd:lpbwdensity} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(p)}}option {cmd:p(}{it:#}{cmd:)}{p_end}
{synopt:{cmd:e(v)}}option {cmd:v(}{it:#}{cmd:)}{p_end}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(bwselect)}}option {cmd:bwselect(}{it:BwMethod}{cmd:)}{p_end}
{synopt:{cmd:e(kernel)}}option {cmd:kernel(}{it:KernelFn}{cmd:)}{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(result)}}estimation result{p_end}

{marker references}{...}
{title:References}

{p 4 8}Cattaneo, M. D., Michael Jansson, and Xinwei Ma. 2020.
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf":Simple Local Polynomial Density Estimators}.{p_end}
{p 8 8}{it:Journal of the American Statistical Association} 115(531): 1449-1455.{p_end}

{p 4 8}Cattaneo, M. D., Michael Jansson, and Xinwei Ma. 2022a.
{browse "https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE.pdf":Local Regression Distribution Estimators}.{p_end}
{p 8 8}{it:Journal of Econometrics}, forthcoming.{p_end}

{p 4 8}Cattaneo, M. D., Michael Jansson, and Xinwei Ma. 2022b.
{browse "https://doi.org/10.18637/jss.v101.i02":lpdensity: Local Polynomial Density Estimation and Inference}.{p_end}
{p 8 8}{it:Journal of Statistical Software} 101(2): 1-25.{p_end}



{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Michael Jansson, University of California Berkeley, Berkeley, CA.
{browse "mailto:mjansson@econ.berkeley.edu":mjansson@econ.berkeley.edu}.{p_end}

{p 4 8}Xinwei Ma, University of California San Diego, La Jolla, CA.
{browse "mailto:x1ma@ucsd.edu":x1ma@ucsd.edu}.{p_end}


