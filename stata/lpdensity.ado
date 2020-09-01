********************************************************************************
* LPDENSITY STATA PACKAGE -- lpdensity
* Authors: Matias D. Cattaneo, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 2.2 2020-09-01

capture program drop lpdensity

program define lpdensity, eclass
syntax 	varlist(max=1) 				///
									///
		[if] [in] [, 				///
									///
		GRId(varname) 				///
		BW(string) 					///
		P(integer 2) 				///
		Q(integer -1) 				///
		V(integer 1) 				///
		BWSelect(string) 			///
		noREGularize 			    ///
		NLOCalmin (integer -1)		///
		NUNIquemin (integer -1)		///
		noMASSpoints				///
		noSTDVar					///
		KERnel(string) 				///
		CWeights(varname) 			///
		PWeights(varname) 			///
		SCAle(real 1) 				///
		GENvars(string) 			///
		RGRId(varname)				///
		RINDex(varname)				///
		Level(real 95) 				///
		CIUniform					///
		CISimul(integer 2000)		///
		SEParator(integer 5)		///
		PLot						///
		ESTYpe(string)				///
		ESLIne_options(string)		///
		ESPOint_options(string)		///
		CITYpe(string)				///
		CIREgion_options(string)	///
		CILIne_options(string)		///
		CIEBar_options(string)		///
		HISTogram		 			///
		HIPLot_options(string)		///
		GRAph_options(string)		///
		]

marksample touse

********************************************************************************
**** error check: main variable
local x "`varlist'"
cap confirm numeric variable `x'
if _rc {
	di as err `"`grid' is not a numeric variable"'
	exit 198
}
else {
	qui count if `x' != . & `touse'
	local n = r(N)
	if (`n' == 0) {
		di as err `"`x' has length 0"'
		exit 198
	}
}

********************************************************************************
**** error check: grid()
if ("`grid'" == "") {
	tempvar temp_grid	
	qui pctile `temp_grid' = `x' if `touse', nq(20)
	local ng = 19
} 
else {
	cap confirm numeric variable `grid'
	if _rc {
		di as err `"grid(): `grid' is not a numeric variable"'
		exit 198
	}
	else {
		if ("`x'" == "`grid'") {
			qui count if `grid' != . & `touse'
		}
		else {
			qui count if `grid' != .
		}
		
		local ng = r(N)
		//disp `ng'
		if (`ng' == 0) {
			di as err `"grid(): `grid' has length 0"'
			exit 198
		}
	}
	local temp_grid "`grid'"
}

********************************************************************************
**** error check: rgrid() rindex()
if ("`rgrid'" != "") {
	cap confirm numeric variable `rgrid'
	if _rc {
		di as err `"rgrid(): `rgrid' is not a numeric variable"'
		exit 198
	}
	else {
		if ("`x'" == "`rgrid'") {
			qui count if `rgrid' != . & `touse'
		}
		else {
			qui count if `rgrid' != .
		}
		
		local nrg = r(N)
		if (`nrg' == 0) {
			di as err `"rgrid(): `rgrid' has length 0"'
			exit 198
		}
	}
}

if ("`rindex'" != "") {
	cap confirm numeric variable `rindex'
	if _rc {
		di as err `"rindex(): `rindex' is not a numeric variable"'
		exit 198
	} 
	else {
		qui count if `rindex' != .
		
		local nrg = r(N)
		if (`nrg' == 0) {
			di as err `"rindex(): `rindex' has length 0"'
			exit 198
		}
	}
}

********************************************************************************
**** error check: p() q() v()
if (`p' < 0 | `p' > 7) {
	di as err `"p(): has to be integer between 0 and 7"'
	exit 198
}
if (`q' == -1) {
	local q = `p' + 1
}
if (`q' < `p' | `q' > 7) {
	di as err `"q(): has to be integer between `p' and 7"'
	exit 198
}
if (`v' < 0 | `v' > `p') {
	di as err `"v(): has to be integer between 0 and `p'"'
	exit 198
}

********************************************************************************
**** error check: bwselect()
if ("`bwselect'" == "") {
	local bwselect = "mse-dpi"
} 
else {
	if ("`bwselect'" != "mse-dpi" & "`bwselect'" != "imse-dpi" & "`bwselect'" != "mse-rot" & "`bwselect'" != "imse-rot") {
		di as err `"bwselect(): incorrectly specified: options(mse-dpi, imse-dpi, mse-rot, imse-rot)"'
		exit 198
	}
}

********************************************************************************
**** error check: kernel()
if ("`kernel'" == "") {
	local kernel = "triangular"
} 
else {
	if ("`kernel'" != "triangular" & "`kernel'" != "uniform" & "`kernel'" != "epanechnikov") {
		di as err `"kernel(): incorrectly specified: options(triangular, uniform, epanechnikov)"'
		exit 198
	}
}

********************************************************************************
**** error check: scale()
if (`scale' < 0) {
	di as err `"scale(): incorrectly specified: should be between 0 and 1"'
	exit 198
}

********************************************************************************
**** error check: level()
if (`level' <= 0 | `level' >= 100) {
	di as err `"level(): incorrectly specified: should be between 0 and 100"'
	exit 198
}

********************************************************************************
**** error check: cisimul()
if (`cisimul' <= 0) {
	di as err `"cisimul(): incorrectly specified: should be strictly positive"'
	exit 198
}

********************************************************************************
**** error check: cweights()
if ("`cweights'" == "") {
	tempvar temp_cweights
	qui gen `temp_cweights' = `x' * 0 + 1
}
else {
	cap confirm numeric variable `cweights'
	if _rc {
		di as err `"cweights(): `cweights' is not a numeric variable"'
		exit 198
	}
	else {
		qui count if `cweights' != . & `touse'
		local nCw = r(N)
		if (`nCw' != `n') {
			di as err `"cweights(): `cweights' has different length as `x'"'
			exit 198
		}
	}
	local temp_cweights "`cweights'"
}

********************************************************************************
**** error check: pweights()
if ("`pweights'" == "") {
	tempvar temp_pweights
	qui gen `temp_pweights' = `x' * 0 + 1
}
else {
	cap confirm numeric variable `pweights'
	if _rc {
		di as err `"pweights(): `pweights' is not a numeric variable"'
		exit 198
	}
	else {
		qui count if `pweights' != . & `touse'
		local nPw = r(N)
		if (`nPw' != `n') {
			di as err `"pweights(): `pweights' has different length as `x'"'
			exit 198
		}
	}
	local temp_pweights "`pweights'"
}

********************************************************************************
**** error check: noregularize nomasspoints nostdvar
**** WARNING: this is a very tricky part of STATA, when an option starts with "no"
if ("`regularize'" == "") {
	local regularize = 1
}
else {
	local regularize = 0
}

if ("`masspoints'" == "") {
	local masspoints = 1
}
else {
	local masspoints = 0
}

if ("`stdvar'" == "") {
	local stdvar = 1
}
else {
	local stdvar = 0
}

********************************************************************************
**** error check: nlocalmin() nuniquemin()
if (`nuniquemin' < 0) {
	local nuniquemin = 20 + `p' + 1
}

if (`nuniquemin' < 0) {
	local nuniquemin = 20 + `p' + 1
}
 
********************************************************************************
**** error check: bw()
if ("`bw'" == "") {
	tempvar temp_bw
	tempvar temp_lpbwdensity
	
	if (`regularize' == 0) {
		local flag_regularize = "noregularize"
	}
	else {
		local flag_regularize = "nlocalmin(`nlocalmin') nuniquemin(`nuniquemin')"
	}
	
	if (`masspoints' == 0) {
		local flag_masspoints = "nomasspoints"
	}
	else {
		local flag_masspoints = ""
	}
	
	if (`stdvar' == 0) {
		local flag_stdvar = "nostdvar"
	}
	else {
		local flag_stdvar = ""
	} 
	
	qui lpbwdensity `x' if `touse', grid(`temp_grid') p(`p') v(`v') bwselect(`bwselect') kernel(`kernel') cweights(`temp_cweights') pweights(`temp_pweights') genvars(`temp_lpbwdensity') separator(1) `flag_regularize' `flag_masspoints' `flag_stdvar'
	
	qui gen `temp_bw' = `temp_lpbwdensity'_bw
	qui capture drop `temp_lpbwdensity'_grid `temp_lpbwdensity'_bw `temp_lpbwdensity'_nh
}
else {

	cap confirm numeric variable `bw'
	if _rc {
		// check if variable exists
		capture confirm variable `bw'
		if (!_rc) {
			di as err `"bw(): `bw' is not a numeric variable"'
			exit 198
		} 
		else {
			tempvar temp_bw			
			qui gen `temp_bw' = "`bw'" if `temp_grid' != .
			capture destring `temp_bw', replace
			cap confirm numeric variable `temp_bw'
			if _rc {
				di as err `"bw(): `bw' is not numeric"'
				exit 198
			}
		}
	}
	else {
		if ("`x'" == "`grid'") {
			qui count if `bw' != . & `touse'
		}
		else {
			qui count if `bw' != .
		}
		
		local nb = r(N)
		if (`nb' != `ng') {
			di as err `"bw(): `bw' has different length as grid()"'
			exit 198
		}
		local temp_bw "`bw'"
		local bwselect = "user provided"
	}
}

********************************************************************************
**** error check: separator()
if (`separator' <= 1) {
	local separator = 1
}

********************************************************************************
**** error check: citype()
if ("`citype'" == "") {
	local citype = "region"
}
else if ("`citype'" != "all" & "`citype'" != "region" & "`citype'" != "line" & "`citype'" != "ebar" & "`citype'" != "none") {
	di as err `"citype(): incorrectly specified: options(region, line, ebar, all, none)"'
	exit 198
}

********************************************************************************
**** error check: estype()
if ("`estype'" == "") {
	local estype = "line"
}
else if ("`estype'" != "both" & "`estype'" != "line" & "`estype'" != "point" & "`estype'" != "none") {
	di as err `"estype(): incorrectly specified: options(line, point, both, none)"'
	exit 198
}

********************************************************************************
**** temporaty varaibles for plotting
if ("`plot'" != "" & "`genvars'" == "") {
	tempvar plot_grid
	tempvar plot_f
	tempvar plot_cil
	tempvar plot_cir
}

********************************************************************************
**** MATA
tempvar temp_touse
qui gen `temp_touse' = `touse'
mata{
	x = st_data(., "`x'", "`temp_touse'") //; x
	
	if ("`x'" == "`grid'") {
		grid     = st_data(., "`temp_grid'"	, "`temp_touse'") //; grid
		bw       = st_data(., "`temp_bw'"	, "`temp_touse'") //; bw
	}
	else {	
		grid     = st_data(., "`temp_grid'"	, 0) //; grid
		bw       = st_data(., "`temp_bw'"	, 0) //; bw
	}

	cweights = st_data(., "`temp_cweights'"	, "`temp_touse'") //; cweights
	pweights = st_data(., "`temp_pweights'"	, "`temp_touse'") //; pweights
	
	if ("`ciuniform'" == "") {
		Result = lpdensity_fn(x, grid, bw, `p', `q', `v', "`kernel'", cweights, pweights, `masspoints', `scale', 1-`level'/100, 0, `cisimul')
	} 
	else {
		Result = lpdensity_fn(x, grid, bw, `p', `q', `v', "`kernel'", cweights, pweights, `masspoints', `scale', 1-`level'/100, 1, `cisimul')
	}
	
	
	genvars = st_local("genvars")
	if (genvars != "") {
		(void) 	st_addvar("double", 	invtokens((genvars, "grid"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "grid"), 		"_"), Result[., 1])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "bw"), 			"_"))
				st_store(Result[., 10], invtokens((genvars, "bw"), 			"_"), Result[., 2])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "nh"), 			"_"))
				st_store(Result[., 10], invtokens((genvars, "nh"), 			"_"), Result[., 3])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "f_p"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "f_p"), 		"_"), Result[., 4])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "se_p"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "se_p"), 		"_"), Result[., 6])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "CI_l"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "CI_l"), 		"_"), Result[., 8])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "CI_r"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "CI_r"), 		"_"), Result[., 9])
					
		if (`q' > `p') {
		(void) st_addvar("double", 		invtokens((genvars, "f_q"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "f_q"), 		"_"), Result[., 5])
		
		(void) st_addvar("double", 		invtokens((genvars, "se_q"), 		"_"))
				st_store(Result[., 10], invtokens((genvars, "se_q"), 		"_"), Result[., 7])
		}
	}
	else if ("`plot'" != "") {
		(void) 	st_addvar("double", 	"`plot_grid'")
				st_store(Result[., 10], "`plot_grid'", Result[., 1])
				
		(void) st_addvar("double", 		"`plot_f'")
				st_store(Result[., 10], "`plot_f'", Result[., 4])
				
		(void) st_addvar("double", 		"`plot_cil'")
				st_store(Result[., 10], "`plot_cil'", Result[., 8])
				
		(void) st_addvar("double", 		"`plot_cir'")
				st_store(Result[., 10], "`plot_cir'", Result[., 9])
	}
	
	if ("`rgrid'" != "") {
		if ("`x'" == "`rgrid'") {
			rgrid     = st_data(., "`rgrid'", "`temp_touse'") //; rgrid
		}
		else {	
			rgrid     = st_data(., "`rgrid'", 0) //; rgrid
		}
		selectIndex = J(length(rgrid), 1, .)
		for (j=1; j<=length(rgrid);j++) {
			selectIndex[j] = lpdensity_whichmin(abs(rgrid[j] :- Result[., 1]))
		}
		st_matrix("Result", Result[selectIndex, 1..10])
	} 
	else {
		if ("`rindex'" != "") {
			rindex     = st_data(., "`rindex'", 0) //; rindex
			selectIndex = J(length(rindex), 1, .)
			for (j=1; j<=length(rindex);j++) {
				selectIndex[j] = lpdensity_whichmin(abs(rindex[j] :- Result[., 10]))
			}
			st_matrix("Result", Result[selectIndex, 1..10])
		}
		else {
			st_matrix("Result", Result[., 1..10])
		}
	}
}

********************************************************************************
**** generate variable labels
if ("`genvars'" != "") {
    label variable `genvars'_grid 	"lpdensity: grid point"
	label variable `genvars'_bw 	"lpdensity: bandwidth"
	label variable `genvars'_nh 	"lpdensity: effective sample size"
	label variable `genvars'_f_p 	"lpdensity: point estimate with pol. order `p'"
	label variable `genvars'_se_p 	"lpdensity: standard error for f_p"
	label variable `genvars'_CI_l 	"lpdensity: left level-`level' CI, conventional"
	label variable `genvars'_CI_r 	"lpdensity: right level-`level' CI, conventional"
	if (`q' > `p') {
	label variable `genvars'_f_q 	"lpdensity: point estimate with pol. order `q'"
	label variable `genvars'_se_q 	"lpdensity: standard error for f_q"
	label variable `genvars'_CI_l 	"lpdensity: left level-`level' CI, robust bias-corrected"
	label variable `genvars'_CI_r 	"lpdensity: right level-`level' CI, robust bias-corrected"
	}	
}

********************************************************************************
**** display
disp ""
disp "Local Polynomial Density Estimation and Inference." 
disp ""

disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %15.0f `n'
disp in smcl in gr "{lalign 1: Polynomial order for point estimation    (p=)    }" _col(19) in ye %15.0f `p'
if (`v' == 1) {
disp in smcl in gr "{lalign 1: Density function estimated               (v=)    }" _col(19) in ye %15.0f `v'
}
else if (`v' == 0) {
disp in smcl in gr "{lalign 1: Distribution function estimated          (v=)    }" _col(19) in ye %15.0f `v'
} 
else {
disp in smcl in gr "{lalign 1: Order of derivative estimated            (v=)    }" _col(19) in ye %15.0f `v'
}
disp in smcl in gr "{lalign 1: Polynomial order for confidence interval (q=)    }" _col(19) in ye %15.0f `q'
disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 15: `kernel'}"
if (`scale' < 1) {
disp in smcl in gr "{lalign 1: Scaling factor                                   }" _col(19) in ye %15.2f `scale'
}
disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 15: `bwselect'}"
disp ""

disp in smcl in gr "{hline 72}"
if (`q' > `p') {
disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: }" _col(14) "{ralign 10: }" _col(24) "{ralign 8: }" _col(32) "{ralign 10: Point}" _col(42) "{ralign 10: Std.}" _col(52) "{ralign 20: Robust B.C.}" _col(72)
} 
else {
disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: }" _col(14) "{ralign 10: }" _col(24) "{ralign 8: }" _col(32) "{ralign 10: Point}" _col(42) "{ralign 10: Std.}" _col(52) "{ralign 20: Conventional}" _col(72)
}
if ("`ciuniform'" == "") {
	disp in smcl in gr "{ralign 4: Index}" _col(4) "{ralign 8: Grid}" _col(14)  "{ralign 10: B.W.}" _col(24) "{ralign 8: Eff.n}" _col(32) "{ralign 10: Est.}" _col(42) "{ralign 10: Error}" _col(52) "{ralign 20: `level'% C.I.}" _col(72)
}
else {
	disp in smcl in gr "{ralign 4: Index}" _col(4) "{ralign 8: Grid}" _col(14)  "{ralign 10: B.W.}" _col(24) "{ralign 8: Eff.n}" _col(32) "{ralign 10: Est.}" _col(42) "{ralign 10: Error}" _col(52) "{ralign 20: Unif. `level'% C.I.}" _col(72)
}

disp in smcl in gr "{hline 72}"

local ng = rowsof(Result)

forvalues i = 1(1)`ng' {
disp in smcl in gr %4.0f Result[`i', 10] _col(4) in ye %10.4f Result[`i', 1] _col(14)  %10.4f Result[`i', 2] _col(24) %8.0f Result[`i', 3] _col(32) %10.4f Result[`i', 4] _col(42) %10.4f Result[`i', 6] _col(52) ///
    in ye %10.4f Result[`i', 8] _col(62) in ye %10.4f Result[`i', 9] _col(72)
	if (`i' != `ng' & `separator' > 1 & mod(`i', `separator') == 0) {
		disp in smcl in gr "{hline 72}"
	}
}
disp in smcl in gr "{hline 72}"

********************************************************************************
**** plot
if ("`plot'" != "" & "`genvars'" != "") {
	local plot_grid	"`genvars'_grid"
	local plot_f	"`genvars'_f_p"
	local plot_cil	"`genvars'_CI_l"
	local plot_cir 	"`genvars'_CI_r"
}

if ("`plot'" != "") {
	// graph option check
	if (`"`graph_options'"' == "" ) {
		local graph_options = `"legend(off) title("lpdensity (p=`p', q=`q', v=`v')", color(gs0)) xtitle("`x'") ytitle("")"'
	}
	
	// ci type check
	if ("`citype'" == "region" | "`citype'" == "all") {
		if ("`ciregion_options'" == "") {
			local ci_plot_region = `"(rarea `plot_cil' `plot_cir' `plot_grid', sort color(red%30))"'
		} 
		else {
			local ci_plot_region = `"(rarea `plot_cil' `plot_cir' `plot_grid', sort `ciregion_options')"'
		}
	} 
	else {
		local ci_plot_region = `""'
	}
	if ("`citype'" == "line" | "`citype'" == "all") {
		if ("`ciline_options'" == "") {
			local ci_plot_line = `"(rline `plot_cil' `plot_cir' `plot_grid', sort color(red%70))"'
		} 
		else {
			local ci_plot_line = `"(rline `plot_cil' `plot_cir' `plot_grid', sort `ciline_options')"'
		}
	}
	else {
		local ci_plot_line = `""'
	}
	if ("`citype'" == "ebar" | "`citype'" == "all") {
		if ("`ciebar_options'" == "") {
			local ci_plot_ebar = `"(rcap `plot_cil' `plot_cir' `plot_grid', sort color(red%70))"'
		} 
		else {
			local ci_plot_ebar = `"(rcap `plot_cil' `plot_cir' `plot_grid', sort `ciebar_options')"'
		}
	}
	else {
		local ci_plot_ebar = `""'
	}
	
	// point est type check
	if ("`estype'" == "line" | "`estype'" == "both") {
		if ("`esline_options'" == "") {
			local es_plot_line = `"(line `plot_f' `plot_grid', sort lcolor(red) lwidth("medthin") lpattern(solid))"'
		} 
		else {
			local es_plot_line = `"(line `plot_f' `plot_grid', sort `esline_options')"'
		}
	} 
	else {
		local es_plot_line = `""'
	}
	if ("`estype'" == "point" | "`estype'" == "both") {
		if ("`espoint_options'" == "") {
			local es_plot_point = `"(scatter `plot_f' `plot_grid', sort color(red))"'
		} 
		else {
			local es_plot_point = `"(scatter `plot_f' `plot_grid', sort `espoint_options')"'
		}
	} 
	else {
		local es_plot_point = `""'
	}
	
	// histogram check
	if ("`histogram'" != "") {
		if ("`hiplot_options'" == "") {
			local hist_plot = `"(histogram `x' if `touse', color(blue%30))"'
		}
		else {
			local hist_plot = `"(histogram `x' if `touse', `hiplot_options')"'
		}
	} 
	else {
		local hist_plot = `""'
	}
	
	twoway 	`hist_plot'	 		///
			`ci_plot_region' 	///
			`ci_plot_line'   	///
			`ci_plot_ebar'   	///
			`es_plot_line'   	///
			`es_plot_point'  	///
			, 					///
			`graph_options'
}

********************************************************************************
**** ereturn

ereturn clear
ereturn scalar N = `n'
ereturn scalar p = `p'
ereturn scalar q = `q'
ereturn scalar v = `v'
ereturn local kernel = "`kernel'"
ereturn local bwselect = "`bwselect'"
ereturn scalar scale = `scale'
ereturn scalar level = `level'
matrix colnames Result = grid bw nh f_p f_q se_p se_q CI_l CI_r
ereturn matrix result = Result

********************************************************************************
**** end

mata: mata clear

end
	
