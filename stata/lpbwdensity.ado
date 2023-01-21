********************************************************************************
* LPDENSITY STATA PACKAGE -- lpbwdensity
* Authors: Matias D. Cattaneo, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 2.4 2023-01-21

capture program drop lpbwdensity

program define lpbwdensity, eclass
syntax 	varlist(max=1) 				///
									///
		[if] [in] [, 				///
									///
		GRId(varname) 				///
		P(integer 2) 				///
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
		GENvars(string) 			///
		RGRId(varname)				///
		RINDex(varname)				///
		SEParator(integer 5) 		///
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
**** error check: p() v()
if (`p' < 0 | `p' > 7) {
	di as err `"p(): has to be integer between 0 and 7"'
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
**** error check: separator()
if (`separator' <= 1) {
	local separator = 1
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
if (`nlocalmin' < 0) {
	local nlocalmin = 20 + `p' + 1
}

if (`nuniquemin' < 0) {
	local nuniquemin = 20 + `p' + 1
}

********************************************************************************
**** MATA

tempvar temp_touse

qui gen `temp_touse' = `touse'

mata{
	x = st_data(., "`x'", "`temp_touse'") //; x
	
	if ("`x'" == "`grid'") {
		grid     = st_data(., "`temp_grid'", "`temp_touse'") //; grid
	}
	else {	
		grid     = st_data(., "`temp_grid'", 0) //; grid
	}
	
	cweights = st_data(., "`temp_cweights'", "`temp_touse'") //; cweights
	pweights = st_data(., "`temp_pweights'", "`temp_touse'") //; pweights
	
	Result = J(length(grid), 4, .)
	Result[., 1] = grid
	
	if ("`bwselect'" == "mse-dpi") {
		Result[., 2] = 	lpdensity_bwMSE( x, grid, `p', `v', "`kernel'", cweights, pweights, `masspoints', `stdvar', `regularize', `nlocalmin', `nuniquemin') 
	}
	else if ("`bwselect'" == "imse-dpi") {
		Result[., 2] = J(length(grid), 1, 
						lpdensity_bwIMSE(x, grid, `p', `v', "`kernel'", cweights, pweights, `masspoints', `stdvar', `regularize', `nlocalmin', `nuniquemin')) 
	}
	else if ("`bwselect'" == "mse-rot") {
		Result[., 2] = 	lpdensity_bwROT( x, grid, `p', `v', "`kernel'", cweights, pweights, `masspoints', `stdvar', `regularize', `nlocalmin', `nuniquemin') 
	}
	else {
		Result[., 2] = J(length(grid), 1, 
						lpdensity_bwIROT(x, grid, `p', `v', "`kernel'", cweights, pweights, `masspoints', `stdvar', `regularize', `nlocalmin', `nuniquemin')) 
	}
	
	x = st_data(., "`x'", "`temp_touse'") //; x
	
	for (j=1; j<=length(grid); j++) {
		Result[j, 3] = sum(abs(x :- Result[j, 1]) :<= Result[j, 2])
		Result[j, 4] = j
		
	}

	genvars = st_local("genvars")
	if (genvars != "") {
		(void) 	st_addvar("double", 	invtokens((genvars, "grid"), 	"_"))
				st_store(Result[., 4], 	invtokens((genvars, "grid"), 	"_"), Result[., 1])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "bw"), 		"_"))
				st_store(Result[., 4], 	invtokens((genvars, "bw"), 		"_"), Result[., 2])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "nh"), 		"_"))
				st_store(Result[., 4], 	invtokens((genvars, "nh"), 		"_"), Result[., 3])
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
		st_matrix("Result", Result[selectIndex, 1..4])
	} 
	else {
		if ("`rindex'" != "") {
			rindex     = st_data(., "`rindex'", 0) //; rindex
			selectIndex = J(length(rindex), 1, .)
			for (j=1; j<=length(rindex);j++) {
				selectIndex[j] = lpdensity_whichmin(abs(rindex[j] :- Result[., 4]))
			}
			st_matrix("Result", Result[selectIndex, 1..4])
		}
		else {
			st_matrix("Result", Result[., 1..4])
		}
	}
}

********************************************************************************
**** generate variable labels
if ("`genvars'" != "") {
    label variable `genvars'_grid 	"lpbwdensity: grid point"
	label variable `genvars'_bw 	"lpbwdensity: bandwidth"
	label variable `genvars'_nh 	"lpbwdensity: effective sample size"
}

********************************************************************************
**** display
disp ""
disp "Bandwidth Selection for Local Polynomial Density Estimation." 
disp ""

disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %12.0f `n'
disp in smcl in gr "{lalign 1: Polynomial order for point estimation    (p=)    }" _col(19) in ye %12.0f `p'
if (`v' == 1) {
disp in smcl in gr "{lalign 1: Density function estimated               (v=)    }" _col(19) in ye %12.0f `v'
}
else if (`v' == 0) {
disp in smcl in gr "{lalign 1: Distribution function estimated          (v=)    }" _col(19) in ye %12.0f `v'
} 
else {
disp in smcl in gr "{lalign 1: Order of derivative estimated            (v=)    }" _col(19) in ye %12.0f `v'
}
disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 12: `kernel'}"
disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 12: `bwselect'}"
disp ""

disp in smcl in gr "{hline 32}"
disp in smcl in gr "{lalign 4: Index}" _col(4) "{ralign 8: Grid}" _col(14)  "{ralign 10: B.W.}" _col(24) "{ralign 8: Eff.n}" _col(32)
disp in smcl in gr "{hline 32}"

local ng = rowsof(Result)

forvalues i = 1(1)`ng' {
disp in smcl in gr %4.0f Result[`i', 4] _col(4) in ye %10.4f Result[`i', 1] _col(14)  %10.4f Result[`i', 2] _col(24) %8.0f Result[`i', 3] _col(32) 
	if (`i' != `ng' & `separator' > 1 & mod(`i', `separator') == 0) {
		disp in smcl in gr "{hline 32}"
	}
}
disp in smcl in gr "{hline 32}"

********************************************************************************
**** return

ereturn clear
ereturn scalar N = `n'
ereturn scalar p = `p'
ereturn scalar v = `v'
ereturn local kernel = "`kernel'"
ereturn local bwselect = "`bwselect'"
matrix colnames Result = grid bw nh
ereturn matrix result = Result

********************************************************************************
**** end

mata: mata clear

end
	
