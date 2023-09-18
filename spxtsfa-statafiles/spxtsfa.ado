*! version 1.0.2
*18Sep2023
*! version 1.0.1
*30July2023
* checkupdate bug fixes
* _checkspmat bug fixes

*! version 1.0
*27May2023
/*-------------------------------------------------------------------------------*/
*! STATA command: spxtsfa *! Author: [Your Name] *! Help file created on [Date]
/*
Description:
    This program estimates the spatial panel SFA model.

Syntax:
    spxtsfa varlist, Uhet(varlist) [INItial(name) NOCONstant NORMalize(string) wu(string) wv(string) ///
                                DELmissing MLPLOT NOGraph MLMODELopt(string) level(real 95) COST wxvars(varlist) ///
                                MLSEarch(string) MLMAXopt(string) DELVE CONSTraints(string) wy(string) wx(string)]

Options:
    varlist         - variables to be used in the model
    Uhet(varlist)   - variables for the unobserved heterogeneity
    INItial(name)   - matrix name of the initial values 
    NOCONstant      - suppress the constant in the model
    NORMalize(string)- normalize the spatial weight matrix
    wu(string)      - weight matrix for the inefficiency term
    wv(string)      - weight matrix for the random error term
    DELmissing      - delete units with missing values
    MLPLOT          - use ml plot to search better initail values
    NOGraph         - suppress the graphs of ml plot
    MLMODELopt(string)- options for the ml model command
    level(real 95)  - confidence level for the confidence intervals
    COST            - specificy the frontier is a cost function. The default is production function.
    wxvars(varlist) - variables to be included in the weights matrix
    MLSEarch(string) - options for ml search
    MLMAXopt(string) - options for ml max
    DELVE           - use fontier command for initial estimation
    CONSTraints(string)- constraints on the coefficients
    wy(string)      - weight matrix for the outcome variable
    wx(string)      - weight matrix for the covariates

Examples:
    . spxtsfa y x1 x2, Uhet(id) wu(w) wv(rho) NOCONstant MLPLOT

    . spxtsfa y x1 x2, Uhet(id) NOCONstant NORMalize(1)

Remarks:
    - The program estimates the the spatial panel SFA model.
    - The program uses a maximum likelihood approach for estimation.
*/


capture program drop spxtsfa
program define spxtsfa
version 16

spxtsfa_vparse `0'
local version `r(version)'
if("`version'"!=""){
	dis "The installed version of spxtsfa is `version'"
	checkupdate spxtsfa
	exit
}

	if replay() {
		if      (`"`e(cmd)'"' == "spxtsfayuv")  spxtsfayuv   `0'
		else if (`"`e(cmd)'"' == "spxtsfayuv0") spxtsfayuv0  `0'
		else if (`"`e(cmd)'"' == "spxtsfayu")   spxtsfayu    `0'
		else if (`"`e(cmd)'"' == "spxtsfayv")   spxtsfayv    `0'
		else if (`"`e(cmd)'"' == "spxtsfay")    spxtsfay     `0'
		else if (`"`e(cmd)'"' == "spxtsfau")    spxtsfau     `0'
		else if (`"`e(cmd)'"' == "spxtsfav")    spxtsfav     `0'
		else if (`"`e(cmd)'"' == "spxtsfauv")   spxtsfauv    `0'
		else if (`"`e(cmd)'"' == "spxtsfauv0")  spxtsfauv0   `0'
		else error 301
	}
	else Estimate `0'
end

program Estimate, eclass sortpreserve

cap which lmoremata.mlib

if _rc ssc install moremata

qui mata mata mlib index

syntax varlist, Uhet(varlist) [INItial(name) NOCONstant NORMalize(string) wu(string) wv(string) ///
                              te(name)  mldisplay(string) ///
                              DELmissing MLPLOT NOGraph MLMODELopt(string) level(real 95) COST wxvars(varlist) ///
							  MLSEarch(string) MLMAXopt(string) DELVE CONSTraints(string) wy(string) wx(string) ///
							  lndetfull lndetmc(numlist >0 min=2 max=2) GENWXVARS NOLOG] 

*******************************
** check the options regarding the spatial weight matirxs
if ("`wy'"=="" & "`wu'"=="" & "`wv'"==""){
	di as error "No spatial weight matrixs are specified. At least one of wy(), wu() and wv() is required."
	exit 198
}

if ("`wxvars'"!="" & "`wx'"==""){
	di as error "varlist is specified in wxvars(), but spmatrix is not specified in wx()"
	exit 198
}

if ("`wxvars'"=="" & "`wx'"!=""){
	di  `"varlist is not specified in wxvars(), wx(`wx') is neglected"'
}

if (`"`wy'"'==""){
	if ("`lndetfull'"!=""){
		di "lndetfull is neglected as wy() is not specified"
	}
	if ("`lndetmc'"!=""){
		di "lndetmc(`lndetmc') is neglected as wy() is not specified"
	}	
}

if("`lndetfull'"!="" & "`lndetmc'"!=""){
	di as error "lndetfull and lndetmc() cannot be combined."
	error 198
}

if (`"`te'"'!="") confirm new var `te'
*******************************
if(`"`wy'"'!=""){
	tightwmat `wy'
	local wy `tightwmat'
} 
if(`"`wx'"'!=""){
	tightwmat `wx'
	local wx `tightwmat'
} 
if(`"`wu'"'!=""){
	tightwmat `wu'
	local wu `tightwmat'
} 
if(`"`wv'"'!=""){
	tightwmat `wv'
	local wv `tightwmat'
} 
******************************
* Restricted cases
** case 1, wy
if ("`wy'"!="" & "`wu'"=="" & "`wv'"==""){
	spxtsfay `0'
	exit
}
** case 2, wy wu
if ("`wy'"!="" & "`wu'"!="" & "`wv'"==""){
	spxtsfayu `0'
	exit
}
** case 3, wy wu wv
if ("`wy'"!="" & "`wu'"=="" & "`wv'"!=""){
	spxtsfayv `0'
	exit
}
** case 4, wu=wv
if ("`wy'"=="" & ("`wu'"!="" & "`wv'"=="`wu'")){
	spxtsfauv0 `0'
	exit
}
** case 5, wu,wv
if ("`wy'"=="" & ("`wu'"!="" & "`wv'"!="")){
	spxtsfauv `0'
	exit
}
** case 6, wu
if ("`wy'"=="" & ("`wu'"!="" & "`wv'"=="")){
	spxtsfau `0'
	exit
}
** case 7, wv
if ("`wy'"=="" & ("`wu'"=="" & "`wv'"!="")){
	spxtsfav `0'
	exit
}

** case 8, wy==wu==wv
if ("`wy'"=="`wu'" & "`wv'"=="`wu'"){
	spxtsfayuv0 `0'
	exit
}
******************************
** case 9, general case: wy wu wv
    spxtsfayuv `0'

end

//include spxtsfa_aug.ado

//////utility comands and function for spxtsfa////

cap program drop genwxvars
program define genwxvars,rclass

version 16

syntax varlist, aname(name) [tvar(varname)]

if `"`tvar'"'==""{
	tempvar tvar 
	qui gen  byte `tavr'=1
}


mata: _genwxvars("`varlist'",`aname',"`tvar'")
return local wxnames `wxnames'

end

////////////////////////////////////////////
capture program drop checkspmat
program define checkspmat, rclass

syntax namelist(name=wnames), time(varname) touse(varname) [DELMissing NORMalize(string)]

//preserve

qui count if `touse'==0
local n0 = r(N)

if `n0'>0 & "`delmissing'"==""{
	di as red "missing values found. use delmissing to remove the units from the spmatrix"
	error 198
}


if `n0'>0 & "`delmissing'"!=""{
	di  "missing values found. The corresponding units are deleted from the spmatrix" _n
}

if "`normalize'"==""{
	local ntype=0
}
else if "`normalize'"=="row"{
	local ntype=1
}
else if "`normalize'"=="col"{
	local ntype=2
}
else if "`normalize'"=="spectral"{
	local ntype=3
}
else if "`normalize'"=="minmax"{
	local ntype=4
}
else{
	di as error "errors in normalize(), one of {row,col,spectral,minmax} should be specified. "
	error 198
}


//mata mata describe
foreach w in `wnames'{
    mata: _checkspmat("`time' `touse'",`w',`ntype')
    return scalar rmin_`w' = r(rmin)
    return scalar rmax_`w' = r(rmax)
}


end


///////////////////////////////////

capture program drop parsespmat0
program define parsespmat0,rclass
syntax namelist(name=wnames),[MATA ARRAY] 
if "`mata'"=="" & "`array'"==""{
    return local ldot=","
}

end

///////////////////////////////////
capture program drop parsespmat1
program define parsespmat1, rclass
syntax namelist(name=wnames),aname(name) [MATA ARRAY] 
local nw: word count `wnames'

local i=1
if "`mata'"=="" & "`array'"==""{
	mata: `aname' = asarray_create("real")
	foreach w in `wnames'{
		tempname w`i' w`i'_id
		spmatrix matafromsp `w`i'' `w`i'_id' = `w'
		local matanames `matanames' `w`i''
		mata: asarray(`aname',`i',`w`i'')
		local i=`i'+1
	}
	cap mata mata drop `matanames'
	
}
else if "`mata'"!=""{
	mata: `aname' = asarray_create("real")
	local matanames `wnames'
	local i=1
	foreach w in `matanames'{
		mata: asarray(`aname',`i',`w')
		local i=`i'+1
	}

}
else{
	mata: _temparray = asarray_create("real")
	mata: keys = asarray_keys(`wnames')
    mata: keys = sort(keys,1) // sort w in time order
	mata: st_local("keytypes",eltype(keys))
	if ("`keytypes'"!="real"){
		di as error "keys in array `wnames' is not real"
		exit 198
	}
    mata: st_numscalar("r(nw)",length(keys))
	local nw = r(nw)
    forv j=1/`nw'{
		mata: asarray(_temparray,`j',asarray(`wnames',keys[`j']))
	}
	mata: `aname' = asarray_create("real")
	forv j=1/`nw'{
		mata: asarray(`aname',`j',asarray(_temparray,`j'))
	}
	cap mata mata drop _temparray

 }

return scalar nw=`nw'

end 


//////////////////////////////////////

capture program drop issorted
program define issorted
	syntax	varlist 
	
	local sorted : dis "`:sortedby'"
	if "`sorted'" != "`varlist'" {
	    noi disp `"sort data by `varlist'"'
		noi disp "make sure that each spmatrix is the same order" _n
	    sort `varlist'
	}

end

// -------------------------------------------------------------------------------------------------



cap program drop spxtsfa_vparse
program define spxtsfa_vparse,rclass

syntax [varlist], [VERSION Gitee *]

if("`version'"!=""){
  return local version 1.0.2
}
if "`gitee'"!="" global gitee gitee

end

cap program drop tightwmat
program define tightwmat
syntax [anything], [Mata Array]

gettoken i anything:anything
local tightmat `i'
while(`"`i'"'!=""){
	gettoken i anything:anything
	local tightmat `tightmat' `i'
}

if ("`mata'"!="") local tightmat `tightmat',mata
if ("`array'"!="") local tightmat `tightmat',array

c_local tightwmat `tightmat'
end

/////////
cap program drop checkupdate
program define checkupdate
local url1 https://raw.githubusercontent.com/kerrydu/spxtsfa/main/spxtsfa-statafiles
local url2 https://gitee.com/kerrydu/spxtsfa/raw/main/spxtsfa-statafiles
if `"$gitee"' != ""{
	cap mata: vfile = cat(`"`url2'/`0'.ado"')
	if _rc exit
}
else{
	cap mata: vfile = cat(`"`url1'/`0'.ado"')
	if _rc{
		cap mata: vfile = cat(`"`url2'/`0'.ado"')
	}
	if _rc exit
}


mata: vfile = select(vfile,vfile:!="")
mata: vfile = usubinstr(vfile,char(9)," ",.)
mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
mata: st_local("versiongit",vfile[1])
local versiongit = ustrregexrf("`versiongit'","^[\D]+","")
gettoken vers versiongit:versiongit, p(", ")
local versiongit `vers'

qui findfile `0'.ado
mata: vfile = cat("`r(fn)'")
mata: vfile = select(vfile,vfile:!="")
mata: vfile = usubinstr(vfile,char(9)," ",.)
mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
mata: st_local("versionuse",vfile[1])
local versionuse = ustrregexrf("`versionuse'","^[\D]+","")
gettoken vers versionuse:versionuse, p(", ")
local versionuse `vers'

if("`versionuse'"!="`versiongit'"){
	di "New version available, `versionuse' =>`versiongit'"
	di "It can be updated by:"
	di "  net install `0',from(`url1') replace force"
	di "or,"
	di "  net install `0',from(`url2') replace force"
}


end