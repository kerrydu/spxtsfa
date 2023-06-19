//from pointer to array

capture program drop spxtsfauv
program define spxtsfauv, eclass sortpreserve
version 16

	if replay() {
		if (`"`e(cmd)'"' != "spxtsfauv") error 301
		Replay `0'
	}
	else	Estimate `0'
end

program Estimate, eclass sortpreserve

/*
syntax varlist(fv ts min=2),  Uhet(varlist fv ts)  [INItial(name) NOCONstant NORMalize(string) ///
                              wu(string) wv(string) ///
                              DELmissing MLPLOT NOGraph MLMODELopt(string) level(real 95) COST ///
							  MLSEarch(string) MLMAXopt(string) DELVE CONSTraints(string) te(name) ///
							  lndetmc(numlist >0 min=2 max=2) GENWXVARS] 
*/
syntax varlist,Uhet(varlist) [INItial(name) NOCONstant NORMalize(string) ///
                              wu(string) wv(string) mldisplay(string) wx(string) wxvars(varlist) ///
                              DELmissing MLPLOT NOGraph MLMODELopt(string) level(real 95) COST ///
							  MLSEarch(string) MLMAXopt(string) DELVE CONSTraints(string) te(name) ///
							  lndetmc(numlist >0 min=2 max=2) GENWXVARS NOLOG] 

//marksample touse 
if ("`nolog'"!="") local nolog qui
local cmdline spxtsfa `0'
if ("`wu'"=="" & "`wv'"==""){
	di as red "Spatial weight matrixs should be specified in wu() or/and wv()"
	error 198
}

/*
if("`wu'"!="" & "`wv'"==""){
    spxtsfau `0'
    exit
}
if("`wu'"=="" & "`wv'"!=""){
    spxtsfav `0'
    exit
}
*/
/*
if("`wu'"== "`wv'"){
    spxtsfauv0 `0'
    exit
}
*/

if ("`wxvars'"!="" & "`wx'"==""){
	di as error "varlist is specified in wxvars(), but spmatrix is not specified in wx()"
	exit 198
}

if ("`wxvars'"=="" & "`wx'"!=""){
	di  `"varlist is not specified in wxvars(), wx(`wx') is neglected"'
}
preserve

marksample touse 

local diopts level(`level') `mldisplay'
mlopts std, `mlmodelopt' `mlmaxopt' `constraints'
local cns constraints(`constraints')
gettoken yvar xvars: varlist 
//_fv_check_depvar `yvar'
if ("`initial'"!="" & "`delve'"!=""){
	di "Warning: initial(`initial') overrides delve"
}
if ("`initial'"!="" & "`mlsearch'"!=""){
	di "Warning: initial(`initial') overrides mlsearch(`mlsearch')"
}
if ("`delve'"!="" & "`mlsearch'"!=""){
	di "Warning: delve overrides mlsearch(`mlsearch')"
}

parsespmat0 `wu' 
parsespmat1 `wu' `r(ldot)' aname(wu_ina)
local nwu = r(nw)
parsespmat0 `wv' 
parsespmat1 `wv' `r(ldot)' aname(wv_ina)
local nwv = r(nw)

// 检查权重矩阵与数据是否匹配
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
    qui keep `varlist' `wxvars' `id' `time' `uhet' `touse'
    tempvar order0
    qui gen int `order0' =_n	
// sort 数据	
	qui issorted `time' `id'	
	//tempvar time2

	qui distinct2 `time'
	local T = r(ndistinct)	

	
	if (`nwv'!=1 & `nwv'!=`T') {
		di as error "Spatial weight matrixs in wv() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}
	if (`nwu'!=1 & `nwu'!=`T') {
		di as error "Spatial weight matrixs in wu() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}    
	tempvar time2
	qui egen `time2' = group(`time')
	global paneltvar `time2'
	//mata mata describe
	checkspmat wu_ina wv_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')

    scalar rumin = max(-0.9999,r(min_wu_ina))
	scalar rumax = min(0.9999,r(max_wu_ina))
	global rumin = rumin
	global rumax = rumax

    scalar rvmin = max(-0.9999,r(min_wv_ina))
	scalar rvmax = min(0.9999,r(max_wv_ina))
	global rvmin = rvmin
	global rvmax = rvmax   

**************
    if ("`wx'"!="" & `"`wxvars'"'!=""){
        parsespmat0 `wx' 
        parsespmat1 `wx' `r(ldot)' aname(wx_ina)
        local nwx = r(nw)
        if (`nwx'!=1 & `nwx'!=`T') {
            di as error "Spatial weight matrixs in wx() are specified as time-varying, but # of spmatrix != # of periods"
            exit 198
        }  
        checkspmat wx_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
    }

   * generating Wx
	if(`"`wxvars'"'!=""){
      qui genwxvars `wxvars', aname(wx_ina) tvar(`time2')
      local wxvars2  `r(wxnames)'
      mata: _order_wx = st_data(.,"`wxvars2'","`touse'")
	  cap mata mata drop wx_ina
	}	
************
    mata: _order_0   = st_data(.,"`order0'","`touse'") // record the row# in the original data
	//mata: marksuse = st_data(.,"`touse'")
	qui count if `touse'==0
	local nummissing = r(N)
	if(`nummissing'>0){
		mata: marksuse = st_data(.,"`touse'")
	}

    qui keep if `touse'
    mata: _pan_tvar =st_data( .,"`time2'")	   
	if("`initial'"=="" & "`delve'"!="") { 
		qui frontier `yvar' `xvars' `wxvars2',`noconstant' uhet(`uhet') iterate(50) `cns'
	    mat b0 =e(b)
        local rhomin = max(rumin,rvmin)
        local rhomax = min(rumax,rvmax)
		local r0 = runiform(`rhomin',`rhomax')
		//local r0 = 0.3
		local r1 = (`r0'-`rhomin')/(`rhomax'-`rhomin')
		mat b0=b0, ln(`r1'/(1-`r1')),ln(`r1'/(1-`r1'))
	}
    if ("`cost'"!=""){
        mata: _cost = -1
    } 
    else{
        mata: _cost = 1
    }
    ereturn clear
	local modeltype = cond("`wxvars'"=="","uv-SAR","xuv-SAR")
	local title Spatial frontier model(`modeltype')
	ml model d0 spxtsfauv() (frontier:`yvar' = `xvars' `wxvars2',`noconstant')(uhet: `uhet',noconstant) ///
	                         /lnsigma2_u  /lnsigma2_v  (Wu:)  (Wv:), nopreserve `cns' `mlmodelopt' title(`title')
	
	
	if("`initial'"=="" & "`delve'"!="") { 
		ml init b0,copy
		//ml init Wu:_cons=0.2
		//ml init Wv:_cons=0.2
	}
	if ("`initial'"=="" & "`delve'"=="") `nolog' ml search, `mlsearch'
	if ("`initial'"!="") ml init `initial', copy
	if ("`mlplot'"!=""){
		if "`'nograph'"!="" set graphics off
		`nolog' ml plot Wu:_cons
		`nolog' ml plot Wv:_cons
		if "`'nograph'"!="" set graphics on
	}
   
   `nolog' ml max, `mlmaxopt' noout 
   ereturn local cmd spxtsfauv
   ereturn local cmdbase ml
   ereturn local cmdline `cmdline'  
   Replay , `diopts'

   if(`"`te'"'!=""){
		tempname bml
		mat `bml' = e(b)
		mata: _b_ml = st_matrix("`bml'")	
	    local nx: word count `xvars' `wxvars2'
		local nz: word count `uhet'
		if("`noconstant'"=="") local noconstant constant
		mata:_te_order=spxtsfauv_te(_b_ml,`nx',`nz',"`yvar'","`xvars' `wxvars2'","`uhet'","`noconstant'")
   }	

  restore
  
  	if(`"`wxvars'"'!=""&"`genwxvars'"!=""){
      foreach v in `wxvars'{
        qui gen double W_`v' = .
        label var W_`v' `"W*`v'"'
      }
	  mata: getdatafmata(_order_wx,_order_0,"`wxvars2'")
      cap mata mata drop  _order_wx

	}

   if(`"`te'"'!=""){
		qui gen double `te' = .
		label var `te' "technical efficiency"
		mata: getdatafmata(_te_order,_order_0,"`te'")
		cap mata mata drop  _te_order		
   }  

   	if(`nummissing'>0){
		di "Missing values found"
		di "The regression sample recorded by variable __e_sample__"
		cap drop __e_sample__
		qui cap gen byte __e_sample__ = 0
		label var __e_sample__ "e(sample)"
		mata: getdatafmata(J(length(_order_0),1,1),_order_0,"__e_sample__")
		//cap mata mata drop  _touse		
	}	    	

   cap mata mata drop _order_0


end


///////////////////////subprograms//////////////////
//include spfrontier_aug.ado
cap program drop Replay
program Replay
	syntax [, Level(cilevel) * ]
	ml display , level(`level')	`options'		///
		diparm(lnsigma2_u , exp prob label("sigma2_u"))	///
		diparm(lnsigma2_v , exp prob label("sigma2_v"))  ///
		diparm(Wu, label(tau) prob function($rumin/(1+exp(@))+$rumax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rumax-$rumin)/(1+exp(@))^2))) ///
	    diparm(Wv, label(gamma) prob function($rvmin/(1+exp(@))+$rvmax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rvmax-$rvmin)/(1+exp(@))^2)))
end


/////////////////////////////
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



