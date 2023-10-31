* version 1.0, 2023-10-10
* subprogram for xtsfsp 
* Case3: y=rhoWy+x+v-u; v=[I-tauW]^{-1}v*

program define xtsfspyv, eclass sortpreserve
version 16

	if replay() {
		if (`"`e(cmd)'"' != "xtsfayv") error 301
		Replay `0'
	}
	else	Estimate `0'
end

program Estimate, eclass sortpreserve

syntax syntax varlist,  Uhet(string) wv(string) 
                             wy(string) [INItial(name) NOCONstant ///
                              NORMalize(string)  te(name) mldisplay(string) ///
                              DELmissing MLPLOT NOGraph MLMODELopt(string) ///
							  level(real 95) COST wxvars(varlist) ///
							  MLSEarch(string) MLMAXopt(string) DELVE ///
							  CONSTraints(string)  wx(string) ///
							  LNDETFULL lndetmc(numlist >0 min=2 max=2) GENWXVARS NOLOG Vhet(string)] 

//marksample touse 
if ("`nolog'"!="") local nolog qui
local cmdline xtsfsp `0'
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



parsespmat0 `wv' 
parsespmat1 `wv' `r(ldot)' aname(wv_ina)
local nwv = r(nw)

// check wether data macthes with sp matrix
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
	gettoken uhet0 : uhet, p(,)
	gettoken vhet0 : vhet, p(,)
    qui keep `varlist' `wxvars' `id' `time' `uhet0' `touse' `vhet0'
    tempvar order0
    qui gen int `order0' =_n
// sort data	
	qui issorted `time' `id'	
	//tempvar time2

	qui distinct2 `time'
	local T = r(ndistinct)	
    //mata mata describe
	//mata: marksuse = st_data(.,"`touse'")
	qui count if `touse'==0
	local nummissing = r(N)
	if(`nummissing'>0){
		mata: marksuse = st_data(.,"`touse'")
	}	
    mata: _order_0   = st_data(.,"`order0'","`touse'") // record the row# in the original data	

	if (`nwv'!=1 & `nwv'!=`T') {
		di as error "Spatial weight matrixs in wu() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}    
	tempvar time2
	qui egen `time2' = group(`time')
	//global paneltvar `time2'
	//mata mata describe
	checkspmat wv_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
    scalar rvmin = max(-0.9999,r(min_wv_ina))
	scalar rvmax = min(0.9999,r(max_wv_ina))
	global rvmin = rvmin
	global rvmax = rvmax
   if(`"`wy'"'==`"`wv'"'){
	 mata: wy_ina = &(wv_ina)
	 mata: wyeq = 2
	 local nwy = `nwv'
     scalar rymin = rvmin
	 scalar rymax = rvmax
	 global rymin = rymin
	 global rymax = rymax

	 	if ("`lndetfull'"!=""){
			local bp bp 
			mata: _rho_lndet_ = panlndetfull(wv_ina,$rymin,$rymax,`T')
		}
		if ("`lndetmc'"!=""){
			local bp bp 
			tokenize `lndetmc'
			mata: _rho_lndet_ = panlndetmc(`1',`2',wv_ina,$rymin,$rymax,`T')
		}
	 
	 
   } 
   else{
	mata: wyeq = 0
	parsespmat0 `wy' 
	parsespmat1 `wy' `r(ldot)' aname(wy_ina)
	local nwy = r(nw)
	if (`nwy'!=1 & `nwy'!=`T') {
		di as error "Spatial weight matrixs in wy() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}  	
    checkspmat wy_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')       
	scalar rymin = max(-0.9999,r(min_wy_ina))
	scalar rymax = min(0.9999,r(max_wy_ina))
	global rymin = rymin
	global rymax = rymax
	
		if ("`lndetfull'"!=""){
			local bp bp 
			mata: _rho_lndet_ = panlndetfull(wy_ina,$rymin,$rymax,`T')
		}
		if ("`lndetmc'"!=""){
			local bp bp 
			tokenize `lndetmc'
			mata: _rho_lndet_ = panlndetmc(`1',`2',wy_ina,$rymin,$rymax,`T')
		}
	 
   }
**************
	if ("`wx'"=="`wy'"){
		local wxwx wy_ina
	}
	else if ("`wx'"=="`wv'"){
		local wxwx wv_ina
	}
	else{
		local wxwx wx_ina
	}

    if ("`wx'"!="" & `"`wxvars'"'!=""& "`wxwx'"=="wx_ina" ){
        parsespmat0 `wx' 
        parsespmat1 `wx' `r(ldot)' aname(wx_ina)
        local nwx = r(nw)
        if (`nwx'!=1 & `nwx'!=`T') {
            di as error "Spatial weight matrixs in wx() are specified as time-varying, but # of spmatrix != # of periods"
            exit 198
        }  
        qui checkspmat wx_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
    }

****************
    qui keep if `touse'
    mata: _pan_tvar =st_data( .,"`time2'")
   * generating Wx
	if(`"`wxvars'"'!=""){
      qui genwxvars `wxvars', aname(`wxwx') tvar(`time2')
      local wxvars2  `r(wxnames)'
      mata: _order_wx = st_data(.,"`wxvars2'","`touse'")
	  //cap mata mata drop wx_ina
	}	
************

	if("`initial'"=="" & "`delve'"!="") {
		qui genwxvars `yvar', aname(wy_ina) tvar(`time2')
	    local wyvar `r(wxnames)'
		qui corr `yvar' `wyvar'
		local rhoy = r(rho)			 
		qui frontier `yvar' `xvars' `wxvars2',`noconstant' uhet(`uhet') iterate(50) `cns'  vhet(`vhet')
	    mat b0 =e(b)
		qui predict double xbfron
		qui gen double ehfron = `yvar'-xbfron
		qui genwxvars ehfron, aname(wv_ina) tvar(`time2')
	    local wehfron `r(wxnames)'
		qui corr ehfron `wehfron'
		local rhouv = r(rho)
		mat b0=b0, ln(`rhoy'/(1-`rhoy')),ln(`rhouv'/(1-`rhouv'))		
		//local r0 = runiform(rymin,rymax)
		//local r0 = 0.3
		//local r1 = (`r0'-rymin)/(rymax-rymin)
		//local r0 = runiform(rymin,rymax)
		//local r0 = 0.3
		//local r2 = (`r0'-rvmin)/(rvmax-rvmin)		
		//mat b0=b0, ln(`r1'/(1-`r1')),ln(`r2'/(1-`r2'))
	}

	if ("`cost'"!=""){
		mata: _cost = -1
	} 
	else{
		mata: _cost = 1
	}	

	local modeltype = cond("`wxvars'"=="","yv-SAR","yxv-SAR")
	local title Spatial frontier model(`modeltype')	
    if `"`vhet'"'=="" {
		local vterm /lnsigv2
	}
	else{
		local vterm (lnsigv2: `vhet') 
	}	
	
//mata mata describe
	ml model d0 xtsfspyv`bp'() (frontier:`yvar' = `xvars' `wxvars2',`noconstant') ///
	                          `vterm'   (uhet: `uhet') ///
	                        (Wy:) (Wv:), nopreserve `cns' `mlmodelopt' title(`title')
	
	
	if("`initial'"=="" & "`delve'"!="") { 
		ml init b0,copy
	}
	if ("`initial'"=="" & "`delve'"==""){
		`nolog' display "searching initial values..."
		`nolog' ml search, `mlsearch'
	} 
	if ("`initial'"!="") ml init `initial', copy
	if ("`mlplot'"!=""){
		`nolog' display "mlplot initial values..."
		if "`nograph'"!="" set graphics off
		`nolog' cap ml plot /lnsigv2
		`nolog' cap ml plot lnsigv2:_cons
		`nolog' ml plot uhet:_cons		
		`nolog' ml plot Wy:_cons
		`nolog' ml plot Wv:_cons
		if "`nograph'"!="" set graphics on
	}
 	local mlmaxopt `mlmaxopt' noout difficult
	local mlmaxopt: list uniq mlmaxopt  
	di "ml maximize..."
   `nolog' ml max, `mlmaxopt' 

   ereturn local cmd xtsfspyv
   ereturn local cmdbase ml
   ereturn local cmdline `cmdline'
   ereturn local wy wy_ina 
   ereturn local wv wv_ina 

   Replay , `diopts' 
   if `"`wxvars'"'!="" di "     W_(`wxvars') represent Spatial Durbin terms W(`wxvars')"
   
    if(`"`te'"'!=""){
		tempname bml
		mat `bml' = e(b)
		mata: _b_ml = st_matrix("`bml'")	
	    local nx: word count `xvars' `wxvars2'
		//local nz: word count `uhet'
		if("`noconstant'"=="") local noconstant constant
		mata:_te_order=xtsfspyv_te(_b_ml,`nx',"`yvar'","`xvars' `wxvars2'","`uhet'","`vhet'","`noconstant'")
   }	

  restore
  
  	if(`"`wxvars'"'!=""&"`genwxvars'"!=""){
      foreach v in `wxvars'{
        qui gen double W_`v' = .
        label var W_`v' `"W*`v'"'
      }
	  mata: getdatafmata(_order_wx,_order_0,"`wxvars2'")
      cap mata mata drop  _order_wx
	  ereturn local wx `wxwx'

	}

   if(`"`te'"'!=""){
		qui gen double `te' = .
		label var `te' "technical efficiency"
		mata: getdatafmata(_te_order,_order_0,"`te'")
		cap mata mata drop  _te_order		
   }  	
   	if(`nummissing'>0){
		di "      Missing values found"
		di "      The regression sample recorded by variable __e_sample__"
		cap drop __e_sample__
		qui cap gen byte __e_sample__ = 0
		label var __e_sample__ "e(sample)"
		mata: getdatafmata(J(length(_order_0),1,1),_order_0,"__e_sample__")
		//cap mata mata drop  _touse		
	}	 
   cap mata mata drop _order_0
end


///////////////////////subprograms//////////////////
//include spxtsfa_aug.ado
cap program drop Replay
program Replay
	syntax [, Level(cilevel) *]
	ml display , level(`level')	`options'		///
		diparm(Wy, label(rho) prob function($rymin/(1+exp(@))+$rymax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rymax-$rymin)/(1+exp(@))^2))) ///  
	   	diparm(Wv, label(gamma) prob function($rvmin/(1+exp(@))+$rvmax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rvmax-$rvmin)/(1+exp(@))^2)))  
	di "Note: Wy:_cons and Wv:_cons are the transformed parameters;"
	di "      rho and gamma  are their origin metrics in spatial components, respectively."	   
	   
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
* This command is browed from Morad Zekhnini's nwxtregress package
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




