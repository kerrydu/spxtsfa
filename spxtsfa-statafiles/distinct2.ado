
cap program drop distinct2
program define distinct2, rclass

syntax varname [if] [in],[CMissing]

tempname mv 
cap confirm string var `varlist'

if _rc {
	mata: `mv'=st_data(.,"`varlist'")
	if "`cmissing'"==""{
		mata: `mv' =select(`mv',`mv':!=.)
	}
}
else{
	mata: `mv'=st_sdata(.,"`varlist'")
	if "`cmissing'"==""{
		mata: `mv' =select(`mv',`mv':!="")
	}	
}


mata: `mv' = uniqrows(`mv')

mata: st_numscalar("r(ndistinct)",length(`mv'))
mata mata drop `mv'

return scalar ndistinct = r(ndistinct)

di
di "# of distinct values in {`varlist'} is `r(ndistinct)'"
if "`cmissing'"==""{
di "Note: missing values are not counted."	
}
else{
	di "Note: missing values are counted as one distinct value"
}


end