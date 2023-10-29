* translog and UCLA Statistical Consulting Group's "graph2tex" packge should be installed in advance
capture log close
cap mkdir output 
cap cd ./data
adopath ++ ../ado
clear all
sjlog using ../output/xtsfsp_DGP1, replace
use xtsfsp_ex1.dta
xtset id t 
* importing spatial weight matrix from xtsfsp_wmat1.mmat
mata mata matuse xtsfsp_w1.mmat,replace
* fitting the model
xtsfsp y x, uhet(z) wu(w1,mata) wy(w1,mata) wv(w1,mata) wx(w1,mata) wxvars(x) te(te) nolog
sjlog close, replace


clear all
sjlog using ../output/xtsfsp_DGP2, replace
* importing spatial weight matrices from xtsfsp_w2.mmat
mata mata matuse xtsfsp_w2.mmat,replace
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
use xtsfsp_ex2.dta
xtset id t 
* initial values for estimated parameters
mat b=(2,0.5,1,-1.5,4,-1.5,0.6,0.6)
* fitting the model
xtsfsp y x, cost uhet(z) wu(w2,mata) wv(w1,mata) wxvars(x) wx(`w',mata) init(b) genwxvars nolog
sjlog close, replace


sjlog using ../output/xtsfsp_DGP2a, replace
* replace some observations of y to be missing
replace y=. if _n==1 | _n==100
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
* estimation is aborted
cap noi xtsfsp y x, cost uhet(z) wu(w2,mata) wv(w1,mata) wxvars(x) wx(`w',mata) init(b) nolog
sjlog close, replace

sjlog using ../output/xtsfsp_DGP2b, replace
* re-estimation with delmissing option 
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
xtsfsp y x, cost uhet(z) wu(w2,mata) wv(w1,mata) wxvars(x) wx(`w',mata) init(b) delmissing nolog
sjlog close, replace


clear all
sjlog using ../output/xtsfsp_DGP3, replace
* importing spatial weight matrices from xtsfsp_w3.mmat
mata mata matuse xtsfsp_w3,replace
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
use xtsfsp_ex3.dta
xtset id t
xtsfsp y x, wu(`w',mata) wv(`w',mata) uhet(z) vhet(d) nolog
sjlog close, replace


* translog command should be installed; ssc install translog
capture log close
clear all
sjlog using ../output/xtsfsp_ex4, replace
* Translate shapefile to Stata format
cap spshape2dta province

use province 
drop if _ID == 26 | _ID>31
spset 
* Create spatial contiguity matrix
spmatrix create contiguity w_con, normalize(none) 
* Obtain spatial matrix as Mata matrix wm from w_con 
spmatrix matafromsp wm id = w_con
* Match the iland (_ID = 21) with the nearest province(_ID =19)
mata: wm[19,21]=1
mata: wm[21,19]=1
* Create spmatrix w_con from Mata matrix wm and rwo-normalized the matrix 
spmatrix spfrommata w_con= wm id, normalize(row) replace
* Obtain the new spatial matrix as Mata matrix wm from w_con 
spmatrix matafromsp wm id = w_con

use chnempirical.dta,clear
* Generate varables for the translog function
qui translog Y K L , time(year) norm 
global x  lnK lnL _t lnK_lnL _t_lnK _t_lnL _t_2 lnK_2 lnL_2
global z fiscal trade fdi
* Fit the model with frontier command
frontier lnY $x,uhet($z) nolog
* Predict the inefficiency term and efficiency scores
predict double uhat, u
gen double te0 = exp(-u) 
* Store the estimated parameters
mat b0=e(b)

* Fit the model with xtsfsp command
xtset _ID year
mat b1 = b0,0.6,0.6,0.6
xtsfsp lnY $x, uhet($z) wy(wm,mata) wu(wm,mata) wv(wm,mata) init(b1) te(tesp1) nolog 
scalar loglikehood1 =  e(ll)
mat b1 = b0,0.6
xtsfsp lnY $x, uhet($z)  wu(wm,mata)  init(b1) te(tesp2)  nolog
scalar loglikehood2 =  e(ll)
local lrtest = -2*(loglikehood2-loglikehood1)
local pvalue = 1- chi2(2,`lrtest')
display "Likelihood-ratio test: LR chi2(2) = `lrtest', Prob > chi2 = `pvalue'"

* Plot the density of estimates of technical efficiency from different models
twoway (kdensity te0, color(black) lpattern(solid)) ///
       (kdensity tesp1,color(red) lpattern(dash))   ///
	   (kdensity tesp2,color(blue) lpattern(longdash)),  ///
	legend(pos(10) ring(0) label(1 Non-spatial Stoc. Frontier) ///
	label(2 Spatial Stoc. Frontier:yuv) label(3 Spatial Stoc. Frontier:u)) ///
	xtitle("Technical Efficiency") ytitle("Density")
graph2tex, epsfile("../output/fig1") caption(Distribution of efficency scores) label(fig1)

sjlog close, replace 
