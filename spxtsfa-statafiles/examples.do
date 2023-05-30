capture log close
clear all
sjlog using spsfa_DGP1, replace
use spxtsfa_DGP1.dta
xtset id t 
* importing spatial weight matrix from spxtsfa_wmat1.mmat
mata mata matuse spxtsfa_wmat1.mmat,replace
* fitting the model
spxtsfa y x, uhet(z) noconstant  wy(w1,mata) wx(w1,mata) wu(w1,mata) wv(w1,mata) wxvars(x) nolog
sjlog close, replace


clear all
sjlog using spsfa_DGP2, replace
use spxtsfa_DGP2.dta
xtset id t 
* importing spatial weight matrices from spxtsfa_wmat2.mmat
mata mata matuse spxtsfa_wmat2.mmat,replace
* fitting the model
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
spxtsfa y x, uhet(z) wu(`w',mata) wv(`w',mata) te(efficiency) nolog
sjlog close, replace



clear all
sjlog using spsfa_DGP3, replace
use spxtsfa_DGP3.dta
xtset id t 
* importing spatial weight matrices from spxtsfa_wmat2.mmat
mata mata matuse spxtsfa_wmat2.mmat,replace
* fitting the model
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
mat b = (1,1,1,1,-1,-1,0.5,0.5)
spxtsfa y x, uhet(z) wu(w2,mata) wv(w1,mata) wxvars(x) ///
             wx(`w',mata) cost init(b) genwxvars nolog
sjlog close, replace
