// 25 Oct 2023 considering vhet function
//Mrho = matinv(I(N)-rho*spwi) search and replace
//version 1.0 
// 14 Oct 2023
// change function name and parameter order
// important: wy-> gamma  wv->rho  wu->tau 
// order: frontier lnsigmav^2  uhet lnsigmau^2 wy wv wu
//version 0.2
//2023-09-08
////////////////////////////////////////////////////////
capture mata mata drop _genwxvars() 
capture mata mata drop loaddata()
capture mata mata drop extrpoi()
capture mata mata drop _checkspmat()
capture mata mata drop MaxEV()
capture mata mata drop MinMaxRC()
capture mata mata drop Wnorm()
capture mata mata drop Wnormi()
capture mata mata drop studentize()
capture mata mata drop _studentizeT()
capture mata mata drop getdatafmata()
capture mata mata drop xtsfspyuv()
capture mata mata drop xtsfspyuv_te()
capture mata mata drop xtsfspyuv0()
capture mata mata drop xtsfspyuv0_te()
capture mata mata drop xtsfspyu()
capture mata mata drop xtsfspyu_te()
capture mata mata drop xtsfspyv()
capture mata mata drop xtsfspyv_te()
capture mata mata drop xtsfspy()
capture mata mata drop xtsfspy_te()
capture mata mata drop xtsfspu()
capture mata mata drop xtsfspu_te()
capture mata mata drop xtsfspv()
capture mata mata drop xtsfspv_te()
capture mata mata drop xtsfspuv()
capture mata mata drop xtsfspuv_te()
capture mata mata drop xtsfspuv0()
capture mata mata drop xtsfspuv0_te()
capture mata mata drop myrowsum()
capture mata mata drop mycolsum()
capture mata mata drop myquadsum()
capture mata mata drop myinvsym()
capture mata mata drop matinv()
capture mata mata drop cholqrpinv()
capture mata mata drop lndetmc()
capture mata mata drop panlndetmc()
capture mata mata drop panlndetfull()
capture mata mata drop lndetfull()
capture mata mata drop xtsfspyuvbp()
capture mata mata drop xtsfspyuv0bp()
capture mata mata drop xtsfspyubp()
capture mata mata drop xtsfspyvbp()
capture mata mata drop xtsfspybp()
capture mata mata drop xtsfspubp()
capture mata mata drop xtsfspvbp()
capture mata mata drop xtsfspuvbp()
capture mata mata drop xtsfspuv0bp()
capture mata mata drop pifun()
capture mata mata drop tvpifun()

mata:

void function pifun(real vector s2_v,
                    real matrix Mr, 
					real matrix iMr, 
					real matrix lndetPi, 
					real matrix invPi)
{
	if (length(s2_v)==1){
		Pi = s2_v*(quadcross(iMr,iMr))' //cross
		invPi = 1/s2_v*(quadcross(Mr,Mr))
	    lndetPi = ln(det(Pi))
	}
}

void function tvpifun(real vector s2_v,
                    real matrix Mr, 
					real matrix iMr, 
					real matrix lndetPi, 
					real matrix invPi,
					real scalar i,
					real matrix info)
{
	if (length(s2_v)>1){
		invPi = quadcross(Mr,1:/panelsubmatrix(s2_v, i, info),Mr)
		Pi = quadcross(iMr,panelsubmatrix(s2_v, i, info),iMr)' //cross
		lndetPi = ln(det(Pi))    					
	}
}



real matrix function loaddata(real colvector r, string scalar v,real scalar len)
{
	if(v==""){
		return(J(len,1,1))
		exit()
	}
	
	vv = tokens(v,",")
	//vv
	if(vv[length(vv)]=="noc" | vv[length(vv)]=="noconstant" | vv[length(vv)]=="nocon" | vv[length(vv)]=="noconst"){
		return(st_data(r,vv[1]))
	}
	else{
		if(vv[1]!=","){
		   data = st_data(r,vv[1])
		   data = data, J(rows(data),1,1)
		}
		else{
			data = J(len,1,1)
		}

		return(data)
	}
}


function myinvsym(  numeric matrix A)
{
	return(invsym(A))
}

function matinv(numeric matrix A, | real scalar issym)
{
    if (args()==1) {
        issym = issymmetric(A)  
    }
    if (issym==1){
        return(invsym(A))
    }
    else{
        return(cholqrpinv(A))
    }
    
}

	function cholqrpinv (  numeric matrix A,
						  | real scalar useqr)
	{
			if (args()==2) useqr = 0
			real matrix C
			if (!useqr) {
					C = cholinv(A)
					if (C[1,1]==.) {
							C = qrinv(A)
						if (C[1,1]==.){
							C = pinv(A)
						}
					}
			}
			else {
					C = qrinv(A)
			}
			return(C)
	}

//////////////////////////////////////////////////////////////////
function lndetmc(real scalar order, real scalar iter, transmorphic matrix wsw, real scalar rmin, real scalar rmax)
 { 
 
 n = rows(wsw)
 //quadsum(wsw)
 td =0 \ quadsum(wsw:^2)/2
 oexact = length(td)

   mavmomi = J(order, iter, 0)
   for (j = 1; j <= iter; j++) { 
    u = rnormal(n, 1, 0, 1)
    v = u
    utu = u' * u 
    for (i = 1; i <= order; i++) { 
        v = wsw * v
        mavmomi[i, j] = n * ((u' * v) / (i * utu))
        } 
    }
	
     mavmomi[1..oexact,.] = J(1,cols(mavmomi),td) 
	 avmomi = (rowsum(mavmomi)/cols(mavmomi))
     alpha = range(rmin,rmax,0.01)
     valphaf = Vandermonde(alpha)
     alomat = -valphaf[., 2..order+1]
     lndetmat = alomat * avmomi

     out = alpha,lndetmat
     return(out)
    
 }


///////////////////////////////////////////////////////////////
function panlndetmc(real scalar order, 
	                real scalar iter, 
	                transmorphic matrix wsw, 
	                real scalar rmin, 
	                real scalar rmax, 
	                real scalar T)
{
	if(eltype(wsw)=="real"){
		out = lndetmc(order,iter,wsw,rmin,rmax)
		out=out[.,1],out[.,2]*T
	}
	else{
		keys = sort(asarray_keys(wsw),1)
		//keys = sort(keys,1)
		if(length(keys)==1){
			wi = extrpoi(asarray(wsw,1))
			out = lndetmc(order,iter,wi,rmin,rmax)
			out=out[.,1],out[.,2]*T			
		}
		else{
			wi = extrpoi(asarray(wsw,1))
			out = lndetmc(order,iter,wi,rmin,rmax)
			out2 = out
			for(t=2;t<=T;t++){
				wi = extrpoi(asarray(wsw,1))
				out = lndetmc(order,iter,wi,rmin,rmax)	
				out2[.,2] =out2[.,2]+out[.,2]
			}
			out=out2
		}
	}
	
	fgrid = range(rmin,rmax,0.001)
	lndet = spline3eval(spline3(out[.,1],out[.,2]),fgrid)
	out = fgrid,lndet
	return(out)
	
}
/////////////////////////////////////////////////////////////////

function panlndetfull(transmorphic matrix wsw, real scalar rmin, real scalar rmax, real scalar T)
{
	if(eltype(wsw)=="real"){
		out = lndetfull(wsw,rmin,rmax)
		out=out[.,1],out[.,2]*T
	}
	else{
		keys =sort( asarray_keys(wsw),1)
		//keys = sort(keys,1)
		if(length(keys)==1){
			wi = extrpoi(asarray(wsw,1))
			out = lndetfull(wi,rmin,rmax)
			out=out[.,1],out[.,2]*T			
		}
		else{
			wi = extrpoi(asarray(wsw,1))
			out = lndetfull(wi,rmin,rmax)
			out2 = out
			for(t=2;t<=T;t++){
				wi = extrpoi(asarray(wsw,1))
				out = lndetfull(wi,rmin,rmax)	
				out2[.,2] =out2[.,2]+out[.,2]
			}
			out=out2
		}
	}
	
	fgrid = range(rmin,rmax,0.001)
	lndet = spline3eval(spline3(out[.,1],out[.,2]),fgrid)
	out = fgrid,lndet
	return(out)
	
}

///////////////////////////////////////////////

function lndetfull(transmorphic matrix wsw, real scalar rmin, real scalar rmax)
 { 
 
 n = rows(wsw)
 
 rho = range(rmin,rmax,0.01)
 
 lndet =J(length(rho),1,.)
 
 for(i=1;i<=length(rho);i++){
 	
	 B = I(n) - rho[i]*wsw 
	 lud(B,L=0,U=0,p=0)
	lu1 = L+U-I(rows(B))
	lu1 = diagonal(lu1)
	lndet[i]= sum(log(lu1))
 }

     out = rho,lndet
     return(out)
    
 }



////////////////////////////_genwxvars////////////////////////////
void function _genwxvars(string scalar vars,transmorphic matrix arr, string scalar tvar)
{
	xnames = tokens(vars)
	data = st_data(.,vars)
	t = st_data(.,tvar)
	uniqt = uniqrows(t)
	keys = sort(asarray_keys(arr),1)
	for(i=1;i<=cols(data);i++){
		for(j=1;j<=length(uniqt);j++){
			ind=mm_which(t:==j)
			if(length(keys)>1){
				data[ind,i]=extrpoi(asarray(arr,keys[j]))*data[ind,i]
			}
			else{
				data[ind,i]=extrpoi(asarray(arr,keys[1]))*data[ind,i]
			}
			
		}
		xnames[i] = "W_"+xnames[i] 
		xxx=st_addvar("double",xnames[i] )
		st_store(.,xxx,data[.,i])
	}
	st_local("wxnames",invtokens(xnames))
	
	//st_view(xxx=.,.,invtokens(xnames))
	//xxx[.,.]=data
	
}


////////////////////////////extrpoi////////////////////////////


real matrix function extrpoi (transmorphic matrix w)
{
	if(eltype(w)=="pointer"){
		return(*w)
	}
	else{
		return(w)
	}
}


////////////////////////////_checkspmat////////////////////////////
/*
function _checkspmat(string scalar ttouse, transmorphic matrix wina,real scalar normalize)
{
	wkeys = sort(asarray_keys(wina),1)
	//wkeys = sort(wkeys,1)
	nw = length(wkeys)
	data = st_data(.,ttouse)
	T = max(data[.,1])
	n0 = sum(data[.,2]:==0)
	//normalize = st_numscalar("ntype")
	lmax = J(nw,1,.)
	lmin = J(nw,1,.)	
    if(n0==0){  // no missing values
		if(nw==1){
			if(normalize!=0){
				asarray(wina,wkeys[1],Wnorm(asarray(wina,wkeys[1]),normalize,0))
			} 
			lamb = Re( eigenvalues(asarray(wina,wkeys[1])))
			//st_numscalar("r(max)",max(1:/lamb))
			//st_numscalar("r(min)",min(1:/lamb))
			st_numscalar("r(max)",1/max(lamb))
			st_numscalar("r(min)",1/min(lamb))
		}
		else{
			for(t=1;t<=T;t++){
				if(normalize!=0){
					asarray(wina,wkeys[t], Wnorm(asarray(wina,wkeys[t]),normalize,0))
				} 
				lamb = Re( eigenvalues(asarray(wina,wkeys[t])))
				lmax[t] = 1/max(lamb)
				lmin[t] = 1/min(lamb)
			}
			st_numscalar("r(max)",min(lmax))
			st_numscalar("r(min)",max(lmin))	
		}
	}
    else{      
        if(nw==1){
            pointer matrix pv 
            pv = J(T,1,NULL)
            wtemp1 = asarray(wina,wkeys[1])
            for(t=1;t<=T;t++){
                pv[t]=&wtemp1
            }
                for(t=1;t<=T;t++){
                    touse = select(data[.,2],data[.,1]:==t)
                    if(sum(touse:!=0)<length(touse)){ // delete units with missing values
                        ind = mm_which(touse)
                        pv[t]=&(Wnorm((*pv[t])[ind,ind],normalize,0))
                        //pv[t] = &(Wnorm(*pv[t],normalize,0))    
                        lamb = Re( eigenvalues(*pv[t]))
                        lmax[t]=1/max(lamb)
                        lmin[t]=1/min(lamb)                             
                    }
                }
                wtemp1 = Wnorm(wtemp1,normalize,0)
                lamb = Re( eigenvalues(wtemp1))
                ind = mm_which(lmin:==.)
                lmax[ind] = J(length(ind),1,1/max(lamb))
                lmin[ind] = J(length(ind),1,1/min(lamb))
                st_numscalar("r(max)",min(lmax))
                st_numscalar("r(min)",max(lmin))
                //wva = asarray_create("real")
                for(t=1;t<=T;t++){
                    asarray(wina,wkeys[t],pv[t])
                }
                //wina = wva 

        }
        else{
            for(t=1;t<=T;t++){
                touse = select(data[.,2],data[.,1]:==t)
                if(sum(touse:!=0)<length(touse)){ // delete units with missing values
                    ind = mm_which(touse)
                    wt = asarray(wina,wkeys[t])[ind,ind]
                    if(normalize!=0) wt = (wt,normalize,0)
                    asarray(wina,wkeys[t],wt) 
                    lamb = Re( eigenvalues(wt))
                    lmax[t]=1/max(lamb)
                    lmin[t]=1/min(lamb) 
                }
                else{
                    wt = asarray(wina,wkeys[t])
                    if(normalize!=0){
                        wt = (wt,normalize,0)
                        asarray(wina,wkeys[t],wt) 
                    } 
                    lamb = Re( eigenvalues(wt))
                    lmax[t]=1/max(lamb)
                    lmin[t]=1/min(lamb)              
                }
            }
            st_numscalar("r(max)",min(lmax))
            st_numscalar("r(min)",max(lmin))
        }
    }

}

*/

////////////////////////////////////////////////////////////////


function _checkspmat(string scalar ttouse, transmorphic matrix wina,real scalar normalize)
{
	wkeys = sort(asarray_keys(wina),1)
	//wkeys = sort(wkeys,1)
	nw = length(wkeys)
	data = st_data(.,ttouse)
	T = max(data[.,1])
	n0 = sum(data[.,2]:==0)
	//normalize = st_numscalar("ntype")
	lmax = J(nw,1,.)
	lmin = J(nw,1,.)	
    if(n0==0){  // no missing values
		if(nw==1){
			wmcol = cols(asarray(wina,wkeys[1]))
			maxnid = ceil(rows(data)/T)
			if(wmcol<maxnid){
				_error(3200,"The dimensions of spmatrix < # of units")
			}
			if(normalize!=0){
				asarray(wina,wkeys[1],Wnorm(asarray(wina,wkeys[1]),normalize,0))
			} 
			lamb = Re( eigenvalues(asarray(wina,wkeys[1])))
			//st_numscalar("r(max)",max(1:/lamb))
			//st_numscalar("r(min)",min(1:/lamb))
			st_numscalar("r(max)",1/max(lamb))
			st_numscalar("r(min)",1/min(lamb))
		}
		else{
			for(t=1;t<=T;t++){
				wmcol = cols(asarray(wina,wkeys[t]))
				maxnid = length(select(data[.,1],data[.,1]:==t))
				if(wmcol<maxnid){
					_error(3200,"The dimensions of spmatrix < # of units")
				}				
				if(normalize!=0){
					asarray(wina,wkeys[t], Wnorm(asarray(wina,wkeys[t]),normalize,0))
				} 
				lamb = Re( eigenvalues(asarray(wina,wkeys[t])))
				lmax[t] = 1/max(lamb)
				lmin[t] = 1/min(lamb)
			}
			st_numscalar("r(max)",min(lmax))
			st_numscalar("r(min)",max(lmin))	
		}
	}
    else{      
        if(nw==1){			
            pointer matrix pv 
            pv = J(T,1,NULL)
            wtemp1 = asarray(wina,wkeys[1])
			lmin = J(T,1,.)
			lmax = J(T,1,.)
            for(t=1;t<=T;t++){
                pv[t]=&wtemp1
            }
                for(t=1;t<=T;t++){
					wmcol = cols(*pv[t])
					maxnid = length(select(data[.,1],data[.,1]:==t))
					if(wmcol<maxnid){
						_error(3200,"The dimensions of spmatrix < # of units")
					}					
                    touse = select(data[.,2],data[.,1]:==t)
					//sum(touse:!=0)
					//length(touse)
                    if(sum(touse:!=0)<length(touse)){ // delete units with missing values
                        ind = mm_which(touse)
                        pv[t]=&(Wnorm((*pv[t])[ind,ind],normalize,0))
						//rows(Wnorm((*pv[t])[ind,ind],normalize,0))
                        //pv[t] = &(Wnorm(*pv[t],normalize,0))    
                        lamb = Re( eigenvalues(*pv[t]))
                        lmax[t]=1/max(lamb)
                        lmin[t]=1/min(lamb)                             
                    }
                }
				
                wtemp1 = Wnorm(wtemp1,normalize,0)
                lamb = Re( eigenvalues(wtemp1))
                ind = mm_which(lmin:==.)
				
                lmax[ind] = J(length(ind),1,1/max(lamb))
                lmin[ind] = J(length(ind),1,1/min(lamb))

                st_numscalar("r(max)",min(lmax))
                st_numscalar("r(min)",max(lmin))
                //wva = asarray_create("real")
				asarray(wina,wkeys[1],pv[1])
                for(t=2;t<=T;t++){
                    asarray(wina,t,pv[t])
                }
                //wina = wva 

        }
        else{
            for(t=1;t<=T;t++){
					wmcol = cols(asarray(wina,wkeys[t]))
					maxnid = length(select(data[.,1],data[.,1]:==t))
					if(wmcol<maxnid){
						_error(3200,"The dimensions of spmatrix < # of units")
					}				
                touse = select(data[.,2],data[.,1]:==t)
                if(sum(touse:!=0)<length(touse)){ // delete units with missing values
                    ind = mm_which(touse)
                    wt = asarray(wina,wkeys[t])[ind,ind]
                    if(normalize!=0) wt = (wt,normalize,0)
                    asarray(wina,wkeys[t],wt) 
                    lamb = Re( eigenvalues(wt))
                    lmax[t]=1/max(lamb)
                    lmin[t]=1/min(lamb) 
                }
                else{
                    wt = asarray(wina,wkeys[t])
                    if(normalize!=0){
                        wt = (wt,normalize,0)
                        asarray(wina,wkeys[t],wt) 
                    } 
                    lamb = Re( eigenvalues(wt))
                    lmax[t]=1/max(lamb)
                    lmin[t]=1/min(lamb)              
                }
            }
            st_numscalar("r(max)",min(lmax))
            st_numscalar("r(min)",max(lmin))
        }
    }

}


////////////////////////////MaxEV////////////////////////////


	function MaxEV(real matrix mat)
	{
		eigensystem(mat,tmp=.,ev=.)
		return(max(Re(ev)))
	}


////////////////////////////MinMaxRC////////////////////////////

	function MinMaxRC(real matrix mat)
	{
		cols = max(quadcolsum(mat))
		rows = max(quadrowsum(mat))
		return(min((cols,rows)))
	}


////////////////////////////Wnorm////////////////////////////


	function Wnorm(	 transmorphic matrix mat, 	/// 	
					real scalar type,			/// 0: nothing, 1: for row normalisation, 2: column, 3: spectral
					real scalar issparse		/// 1: sparse matrix 
					)


	{
		
		if (type > 0 ) {
			if (issparse == 1) {
				"normalise sparse matrix"
				if (type == 1)	{
					"row norm"
					loc = (1,2)
					spec = 0					
				}
				else if (type == 2) {
					"col norm"
					loc = (2,1)
					spec = 0
				}
				else if (type == 3) {
					"spectral"
					loc = (1,2)
					spec = 1
				}
				else if (type == 4) {
					"min max"
					loc = (1,2)
					spec = 2
				}

				if (eltype(mat) == "real") {
					"is real"
					mat = Wnormi(mat,loc,spec)
				} 
				else {
					"non real"
					keys = sort(asarray_keys(mat),1)
					nkeys = rows(keys)
					i = 1
					while (i<= nkeys) {
						Wi = asarray(mat,keys[i])
						
						/// remove diagonals
						Wi = Wi[selectindex(Wi[.,1]:!=Wi[.,2]),.]
						Wi = Wnormi(Wi,loc,spec)
						asarray(mat,keys[i],Wi)						
						
						i++
					}
				}

			}
			else {
				pointer(function) fn

				if (type == 1) {
					fn = &myrowsum()
				}
				else if (type == 2) {
					fn = &mycolsum()
				}
				else  if (type == 3) {
					fn = &MaxEV()
				}
				else  if (type == 4) {
					fn = &MinMaxRC()
				}

				if (eltype(mat) == "real") {
					mat = mat:/(*fn)(mat)
				}
				else {
					keys = asarray_keys(mat)
					keys = sort(keys,1)
					nkeys = rows(keys)
					i = 1
					while (i<= nkeys) {
						Wi = asarray(mat,keys[i])
						Wi = Wi :/ (*fn)(Wi)
						asarray(mat,keys[i],Wi)
						i++
					}
				}
			}			
		}
		return(mat)
	}


////////////////////////////Wnormi////////////////////////////


	function Wnormi(real matrix mat, real matrix loc,real matrix spec)
	{
		
		if (spec:==0) {
			mat = sort(mat,(loc))
			index = panelsetup(mat,loc[1])
			i = panelstats(index)[1]

			while (i>0) {			
				wi = panelsubmatrix(mat,i,index)
				ws = quadsum(wi[.,3])
				
				mat[|index[i,1],3 \ index[i,2],3|] = wi[.,3] :/ ws
				///
				i--
			}
		}
		else if (spec:==1) {
			idt = mat[,(1,2)]
			mat_nsp = SparseUndefine(mat,idt,1)
			maxev = MaxEV(mat_nsp)
			mat[.,3] = mat[.,3] :/ maxev
		}
		else if (spec:==2) {
			mat[.,3] = mat[.,3] :/ MinMaxRC(mat[.,3])
		}

		return(mat)
	}


////////////////////////////studentize////////////////////////////


	function studentize(real matrix mat)
	{
		nvar = cols(mat)
		nobs = rows(mat)

		if (nvar == 1) {
			var = quadmeanvariance(mat)
			mean = var[1]
			sd = sqrt(var[2])
	
			if (sd == 0) {
				res = J(nobs,1,0)
			}
			else {
				res = (mat :- mean) :/ sd
			}
		}
		else {
			var = quadmeanvariance(mat)
			mean = var[1,.]
			sd = sqrt(diagonal(var[(2..nvar+1),.]))
			res = (mat :- mean) :/ sd'
			///_editmissing(res,0)
		}
		return(res)
	}

////////////////////////////_studentizeT////////////////////////////


	function _studentizeT(real matrix mat, real matrix T, real matrix idt,  real matrix uniqueT)
	{
		
		t = 1
		while (t<=T) {
			index = selectindex(idt[.,2]:==uniqueT[t])
			mat[index,.] = studentize((mat[index,.]))
			t++
		}
	}


////////////////////////////getdatafmata////////////////////////////


	void function getdatafmata(real matrix x, real colvector ind, string scalar var)
	{
		st_view(y=.,.,var)
		y[ind,.]=x
	}




////////////////////////////xtsfspyuv////////////////////////////


void function xtsfspyuv(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wu_ina
    external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
    external real colvector _pan_tvar
	external real scalar _cost
    //external real scalar wueq
    external real scalar wveq
    external real scalar wyeq
	xb = moptimize_util_xb(M,b,1)
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_u=1
	// wy:->gamma ; wv:->rho ; wu:->tau
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= moptimize_util_xb(M,b,6)	// scalar
	taup = exp(taup)/(1+exp(taup))
	rmin = st_numscalar("rumin")
	rmax = st_numscalar("rumax")	
	tau  = rmin*(1-taup)+taup*rmax

	rhop	= moptimize_util_xb(M,b,5)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rvmin")
	rmax = st_numscalar("rvmax")		
	rho  = rmin*(1-rhop)+rhop*rmax	


	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	lndetPi=.
	invPi=.
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(wu_ina),1)
   // ukeys = sort(ukeys,1)
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
        if(wveq==1){
           //Mrho = matinv(I(N)-rho*spwi)
		   Mr = I(N)-rho*spwi
		   iMr = matinv( Mr)
		   pifun(s2_v,Mr,iMr,lndetPi,invPi)
        }
        if(wyeq==1){
           lndetIrhoW = ln(det(I(N)-gamma*spwi))	
        }
	}
    if(wveq==0){
        vkeys = sort(asarray_keys(wv_ina),1)
        if(length(vkeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wv_ina,vkeys[1]))		
            //Mrho = matinv(I(N)-rho*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv(Mr )
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
            if(wyeq==2){
              lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }
        }        
    }

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,vkeys[1]))
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	   
        }  
    }
      

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
            if(wveq==1){
            	Mr = I(N)-rho*spwi
				iMr = matinv(Mr)
				pifun(s2_v,Mr,iMr,lndetPi,invPi)			
            }
            if(wyeq==1){
               lndetIrhoW = ln(det(I(N)-gamma*spwi))
			//esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)	
            }            
		}
        if(wveq==0){
            if(length(vkeys)>1){
                spwi = extrpoi(asarray(wv_ina,i))
                Mr = I(N)-rho*spwi
				iMr = matinv(Mr)
				pifun(s2_v,Mr,iMr,lndetPi,invPi)              
                if(wyeq==2){
                lndetIrhoW = ln(det(I(N)-gamma*spwi))
				//esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)	
                }
            }            
        }
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))
                lndetIrhoW = ln(det(I(N)-gamma*spwi))	
				//esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)
            } 
			//else{
			//	spwi = extrpoi(asarray(wy_ina,1))
			//	esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)
			//}           
        }
       
		
		if(wyeq==1 & length(ukeys)==1){
			spwi = extrpoi(asarray(wu_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}
		else if(wyeq==1 & length(ukeys)>1){
			spwi = extrpoi(asarray(wu_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}
		else if(wyeq==2 & length(vkeys)==1){
			spwi = extrpoi(asarray(wv_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}		
		else if(wyeq==2 & length(vkeys)>1){
			spwi = extrpoi(asarray(wv_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}	
		else if(wyeq==0 & length(ykeys)==1){
			spwi = extrpoi(asarray(wy_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}	
		else {
			spwi = extrpoi(asarray(wy_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}							
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info) 

		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}

////////////////////////////////////////////////////////////////////


void function xtsfspyuvbp(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wu_ina
    external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
    external real colvector _pan_tvar
	external real scalar _cost
    //external real scalar wueq
    external real scalar wveq
    external real scalar wyeq
    external real matrix _rho_lndet_
	xb = moptimize_util_xb(M,b,1)
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	
	
	// wy:->gamma ; wv:->rho ; wu:->tau	
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= moptimize_util_xb(M,b,6)	// scalar
	taup = exp(taup)/(1+exp(taup))
	rmin = st_numscalar("rumin")
	rmax = st_numscalar("rumax")	
	tau  = rmin*(1-taup)+taup*rmax

	rhop	= moptimize_util_xb(M,b,5)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rvmin")
	rmax = st_numscalar("rvmax")		
	rho  = rmin*(1-rhop)+rhop*rmax	

	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	lndetPi=.
	invPi=.
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(wu_ina),1)
   // ukeys = sort(ukeys,1)
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
        if(wveq==1){
           //Mrho = matinv(I(N)-rho*spwi)
		   Mr = I(N)-rho*spwi
		   iMr = matinv( Mr)
		   pifun(s2_v,Mr,iMr,lndetPi,invPi)
        }
        if(wyeq==1){
           lndetIrhoW = ln(det(I(N)-gamma*spwi))	
        }
	}
    if(wveq==0){
        vkeys = sort(asarray_keys(wv_ina),1)
        if(length(vkeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wv_ina,vkeys[1]))		
            //Mrho = matinv(I(N)-rho*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv(Mr )
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
            if(wyeq==2){
              lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }
        }        
    }

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,vkeys[1]))
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	   
        }  
    }
      

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
            if(wveq==1){
            	Mr = I(N)-rho*spwi
				iMr = matinv(Mr)
				pifun(s2_v,Mr,iMr,lndetPi,invPi)			
            }
            if(wyeq==1){
               lndetIrhoW = ln(det(I(N)-gamma*spwi))
			//esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)	
            }            
		}
        if(wveq==0){
            if(length(vkeys)>1){
                spwi = extrpoi(asarray(wv_ina,i))
                Mr = I(N)-rho*spwi
				iMr = matinv(Mr)
				pifun(s2_v,Mr,iMr,lndetPi,invPi)              
                if(wyeq==2){
                lndetIrhoW = ln(det(I(N)-gamma*spwi))
				//esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)	
                }
            }            
        }
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))
                lndetIrhoW = ln(det(I(N)-gamma*spwi))	
				//esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)
            } 
			//else{
			//	spwi = extrpoi(asarray(wy_ina,1))
			//	esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info)
			//}           
        }
       
		
		if(wyeq==1 & length(ukeys)==1){
			spwi = extrpoi(asarray(wu_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}
		else if(wyeq==1 & length(ukeys)>1){
			spwi = extrpoi(asarray(wu_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}
		else if(wyeq==2 & length(vkeys)==1){
			spwi = extrpoi(asarray(wv_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}		
		else if(wyeq==2 & length(vkeys)>1){
			spwi = extrpoi(asarray(wv_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}	
		else if(wyeq==0 & length(ykeys)==1){
			spwi = extrpoi(asarray(wy_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}	
		else {
			spwi = extrpoi(asarray(wy_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		}							
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info) 

		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}
	lnf = moptimize_util_sum(M, lnfj)
	
}




////////////////////////////xtsfspyuv_te////////////////////////////




    function xtsfspyuv_te(   real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix wu_ina
    external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
    external real colvector _pan_tvar
	external real scalar _cost
    //external real scalar wueq
    external real scalar wveq
    external real scalar wyeq
	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1

	
   // wy:->gamma    wv:->rho  wu:->tau
	gammap	= b[nx+nv+nz+1]	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= b[nx+nv+nz+3]	// scalar
	taup = exp(taup)/(1+exp(taup))
	rmin = st_numscalar("rumin")
	rmax = st_numscalar("rumax")	
	tau  = rmin*(1-taup)+taup*rmax

	rhop	= b[nx+nv+nz+2]	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rvmin")
	rmax = st_numscalar("rvmax")		
	rho  = rmin*(1-rhop)+rhop*rmax

	
	es=_cost*(yy-xb)
	Eie = J(0,1,.)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	//lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(wu_ina),1)
    lndetPi=.
	invPi=.
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
        if(wveq==1){
           Mr = I(N)-rho*spwi
		   iMr = matinv(Mr)
		   pifun(s2_v,Mr,iMr,lndetPi,invPi)
        }
	}
    if(wveq==0){
        vkeys = sort(asarray_keys(wv_ina),1)
        if(length(vkeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wv_ina,vkeys[1]))		
			Mr = I(N)-rho*spwi
			iMr = matinv(Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
			
        }        
    }

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1) 
    }

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
            if(wveq==1){
				Mr = I(N)-rho*spwi
				iMr = matinv(Mr)
				pifun(s2_v,Mr,iMr,lndetPi,invPi)		
            }           
		}
        if(wveq==0){
            if(length(vkeys)>1){
                spwi = extrpoi(asarray(wv_ina,i))
				Mr = I(N)-rho*spwi
				iMr = matinv(Mr)
				pifun(s2_v,Mr,iMr,lndetPi,invPi)
            }            
        }
   
       
		
		if(wyeq==1 & length(ukeys)==1){
			spwi = extrpoi(asarray(wu_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		}
		else if(wyeq==1 & length(ukeys)>1){
			spwi = extrpoi(asarray(wu_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		}
		else if(wyeq==2 & length(vkeys)==1){
			spwi = extrpoi(asarray(wv_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		}		
		else if(wyeq==2 & length(vkeys)>1){
			spwi = extrpoi(asarray(wv_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		}	
		else if(wyeq==0 & length(ykeys)==1){
			spwi = extrpoi(asarray(wy_ina,1))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		}	
		else {
			spwi = extrpoi(asarray(wy_ina,i))
			esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		}
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi)

		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
	}
	return(exp(-Eie))

}




////////////////////////////xtsfspyuv0////////////////////////////




void function xtsfspyuv0(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
                        //real colvector lnfj,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar

	xb = moptimize_util_xb(M,b,1)
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_u=1
	
   // wy:->gamma    wv:->rho  wu:->tau	
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= moptimize_util_xb(M,b,6)	// scalar
	taup = exp(taup)/(1+exp(taup))
	tau  = rmin*(1-taup)+taup*rmax


	rhop	= moptimize_util_xb(M,b,5)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rho  = rmin*(1-rhop)+rhop*rmax

	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(w_ina),1)
    //vkeys = asarray_keys(wv_ina)
    lndetPi=.
	invPi=.
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
		//Mrho = matinv(I(N)-rho*spwi)
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
		lndetIrhoW = ln(det(I(N)-gamma*spwi))	
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		    lndetIrhoW = ln(det(I(N)-gamma*spwi))

		}
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)


		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}

//////////////////////////////////////////////////////////////////////



void function xtsfspyuv0bp(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
                        //real colvector lnfj,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar
    external real matrix _rho_lndet_

	xb = moptimize_util_xb(M,b,1)
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_u=1
	
  // wy:->gamma    wv:->rho  wu:->tau		
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= moptimize_util_xb(M,b,6)	// scalar
	taup = exp(taup)/(1+exp(taup))
	tau  = rmin*(1-taup)+taup*rmax


	rhop	= moptimize_util_xb(M,b,5)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rho  = rmin*(1-rhop)+rhop*rmax

	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

	lmax = max(mm_which(_rho_lndet_[.,1]:<=gamma))
	if(length(lmax)==0) lmax=1
	gmin = min(mm_which(_rho_lndet_[.,1]:>=gamma))
	if(length(gmin)==0) gmin=length(_rho_lndet_[.,1])
	lndetIrhoW=_rho_lndet_[round(0.5*lmax+0.5*gmin),2]/nt

    ukeys = sort(asarray_keys(w_ina),1)
    //vkeys = asarray_keys(wv_ina)
    lndetPi=.
	invPi=.
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
		//Mrho = matinv(I(N)-rho*spwi)
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
		lndetIrhoW = ln(det(I(N)-gamma*spwi))	
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		    lndetIrhoW = ln(det(I(N)-gamma*spwi))

		}
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)


		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}


////////////////////////////xtsfspyuv0_te////////////////////////////


     function xtsfspyuv0_te( real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar

	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1

 // wy:->gamma    wv:->rho  wu:->tau
	gammap	= b[nx+nv+nz+1]	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= b[nx+nv+nz+3]	// scalar
	taup = exp(taup)/(1+exp(taup))	
	tau  = rmin*(1-taup)+taup*rmax

	rhop	= b[nx+nv+nz+2]	// scalar
	rhop = exp(rhop)/(1+exp(rhop))		
	rho  = rmin*(1-rhop)+rhop*rmax	
	
	es=_cost*(yy-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
    lndetPi=.
	invPi=.
    ukeys = sort(asarray_keys(w_ina),1)
    Eie = J(0,1,.)

	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)

	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)		
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
	}
	return(exp(-Eie))
	
}





////////////////////////////xtsfspyu////////////////////////////



void function xtsfspyu(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wu_ina
    //external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
	external real scalar _cost
    external real scalar wyeq
	external real colvector _pan_tvar
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
    hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	s2_u=1
	
// wy:->gamma    wv:->rho  wu:->tau	
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= moptimize_util_xb(M,b,5)	// scalar
	taup = exp(taup)/(1+exp(taup))
	rmin = st_numscalar("rumin")
	rmax = st_numscalar("rumax")		
	tau  = rmin*(1-taup)+taup*rmax	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort( asarray_keys(wu_ina),1)
    lndetPi=.
	invPi=.
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
        if(wyeq==1){
           lndetIrhoW = ln(det(I(N)-gamma*spwi))	
        }
	}

    if(wyeq==0){
        ykeys =sort( asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,ykeys[1]))
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	  
        }  
    }

      

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
            if(wyeq==1){
              lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }            
		}
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))
                lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }
            else{
                spwi = extrpoi(asarray(wy_ina,ykeys[1]))
            }            
        } 

		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}


		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
        hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}

//////////////////////////////////////////////////////////////////

void function xtsfspyubp(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wu_ina
    //external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
	external real scalar _cost
    external real scalar wyeq
	external real colvector _pan_tvar
	external real matrix _rho_lndet_
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
    hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	s2_u =1
	// wy:->gamma    wv:->rho  wu:->tau
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax


	taup	= moptimize_util_xb(M,b,5)	// scalar
	taup = exp(taup)/(1+exp(taup))
	rmin = st_numscalar("rumin")
	rmax = st_numscalar("rumax")		
	tau  = rmin*(1-taup)+taup*rmax	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)
	lmax = max(mm_which(_rho_lndet_[.,1]:<=gamma))
	if(length(lmax)==0) lmax=1
	gmin = min(mm_which(_rho_lndet_[.,1]:>=gamma))
	if(length(gmin)==0) gmin=length(_rho_lndet_[.,1])
	lndetIrhoW= lndetIrhoW=_rho_lndet_[round(0.5*lmax+0.5*gmin),2]/nt

    ukeys = sort(asarray_keys(wu_ina),1)
    lndetPi=.
	invPi=.
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
        //if(wyeq==1){
        //   lndetIrhoW = ln(det(I(N)-gamma*spwi))	
        //}
	}

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,ykeys[1]))
            //lndetIrhoW = ln(det(I(N)-gamma*spwi))	  
        }  
    }

      

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,i))
            Mtau = matinv(I(N)-tau*spwi)
            //if(wyeq==1){
            //lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            //}            
		}
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))
                //lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }
            else{
                spwi = extrpoi(asarray(wy_ina,ykeys[1]))
            }            
        } 

		
		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}


		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
        hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}


////////////////////////////xtsfspyu_te////////////////////////////


 function xtsfspyu_te(       real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix wu_ina
    //external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
	external real scalar _cost
    external real scalar wyeq
	external real colvector _pan_tvar
 	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1
	// wy:->gamma    wv:->rho  wu:->tau
	gammap	= b[nx+nv+nz+1]	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax

	taup	= b[nx+nv+nz+2]	// scalar
	taup = exp(taup)/(1+exp(taup))
	rmin = st_numscalar("rumin")
	rmax = st_numscalar("rumax")		
	tau  = rmin*(1-taup)+taup*rmax	
	
	es=_cost*(yy-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	
    lndetPi=.
	invPi=.
    ukeys = sort(asarray_keys(wu_ina),1)
 
	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
	}

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,ykeys[1]))	  
        }  
    }

      
	Eie = J(0,1,.)
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,i))
            Mtau = matinv(I(N)-tau*spwi)           
		}
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))	
            }
            else{
                spwi = extrpoi(asarray(wy_ina,ykeys[1]))
            }            
        } 

		
		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    //lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			//lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}


		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
        hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))

	
}


////////////////////////////xtsfspyv////////////////////////////



void function xtsfspyv(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
	external real scalar _cost
    external real scalar wyeq
	external real colvector _pan_tvar
	//external string scalar wxvars
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	s2_u=1
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
// wy:->gamma    wv:->rho  wu:->tau
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax

	rhop	= moptimize_util_xb(M,b,5)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rvmin")
	rmax = st_numscalar("rvmax")	
	rho  = rmin*(1-rhop)+rhop*rmax
	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)
    lndetPi=.
	invPi=.
    vkeys = sort(asarray_keys(wv_ina),1)
    //vkeys = asarray_keys(wv_ina)
	
	if(length(vkeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wv_ina,vkeys[1]))
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
        if(wyeq!=0){
           lndetIrhoW = ln(det(I(N)-gamma*spwi))	
        }
	}

    if(wyeq==0){
        ykeys =sort( asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,vkeys[1]))
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	  
        }  
    }

    Mtau =1

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(vkeys)>1){
            spwi = extrpoi(asarray(wv_ina,i))
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
            if(wyeq==2){
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }            
		}
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))
                lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }
			else{
				spwi = extrpoi(asarray(wy_ina,1))
			}            
        } 

		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)
						
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}


/////////////////////////////////////////////////////////////


void function xtsfspyvbp(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
	external real scalar _cost
    external real scalar wyeq
	external real colvector _pan_tvar
	external real matrix  _rho_lndet_
	//external string scalar wxvars
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	s2_u =1
	// wy:->gamma    wv:->rho  wu:->tau
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax

	rhop	= moptimize_util_xb(M,b,5)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rvmin")
	rmax = st_numscalar("rvmax")	
	rho  = rmin*(1-rhop)+rhop*rmax
	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
    lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)
	lmax = max(mm_which(_rho_lndet_[.,1]:<=gamma))
	if(length(lmax)==0) lmax=1
	gmin = min(mm_which(_rho_lndet_[.,1]:>=gamma))
	if(length(gmin)==0) gmin=length(_rho_lndet_[.,1])
	lndetIrhoW=_rho_lndet_[round(0.5*lmax+0.5*gmin),2]/nt

    vkeys =sort( asarray_keys(wv_ina),1)
    //vkeys = asarray_keys(wv_ina)
	
	if(length(vkeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wv_ina,vkeys[1]))
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
		/*	
        if(wyeq!=0){
           lndetIrhoW = ln(det(I(N)-gamma*spwi))	
        }
        */
	}

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,vkeys[1]))
            //lndetIrhoW = ln(det(I(N)-gamma*spwi))	  
        }  
    }

    Mtau =1

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(vkeys)>1){
            spwi = extrpoi(asarray(wv_ina,i))
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
			/*			
            if(wyeq==2){
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            } 
            */           
		}
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))
                //lndetIrhoW = ln(det(I(N)-gamma*spwi))	
            }
			else{
				spwi = extrpoi(asarray(wy_ina,1))
			}            
        } 

		
		if (length(s2_v)>1){
			Pi = quadcross(Mrho',panelsubmatrix(s2_v, i, info),Mrho')
			invPi = invsym(Pi)
			lndetPi = ln(det(Pi))    					
		}
						
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}


////////////////////////////xtsfspyv_te////////////////////////////




real colvector function xtsfspyv_te(   real rowvector b,
										real scalar nx,
										string scalar yvar,
										string scalar xvars,
										string scalar zvars,
										string scalar vhvars,
										string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix wv_ina
    external transmorphic  matrix wy_ina
	external real scalar _cost
    external real scalar wyeq
	external real colvector _pan_tvar
	//external string scalar wxvars
 	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	//nx = cols(xx)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1

// wy:->gamma    wv:->rho  wu:->tau
	gammap	= b[nx+nv+nz+1]	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rymin")
	rmax = st_numscalar("rymax")
	gamma  = rmin*(1-gammap)+gammap*rmax

	rhop	= b[nx+nv+nz+2]	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rvmin")
	rmax = st_numscalar("rvmax")		
	rho  = rmin*(1-rhop)+rhop*rmax
	
	es=_cost*(yy-xb)
   lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	Eie = J(0,1,.)
    vkeys = sort(asarray_keys(wv_ina),1)
    //vkeys = asarray_keys(wv_ina)
	
	if(length(vkeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wv_ina,vkeys[1]))
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
	}

    if(wyeq==0){
        ykeys = sort(asarray_keys(wy_ina),1)
        if(length(ykeys)==1){
            N = info[1,2]-info[1,1]+1
            spwi = extrpoi(asarray(wy_ina,vkeys[1]))	  
        }  
    }

    Mtau =1

	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(vkeys)>1){
            spwi = extrpoi(asarray(wv_ina,i))
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}
        if(wyeq==0){
            if(length(ykeys)>1){
                spwi = extrpoi(asarray(wy_ina,i))	
            }
			else{
				spwi = extrpoi(asarray(wy_ina,1))
			}            
        } 

        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)		
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))
	
}


////////////////////////////xtsfspy////////////////////////////



void function xtsfspy(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    external transmorphic  matrix w_ina
    external real colvector _pan_tvar
	external real scalar _cost
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
    hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	s2_u =1
	
	// wy:->gamma    wv:->rho  wu:->tau
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	gamma  = rmin*(1-gammap)+gammap*rmax
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)
    lndetPi=.
	invPi=.
    ykeys = sort(asarray_keys(w_ina),1)

	if(length(ykeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ykeys[1]))
        lndetIrhoW = ln(det(I(N)-gamma*spwi))	
   
	}
    
    Mtau=1
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ykeys)>1){
            spwi = extrpoi(asarray(w_ina,i))
            lndetIrhoW = ln(det(I(N)-gamma*spwi))	
          
		}
		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}          
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
        hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}

/////////////////////////////////////////////////////////////////


void function xtsfspybp(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    external transmorphic  matrix w_ina
    external real colvector _pan_tvar
	external real scalar _cost
	external real matrix _rho_lndet_
	
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
    hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	s2_u = 1
	gammap	= moptimize_util_xb(M,b,4)	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	gamma  = rmin*(1-gammap)+gammap*rmax
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
    lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)
	lmax = max(mm_which(_rho_lndet_[.,1]:<=gamma))
	if(length(lmax)==0) lmax=1
	gmin = min(mm_which(_rho_lndet_[.,1]:>=gamma))
	if(length(gmin)==0) gmin=length(_rho_lndet_[.,1])
	lndetIrhoW=_rho_lndet_[round(0.5*lmax+0.5*gmin),2]/nt
	

    ykeys = sort(asarray_keys(w_ina),1)

	if(length(ykeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ykeys[1]))
        //lndetIrhoW = ln(det(I(N)-gamma*spwi))

	}
    
    Mtau=1
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ykeys)>1){
            spwi = extrpoi(asarray(w_ina,i))
            //lndetIrhoW = ln(det(I(N)-gamma*spwi))	
       
		}
 		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}         
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(y, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
        hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = lndetIrhoW-0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) ///
		        + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}



////////////////////////////xtsfspy_te////////////////////////////



real function xtsfspy_te(   real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    external transmorphic  matrix w_ina
    external real colvector _pan_tvar
	external real scalar _cost
	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1	

	gammap	= b[nx+nv+nz+1]	// scalar
	gammap = exp(gammap)/(1+exp(gammap))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	gamma  = rmin*(1-gammap)+gammap*rmax
	
	es=_cost*(yy-xb)
   lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	Eie = J(0,1,.)

    ykeys = sort(asarray_keys(w_ina),1)

	if(length(ykeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ykeys[1]))	
	}
    
    Mtau=1
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ykeys)>1){
            spwi = extrpoi(asarray(w_ina,i))     
		}
		
		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    //lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			//lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}		
          
		esi = panelsubmatrix(es, i, info)-_cost*gamma*spwi*panelsubmatrix(yy, i, info) //
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
        hbi	= Mtau*panelsubmatrix(hb, i, info)
		//rows(invPi),cols(invPi)
		//rows(hbi),cols(hbi)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))
	
}


////////////////////////////xtsfspu////////////////////////////



void function xtsfspu(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
                        //real colvector lnfj,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar
    //external transmorphic  matrix wu_ina
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar
	s2_u = 1
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	//rho  = rmin*(1-rhop)+rhop*rmax
	taup	= moptimize_util_xb(M,b,4)	// scalar
	taup = exp(taup)/(1+exp(taup))
	tau  = rmin*(1-taup)+taup*rmax	
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
    lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(w_ina),1)
	if(length(ukeys)==1 ){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)			
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,ukeys[i]))
            Mtau = matinv(I(N)-tau*spwi)
		}     
		
		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}		

		
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = -0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}



////////////////////////////xtsfspu_te////////////////////////////



  function xtsfspu_te(   real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar
    //external transmorphic  matrix wu_ina
	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	//rho  = rmin*(1-rhop)+rhop*rmax
	taup	= b[nx+nv+nz+1]	// scalar
	taup = exp(taup)/(1+exp(taup))
	tau  = rmin*(1-taup)+taup*rmax	

	es=_cost*(yy-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
   lndetPi=.
	invPi=.	
	Eie = J(0,1,.)
    ukeys = sort(asarray_keys(w_ina),1)
	if(length(ukeys)==1 ){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)			

	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,ukeys[i]))
            Mtau = matinv(I(N)-tau*spwi)
		}
		
		if (length(s2_v)==1){
			invPi = 1/s2_v*I(N)
		    //lndetPi =N*ln(s2_v) //
		}
		else{
			invPi = diag(1:/panelsubmatrix(s2_v, i, info))
			//lndetPi = sum(ln(panelsubmatrix(s2_v, i, info)))
		}		

		
		esi = panelsubmatrix(es, i, info)
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))
	
}


////////////////////////////xtsfspv////////////////////////////



void function xtsfspv(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
                        //real colvector lnfj,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	s2_u = 1
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	rhop	= moptimize_util_xb(M,b,4)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	rho  = rmin*(1-rhop)+rhop*rmax
	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)
    lndetPi=.
	invPi=.
    ukeys =sort(asarray_keys(w_ina),1)
    Mtau=1
	if(length(ukeys)==1 ){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))			
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,ukeys[i]))			
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}
		
        tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)	
	
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = -0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}


////////////////////////////xtsfspv_te////////////////////////////



real function xtsfspv_te(   real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
	external transmorphic  matrix w_ina
	external real scalar _cost
	external real colvector _pan_tvar
	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1
	rhop	= b[nx+nv+nz+1]	// scalar
	rhop = exp(rhop)/(1+exp(rhop))
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")
	rho  = rmin*(1-rhop)+rhop*rmax

	es=_cost*(yy-xb)
   lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	Eie = J(0,1,.)

    ukeys = sort(asarray_keys(w_ina),1)
    Mtau=1
	if(length(ukeys)==1 ){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))			
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)
        
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,ukeys[i]))			
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}
		
		tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)	
	
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))
	
}


////////////////////////////xtsfspuv////////////////////////////



void function xtsfspuv(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
                        //real colvector lnfj,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix wv_ina
    external transmorphic  matrix wu_ina
	external real colvector _pan_tvar
	external real scalar _cost
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	s2_u =1 
	taup	= moptimize_util_xb(M,b,5)	// scalar
	taup = exp(taup)/(1+exp(taup))
    rumin = st_numscalar("rumin")
	rumax = st_numscalar("rumax")
	tau  = rumin*(1-taup)+taup*rumax

	rvmin = st_numscalar("rvmin")
	rvmax = st_numscalar("rvmax")
	rhop	= moptimize_util_xb(M,b,4)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))    
	rho  = rvmin*(1-rhop)+rhop*rvmax
	
	y =moptimize_util_depvar(M,1)
	//if (st_global("cost")=="cost") es=-(y-xb)
	//else es=y-xb
	es = _cost*(y-xb)
    lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(wu_ina),1)
    vkeys = sort(asarray_keys(wv_ina),1)

	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
	}
    if(length(vkeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wv_ina,vkeys[1]))
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)	
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,ukeys[i]))
            Mtau = matinv(I(N)-tau*spwi)
		}
 		if(length(vkeys)>1){
            spwi = extrpoi(asarray(wv_ina,vkeys[i]))
            Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}       

		tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = -0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}


////////////////////////////xtsfspuv_te////////////////////////////


real function xtsfspuv_te(   real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix wv_ina
    external transmorphic  matrix wu_ina
	external real colvector _pan_tvar
	external real scalar _cost
	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	//nv
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1

	taup	= b[nx+nv+nz+2]	// scalar
	taup = exp(taup)/(1+exp(taup))
    rumin = st_numscalar("rumin")
	rumax = st_numscalar("rumax")
	tau  = rumin*(1-taup)+taup*rumax

	rvmin = st_numscalar("rvmin")
	rvmax = st_numscalar("rvmax")
	rhop	=  b[nx+nv+nz+1]	// scalar
	rhop = exp(rhop)/(1+exp(rhop))    
	rho  = rvmin*(1-rhop)+rhop*rvmax
//	printf("test1..")
	es = _cost*(yy-xb)
    lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
 	Eie = J(0,1,.)   

    ukeys = sort(asarray_keys(wu_ina),1)
    vkeys = sort(asarray_keys(wv_ina),1)

	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wu_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
	}
    if(length(vkeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(wv_ina,vkeys[1]))
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)		
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(wu_ina,ukeys[i]))
            Mtau = matinv(I(N)-tau*spwi)
		}
 		if(length(vkeys)>1){
            spwi = extrpoi(asarray(wv_ina,vkeys[i]))
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}       

		tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
//		printf("test3..")
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))
	
}


////////////////////////////xtsfspuv0////////////////////////////



void function xtsfspuv0(transmorphic scalar M,
					    real scalar todo,
                        real rowvector b,
                        //real colvector lnfj,
					    real scalar lnf,
					    real rowvector g, 
					    real matrix H)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real colvector _pan_tvar	
    //external transmorphic  matrix wu_ina
	external real scalar _cost
	xb = moptimize_util_xb(M,b,1)
	hb = moptimize_util_xb(M,b,3)
	hb = exp(0.5*hb)
	s2_u =1
	//s2_u	= exp(moptimize_util_xb(M,b,4))	// scalar
	s2_v	= exp(moptimize_util_xb(M,b,2))	// scalar	
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")

	taup	= moptimize_util_xb(M,b,5)	// scalar
	taup = exp(taup)/(1+exp(taup))
	tau  = rmin*(1-taup)+taup*rmax	
		
	rhop	= moptimize_util_xb(M,b,4)	// scalar
	rhop = exp(rhop)/(1+exp(rhop))	
	rho  = rmin*(1-rhop)+rhop*rmax

	y =moptimize_util_depvar(M,1)
	es=_cost*(y-xb)
    lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	lnfj = J(nt,1,.)

    ukeys = sort(asarray_keys(w_ina),1)
    //vkeys = asarray_keys(wv_ina)

	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)		
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,ukeys[i]))
            Mtau = matinv(I(N)-tau*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}       

		tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		lnfj[i] = -0.5*N*ln(2*pi())-0.5*lndetPi-0.5*quadcross(esinvPi',esi) + 0.5*(us^2/s2_s-mu^2/s2_u)+ln(sqrt(s2_s)*normal(us/sqrt(s2_s)))-ln(sqrt(s2_u)*normal(mu/sqrt(s2_u)))
		
	}

	lnf = moptimize_util_sum(M, lnfj)
	
}



////////////////////////////xtsfspuv0_te////////////////////////////



real function xtsfspuv0_te(   real rowvector b,
							 real scalar nx,
							 string scalar yvar,
							 string scalar xvars,
							 string scalar zvars,
							 string scalar vhvars,
							 string scalar noconstant)
{
    //external real matrix spw
	external transmorphic  matrix w_ina
	external real colvector _pan_tvar	
    //external transmorphic  matrix wu_ina
	external real scalar _cost
	xx = st_data(.,xvars)
	zz = loaddata(.,zvars,rows(xx))
	vh = loaddata(.,vhvars,rows(xx))
	yy = st_data(.,yvar)
	if(noconstant=="constant"){
		nx= nx+1
		xx = xx, J(rows(xx),1,1)
	}
	nz = cols(zz)
	nv = cols(vh)
	
	xb = xx*(b[1..nx])'	
	s2_v	= exp(vh*(b[(nx+1)..(nx+nv)])')	// scalar	
	hb = exp(0.5*zz*(b[(nx+nv+1)..(nx+nv+nz)])')
	//s2_u	= exp(b[nx+nz+2])	// scalar
	s2_u=1	
	rmin = st_numscalar("rmin")
	rmax = st_numscalar("rmax")

	taup	= b[nx+nv+nz+2]	// scalar
	taup = exp(taup)/(1+exp(taup))
	tau  = rmin*(1-taup)+taup*rmax	
		
	rhop	= b[nx+nv+nz+1]	// scalar
	rhop = exp(rhop)/(1+exp(rhop))	
	rho  = rmin*(1-rhop)+rhop*rmax

	es=_cost*(yy-xb)
   lndetPi=.
	invPi=.	
	//create view of panels
	//st_view(_pan_tvar = .,.,st_global("paneltvar"))
	info = panelsetup(_pan_tvar,1)
	nt = panelstats(info)[1]
	mu = 0
	Eie = J(0,1,.)   

    ukeys = sort(asarray_keys(w_ina),1)
    //vkeys = asarray_keys(wv_ina)

	if(length(ukeys)==1){
		N = info[1,2]-info[1,1]+1
		spwi = extrpoi(asarray(w_ina,ukeys[1]))
		Mtau = matinv(I(N)-tau*spwi)		
		Mr = I(N)-rho*spwi
		iMr = matinv( Mr)
		pifun(s2_v,Mr,iMr,lndetPi,invPi)		
	}
	for(i=1;i<=nt;i++){
		N = info[i,2]-info[i,1]+1
		if(length(ukeys)>1){
            spwi = extrpoi(asarray(w_ina,ukeys[i]))
            Mtau = matinv(I(N)-tau*spwi)
			Mr = I(N)-rho*spwi
			iMr = matinv( Mr)
			pifun(s2_v,Mr,iMr,lndetPi,invPi)
		}       

		tvpifun(s2_v,Mr,iMr,lndetPi,invPi,i,info)
		//lndetPi = ln(det(Pi))
		esi = panelsubmatrix(es, i, info)
		//hbi	= Mtau*exp(panelsubmatrix(hb, i, info))
		hbi	= Mtau*panelsubmatrix(hb, i, info)
		hinvPi = quadcross(hbi,invPi)
		hinvPihi = quadcross(hinvPi',hbi)
		esinvPi = quadcross(esi,invPi)
		//us = (mu/s2_u - quadcross(esinvPi',hbi))/(hinvPihi+1/s2_u)
		//s2_s=1/(hbi'*invPi*hbi+1/s2_u)
		s2_s = 1/(hinvPihi+1/s2_u)
		us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
		Eie = Eie \ hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
		
	}
	return(exp(-Eie))
	
}



function myrowsum(x) return(rowsum(x))
function mycolsum(x) return(colsum(x))
function myquadsum(x) return(quadsum(x))


 mata mlib create lxtsfsp, replace 
 mata mlib add lxtsfsp *()
 mata mlib index

end


