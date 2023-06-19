{smcl}
{cmd:help spxtsfa}
{hline}

{title:Title}

{p2colset 5 13 15 2}{...}
{p2col :{hi:spxtsfa} {hline 2}}Spatial panel stochastic frontier models in the style of {help spxtsfa##OA2019:{bind:Orea and Álvarez (2019)}} and {help spxtsfa##Galli2022:{bind:Galli (2022)}} {p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Estimation syntax

{p 8 17 2}
{cmd:spxtsfa} {depvar} [{indepvars}] [{cmd:,} {it:options}]

{pstd}
Version syntax

{p 8 17 2}
{cmd:spxtsfa}{cmd:,} {opt ver:sion}

{pstd}
Replay syntax

{p 8 17 2}
{cmd:spxtsfa} [{cmd:,} {cmdab:l:evel(}{help level##remarks:{it:#}}{cmd:)}]


{synoptset 31 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Frontier}
{synopt :{opt nocons:tant}}suppress constant term{p_end}
{synopt :{opt cost}}fit cost frontier model; default is {cmd:production}{p_end}
{synopt :{cmd:wxvars({it:varlist})}}spatially lagged independent variables{p_end}

{syntab :Inefficiency variance function}
{synopt :{cmd:uhet({it:varlist})}}explanatory variables for technical inefficiency variance function{p_end}

{syntab :Spatial weight matrix}
{synopt :{cmd:wy(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix for lagged dependent variable{p_end}
{synopt :{cmd:wx(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix for lagged independent variables{p_end}
{synopt :{cmd:wu(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix for technical inefficiency{p_end}
{synopt :{cmd:wv(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix for the error term{p_end}
{synopt :{cmd:normalize(}{it:string}{cmd:)}}specify the normalized method of spatial weight matrixs{p_end}

{syntab :Regression}
{synopt :{cmdab:init:ial(}{it:{help matrix:matname}}{cmd:)}}specify initial values matrix{p_end}
{synopt :{cmd:mlsearch(}{it:{help ml##model_options:search_options}}{cmd:)}}specify options for searching initial values{p_end}
{synopt :{opt delve}}delve into maximization problem to find initial values{p_end}
{synopt :{opt mlplot}}use ml plot to find better initial values{p_end}
{synopt :{cmd:mlmodel(}{it:{help ml##model_options:model_options}}{cmd:)}}control {cmd:ml model} options{p_end}
{synopt :{cmd:mlmax({it:{help ml##ml_maximize_options:maximize_options}})}}control {cmd:ml maximize} options{p_end}
{synopt :{cmd:lndetmc({it:numlist})}}settings for BarryPace Trick. Order is iterations, maxorder.{p_end}

{syntab :Reporting}
{synopt :{cmd:nolog}}omit the display of the criterion function iteration log{p_end}
{synopt :{cmdab:mldis:play(}{it:{help ml##display_options:display_options}}{cmd:)}}control {cmd:ml display} options; seldom used{p_end}
{synopt :{cmd:te(}{it:{help newvar:effvar}}{cmd:)}}create efficiency variables{p_end}
{synopt :{opt genwxvars}}generate the spatial Durbin terms{p_end}

{syntab :Other}
{synopt :{cmdab:constraints(}{it:{help estimation options##constraints():constraints}}{cmd:)}}apply specified linear constraints{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
uhet() is required. At least one of wy(),wu() and wv() should be specified. 
A panel variable and a time variable must be specified; use xtset. {p_end}


{marker description}{...}
{title:Description}

{pstd}
{opt spxtsfa} fits spatial panel stochastic production or cost frontier models
following the methodology of 
{help spxtsfa##OA2019:{bind:Orea and Álvarez (2019)}} and {help spxtsfa##Galli2022:{bind:Galli (2022)}}. See 
{help spxtsfa##OA2019:{bind:Orea and Álvarez (2019)}} and {help spxtsfa##Galli2022:{bind:Galli (2022)}} for a detailed
explanation of their methodology and empirical analyses.


{marker options}{...}
{title:Options for the estimation syntax}

{dlgtab:Frontier}

{phang}
{cmd:noconstant} suppresses the constant term (intercept) in the frontier.

{phang}
{cmd:cost} specifies that the model to be fit is a cost frontier model.  The
default is {cmd:production}.

{phang}
{cmd: wxvars({it:varlist})} specifies spatially lagged independent variables.

{dlgtab:Ineffciency}

{phang}
{cmd:uhet({it:varlist})} specifies that the variables
in the technical inefficiency variance function.

{dlgtab:Spatial weught matrix}

{phang}
{opt wy(W1 [W2 ... WT][,mata array])} specifies that the spatial weight matrix for lagged dependent variable. 

{phang}
{opt wx(W1 [W2 ... WT][,mata array])} specifies that the spatial weight matrix for lagged independent variable. 

{phang}
{opt wu(W1 [W2 ... WT][,mata array])} specifies spatial weight matrix for technical inefficiency{p_end}

{phang}
{opt wv(W1 [W2 ... WT][,mata array])} specifies spatial weight matrix for the error term{p_end}

{phang}
By default, the weight matrices are Sp object. mata declares weight matrices are mata matrix. 
If one weight matrix is specified, it assumes time-constant weight matrix. 
For time-varying cases, T weight matrices should be specified in time order. 
Alternatively, using array to declares weight matrices are store in a array.  
If only one matrix is stored in the specified array, the time-constant weight matrix is assumed.  
Otherwise, the keys of the array specifies time information and 
the values store time-specific weight matrices.  

{phang}
{opt normalize(row|col|spectral|minmax)} specifies the normalized method of spatial weight matrixs. 
By default, the command would not normalization the spatial weight matrixs. normalize(row) is row normalisation;
normalize(col) is collumn normalisation; normalize(spectral) is spectral normalisation;
normalize(minmax) is minmax normalisation.


{dlgtab:Regression}

{phang}
{cmd:initial(}{it:{help matrix:matname}}{cmd:)} specifies that {it:matname} is
the initial value matrix.

{phang}
{cmd:mlsearch(}{it:{help ml##model_options:search_options}}{cmd:)} specifies ml search options for searching initial values.

{phang}
{cmd:delve} provides a regression-based methodology to search for 
initial values. The default is to use {helpb ml search:ml search} with default options.

{phang}
{opt mlplot} specifes using ml plot to find better initial values.

{phang}
{cmd:mlmodel({it:{help ml##mlmode:model_options}})} controls the {cmd:ml}
{cmd:model} options; it is seldom used.

{phang}
{cmd:mlmax({it:{help ml##ml_max_descript:maximize_options}})} controls the
{cmd:ml max} options; it is seldom used.

{phang}
{cmd:lndetmc({it:numlist})} settings for BarryPace Trick to solve the inverse of (I−ρW). 
Order is iterations, maxorder. lndetmc(50 100) is recommended.

{phang}
{opt delmissing} deletes the units with missing oservations from spmatrix.

{dlgtab:Reporting}

{phang}
{cmd:nolog} suppresses the display of the criterion function iteration log.

{phang}
{cmd:mldisplay({it:{help ml##mldisplay:display_options}})} controls the
{cmd:ml display} options; it is seldom used.

{phang}
{cmd:te(}{it:{help newvar:effvar}}{cmd:)} generates
the production or cost efficiency variable.

{phang}
{cmd:genwxvars} generates the spatial Durbin terms. 
The option automatically extends any specified variable name in wxvars() 
with W_ prefix.


{marker optionsversion}{...}
{title:Options for the version and replay syntax}

{phang}
{cmd:version} displays the version of {cmd:spxtsfa} installed on Stata and the
program author information.  This option can be used only in version syntax.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals.  The default is {cmd:level(95)} or as set by {helpb set level}.
This option can be used in the replay syntax or in 
{cmd:mldisplay(}{it:{help ml##mldisplay:display_options}}{cmd:)}.

{marker optionsversion}{...}
{title:Other}

{phang}
{cmdab:constraints(}{it:{help estimation options##constraints():constraints}}{cmd:)}
specifies linear constraints for the estimated model.

{marker examples}{...}
{title:Examples}

    {title:yuv-SAR model with time-varying spatial weight matrix}

{pstd}
Setup{p_end}
{phang2}{bf:. {stata "mata mata matuse spxtsfa_wmat,replace"}}{p_end}
{phang2}{bf:. {stata "use yuvexample.dta"}}{p_end}
{phang2}{bf:. {stata "xtset id t"}}{p_end}

{pstd}
Stochastic production model with three different sources of spatial cross-sectional dependence {p_end}
{phang2}{bf:. {stata "spxtsfa y x, noconstant uhet(z) wu(w1 w2 w3 w4 w5,mata) wy(w1 w2 w3 w4 w5,mata) wv(w1 w2 w3 w4 w5,mata) te(te)"}}{p_end}


    {title:uv-SAR model with time-invariant spatial weight matrixs}

{pstd}
Setup{p_end}
{phang2}{bf:. {stata "mata mata matuse spxtsfa_wmat,replace"}}{p_end}
{phang2}{bf:. {stata "use uvexample.dta"}}{p_end}
{phang2}{bf:. {stata "xtset id t"}}{p_end}

{pstd}
Stochastic cost model with spatial cross-sectional dependence in inefficiency and random error terms. {p_end}
{phang2}{bf:. {stata "spxtsfa y x, cost noconstant wu(w2,mata) wv(w1,mata) uhet(z) te(te)"}}{p_end}
{phang2}{bf:. {stata "xtset id t"}}{p_end}
  
    {title:yx-SAR model with time-invariant spatial weight matrixs}

{pstd}
Setup{p_end}
{phang2}{bf:. {stata "mata mata matuse spxtsfa_wmat,replace"}}{p_end}
{phang2}{bf:. {stata "use yxexample.dta"}}{p_end}

{pstd}
Stochastic cost model with spatial cross-sectional dependence in the dependent variable. {p_end}
{phang2}{bf:. {stata "spxtsfa y x, cost noconstant wy(w2,mata) wx(w1,mata) wxvars(x) uhet(z) te(te) genwxvars"}}{p_end}  


{marker results}{...} 
{title:Stored results}

{pstd}
{cmd:spxtsfa} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(p)}}significance{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}command used for estimation{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared test{p_end}
{synopt:{cmd:e(vce)}}{it:oim}{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(which)}}{cmd:max} or {cmd:min}; whether optimizer is to perform maximization or minimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(cmdbase)}}base command{p_end}
{synopt:{cmd:e(function)}}production or cost{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}
We are grateful to Federica Galli for his Matlab codes, Belotti, Daidone, Ilardi and Atella for the sfpanel package,  Karakaplan for the sfkk package, 
and Jan Ditzen, William Grieser and Morad Zekhnini for the nwxtregress package. All remaining errors are our own.

{pstd}
Kerui Du acknowledges financial support from the National Natural Science Foundation of China(Grant no. 72074184).



{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
{cmd:spxtsfa} is not an official Stata command.  It is a third-party command
programmed as a free contribution
to the research society.  By choosing to download, install, and use the
{cmd:spxtsfa} package, users assume all the liability for any
{cmd:spxtsfa}-package-related risk.  If you encounter any problems with the
{cmd:spxtsfa} package, or if you have comments, suggestions, or questions, please
send an email to Kerui Du at 
{browse "mailto:kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}.


{marker citation}{...}
{title:References}

{marker OA2019}{...}
{phang}
Orea, Luis, and Inmaculada C. Álvarez. “A New Stochastic Frontier Model with Cross-Sectional Effects in Both Noise and Inefficiency Terms.” 
Journal of Econometrics 213, 
no. 2 (December 1, 2019): 556–77. https://doi.org/10.1016/j.jeconom.2019.07.004.

{marker Galli2022}{...}
{phang}
Galli, Federica. “A Spatial Stochastic Frontier Model Including Both Frontier and Error-Based Spatial 
Cross-Sectional Dependence.” 
Spatial Economic Analysis 0, no. 0 (July 26, 2022): 1–20. 
https://doi.org/10.1080/17421772.2022.2097729.

{phang}
Belotti, F., Daidone, S., Ilardi, G.,  Atella, V. "Stochastic frontier analysis using Stata", Stata Journal, 2013, 13(4):719-758.

{phang}
 Karakaplan, M. U. "Estimating endogenous stochastic frontier models in Stata", Stata Journal, 2017,17: 39-55. 

{phang}
Ditzen, Grieser, Zekhnini. "nwxtregress - network regression in Stata", 2023. https://github.com/JanDitzen/nwxtregress



{marker author}{...}
{title:Author}

{pstd}
Kerui Du{break}
Xiamen University{break}
School of Management{break}
China{break}
{browse "kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}{break}


{pstd}
Luis Orea{break}
University of Oviedo{break}
Department of Economics{break}
Spain{break}
{browse "lorea@uniovi.es":lorea@uniovi.es}{break}


{pstd}
Inmaculada C. Álvarez{break}
Universidad Autónoma de Madrid{break}
Department of Economics{break}
Spain{break}
{browse "inmaculada.alvarez@uam.es":inmaculada.alvarez@uam.es}{break}


{marker see}{...}
{title:Also see}

{p 7 14 2}{manhelp frontier R}, 
{manhelp xtfrontier XT}{p_end} 
{p 7 14 2}{helpb sfpanel}, {helpb sfkk},  {helpb nwxtregress} (if installed)
