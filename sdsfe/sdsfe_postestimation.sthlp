{smcl}
{* *! version 1.0.0  10mar2015}{...}
{* revised: }{...}
{cmd:help sdsfe postestimation}{right:also see: {help sdsfe}}
{hline}

{title:Title}

{p2colset 5 32 38 2}{...}
{p2col :{hi:sdsfe postestimation} {hline 2}}Postestimation tools for sdsfe{p_end}
{p2colreset}{...}


{title:Description}

{pstd}
The following postestimation commands are available for {cmd:sdsfe}.

{synoptset 13}{...}
{p2coldent :command}description{p_end}
{synoptline}
{synopt:{bf:{help estat}}}AIC, BIC, VCE, and estimation sample summary{p_end}
INCLUDE help post_estimates
INCLUDE help post_lincom
INCLUDE help post_lrtest
INCLUDE help post_nlcom
{synopt :{helpb sdsfe postestimation##predict:predict}}predictions, residuals, influence statistics, and other diagnostic measures{p_end}
INCLUDE help post_predictnl
INCLUDE help post_test
INCLUDE help post_testnl
{synoptline}
{p2colreset}{...}


{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}{cmd:predict} {dtype} {newvar} {ifin} [{cmd:,} {it:statistic}]

{synoptset 15 tabbed}{...}
{synopthdr :statistic}
{synoptline}
{syntab:Main}
{synopt :{opt xb}}linear prediction; the default{p_end}
{synopt :{opt residuals}} estimates of residuals (depvar - xb) {p_end}
{synopt :{opt u}}estimates of (technical or cost) inefficiency via {it:E}(u|e) (Jondrow et al., 1982){p_end}
{synopt :{opt jlms}}estimates of (technical or cost) efficiency via exp[-E(u|e)]{p_end}
{synopt :{opt bc}}estimates of (technical or cost) efficiency via E[exp(-u)|e] (Battese and Coelli 1988){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
These statistics are only available for the estimation sample.


INCLUDE help menu_predict


{title:Options for predict}

{dlgtab:Main}

{phang}
{opt xb}, the default, calculates the linear prediction.

{phang}
{opt residuals} calculates the composite residuals (depvar - xb).

 
{phang}
{opt u} produces estimates of (technical or cost) inefficiency via E(u|e) using the Jondrow et al. (1982) estimator.

{phang}
{opt jlms} produces estimates of (technical or cost) efficiency via exp[-E(u|e)]. 


{phang}
{opt bc} produces estimates of (technical or cost) efficiency via E[exp(-u)|e]. 


{title:Remarks}

{pstd}When the {cmd:sdsfe} command is used to estimate production frontiers, {cmd:predict} will provide the post-estimation of technical (in)efficiency, 
while when the {cmd:sdsfe} command is used to estimate cost frontiers, {cmd:predict} will provide the post-estimation of cost (in)efficiency. 
It is worth noting that {cmd:sdsfe} and the related {cmd:predict} command follow the definitions of technical and cost (in)efficiency given in Kumbhakar and Lovell (2000).{p_end}



    {marker examples}{...}
    {title:Examples}
    
        {title:SD-SF model with time-constant spatial weight matrix}
    
    {pstd}
    Setup{p_end}
    {phang2}{bf:. {stata "mata mata matuse w_ex1,replace"}}{p_end}
    {phang2}{bf:. {stata "use sdsfbc_ex1.dta"}}{p_end}
    
    {pstd}
    Stochastic Durbin production model {p_end}
    {phang2}{bf:. {stata "sdsfe y x, id(id) time(t) noconstant wmat(wm,mata) mu(z) wxvars(x) wmuvars(z) genwvars "}}{p_end}
    {phang2}{bf:. {stata "predict uhat,u"}}{p_end}
    {phang2}{bf:. {stata "predict te_bc,bc"}}{p_end}
    {phang2}{bf:. {stata "predict te_jlms,jlms"}}{p_end}
    

{marker author}{...}
{title:Author}

{pstd}
Kerui Du{break}
Xiamen University{break}
School of Management{break}
China{break}
{browse "kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}{break}


{pstd}
Federica Galli{break}
University of Bologna{break}
Department of Statistical Sciences “Paolo Fortunati”{break}
Italy{break}
{browse "federica.galli14@unibo.it":federica.galli14@unibo.it}{break}


{pstd}
Luojia Wang{break}
Xiamen University{break}
School of Management{break}
China{break}



{title:Also see}

{psee}
{space 2}Help:  {help sdsfe}.
{p_end}
