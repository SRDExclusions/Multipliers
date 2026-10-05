
***** Multiplier SVAR analysis *****
// By Siyanda Baduza - 2026 IEJ

**** 3 key variables - GDP, Social Protection spending, Government Revenue
// Do file assumes that the variables are already in real terms and already seasonally adjusted. Data is already set to time series. 

*** exogenous (dummies) to be used - user determined, example:

gen covid_dip = (t == tq(2020q2))
gen covid_peak = (t == tq(2020q3) | t == tq(2020q4))
gen dum06change = (t == tq(2006q2))
gen dum0809 = (t >= tq(2008q3) & t < tq(2009q2) )


**** Setting up the key variables *****

* Take logs of the variables
gen ln_g = ln()
gen ln_t = ln()
gen ln_y = ln()

* Difference the logged terms
gen dlng = d.ln_g
gen dlnt = d.ln_t
gen dlny = d.ln_y


/// Check for stationarity - ADF, KPSS, Phillips-Perron (analyis is already done for logged-differences though)

* on levels

dfuller ln_g, lags(4) trend
dfuller ln_t, lags(4) trend
dfuller ln_y, lags(4) trend

pperron ln_g, lags(4) trend 
pperron ln_t, lags(4) trend
pperron ln_y, lags(4) trend

kpss ln_g
kpss ln_t
kpss ln_y

* on differences

dfuller dlng, lags(4) trend
dfuller dlnt, lags(4) trend
dfuller dlny, lags(4) trend

pperron dlng, lags(4) trend 
pperron dlnt, lags(4) trend
pperron dlny, lags(4) trend
 
kpss dlng
kpss dlnt
kpss dlny


* Lag selection criteria


varsoc dlng dlnt dlny, maxlag(8) exog(covid_dip covid_peak dum0809 dum06change) 
local optimal_lag =  // selected lag, add based on above results

** Reduced VAR estimation

// Country specific dummies

var dlng dlnt dlny, lags(1/`optimal_lag') exog(covid_dip covid_peak dum06change dum0809)
varlmar, mlag(8) /// check autocorrelation
varstable // check stability
varstable, graph

*** Extract VAR residuals (reduced form residuals u)

predict u_g, equation(dlng) residuals
predict u_t, equation(dlnt) residuals
predict u_y, equation(dlny) residuals

corr u_g u_t u_y // should be correlated

*** Structural shock decomposition *****

// IMF method - alpha_ty

gen trend = _n


ardl ln_t ln_y, ec exog(trend)
predict prelim_resid, residuals
sum prelim_resid, d
gen outlier = abs(prelim_resid) >2.5*r(sd)


*** ardl with a trend and outlier variables 

ardl ln_t ln_y if t >= tq(2000q1), ec exog(trend outlier)
local alpha_ty = _b[SR: D1.ln_y] // estimated alpha_ty

// alpha_gy and e_g

local alpha_gy = 0 // restriction - set it to 0.

gen u_g_ca = u_g - `alpha_gy' * u_y
gen u_t_ca = u_t - `alpha_ty' * u_y

gen e_g = u_g_ca // restriction leads to u_g_ca = u_g = e_g

// beta_tg and e_t

reg u_t_ca e_g, noconstant
local beta_tg = _b[e_g]
predict e_t, residuals

// gamma yg, yt and e_y

ivregress 2sls u_y (u_g u_t = e_g e_t), noconstant perfect
predict e_y, residuals
local gamma_yg = _b[u_g]
local gamma_yt = _b[u_t]

/// check correlation, should not be correlated
corr e_g e_t e_y

*** SVAR estimation (can also let the diagonal elements be freely estimated, leads to the same results) ***

matrix A_restrict = (1, 0, 0 \ 0, 1, -`alpha_ty' \ -`gamma_yg', -`gamma_yt', 1)
matrix B_restrict = (1, 0, 0\ `beta_tg', 1, 0 \ 0, 0, 1)



*** IRF and multiplier calculations***

svar dlng dlnt dlny if t >= tq(2000q1),aeq(A_restrict) beq(B_restrict) lags(1/`optimal_lag') exog(covid_dip covid_peak dum06change dum0809)
irf create svar_1, set(irf_general_output) step(10) bsp reps(1000) replace

irf graph cirf, impulse(dlng) response(dlny)

irf table cirf, impulse(dlng) response(dlny)
irf table cirf, impulse(dlng) response(dlng)


// Impact multiplier = cirf(at 1) for g on y / cirf(at 1) for g on g all divided by the ratio of social spending on gdp
// Cumulative multiplier = cirf(at 10) for g on y/ cirf(at 10) for g on g all divided by the ratio of social spending on gdp

// GDP components, switch y for either consumption or investment