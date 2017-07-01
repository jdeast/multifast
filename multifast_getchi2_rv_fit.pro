;+
; NAME:
;   EXOFAST_GETCHI2_RV_FIT
;
; PURPOSE: 
;   Computes the chi^2 of a single planet decribed by PARS, while
;   analytically fitting for K, gamma, slope
;
; CALLING SEQUENCE:
;    chi2 = exofast_getchi2_rv_fit(pars)
;
; INPUTS:
;    The best-fit parameters for the RV fit of a single planet.
;
;     pars[0] = time of transit center
;     pars[1] = period
;     pars[2] = e*cos(omega)
;     pars[3] = e*sin(omega)
;     pars[4] = K            (will be overwritten)
;     pars[5] = gamma        (will be overwritten)
;     pars[6] = slope        (will be overwritten)
;
; RESULT:
;   The chi^2 of the parameters.
;
; COMMON BLOCKS:
;   RV_BLOCK - See exofast.pro for definition
;
; MODIFICATION HISTORY 
;  2012/06 -- Public release -- Jason Eastman (LCOGT)
;-

function multifast_getchi2_rv_fit, pars

COMMON rv_block, data

; data is a structure containing bjd, rv, err
; pars[0] = time of transit center
; pars[1] = period
; pars[2] = e*cos(omega)
; pars[3] = e*sin(omega)
; pars[4] = K
; pars[5] = gamma
; pars[6] = slope

COMMON boundary_block, be, ba, bteff, bfeh, blogg, bdp, bds, bar, bK, bTTV, bF0 
if be    EQ !null then be    = [0.0d0, 1.0d0]  ;eccentrcity bounds  
if bK    EQ !null then bK    = [0.0d0, !values.d_infinity] ;velocity semi-amplitude bounds (m/s)

;; recalculate the relevant parameters
e = sqrt(pars[2]^2 + pars[3]^2)

;; not an allowed planet
if e lt be[0] or e gt be[1] then return, !values.d_infinity 
    
;; compute the model inputs
if e eq 0 then omega = !dpi/2.d0 $
else omega = atan(pars[3]/pars[2])
if pars[2] lt 0 then omega += !dpi

mintime = min(data.bjd,max=maxtime)
t0 = (mintime + maxtime)/2.d0

;; calculate the time of periastron
phase = exofast_getphase(e,omega,/primary)
t_periastron = pars[0] - pars[1]*phase

;; calculate the true anomaly for each time
meananom = 2.d0*!dpi*(1.d0 + (data.bjd - t_periastron)/pars[1] mod 1)
eccanom = exofast_keplereq(meananom, e)
trueanom = 2.d0*atan(sqrt((1.d0 + e)/(1.d0 - e))*tan(eccanom/2.d0))

;; fit the amplitude, offset, and slope analytically
derivs = transpose([[(cos(trueanom+omega) + e*cos(omega))/data.err],$
                    [(data.bjd-t0)/data.err],$
                    [1.d0/data.err]])

npars = n_elements(derivs[*,0])
datarr = replicate(1,npars)#data.rv
errarr = replicate(1,npars)#data.err
b = matrix_multiply(derivs,derivs,/btranspose)
d = total(derivs*datarr/errarr,2)
a = invert(b)#d

pars[4] = a[0]
pars[5] = a[2]
pars[6] = a[1]

if pars[4] lt bK[0] or pars[4] gt bK[1] then return, !values.d_infinity

;; calculate the model
model = pars[4]*(cos(trueanom+omega) + e*cos(omega)) + $
  pars[5] + pars[6]*(data.bjd - t0)




;plot, data.bjd, model-data.rv, yrange=[-500,500],psym=1
;wait, 0.05


;; compute the chi2
chi2 = total(((data.rv - model)/data.err)^2)

return, chi2

end

