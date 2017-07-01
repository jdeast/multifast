;+
; NAME:
;   EXOFAST_GETCHI2_RV_FITCIRNOM
;
; PURPOSE: 
;   Computes the chi^2 of a single, circular planet decribed by PARS, while
;   analytically fitting for Tc, K, gamma
;
; CALLING SEQUENCE:
;    chi2 = exofast_getchi2_rv_fitcirnom(pars)
;
; INPUTS:
;    The best-fit parameters for the RV fit of a single planet.
;
;     pars[0] = time of transit center (will be overwritten)
;     pars[1] = period
;     pars[2] = e*cos(omega) (assumed 0)
;     pars[3] = e*sin(omega) (assumed 0)
;     pars[4] = K            (will be overwritten)
;     pars[5] = gamma        (will be overwritten)
;     pars[6] = slope        (assumed 0)
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

function multifast_getchi2_rv_fitcirnom, pars

COMMON rv_block, data

; data is a structure containing bjd, rv, err
; pars[0] = time of transit center
; pars[1] = period
; pars[2] = e*cos(omega)
; pars[3] = e*sin(omega)
; pars[4] = K
; pars[5] = gamma
; pars[6] = slope

;; as defined by circular orbit
tc = pars[0]
period = pars[1]
pars[2] = 0.d0
pars[3] = 0.d0
e = 0
omega = !dpi/2.d0

phase = exofast_getphase(e,omega,/primary)
TPeriastron = tc - period*phase


mintime = min(data.bjd,max=maxtime)
t0 = (mintime + maxtime)/2.d0

;; fit the amplitude, phase, offset, and slope analytically 
derivs = transpose([[cos(2.d0*!dpi*data.bjd/pars[1])/data.err],$
                    [sin(2.d0*!dpi*data.bjd/pars[1])/data.err],$
                    [1.d0/data.err]])

;; the magic
npars = n_elements(derivs[*,0])
datarr = replicate(1,npars)#data.rv
errarr = replicate(1,npars)#data.err
b = matrix_multiply(derivs,derivs,/btranspose)
d = total(derivs*datarr/errarr,2)
a = invert(b)#d

K = sqrt(a[0]^2 + a[1]^2)
phi = -atan(a[1]/a[0])
if a[0] lt 0 then phi += !dpi

pars[0] = pars[1]*(omega-phi)/(2.d0*!dpi)
pars[4] = K
pars[5] = a[2]

;; calculate the model
model = K*cos(2*!dpi*data.bjd/pars[1] + phi) + pars[6]*(data.bjd-t0) + pars[5]

;; compute the chi2
chi2 = total(((data.rv - model)/data.err)^2)

return, chi2

end

