;+
; NAME:
;   EXOFAST_PREFITRV
;
; PURPOSE: 
;   Determines the global fit of RV data for a single planet.
;
; CALLING SEQUENCE:
;    pars = multifast_prefitrv([/CHI2PERDOF,/CIRCULAR,CHI2=,])
;
; OPTIONAL INPUTS:
;    MINPERIOD - The minimum period to consider (default=1 day)
;    MAXPERIOD - The maxmimum period to consider (default=range of
;                input times)
;    NMIN      - The number of local minimums (from the Lomb-Scargle
;                periodogram) to fit fully. This should be increased
;                for highly eccentric planets. Default 5. Only the
;                parameters for the best fit (lowest chi2) are returned.
;    
; OPTIONAL OUTPUTS:
;    CHI2       - The chi2 of the fit (unaffected by /SCALEERR)
;    CHI2PERDOF - The chi2/dof of the fit (unaffected by /SCALEERR)
;
; OPTIONAL KEYWORDS:
;    CIRCULAR - If set, the RV orbit will be assumed circulr (e=0)
;    NOSLOPE  - If set, no slope will be fit to the RV data
;    SCALEERR - If set, the errors will be scaled such that the
;               probability of the chi^2 is 0.5.
;    PLOT     - If set, plots of the best fit model and Lomb Scargle
;               periodograms will be displayed (useful for debgging)
;
; RESULT:
;    The best-fit parameters for the RV fit of a single planet.
;
;     pars[0] = time of transit center
;     pars[1] = period
;     pars[2] = e*cos(omega)
;     pars[3] = e*sin(omega)
;     pars[4] = K
;     pars[5] = gamma
;     pars[6] = slope
;
; COMMON BLOCKS:
;   RV_BLOCK - See exofast.pro for definition
;
; MODIFICATION HISTORY 
;  2012/06 -- Public release -- Jason Eastman (LCOGT)
;  2012/12 -- Add period, perr keywords to fix period
;  2013/01 -- Add error checking if no best fit is found
;  2013/02 -- Add degree of freedom when period is given
;-
function multifast_prefitrv, circular=circular, noslope=noslope, plot=plot,$
                           chisq=chi2, chi2perdof=chi2perdof, scaleerr=scaleerr,$
                           nmin=nmin, minperiod=minperiod,maxperiod=maxperiod, period=period, perr=perr, log=log

;; common rv data (structure with tags, bjd, rv, err)
COMMON rv_block, data

if n_elements(nmin) eq 0 then nmin=5


;; use the lomb-scargle method to get the 5 most promising periods
;; (assumes circular orbit)

if n_elements(period) ne 0 then begin
   periods = period 
   if n_elements(perr) ne 0 then perscale = perr $
   else perscale = 0d0
   nsteps = 1
   dof = 1
endif else begin
   periods = exofast_lombscargle(data.bjd, data.rv, bestpars=lspars, $
                                 scale=perscale,nmin=nmin,$
                                 minperiod=minperiod,maxperiod=maxperiod,$
                                 plot=plot,noslope=noslope)
   dof = 0
endelse

format='(f14.6,x,f12.6,x,f12.6,x,f12.6,x,f12.6,x,f12.6,x,f12.6,x,f12.6,x,f12.6)'

;resolve_all,/continue_on_error ;; makes the type pretty
printf, log, "     T_C            Period        ecosw        esinw        K          gamma          slope      chi^2      chi^2/dof"
print, "     T_C            Period        ecosw        esinw        K          gamma          slope      chi^2      chi^2/dof"

;; now refine with amoeba
nsteps = n_elements(periods)
bestchi2 = !values.d_infinity
amobtol = 1d-8
if keyword_set(circular) then begin
    if keyword_set(noslope) then begin
        dof += n_elements(data.bjd) - 4
        ;; do an ameoba fit for each best period
        ;; fit non-linearly for period
        ;; fit linearly for phi, K, gamma, and slope
        ;; e=0, omega=pi/2 (by definition)
        omega = !dpi/2.d0
        e = 0.d0
        for i=0, nsteps-1 do begin

            ;; define the parameters and scale
            pars = [0.d0,periods[i],0.d0,0.d0, 0.d0, 0.d0, 0.d0]
            scale = [0.d0,perscale[i],0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
            
            ;; amoeba fit
            initpars = multifast_amoeba(amobtol, function_name='exofast_getchi2_rv_fitcirnom',$
                                      p0=pars,scale=scale, log=log)    
            
            ;; if the fit converged, see if it's better than the previous one
            if initpars[0] ne -1 then begin
                chi2 = exofast_getchi2_rv_fitcirnom(initpars)
                if chi2 lt bestchi2 then begin
                    bestchi2 = chi2
                    bestpars = initpars
                endif
                printf, log, initpars, chi2, chi2/dof, format=format
                print, initpars, chi2, chi2/dof, format=format
            endif
        endfor
    endif else begin

        ;; number of degrees of freedom
        dof += n_elements(data.bjd) - 5
        
        ;; do an ameoba fit for each best period
        ;; fit non-linearly for period
        ;; fit linearly for phi, K, gamma, and slope
        ;; e=0, omega=pi/2 (by definition)
        omega = !dpi/2.d0
        e = 0.d0
        for i=0, nsteps-1 do begin
            
            ;; define the parameters and scale
            pars = [0.d0,periods[i],0.d0,0.d0, 0.d0, 0.d0, 0.d0]
            scale = [0.d0,perscale[i],0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
            
            ;; amoeba fit
            initpars = multifast_amoeba(amobtol, function_name='exofast_getchi2_rv_fitcir',$
                                      p0=pars,scale=scale,log=log)    
            
            ;; if the fit converged, see if it's better than the previous one
            if initpars[0] ne -1 then begin
                chi2 = exofast_getchi2_rv_fitcir(initpars)
                if chi2 lt bestchi2 then begin
                    bestchi2 = chi2
                    bestpars = initpars
                endif
                printf, log, initpars, chi2, chi2/dof, format=format
                print, initpars, chi2, chi2/dof, format=format
            endif
        endfor
    endelse

    jd_0 = total(data.bjd/data.err^2)/total(1.d0/data.err^2)
    nper = floor((jd_0 - bestpars[0])/bestpars[1])
    bestpars[0] += bestpars[1]*nper

endif else begin
    if keyword_set(noslope) then begin
        ;; number of degrees of freedom
        dof += n_elements(data.bjd) - 6 
        
        ;; do an ameoba fit for each best period
        ;; fit nonlinearly for t_periastron, period, e, and omega
        ;; fit linearly for K, gamma, and slope
        for i=0, nsteps-1 do begin
            
            jd_0 = total(data.bjd/data.err^2)/total(1.d0/data.err^2)
            
            pars = [jd_0,periods[i],0.d0,0.d0, 0.d0, 0.d0, 0.d0]
            scale = [periods[i]/2.d0,perscale[i],0.5d0,0.5d0, 0.d0, 0.d0, 0.d0]
            initpars=multifast_amoeba(amobtol,function_name='exofast_getchi2_rv_fitnom',p0=pars,scale=scale,log=log)
            
            ;; if the fit is ok, see if it's better than the previous one
            if initpars[0] ne -1 then begin
                chi2 = exofast_getchi2_rv_fitnom(initpars)
                if chi2 lt bestchi2 then begin
                    bestchi2 = chi2
                    bestpars = initpars
                endif
                print, initpars, chi2, chi2/dof, format=format
            endif
            


        endfor     
    
    endif else begin
        ;; number of degrees of freedom
        dof += n_elements(data.bjd) - 7 
        
        ;; do an ameoba fit for each best period
        ;; fit nonlinearly for t_periastron, period, e, and omega
        ;; fit linearly for K, gamma, and slope
        for i=0, nsteps-1 do begin
            
            jd_0 = total(data.bjd/data.err^2)/total(1.d0/data.err^2)
            
            pars = [jd_0,periods[i],0.d0,0.d0, 0.d0, 0.d0, 0.d0]
            scale = [periods[i]/2.75d0,perscale[i],0.5d0,0.5d0, 0.d0, 0.d0, 0.d0]
            initpars=multifast_amoeba(amobtol,function_name='exofast_getchi2_rv_fit',p0=pars,scale=scale,log=log)

            ;; if the fit is ok, see if it's better than the previous one
            if initpars[0] ne -1 then begin
                chi2 = exofast_getchi2_rv_fit(initpars)
                if chi2 lt bestchi2 then begin
                    bestchi2 = chi2
                    bestpars = initpars
                endif
                printf, log, initpars, chi2, chi2/dof, format=format
                print, initpars, chi2, chi2/dof, format=format
            endif           

        endfor
    endelse
    
endelse

if n_elements(bestpars) eq 0 then return, -1

nsteps = 1000
minjd = min(data.bjd,max=maxjd)
prettyjd = minjd + (maxjd - minjd)*dindgen(nsteps)/(nsteps-1.d0)

e = sqrt(bestpars[2]^2 + bestpars[3]^2)
;; compute the model inputs
if e eq 0 then omega = !dpi/2.d0 $
else omega = atan(bestpars[3]/bestpars[2])
if bestpars[2] lt 0 then omega += !dpi
phase = exofast_getphase(e,omega,/primary)
t_periastron = bestpars[0] - bestpars[1]*phase
model = exofast_rv(prettyjd,t_periastron,bestpars[1],bestpars[5],bestpars[4],e,omega,slope=bestpars[6])

if keyword_set(plot) then begin
    plot, prettyjd, model
    oplot, data.bjd, data.rv, psym=1
endif


chi2 = bestchi2
chi2perdof = chi2/dof

if keyword_set(scaleerr) then begin

   errscale = sqrt(chi2/chisqr_cvf(0.5,dof)) ;; scale such that P(chi^2) = 0.5
   printf, log, 'RV Chi^2/dof = ' + strtrim(chi2perdof,2)
   printf, log, 'Scaling RV errors by ' + strtrim(errscale,2)
   print, 'RV Chi^2/dof = ' + strtrim(chi2perdof,2)
   print, 'Scaling RV errors by ' + strtrim(errscale,2)
   data.err *= errscale
   
endif

return, bestpars

end
