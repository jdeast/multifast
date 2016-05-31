function multifast_chi2, pars, determinant=determinant, $
                    modelrv=modelrv, modelflux=modelflux, psname=psname, derived=derived

COMMON chi2_block, rvptrs, transitptrs, priors, debug, rvfit, tranfit, log
COMMON dopptom_block, useDoppTom, dopptomptrs
COMMON sinusoid_block, nSinusoids, frequencies, isDEMC, nAmoebaIterations

if lambdaRange EQ !null then lambdaRange = 0     ;; set to 0 for full lambda range -pi to pi, set to 1 to limit to 0 to pi, set to -1 to limit to -pi to 0

if not isDEMC then nAmoebaIterations++
if not isDEMC and (nAmoebaIterations mod 10000 eq 0) then print, nAmoebaIterations, ' amoeba iterations'

slope     = pars[0]       ;; slope in RV
tc        = pars[1]       ;; transit center time
Period    = 10^(pars[2])  ;; Period of orbit
sqrtecosw = pars[3]       ;; sqrtecosw
sqrtesinw = pars[4]       ;; sqrtesinw
K         = 10^(pars[5])  ;; velocity semi-amplitude
inc       = acos(pars[6]) ;; inclination of the orbit
p         = pars[7]       ;; rp/rstar
ar        = 10^(pars[8])
logg      = pars[9]       ;; density of the Primary (log(rho_star)/(rho_sun))
teff      = pars[10]      ;; effective temperature of the Primary
feh       = pars[11]      ;; metallicity of the Primary
depth2    = pars[12]      ;; secondary depth
vsini     = pars[13]      ;; vsini of the star
lambda    = pars[14]      ;; sky-projected spin-orbit angle
macturb   = pars[15]      ;; non-rotating line width    ;;###
;; pars[15:14+nrv]        ;; RV zero points (gamma)
;; pars[15+nrv:??]        ;; TTV, F0, ndetrend...

;debug = 0
determinant = 1d0

COMMON control, useYY, nInterpModel, exposuretime, transitSpacing, noresiduals, residualOffsetFromTransitBaseline, labelOffset, leftPlotWidthHrs, rightPlotWidthHrs,  $
                transitResidualSpacing, transitResidualLabelOffset, binnedLcYMax, binnedLcYMin, binnedLcResyRange, $
                secondarySpacing, nosecondaryresiduals, residualOffsetFromSecondaryBaseline, secondarylabelOffset, $
                secondaryResidualSpacing, secondaryResidualLabelOffset, binnedSecondaryLcYMax, binnedSecondaryLcYMin, binnedSecondaryLcResyRange, $
                RVresyrange,RVresytick, RMxrange, RMyrange, RMresyrange, RMresytick, BSresyrange, BSresytick, useStellarRadiusPrior, stellarRadiusPriorCenter, stellarRadiusPriorWidth, $
                useImpactParamPrior, impactParamPriorCenter, impactParamPriorWidth, useBlend, blendVals, useBisectorPlot, nTransitsPerPlot, plotRawLightCurves, $
                leftSecondaryPlotWidthHrs, rightSecondaryPlotWidthHrs, DTrangelow, DTrangehigh
                

;;These variable should be set in the main control program, but if not, they will be set here
;
;
;;select method to calculate "priors" on stellar mass and radius from spectroscopic parameters Teff, [Fe/H], and logg
if useYY EQ !null then useYY = 0     ;; set to 1 to use YY-ISOCHRONES or 0 to use Torres relations

options = {ninterp: 1, exptime:1.0d0}
if (nInterpModel EQ !null) or (nInterpModel lt 1) then options.ninterp = 1 else options.ninterp = nInterpModel

;;convert exposure time in seconds to minutes (used only if nInterpModel is greater than 1)
if (exposuretime EQ !null) or (exposuretime lt 0) then options.exptime = 1.0d0 else options.exptime = exposuretime/60.d0;

;;transit plotting parameters
if plotRawLightCurves EQ !null then plotRawLightCurves = 0
if nTransitsPerPlot EQ !null then nTransitsPerPlot = 999999
if nTransitsPerPlot LT 1 then nTransitsPerPlot = 1
if nSinusoids EQ !null then nSinusoids = 0
if transitSpacing EQ !null then transitSpacing = 0.030d0
if noresiduals EQ !null then noresiduals = 1   ;;  0 = show residuals in multi-lightcurve plot, 1 = don't show residuals in multi-lightcurve plot
if residualOffsetFromTransitBaseline EQ !null then residualOffsetFromTransitBaseline = -0.020  ;;used if noresduals = 0 above, negative is below transit baseline, positive above
if labelOffset EQ !null then labelOffset = 0.000d0   ;; negative is below transit baseline, positive above
if leftPlotWidthHrs EQ !null then leftPlotWidthHrs =  4.0d0  ;;positive value
if rightPlotWidthHrs EQ !null then rightPlotWidthHrs = leftPlotWidthHrs  ;;positive value

;;transit residuals plotting parameters
if transitResidualSpacing EQ !null then transitResidualSpacing = transitSpacing
if transitResidualLabelOffset EQ !null then transitResidualLabelOffset = transitResidualSpacing * 0.3d0   ;; negative is below transit baseline, positive above

;;binned transit plotting parameters
if binnedLcYMax EQ !null then binnedLcYMax = 1.005d0
if binnedLcYMin EQ !null then binnedLcYMin = 0.965d0
if binnedLcResyRange EQ !null then binnedLcResyRange = 0.005d0

;;secondary plotting parameters
if leftSecondaryPlotWidthHrs EQ !null then leftSecondaryPlotWidthHrs =  leftPlotWidthHrs  ;;positive value
if rightSecondaryPlotWidthHrs EQ !null then rightSecondaryPlotWidthHrs = rightPlotWidthHrs  ;;positive value
if secondarySpacing EQ !null then secondarySpacing = 0.010d0
if nosecondaryresiduals EQ !null then nosecondaryresiduals = 1   ;;  0 = show residuals in multi-lightcurve plot, 1 = don't show residuals in multi-lightcurve plot
if residualOffsetFromSecondaryBaseline EQ !null then residualOffsetFromSecondaryBaseline = -0.010d0  ;;used if noresduals = 0 above, negative is below transit baseline, positive above
if secondarylabelOffset EQ !null then secondarylabelOffset = 0.005d0   ;; negative is below transit baseline, positive above

;;transit residuals plotting parameters
if secondaryResidualSpacing EQ !null then secondaryResidualSpacing = secondarySpacing
if secondaryResidualLabelOffset EQ !null then secondaryResidualLabelOffset = secondaryResidualSpacing * 0.3d0   ;; negative is below transit baseline, positive above

;;binned transit plotting parameters
if binnedSecondaryLcYMax EQ !null then binnedSecondaryLcYMax = 1.010d0
if binnedSecondaryLcYMin EQ !null then binnedSecondaryLcYMin = 0.990d0
if binnedSecondaryLcResyRange EQ !null then binnedSecondaryLcResyRange = 0.005d0

;;RV plotting parameters 
if RVresyrange EQ !null then RVresyrange = [-100d0,100d0]
if BSresyrange EQ !null then BSresyrange = [-100d0,100d0]
if useBisectorPlot EQ !null then useBisectorPlot = 0   


;;RM plotting parameters
if RMxrange EQ !null then RMxrange = [-leftPlotWidthHrs,rightPlotWidthHrs]
if RMyrange EQ !null then RMyrange = [-100d0,100d0]
if RMresyrange EQ !null then RMresyrange = [-100d0, 100d0]
if RMresytick EQ !null then RMresytick = 50d0

;;Stellar Radius Prior settings
if useStellarRadiusPrior EQ !null then useStellarRadiusPrior = 0
if stellarRadiusPriorCenter EQ !null then stellarRadiusPriorCenter = 1.0d0
if stellarRadiusPriorWidth EQ !null then stellarRadiusPriorWidth = !values.d_infinity

;;Impact Parameter Prior settings
if useImpactParamPrior EQ !null then useImpactParamPrior = 0
if impactParamPriorCenter EQ !null then impactParamPriorCenter = 0.0d0
if impactParamPriorWidth EQ !null then impactParamPriorWidth = !values.d_infinity

if useBlend EQ !null then useBlend = 0 
if blendVals EQ !null then blendVals = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, $
                                        0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]

if useDoppTom EQ !null then useDoppTom = 0 

if useDoppTom then begin
  DTsize = n_elements(dopptomptrs)
  if DTrangelow EQ !null or DTsize ne n_elements(DTrangelow) then DTrangelow = MAKE_ARRAY(DTsize, /DOUBLE, VALUE = -0.06d0)
  if DTrangehigh EQ !null or DTsize ne n_elements(DTrangehigh) then DTrangehigh = MAKE_ARRAY(DTsize, /DOUBLE, VALUE = 0.06d0)
endif

;; initialize the chi^2
chi2 = 0.d0

lambdaRange = 0
if useDoppTom then lambdaRange = (*(dopptomptrs[0])).lambdaRange

e = sqrtecosw^2 + sqrtesinw^2
if e eq 0d0 then omega = !dpi/2d0 $
else omega = atan(sqrtesinw, sqrtecosw)

;; boundary checking
if e ge 1d0 then return, !values.d_infinity
if pars[6] gt 1d0 or pars[6] lt 0d0 then return, !values.d_infinity

;; negative (or zero) macturb is unphysical
if usedopptom and (macturb le 0.d0) then return, !values.d_infinity

;; negative vsini implies lambda is 180 degrees off
if vsini lt 0 then begin
   if isDEMC then begin
      vsini = -vsini
      lambda = lambda + !dpi
      pars[13] = vsini
   endif else begin
      return, !values.d_infinity
   endelse
endif


;; rescale lambda
if isDEMC then begin
    lambda = lambda mod (2.d0*!dpi)
    if lambda lt -!dpi then lambda += 2.d0*!dpi
    if lambda ge  !dpi then lambda -= 2.d0*!dpi
    pars[14] = lambda
    if useDoppTom and (lambdaRange gt 0) and (lambda lt 0.0d0) then return, !values.d_infinity
    if useDoppTom and (lambdaRange lt 0) and (lambda gt 0.0d0) then return, !values.d_infinity
endif

ntransits = n_elements(tranfit)

;; negative frequency amplitude implies phase is 180 degrees off
;; also rescale phase
if tranfit[0] ne -1 and ntransits gt 0 then begin
  if nSinusoids gt 0 then begin
     for i=0, ntransits-1 do begin
       for f=0, nSinusoids-1 do begin
         if isDEMC then begin
           if pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + f*2] lt 0d0 then begin
              pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + f*2] = $
                  -pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + f*2]
              pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] = $
                pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] + !dpi
           endif
           pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] = $
             pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] mod (2.d0*!dpi)
           if pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] lt -!dpi then $
             pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] += 2.d0*!dpi
           if pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] ge  !dpi then $
             pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + 1 + f*2] -= 2.d0*!dpi        
         endif else begin
           if pars[(*(transitptrs[tranfit[i]])).ndx+2+(*(transitptrs[tranfit[i]])).ndetrend + f*2] lt 0d0 then begin
              return, !values.d_infinity
           endif       
         endelse
       endfor
     endfor  
  endif
endif

nsecondaries = 0;
for i=0, ntransits-1 do begin
   if (((*(transitptrs[tranfit[i]])).secondary) eq 1) then nsecondaries += 1;
endfor


;; incorportate "priors" on Mstar, Rstar from YY Isochrones or the Torres relations
;; (only possible with transit data)
if tranfit[0] ne -1 then begin

   ;; if orbiting inside star, that's a problem 
   ;; probably a problem well before that, since tides aren't considered...
   if ar lt (1d0+p) then return, !values.d_infinity 
   
   b = ar*cos(inc)*(1d0-e^2)/(1d0+e*sin(omega))
   if b ge 1d0+p then return, !values.d_infinity
   
   if nsecondaries gt 0 then begin
     bs = ar*cos(inc)*(1d0-e^2)/(1d0-e*sin(omega))
     if bs ge 1d0+p then return, !values.d_infinity
   endif
   
   AU = 215.094177d0 ;; r_sun
   sini = sin(inc)
   G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
   a = (10^logg*(period*86400d0)^2/(4d0*!dpi^2*ar^2*100d0) + $
        K*(Period*86400d0)*sqrt(1d0-e^2)/(2d0*!dpi*sini))/6.9566d8     ;; R_sun
   rstar = a/ar                                                        ;; R_sun
   mp = 2d0*!dpi*(K*86400d0/6.9566d8)*a^2*sqrt(1d0-e^2)/(G*period*sini);; M_sun
   mstar = 4d0*!dpi^2*a^3/(G*period^2) - mp                            ;; M_sun

   ;;  with an estimate of Ma or distance, we can provide another
   ;;  constraint on the stellar radius 
   sigmab = 5.670373d-5/3.839d33*6.9566d10^2 ;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
   lstar = 4d0*!dpi*rstar^2*teff^4*sigmaB    ;; L_sun
;   BC = -5.05244d0 + $           ; +/- 1.89092
;         1.965909d-3*Teff - $    ; +/- 8.869243d-4
;         2.677702d-7*Teff^2 + $  ; +/- 1.382836d-7
;         1.256253d-11*Teff^3     ; +/- 7.166632d-12
;   Mv = -2.5d0*alog10(lstar)+4.60-BC-Av ;; Absolute V-band Magnitude
;
;   flowercoeffs = [-0.370510203809015d5,$
;                   0.385672629965804d5,$
;                   -0.150651486316025d5,$
;                   0.261724637119416d4,$
;                   -0.170623810323864d3]
;   logteff = alog10(teff)
;   BC = total(flowercoeffs*[1d0,logteff,logteff^2,logteff^3,logteff^4])
;  
;   Av = 0.06
;   Mv = -2.5d0*alog10(lstar)+4.732-BC ;; Absolute V-band Magnitude
;   Ma = 8.77d0                           ;; Apparent V-band Magnitude
;   distance = 10d0^((Ma-Mv-Av)/5d0 + 1d0)
;
;   chi2 += ((distance - 110d0)/15d0)^2
;

 
   
  if keyword_set(useYY) then begin  ;;USE YY-ISOCHRONES RATHER THAN TORRES RELATIONS 
     chi2 += massradius_yy2(mstar, feh, teff, rstar, priors[1,11], age=age,debug=debug)
     if ~finite(chi2) then return, !values.d_infinity
     derived = age
  endif else begin  ;;USE TORRES RELATIONS RATHER THAN YY-ISOCHRONES 
     massradius_torres, logg, teff, feh, mstar_torres, rstar_torres
  
     ;; propagate errors in logumstar/logrstar from Torres relation
     umstar = 0.027d0*(mstar_torres*2.3025850929940459d0) 
     urstar = 0.014d0*(rstar_torres*2.3025850929940459d0)
  
     ;; add "prior" penalty
     chi2 += ((mstar-mstar_torres)/umstar)^2
     chi2 += ((rstar-rstar_torres)/urstar)^2     
  endelse
  if not finite(chi2) then stop  
  if keyword_set(useStellarRadiusPrior) then chi2 += ((rstar-stellarRadiusPriorCenter) / stellarRadiusPriorWidth)^2
  if not finite(chi2) then stop  
  if keyword_set(useImpactParamPrior) then chi2 += ((pars[6]*ar-impactParamPriorCenter) / impactParamPriorWidth)^2
  ;;  if keyword_set(useImpactParamPrior) then print, 'b=', pars[6]*ar, 'b_chi2_penalty=', ((pars[6]*ar-impactParamPriorCenter) / impactParamPriorWidth)^2
  if not finite(chi2) then stop  
  ;; Doppler tomography
  if useDoppTom then begin
    chi2 += multifast_dopptom(pars)
    if not finite(chi2) then printf, log, 'ignoring iteration, dopptom chi2=', chi2
    if not finite(chi2) then print, 'ignoring iteration, dopptom chi2=', chi2
    if not finite(chi2) then return, !values.d_infinity
  endif

endif

;; linear limb darkening law for RM effect
;; TRES is more like B+V... does it matter?
u1lin = linld(logg,teff,feh,'V')
if not finite(u1lin) then  return, !values.d_infinity

;; use the best-constrained argument of periastron
phase = exofast_getphase(e,omega,/primary)
tp = tc - period*phase

phase2 = exofast_getphase(e,omega, /secondary)
;; secondary eclipse times

ts = tc + period*(phase2-phase)

;; prepare the plotting device
if keyword_set(debug) or keyword_set(psname) then begin
   aspect_ratio=2
   if keyword_set(psname) then begin
      ;; astrobetter.com tip on making pretty IDL plots
      mydevice=!d.name
      set_plot, 'PS'
      xsize=20;10.5
      ysize=xsize/aspect_ratio
      !p.font=0
      !p.multi=0
      device, filename=psname, /color, bits=24,  /PORTRAIT
      device, xsize=xsize,ysize=ysize,xoffset=2.0,yoffset=2.0
      loadct, 39, /silent
      red = 254
      black = 0
      blue = 75
      green = 159
      gray = 2
      symsize = 0.5
      position1 = [0.05, 0.40, 0.90, 0.95]    ;; data plot
      position2 = [0.05, 0.20, 0.90, 0.40]    ;; residual plot
      
      position1AllLCs =  [0.05, 0.00, 0.90, 0.95]   ;; data and residual plot
      position10AllLCs = [0.02, 0.00, 0.31, 0.95]   ;; data and residual plot
      position11AllLCs = [0.33, 0.00, 0.62, 0.95]   ;; data and residual plot
      position12AllLCs = [0.64, 0.00, 0.93, 0.95]   ;; data and residual plot
      
      position1BinLC = [0.05, 0.40, 0.90, 0.95]    ;; data plot
      position2BinLC = [0.05, 0.20, 0.90, 0.40]    ;; residual plot

      position1rv = [0.05, 0.50, 0.90, 0.95]    ;; data plot
      position2rv = [0.05, 0.35, 0.90, 0.50]    ;; residual plot
      position3rv = [0.05, 0.20, 0.90, 0.35]    ;; residual plot
      charsize=1.3

   endif else begin
      red = '0000ff'x
      black = 'ffffff'x
      blue = 'ff0000'x
      green = '00ff00'x
      gray = '656565'x

      symsize = 1
      screen = GET_SCREEN_SIZE()
      xsize = screen[0]/3d0
      ysize = screen[1]/3d0
      device,window_state=win_state
      if win_state(0) eq 1 then wset, 0 $
      else window, 0, xsize=xsize, ysize=ysize, xpos=0, ypos=screen[1]-ysize
      position1 = [0.10, 0.25, 0.97, 0.95]    ;; data plot
      position2 = [0.10, 0.10, 0.97, 0.25]    ;; residual plot

      position1rv = [0.10, 0.40, 0.97, 0.95]    ;; data plot
      position2rv = [0.10, 0.25, 0.97, 0.40]    ;; residual plot
      position3rv = [0.10, 0.10, 0.97, 0.25]    ;; residual plot
      charsize=1.3
   endelse
;     0 - circle
;     1 - downward arrow (upper limit), base of arrow begins at plot value
;     2 - upward arrow (lower limt)
;     3 - 5 pointed star
;     4 - triangle
;     5 - upside down triangle
;     6 - left pointing arrow
;     7 - right pointing arrow
;     8 - square
;   psyms = [0,8,4]
;   fills = [1,1,1]
;   colors = [black,green,blue]
   psyms = [0,0,8,4,0,8]     ;for special HD189733b plot
   fills = [1,1,1,1,1,1]     ;for special HD189733b plot
   colors = [black,blue,red,red,blue,blue]     ;for special HD189733b plot
;   psyms = [8,0,4]  ;for special TRES plot
;   fills = [1,1,1]  ;for special TRES plot
;   colors = [black,black,blue]   ;for special TRES plot
endif

;; ******** RV ***********
nrv = n_elements(rvptrs)
t0 = tp
ymin1 = !values.d_infinity
ymax1 = -!values.d_infinity
ymin2 = !values.d_infinity
ymax2 = -!values.d_infinity
ymin3 = !values.d_infinity
ymax3 = -!values.d_infinity
mindate = !values.d_infinity
maxdate = -!values.d_infinity

if rvfit[0] ne -1 then begin
   for i=0, n_elements(rvfit)-1 do begin
   
      rvdata = *(rvptrs[rvfit[i]])

      ;; correct all times to target barycenter (?)
      ;; uncomment this line and comment out the next for speed -- see Eastman et al, 2011
      ;; time in target frame (expensive)
      if n_elements(a) eq 0 then rvbjd = rvdata.bjd $
      else rvbjd = bjd2target(rvdata.bjd, inclination=inc, a=a/AU, tp=tp, $
                              period=period, e=e,omega=omega,/primary)
      ;rvbjd = rvdata.bjd  
      rv = rvdata.rv
      rverr = rvdata.err
   
      gamma = pars[16+i]   ;;###
      
      if keyword_set(debug) then begin
         printf, log, 'rvbjd=', rvbjd, '  tp=', tp, '  period=', period, '  gamma=', gamma, '  K=', K, '  e=', e, '  omega=', omega, '  slope=', slope, '  i=', i, '  ar=', ar, '  p=', p,'  vsini=', vsini, '  lambda=', lambda, '  u1=', u1lin, '  t0=', t0 
         print, 'rvbjd=', rvbjd, '  tp=', tp, '  period=', period, '  gamma=', gamma, '  K=', K, '  e=', e, '  omega=', omega, '  slope=', slope, '  i=', i, '  ar=', ar, '  p=', p,'  vsini=', vsini, '  lambda=', lambda, '  u1=', u1lin, '  t0=', t0 
      endif   
         
      modelrv = exofast_rv(rvbjd,tp,period,gamma,K,e,omega,slope=slope,$
                           /rossiter,i=inc,a=ar,p=p,vsini=vsini,lambda=lambda,$
                           u1=u1lin, t0=t0,deltarv=deltarv)
;      if i eq 1 then begin
;        print, 'sqrt(chi^2/dof)=', SQRT(total(((rv - modelrv)/rverr)^2)/(n_elements(rv)-3))
;        stop
;      endif
      
      chi2 += total(((rv - modelrv)/rverr)^2)
      if not finite(chi2) then stop
      (*(rvptrs[rvfit[i]])).residuals = rv - modelrv
      (*(rvptrs[rvfit[i]])).rm = deltarv
      
      if keyword_set(debug) or keyword_set(psname) then begin

         mindate = mindate < min(rvbjd)
         maxdate = maxdate > max(rvbjd)
         ymin1 = ymin1 < min(rv - rverr - gamma - (rvdata.bjd - t0)*slope) 
         ymax1 = ymax1 > max(rv + rverr - gamma - (rvdata.bjd - t0)*slope) 
         ymin2 = ymin2 < min(rv-modelrv)
         ymax2 = ymax2 > max(rv-modelrv)
         ymin3 = ymin3 < min(rvdata.bi)
         ymax3 = ymax3 > max(rvdata.bi)

         if i eq nrv-1 then begin
            
            ;; the pretty model
            nsteps = (maxdate-mindate)*1440d0/15d0 ;; every 15 minutes
            prettytime = mindate + (maxdate-mindate)*dindgen(nsteps)/(nsteps-1.d0)
            prettymodel = exofast_rv(prettytime,tp,period,0,K,e,omega,slope=0d0,$
                                     /rossiter,i=inc,a=ar,p=p,vsini=vsini,$  
                                     lambda=lambda,u1=u1lin,t0=t0,deltarv=deltarvm)
            
            ;; pad the plot to the nearest 5 in the second sig-fig    
            ymin1 = round5(ymin1 < min(prettymodel))
            ymax1 = round5(ymax1 > max(prettymodel))
            
            ;; center the phase curve so transit at phase 0.25, secondary at 0.75
            ;; not standard definition of 'phase' but more informative
            phasejd = (((prettytime - tc) mod period)/period + 1.25d0) mod 1
            sorted = sort(phasejd)
            xrange=[0,1]
            
            ;;  pad the plot to the nearest 5 in the second sig-fig
            ymin2 = round5(ymin2) & ymax2 = round5(ymax2)
            if ymin2 lt -ymax2 then ymax2 = -ymin2
            if ymax2 gt -ymin2 then ymin2 = -ymax2

            ymin3 = round5(ymin3) & ymax3 = round5(ymax3)
            if ymin3 lt -ymax3 then ymax3 = -ymin3
            if ymax3 gt -ymin3 then ymin3 = -ymax3
            


            ;;**** plot the phased model ****
            if useBisectorplot ne 0 then begin
                plot, [0],[0], position=position1rv, xrange=[0,1], $
                      xtickformat='(A1)', ytitle='RV (m/s)', yrange=[ymin1,ymax1],$
                      /ystyle,charsize=charsize            
                oplot, phasejd[sorted], prettymodel[sorted], color=red            

                ;; plot the residuals below
                plot, [0],[0], position=position2rv, /noerase, $
                      xrange=xrange, xtickformat='(A1)', /xstyle, $
                      yrange=RVresyrange, ytitle='O-C!C(m/s)',$
                      /ystyle, yminor=2,yticks=2,charsize=charsize,ytickv=[-RVresytick,0,RVresytick]
                oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  

                ;; and finally, bisectors below that
                plot, [0],[0], position=position3rv, /noerase, xrange=xrange, $
                      xtitle=TeXtoIDL('Orbital Phase + 0.25'), /xstyle, $
                      yrange=BSresyrange, ytitle='Bisec!C(m/s)', /ystyle,$
                      charsize=charsize,yminor=2,yticks=2,ytickv=[-BSresytick,0,BSresytick]
                oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  

                for j=0, n_elements(rvfit)-1 do begin
                   plotsym, psyms[j], symsize, fill=fills[j], color=colors[j]

                   rvdata = *(rvptrs[rvfit[j]])
                   time = (((rvdata.bjd - tc) mod period)/period + 1.25d0) mod 1

                   ;; phased data
                   plot, [0],[0], position=position1rv, xrange=[0,1], $
                         xtickformat='(A1)', ytitle='RV (m/s)', yrange=[ymin1,ymax1],$
                         /ystyle,charsize=charsize,/noerase            
 
                   oploterr, time, rvdata.rv - (rvdata.bjd - t0)*slope - pars[16+j], rvdata.err, 8   ;;###

                   ;; plot the residuals below
                   plot, [0],[0], position=position2rv, /noerase, $
                         xrange=xrange, xtickformat='(A1)', /xstyle, $
                         yrange=RVresyrange, ytitle='O-C!C(m/s)',$
                         /ystyle, yminor=2,yticks=2,charsize=charsize,ytickv=[-RVresytick,0,RVresytick]
                   oploterr, time, rvdata.residuals, rvdata.err, 8 

                   printf, log, 'RMS of RV residuals: ' + strtrim(stddev(rvdata.residuals),2)
                   print, 'RMS of RV residuals: ' + strtrim(stddev(rvdata.residuals),2)
                   bisectorStdDev = stddev(rvdata.bi)
                   ;; and finally, bisectors below that
                   if bisectorStdDev ne 0d0 then begin
                       plot, [0],[0], position=position3rv, /noerase, xrange=xrange, $
                             xtitle=TeXtoIDL('Orbital Phase + 0.25'), /xstyle, $
                             yrange=BSresyrange, ytitle='Bisec!C(m/s)', /ystyle,$
                             charsize=charsize,yminor=2,yticks=2,ytickv=[-BSresytick,0,BSresytick]
                       oploterr, time, rvdata.bi, rvdata.bierr, 8
                       printf, log, 'RMS of ' + rvdata.label + ' bisectors: ' + strtrim(bisectorStdDev,2)
                       print, 'RMS of ' + rvdata.label + ' bisectors: ' + strtrim(bisectorStdDev,2)
                   endif
                endfor
            endif else begin
                position2rv=[0.05, 0.30, 0.90, 0.50]
                ;;**** plot the phased model ****
                plot, [0],[0], position=position1rv, xrange=[0,1], $
                      xtickformat='(A1)', ytitle='RV (m/s)', yrange=[ymin1,ymax1],$
                      /ystyle,charsize=charsize            
                oplot, phasejd[sorted], prettymodel[sorted], color=red            

                ;; plot the residuals below
                plot, [0],[0], position=position2rv, /noerase, $
                      xrange=[0,1],  /xstyle, $
                      yrange=RVresyrange, ytitle='O-C!C(m/s)',$
                      /ystyle, yminor=2,yticks=2,charsize=charsize,ytickv=[-RVresytick,0,RVresytick] ,xtitle=TeXtoIDL('Orbital Phase + 0.25'),xTICKLEN=0.05, yticklen=0.02
                oplot, [-9d9,9d9],[0,0],linestyle=2,color=red

                for j=0, n_elements(rvfit)-1 do begin
                   plotsym, psyms[j], symsize, fill=fills[j], color=colors[j]

                   rvdata = *(rvptrs[rvfit[j]])
                   time = (((rvdata.bjd - tc) mod period)/period + 1.25d0) mod 1

                   ;; phased data
                   plot, [0],[0], position=position1rv, xrange=[0,1], $
                         xtickformat='(A1)', ytitle='RV (m/s)', yrange=[ymin1,ymax1],$
                         /ystyle,charsize=charsize,/noerase            
    ;               if j eq 0 then begin  ;; for special TRES plot
                  oploterr, time, rvdata.rv - (rvdata.bjd - t0)*slope - pars[16+j], rvdata.err, 8   ;;###
    ;               endif else begin ;; for special TRES plot
    ;                oploterr, time, rvdata.rv - (rvdata.bjd - t0)*slope - pars[16+j], rvdata.err/10, 8  ;; for special TRES plot
    ;               endelse  ;; for special TRES plot

                 ;; plot the residuals below
                   plot, [0],[0], position=position2rv, /noerase, $
                         xrange=[0,1], /xstyle, $
                         yrange=RVresyrange, ytitle='O-C!C(m/s)',$
                         /ystyle, yminor=2,yticks=2,charsize=charsize,ytickv=[-RVresytick,0,RVresytick],xtitle=TeXtoIDL('Orbital Phase + 0.25'),xTICKLEN=0.05, yticklen=0.02
    ;;               if j eq 0 then begin  ;;for special TRES plot
                   oploterr, time, rvdata.residuals, rvdata.err, 8 
                endfor
           endelse



            ;;*******************************
;            prettymodel = exofast_rv(prettytime,tp,period,0,K,e,omega,slope=0d0,$
;                                     i=inc,a=ar,p=p,vsini=vsini,$  
;                                     lambda=lambda,u1=u1lin,t0=t0,deltarv=deltarvm)            
            ;; plot the unphased data
            ;;add buffer
            buffer = (maxdate-mindate)*0.02
            mindate = mindate-buffer
            maxdate = maxdate+buffer
            nsteps = (maxdate-mindate)*1440d0/15d0 ;; every 15 minutes
            prettytime = mindate + (maxdate-mindate)*dindgen(nsteps)/(nsteps-1.d0)
            prettymodel = exofast_rv(prettytime,tp,period,0,K,e,omega,slope=slope,$
                                     i=inc,a=ar,p=p,vsini=vsini,$  
                                     lambda=lambda,u1=u1lin,t0=t0,deltarv=deltarvm)
            ;; subtract a zero point to make the plot more legible
            roundto = 10^(strlen(strtrim(floor(maxdate-mindate),2))+1) 
            bjd0 = floor(mindate/roundto)*roundto

            ;; plot the unphased model and data
            if not keyword_set(psname) then begin
               if win_state(1) then wset, 1 $
               else window, 1, xsize=xsize, $
                            ysize=ysize, xpos=0, ypos=screen[1]-2*ysize
            endif
            
            ymin1 = ymin1 < min(rv - rverr - gamma) 
            ymax1 = ymax1 > max(rv + rverr - gamma)
            ymin1 = round5(ymin1 < min(prettymodel)); + (prettytime-t0)*slope))
            ymax1 = round5(ymax1 > max(prettymodel)); + (prettytime-t0)*slope))

            plot, [0],[0], position=position1, $
                  xrange=[mindate,maxdate]-bjd0, xtickformat='(A1)', XTICK_GET = hticks, $
                  ytitle='RV (m/s)', yrange=[ymin1-(ymax1-ymin1)*0.05,ymax1+(ymax1-ymin1)*0.05], /ystyle, charsize = charsize
            maxdate0 = max(hticks)+bjd0
            mindate0 = min(hticks)+bjd0
            buffer0 = (maxdate0-mindate0)*0.02
            mindate0 = mindate0+buffer0
            maxdate0 = maxdate0-buffer0 
            nsteps0 = (maxdate0-mindate0)*1440d0/15d0 ;; every 15 minutes          
            prettytime = mindate0 + (maxdate0-mindate0)*dindgen(nsteps0)/(nsteps0-1.d0)
            prettymodel = exofast_rv(prettytime,tp,period,0,K,e,omega,slope=slope,$
                                     i=inc,a=ar,p=p,vsini=vsini,$  
                                     lambda=lambda,u1=u1lin,t0=t0,deltarv=deltarvm)      
                  
            oplot, prettytime-bjd0, prettymodel, color=red ; + (prettytime-t0)*slope
            ;; plot the residuals below
            plot, [0],[0], position=position2, /noerase, $
                  xrange=[mindate,maxdate]-bjd0, $
                  xtitle=textoidl('BJD_{TDB} - ')+string(bjd0,format='(i7)'), $
                  yrange=RVresyrange,ytitle='O-C (m/s)',/ystyle,yminor=2,yticks=2,ytickv=[-RVresytick,0,RVresytick], charsize = charsize
            oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
            
            for j=0, n_elements(rvfit)-1 do begin
               plotsym, psyms[j], symsize, fill=fills[j], color=colors[j]
               rvdata = *(rvptrs[rvfit[j]])

               ;; unphased data
               plot, [0],[0], position=position1, /noerase, $
                     xrange=[mindate,maxdate]-bjd0, xtickformat='(A1)', $
                     ytitle='RV (m/s)',  yrange=[ymin1-(ymax1-ymin1)*0.05,ymax1+(ymax1-ymin1)*0.05], /ystyle, charsize = charsize
;               oploterr, rvdata.bjd-bjd0, rvdata.rv - (rvdata.bjd - t0)*slope - pars[15+j], rvdata.err, 8  ;;remove slope
;               if j eq 0 then begin  ;; for special TRES plot
               oploterr, rvdata.bjd-bjd0, rvdata.rv - pars[16+j], rvdata.err, 8  ;; retain slope in linear plot  ;;###
;               endif else begin    ;; for special TRES plot
;                 oploterr, rvdata.bjd-bjd0, rvdata.rv - pars[15+j], rvdata.err/10, 8 ;; for special TRES plot
;               endelse  ;; for special TRES plot
               
               ;; plot the residuals below
               plot, [0],[0], position=position2, /noerase, $
                     xrange=[mindate,maxdate]-bjd0, $
                     xtitle=textoidl('BJD_{TDB} - ')+string(bjd0,format='(i7)'), $
                     yrange=RVresyrange,ytitle='O-C (m/s)',/ystyle,yminor=2,yticks=2,ytickv=[-RVresytick,0,RVresytick], charsize = charsize
;               if j eq 0 then begin   ;; for special TRES plot
               oploterr, rvdata.bjd-bjd0, rvdata.residuals, rvdata.err, 8
;               endif else begin   ;; for special TRES plot
;                 oploterr, rvdata.bjd-bjd0, rvdata.residuals, rvdata.err/10, 8   ;; for special TRES plot
;               endelse   ;; for special TRES plot

            endfor

            if not keyword_set(psname) then begin
               if win_state(2) then wset, 2 $
               else window, 2, xsize=xsize, ysize=ysize,$
                            xpos=0, ypos=screen[1]-3*ysize
            endif

            ;; ****plot the RM effect****
            prettymodel = exofast_rv(prettytime,tp,period,0,K,e,omega,slope=0d0,$
                         /rossiter,i=inc,a=ar,p=p,vsini=vsini,$  
                         lambda=lambda,u1=u1lin,t0=t0,deltarv=deltarvm)
            prettyrmtime = (((prettytime - tc) mod period + period) mod period)*24
            toohigh = where(prettyrmtime gt period*12d0)
            if toohigh[0] ne -1 then prettyrmtime[toohigh] -= (period*24d0)
            ;;xrange = [-4,4]
            ;;yrange = [-100,100]
            sorted = sort(prettyrmtime)
            plot, [0],[0],position=position1,$
                  xrange=RMxrange, xtickformat='(A1)', /ystyle, $
                  ytitle='RV (m/s)', yrange=RMyrange,/xs, charsize = charsize
            oplot, prettyrmtime[sorted], deltarvm[sorted], color=red
            ;; plot the residuals below
            plot, [0],[0], position=position2, /noerase,/xs, $
                  xrange=RMxrange,xtitle=textoidl('Time - T_C (hrs)'), $
                  yrange=RMresyrange, ytitle='O-C (m/s)', /ystyle, yminor=2,yticks=2,ytickv=[-RMresytick,0,RMresytick], charsize = charsize
            oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
            
            for j=0, n_elements(rvfit)-1 do begin
               plotsym, psyms[j], symsize, fill=fills[j], color=colors[j]
               rvdata = *(rvptrs[rvfit[j]])

               rmtime = (rvdata.bjd - tc) mod period
               toohigh = where(rmtime gt period/2d0)
               if toohigh[0] ne -1 then rmtime[toohigh] -= period
               toolow = where(rmtime lt -period/2d0)
               if toolow[0] ne -1 then rmtime[toolow] += period
               rmtime *= 24d0
                              
               plot, [0], [0], position=position1, /noerase, $
                     xrange=RMxrange, xtickformat='(A1)', /ystyle, $
                     ytitle='RV (m/s)', yrange=RMyrange,/xs, charsize = charsize
;               if j eq 0 then begin  ;; for special TRES plot
               oploterr, rmtime, rvdata.residuals+rvdata.rm, rvdata.err, 8
;               endif else begin ;; for special TRES plot
;                 oploterr, rmtime, rvdata.residuals+rvdata.rm, rvdata.err/10, 8  ;; for special TRES plot
;               endelse  ;; for special TRES plot
               ;; plot the residuals below
               plot, [0],[0], position=position2, /noerase,/xs, $
                     xrange=RMxrange,xtitle=textoidl('Time - T_C (hrs)'), $
                     yrange=RMresyrange, ytitle='O-C (m/s)', /ystyle, yminor=2,yticks=2,ytickv=[-RMresytick,0,RMresytick], charsize = charsize
;               if j eq 0 then begin  ;; for special TRES plot
               oploterr, rmtime, rvdata.residuals, rvdata.err, 8 
;               endif else begin ;; for special TRES plot
;                oploterr, rmtime, rvdata.residuals, rvdata.err/10, 8 ;; for special TRES plot
;               endelse  ;; for special TRES plot
            endfor





           ;;**** plot the Bisectors vs. RV ****

            ymin = !Values.D_infinity
            ymax = -!Values.D_infinity
            xmin = !Values.D_infinity
            xmax = -!Values.D_infinity    
            hasBisectors = 0            
            for j=0, n_elements(rvfit)-1 do begin
               rvdata = *(rvptrs[rvfit[j]])
               bisectorStdDev = stddev(rvdata.bi)
               if bisectorStdDev eq 0d0 then continue
               hasBisectors = 1
               ymin = min([ymin,rvdata.bi-rvdata.bierr])
               ymax = max([ymax,rvdata.bi+rvdata.bierr])
               xmin = min([xmin,rvdata.rv-rvdata.err])
               xmax = max([xmax,rvdata.rv+rvdata.err])                   
            endfor
            
            if hasBisectors eq 1 then begin
                ymin = ymin - (ymax-ymin)*0.05
                ymax = ymax + (ymax-ymin)*0.05
                xmin = xmin - (xmax-xmin)*0.05
                xmax = xmax + (xmax-xmin)*0.05                
                device, xsize=20, ysize=20, yoffset=2.0, xoffset=2.0, /PORTRAIT
                plot, [0],[0], position=[0.05, 0.05, 0.90, 0.90], xrange=[xmin,xmax], $
                      xtitle=TeXtoIDL('Radial Velocity (m/s)'), /xstyle, $
                      yrange=[ymin,ymax], ytitle='Bisector Span (m/s)', /ystyle, $
                      charsize=charsize;,yminor=2,yticks=2,ytickv=[-BSresytick,0,BSresytick]
                foundfirstbisectors = 0
                for j=0, n_elements(rvfit)-1 do begin
                   plotsym, psyms[j], symsize, fill=fills[j], color=colors[j]
                   rvdata = *(rvptrs[rvfit[j]])
                   bisectorStdDev = stddev(rvdata.bi)
                   if bisectorStdDev ne 0d0 then begin
                       if foundfirstbisectors eq 0 then begin
                          allbisectors = rvdata.bi
                          allrvs = rvdata.rv
                          foundfirstbisectors = 1
                       endif else begin
                          allbisectors = [allbisectors, rvdata.bi]
                          allrvs = [allrvs, rvdata.rv]
                       endelse
                       
                       oploterror, rvdata.rv, rvdata.bi, rvdata.err, rvdata.bierr, psym=8, color=colors[j]
                       ;oplot, rvdata.rv, rvdata.bi, psym=8, color=colors[j]
                   endif
                endfor
                if foundfirstbisectors ne 0 then begin
                   spearman = R_CORRELATE(allrvs, allbisectors, PROBD=prob)
                   speartext = STRING(spearman[0], FORMAT='(F5.2)')
                   spearprob = STRING(prob, FORMAT='(F6.4)')
                   xyouts, (xmax+xmin)/2.0d0, ymax - (ymax-ymin)*0.05, 'Spearman Correlation = '+speartext+', p = '+spearprob, charsize=charsize, alignment=0.5, color = red
                endif
            endif
         endif
      endif
   endfor
endif

;;********************************************************

;;********************* TRANSIT data **********************
;; each transit 

nTransPerPlot = nTransitsPerPlot ;savelocal version so as not to overwrite global version
if nTransPerPlot gt (ntransits-nsecondaries) then nTransPerPlot = ntransits-nsecondaries
nTransitPlots = ceil((double(ntransits-nsecondaries))/double(nTransPerPlot))

;print, 'nTransitPlots=',nTransitPlots
if (tranfit[0] ne -1) and (ntransits - nsecondaries gt 0) then begin
   n  = -1 
   nn = -1
   istart = -1
   pagenum = 0
   plotnum = 0
   nplots = 1
   if plotRawLightCurves and (keyword_set(psname) or keyword_set(debug)) then begin
     !p.multi=[0,3,1]  ; landscape plots on each page showing raw (left) and detrended (center) and residuals (right)
   endif else begin
     !p.multi=0        ;a single column plot for detrended light curves only
   endelse
    
   for i=0, ntransits-1 do begin
      trandata = (*(transitptrs[tranfit[i]])) 
      if (trandata.secondary eq 1) then continue;
      n  += 1
      if plotnum eq 0 then nn += 1
      if n mod nTransPerPlot eq 0 then begin
         n = 0
         if plotnum eq 0 then begin
           pagenum += 1
           lastnn = nn
         endif
         if pagenum eq nTransitPlots then begin
            n = nTransPerPlot - (ntransits - nsecondaries - (plotnum eq 0 ? nn : lastnn))
            ;print, 'lastStartingN=',n
         endif else begin
            n = 0
         endelse
         if keyword_set(psname) or keyword_set(debug) then begin
            if keyword_set(psname) then begin
               !p.font=0
               if plotRawLightCurves then begin
                   xsize = 25
                   ysize = 5 + nTransPerPlot*2
                   maxy = 20
                   if ysize ge maxy then ysize = maxy
                   device, xsize=xsize,ysize=ysize, yoffset=xsize, xoffset = 1.6, /LANDSCAPE
               endif else begin
                   xsize = 20
                   ysize=20/aspect_ratio + (nTransPerPlot)*2
                   maxy = 25
                   if ysize ge maxy then ysize = maxy
                   device, xsize=xsize,ysize=ysize, yoffset=2.0, xoffset=2.0, /PORTRAIT
               endelse
               loadct, 39, /silent
               red = 254
               symsize = 0.33
            endif else if keyword_set(debug) then begin
               !p.multi=0
               xsize = 600
               ysize=(xsize/aspect_ratio + (nTransPerPlot)*150) < screen[1]
               if win_state(3) then wset, 3 $
               else window, 3, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
               red = '0000ff'x
               symsize = 1         
            endif
            spacing = transitSpacing
      
            ;; position keyword required for proper bounding box
            plot, [0],[0],yrange=[1.0d0-(spacing*1.0d0),1.00d0+(spacing*0.5d0)+spacing*((nTransPerPlot-1d0)>0)],$
                  xrange=[-leftPlotWidthHrs,rightPlotWidthHrs],/xstyle,/ystyle, $
                  position=(plotRawLightCurves ? ((plotnum eq 0) ? position10AllLCs : ((plotnum eq 1) ? position11AllLCs : position12AllLCs)) : position1AllLCs), $
                  ytitle=((plotRawLightCurves and (plotnum gt 0)) ? '' : 'Normalized Transit Flux + Constant'), xtitle = TextoIDL('Time - T_C (hrs)'), charsize = charsize*(plotRawLightCurves?2:1), $
                  YTICKFORMAT=((plotRawLightCurves and (plotnum gt 0)) ?"(A1)":""), $
                  title=(plotRawLightCurves ? ((plotnum eq 0) ? 'Undetrended' : ((plotnum eq 1) ? "Detrended" : "Residuals" )) : '')
         endif    
      endif
       ;; If long exposures, create several model points and average 
       npoints = n_elements(trandata.bjd)                                     ;;********
       if options.ninterp gt 1 then begin                                     ;;********
          transitbjd = trandata.bjd#(dblarr(options.ninterp)+1d0) + $          ;;********
                       ((dindgen(options.ninterp)/(options.ninterp-1d0)-0.5d0)/$   ;;********
                        1440d*options.exptime)##(dblarr(npoints)+1d)                ;;********
          modelflux = dblarr(npoints,options.ninterp) + 1d0                       ;;********
          planetvisible = dblarr(npoints,options.ninterp) + 1d0                 ;;********
;;          print, 'interpolating ', npoints, 'points => ', n_elements(transitbjd), ' (exptime=',options.exptime,')'
       endif else begin                                                         ;;********
          transitbjd = trandata.bjd                                             ;;********
          modelflux = dblarr(npoints) + 1d0                                     ;;********
          planetvisible = dblarr(npoints) + 1d0                                 ;;********
;;         print, 'NO interpolation of ', npoints, ' points', ' (exptime=',options.exptime,')'
       endelse                                                                  ;;********

      ;; time of periastron corresponding to transit
      tpn = tp + trandata.epoch*period + pars[trandata.ndx] 

      ;; correct all times to target barycenter (?)
      transitbjd = bjd2target(transitbjd, inclination=inc, $
                              a=a/AU,tp=tpn, period=period, e=e,omega=omega) 
      transitflux = trandata.flux
      transiterr = trandata.err
      
      ;; get the quadratic limb-darkening for each transit 
      ldcoeffs = quadld(logg, teff, feh, trandata.band)
      u1 = ldcoeffs[0]
      u2 = ldcoeffs[1]
      if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity
      
      ;; the impact parameter for each BJD
      z = exofast_getb(transitbjd, i=inc, a=ar, tperiastron=tpn, $
                       period=period, e=e,omega=omega,z=depth)
      
      ;; the flux during transit
;;      npoints = n_elements(transitbjd)                                      ;;********
;;      modelflux = dblarr(npoints) + 1d0                                     ;;********
      primary = where(depth gt 0,complement=secondary)     
      if primary[0] ne -1 then begin
         exofast_occultquad, z[primary], u1, u2, p, mu1
         modelflux[primary] = mu1
      endif
      if n_elements(where(finite(modelflux) eq 0)) gt 1 then begin
         printf, log, "NaNs ahoy!"
         printf, log, u1, u2
         printf, log, logg, teff, feh, trandata.band        
         print, "NaNs ahoy!"
         print, u1, u2
         print, logg, teff, feh, trandata.band
      endif
      ;; add the flux from the planet 
      ;; i.e., everywhere but during secondary eclipse
;;      planetvisible = dblarr(npoints) + 1d0                                 ;;********
      if secondary[0] ne - 1 then begin
         exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
         planetvisible[secondary] = mu1
      endif
      modelflux += depth2*planetvisible - depth2
      
       ;; now integrate the data points (before detrending)                            ;;******
       if options.ninterp gt 1 then modelflux = total(modelflux,2)/options.ninterp     ;;******      
      
      ;; detrending; this could be done analytically 
      ;; but what about covariances with non-linear parameters?
      if trandata.ndetrend gt 0 then begin
         d0 = trandata.detrend - (replicate(1,npoints)##total(trandata.detrend,2))/npoints ;; zero average
         detrend = transpose(pars[trandata.ndx+2:trandata.ndx+2+trandata.ndetrend-1]#d0)
         if nSinusoids gt 0 then begin
            for f=0, nSinusoids-1 do begin
               detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                          pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
            endfor
         endif
         modelflux += detrend ;; add all detrending
      endif else if nSinusoids gt 0 then begin
         detrend = transitbjd*0d0
         for f=0, nSinusoids-1 do begin
            detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                      pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
         endfor   
         modelflux += detrend ;; add all detrending        
      endif else detrend = 0d0
      
      ;; include the blended flux in the filter, blend=(companion flux in filter)/(companion+target flux in filter)
      blend = 0.0d0
      if useBlend eq 1 then begin ;;set to 1 to include all specified filter blends, 0 to NOT include any blends even if specifed
        if trandata.band eq 'U' then blend = blendVals[0] $
        else if trandata.band eq 'B' then blend = blendVals[1] $
        else if trandata.band eq 'V' then blend = blendVals[2] $
        else if trandata.band eq 'R' then blend = blendVals[3] $
        else if trandata.band eq 'I' then blend = blendVals[4] $
        else if trandata.band eq 'J' then blend = blendVals[5] $
        else if trandata.band eq 'H' then blend = blendVals[6] $
        else if trandata.band eq 'K' then blend = blendVals[7] $
        else if trandata.band eq 'Sloanu' then blend = blendVals[8] $
        else if trandata.band eq 'Sloang' then blend = blendVals[9] $
        else if trandata.band eq 'Sloanr' then blend = blendVals[10] $
        else if trandata.band eq 'Sloani' then blend = blendVals[11] $
        else if trandata.band eq 'Sloanz' then blend = blendVals[12] $
        else if trandata.band eq 'b' then blend = blendVals[13] $
        else if trandata.band eq 'v' then blend = blendVals[14] $
        else if trandata.band eq 'y' then blend = blendVals[15] $
        else if trandata.band eq 'Kepler' then blend = blendVals[16] $
        else if trandata.band eq 'CoRoT' then blend = blendVals[17] $
        else if trandata.band eq 'Spit36' then blend = blendVals[18] $
        else if trandata.band eq 'Spit45' then blend = blendVals[19] $
        else if trandata.band eq 'Spit58' then blend = blendVals[20] $
        else if trandata.band eq 'Spit80' then blend = blendVals[21] $
        else begin
           printf, log, 'ERROR: bandname not recognized: ' + trandata.band
           message, 'ERROR: bandname not recognized: ' + trandata.band
        endelse
      endif
      
      
      ;; normalize the model
      modelflux = pars[trandata.ndx+1]*(modelflux*(1d0-blend)+blend)


      if plotnum eq 0 then chi2 += total(((transitflux - modelflux)/transiterr)^2)
      if not finite(chi2) then begin
         printf, log, pars
         print, pars
         stop
      endif

      if keyword_set(debug) or keyword_set(psname) then begin

         ;; *** the fully-sampled model for pretty plots ***
         minbjd = min(trandata.bjd,max=maxbjd)
         minbjd -= 1.0d0 & maxbjd += 1.0d0
         npretty = ceil((maxbjd-minbjd)*1440d0) ;; 1 per minute
         prettybjd = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
         z = exofast_getb(prettybjd, i=inc, a=ar, tperiastron=tpn, $
                          period=period, e=e,omega=omega,z=depth)
         prettyflux = dblarr(npretty) + 1d0
         primary = where(depth gt 0,complement=secondary)
         if primary[0] ne -1 then begin
            exofast_occultquad, z[primary], u1, u2, p, mu1
            prettyflux[primary] = mu1
         endif
      
         ;; add the flux from the planet 
         ;; i.e., everywhere but during secondary eclipse
         planetvisible = dblarr(npretty) + 1d0
         if secondary[0] ne - 1 then begin
            exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
            planetvisible[secondary] = mu1
         endif
         prettyflux += depth2*planetvisible - depth2
         prettytime = (prettybjd - tc - trandata.epoch*period - pars[trandata.ndx])*24.d0
         ;; detrend/f0 are removed from data before plotting

         ;; ************************************

         ;; x axis is time - tc in hours
         time = (trandata.bjd - tc - trandata.epoch*period - pars[trandata.ndx])*24.d0

         ;; subtract blend and normalize the data
         d = (transitflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         modelflux0 = (modelflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         
         ;; subtract the detrending from data
         dd = d - detrend
         modelflux00 = modelflux0 - detrend
                  
;         forprint, trandata.bjd, dd, trandata.err, textout='transit.' + strtrim(i,2) + '.txt',/nocomment, format='(f14.6,x,f7.5,x,f7.5)'

         xmin=min(time,max=xmax)
         plotsym, 0, /fill
         
         ;; transit + model
         if plotRawLightCurves and plotnum eq 0 then begin
            oplot, time, d + spacing*(nTransPerPlot-n-1), psym=8, symsize=symsize
            oplot, time, modelflux0 + spacing*(nTransPerPlot-n-1), thick=2, color=red, linestyle=0    
         endif else if plotRawLightCurves and plotnum eq 2 then begin         
            oplot, time, trandata.flux-modelflux + 1.0d0 + spacing*(nTransPerPlot-n-1), psym=8, symsize=symsize, color=black
            oplot, [-9d9,9d9], [0,0] + 1.0d0 + spacing*(nTransPerPlot-n-1), thick=2, color=red, linestyle=1
         endif else begin
            oplot, time, dd + spacing*(nTransPerPlot-n-1), psym=8, symsize=symsize
            oplot, prettytime, prettyflux + spacing*(nTransPerPlot-n-1), thick=2, color=red, linestyle=0              
         endelse

         ;; residuals
         if (not plotRawLightCurves and not keyword_set(noresiduals)) then begin
            oplot, time, trandata.flux-modelflux + 1.0d0+residualOffsetFromTransitBaseline + spacing*(nTransPerPlot-n-1), $
                   psym=4, symsize=symsize, color=red
            oplot, [-9d9,9d9], [0,0] + 1.0d0+residualOffsetFromTransitBaseline + spacing*(nTransPerPlot-n-1), thick=2, $
                   color=black, linestyle=1
         endif
         ;; centered above transit
         if (not plotRawLightCurves) or (plotnum eq 1) then begin  ;;or (not plotRawLightCurves)
              legendcharsize = 1.0
              if plotRawLightCurves then legendcharsize = 0.6
              xyouts, 0, 1.0d0+labelOffset + spacing*(nTransPerPlot-n-1), trandata.label, charsize=legendcharsize, alignment=0.5, color = red   ;1.0025 *********************
         endif
         
         if plotnum eq 0 then begin
           if nn eq 0 then begin
              alltime = time/24d0 
              allflux = dd
              allerr = trandata.err
              allmodel = modelflux00
           endif else begin
              alltime = [alltime,time/24d0]
              allflux = [allflux,dd]
              allerr = [allerr,trandata.err]
              allmodel = [allmodel,modelflux00]
           endelse
         endif
;         print, 'i=', i, ' ', 'n=', n, ' ', 'nn=', nn, ' ', 'plotnum=', plotnum, ' ', 'pagenum=', pagenum
         if plotRawLightCurves and ((nTransPerPlot eq 1) or ((nn ne 0) and (((n + 1) mod nTransPerPlot) eq 0))) then begin
           plotnum = (plotnum + 1) mod 3
           print, ' '
           if plotnum eq 1 or plotnum eq 2 then begin
              i = istart
           endif else begin
              istart = i
           endelse  
           
         endif         
       endif
      
   endfor
endif      

!p.multi=0   
;;********************* Separate Transit Residual Plot **********************

if (tranfit[0] ne -1) and (ntransits - nsecondaries gt 0) and (keyword_set(psname) or keyword_set(debug)) and (not plotRawLightCurves) then begin  
   if keyword_set(psname) or keyword_set(debug) then begin
      if keyword_set(psname) then begin
         xsize = 20;12.5
         ysize=xsize/aspect_ratio + (ntransits-nsecondaries)*2
         if ysize ge 25 then ysize = 25
         !p.font=0
         !p.multi=0
         device, xsize=xsize,ysize=ysize, yoffset=2.0, xoffset=2.0, /PORTRAIT
         loadct, 39, /silent
         red = 254
         symsize = 0.33
      endif else if keyword_set(debug) then begin
         !p.multi=0
         xsize = 600
         ysize=(xsize/aspect_ratio + (ntransits-nsecondaries)*150) < screen[1]
         if win_state(3) then wset, 3 $
         else window, 3, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
         red = '0000ff'x
         symsize = 1         
      endif
      spacing = transitResidualSpacing

      ;; position keyword required for proper bounding box
      plot, [0],[0],yrange=[1.0d0-(spacing*1.0d0),1.00d0+(spacing*0.5d0)+spacing*((ntransits-nsecondaries-1d0)>0)],$
            xrange=[-leftPlotWidthHrs,rightPlotWidthHrs],/xstyle,/ystyle,position=position1AllLCs,$
            ytitle='Residual Transit Flux + Constant',xtitle = TextoIDL('Time - T_C (hrs)'), charsize = charsize
   endif
   n = -1;
   for i=0, ntransits-1 do begin
       trandata = (*(transitptrs[tranfit[i]])) 
       if (trandata.secondary eq 1) then continue;
       n += 1;
       ;; If long exposures, create several model points and average 
       npoints = n_elements(trandata.bjd)                                     ;;********
       if options.ninterp gt 1 then begin                                     ;;********
          transitbjd = trandata.bjd#(dblarr(options.ninterp)+1d0) + $          ;;********
                       ((dindgen(options.ninterp)/(options.ninterp-1d0)-0.5d0)/$   ;;********
                        1440d*options.exptime)##(dblarr(npoints)+1d)                ;;********
          modelflux = dblarr(npoints,options.ninterp) + 1d0                       ;;********
          planetvisible = dblarr(npoints,options.ninterp) + 1d0                 ;;********
;;          print, 'interpolating ', npoints, 'points => ', n_elements(transitbjd), ' (exptime=',options.exptime,')'
       endif else begin                                                         ;;********
          transitbjd = trandata.bjd                                             ;;********
          modelflux = dblarr(npoints) + 1d0                                     ;;********
          planetvisible = dblarr(npoints) + 1d0                                 ;;********
;;         print, 'NO interpolation of ', npoints, ' points', ' (exptime=',options.exptime,')'
       endelse                                                                  ;;********

      ;; time of periastron corresponding to transit
      tpn = tp + trandata.epoch*period + pars[trandata.ndx] 

      ;; correct all times to target barycenter (?)
      transitbjd = bjd2target(transitbjd, inclination=inc, $
                              a=a/AU,tp=tpn, period=period, e=e,omega=omega) 
      transitflux = trandata.flux
      transiterr = trandata.err
      
      ;; get the quadratic limb-darkening for each transit 
      ldcoeffs = quadld(logg, teff, feh, trandata.band)
      u1 = ldcoeffs[0]
      u2 = ldcoeffs[1]
      if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity
      
      ;; the impact parameter for each BJD
      z = exofast_getb(transitbjd, i=inc, a=ar, tperiastron=tpn, $
                       period=period, e=e,omega=omega,z=depth)
      
      ;; the flux during transit
;;      npoints = n_elements(transitbjd)                                      ;;********
;;      modelflux = dblarr(npoints) + 1d0                                     ;;********
      primary = where(depth gt 0,complement=secondary)     
      if primary[0] ne -1 then begin
         exofast_occultquad, z[primary], u1, u2, p, mu1
         modelflux[primary] = mu1
      endif
      
      ;; add the flux from the planet 
      ;; i.e., everywhere but during secondary eclipse
;;      planetvisible = dblarr(npoints) + 1d0                                 ;;********
      if secondary[0] ne - 1 then begin
         exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
         planetvisible[secondary] = mu1
      endif
      modelflux += depth2*planetvisible - depth2
      
       ;; now integrate the data points (before detrending)                            ;;******
       if options.ninterp gt 1 then modelflux = total(modelflux,2)/options.ninterp     ;;******      
      
      ;; detrending; this could be done analytically 
      ;; but what about covariances with non-linear parameters?
      if trandata.ndetrend gt 0 then begin
         d0 = trandata.detrend - (replicate(1,npoints)##total(trandata.detrend,2))/npoints ;; zero average
         detrend = transpose(pars[trandata.ndx+2:trandata.ndx+2+trandata.ndetrend-1]#d0)
         if nSinusoids gt 0 then begin
            for f=0, nSinusoids-1 do begin
               detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                          pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
            endfor
         endif
         modelflux += detrend ;; add all detrending
      endif else if nSinusoids gt 0 then begin
         detrend = transitbjd*0d0
         for f=0, nSinusoids-1 do begin
            detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                      pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
         endfor   
         modelflux += detrend ;; add all detrending        
      endif else detrend = 0d0
      
      ;; include the blended flux in the filter, blend=(companion flux in filter)/(companion+target flux in filter)
      blend = 0d0
      if useBlend eq 1 then begin ;;set to 1 to include all specified filter blends, 0 to NOT include any blends even if specifed
        if trandata.band eq 'U' then blend = blendVals[0] $
        else if trandata.band eq 'B' then blend = blendVals[1] $
        else if trandata.band eq 'V' then blend = blendVals[2] $
        else if trandata.band eq 'R' then blend = blendVals[3] $
        else if trandata.band eq 'I' then blend = blendVals[4] $
        else if trandata.band eq 'J' then blend = blendVals[5] $
        else if trandata.band eq 'H' then blend = blendVals[6] $
        else if trandata.band eq 'K' then blend = blendVals[7] $
        else if trandata.band eq 'Sloanu' then blend = blendVals[8] $
        else if trandata.band eq 'Sloang' then blend = blendVals[9] $
        else if trandata.band eq 'Sloanr' then blend = blendVals[10] $
        else if trandata.band eq 'Sloani' then blend = blendVals[11] $
        else if trandata.band eq 'Sloanz' then blend = blendVals[12] $
        else if trandata.band eq 'b' then blend = blendVals[13] $
        else if trandata.band eq 'v' then blend = blendVals[14] $
        else if trandata.band eq 'y' then blend = blendVals[15] $
        else if trandata.band eq 'Kepler' then blend = blendVals[16] $
        else if trandata.band eq 'CoRoT' then blend = blendVals[17] $
        else if trandata.band eq 'Spit36' then blend = blendVals[18] $
        else if trandata.band eq 'Spit45' then blend = blendVals[19] $
        else if trandata.band eq 'Spit58' then blend = blendVals[20] $
        else if trandata.band eq 'Spit80' then blend = blendVals[21] $
        else begin
           printf, log, 'ERROR: bandname not recognized: ' + trandata.band
           message, 'ERROR: bandname not recognized: ' + trandata.band
        endelse
      endif
      
      
      ;; normalize the model
      modelflux = pars[trandata.ndx+1]*(modelflux*(1d0-blend)+blend)

;;      chi2 += total(((transitflux - modelflux)/transiterr)^2)

      if keyword_set(debug) or keyword_set(psname) then begin

         ;; *** the fully-sampled model for pretty plots ***
         minbjd = min(trandata.bjd,max=maxbjd)
         minbjd -= 1.0d0 & maxbjd += 1.0d0
         npretty = ceil((maxbjd-minbjd)*1440d0) ;; 1 per minute
         prettybjd = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
         z = exofast_getb(prettybjd, i=inc, a=ar, tperiastron=tpn, $
                          period=period, e=e,omega=omega,z=depth)
         prettyflux = dblarr(npretty) + 1d0
         primary = where(depth gt 0,complement=secondary)
         if primary[0] ne -1 then begin
            exofast_occultquad, z[primary], u1, u2, p, mu1
            prettyflux[primary] = mu1
         endif
      
         ;; add the flux from the planet 
         ;; i.e., everywhere but during secondary eclipse
         planetvisible = dblarr(npretty) + 1d0
         if secondary[0] ne - 1 then begin
            exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
            planetvisible[secondary] = mu1
         endif
         prettyflux += depth2*planetvisible - depth2 
         prettytime = (prettybjd - tc - trandata.epoch*period - pars[trandata.ndx])*24.d0
         ;; detrend/f0 are removed from data before plotting

         ;; ************************************

         ;; x axis is time - tc in hours
         time = (trandata.bjd - tc - trandata.epoch*period - pars[trandata.ndx])*24.d0

         ;; subtract blend and normalize the data
         d = (transitflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         modelflux0 = (modelflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         
         ;; subtract the detrending from data
         d -= detrend
         modelflux0 -= detrend
                  
;         forprint, trandata.bjd, d, trandata.err, textout='transit.' + strtrim(i,2) + '.txt',/nocomment, format='(f14.6,x,f7.5,x,f7.5)'

         xmin=min(time,max=xmax)
         plotsym, 0, /fill
         
         ;; residuals
         oplot, time, trandata.flux-modelflux + 1.0d0 + spacing*(ntransits-nsecondaries-n-1), $
                 psym=8, symsize=symsize, color=black
         oplot, [-9d9,9d9], [0,0] + 1.0d0 + spacing*(ntransits-nsecondaries-n-1), thick=2, $
                 color=red, linestyle=1

         xyouts, 0, 1.0d0+transitResidualLabelOffset + spacing*(ntransits-nsecondaries-n-1), trandata.label,charsize=1.0,alignment=0.5, color = red

      endif
   endfor
endif   


;************* PLOT Combined and Binned Transit Light Curves****************** 

if (tranfit[0] ne -1) and (ntransits - nsecondaries gt 0) and (keyword_set(psname) or keyword_set(debug)) then begin    
   if keyword_set(debug) or keyword_set(psname) then begin
      if keyword_set(psname) then begin
         aspect_ratio=1.5
         xsize=20;;10.5
         ysize=xsize/aspect_ratio
         device, xsize=xsize, ysize=ysize, yoffset=2.0, xoffset=2.0, /PORTRAIT
      endif
      
      bintime, alltime, allflux, allerr, 5d0, binjd, binflux, binerr
      
      xrange = [-leftPlotWidthHrs,rightPlotWidthHrs] ;;[min(binjd)-0.01d0,max(binjd)+0.01d0]*24d0
      xtitle = TextoIDL('Time - T_C (hrs)')

      if not keyword_set(psname) then begin
         if win_state(4) then wset, 4 $
         else window, 4, xpos=2d0*screen[0]/3d0, ypos=2d0*screen[1]/3d0, $
                      xsize=screen[0]/3d0, ysize=screen[1]/3d0
      endif

      ;; binned data
      plot, [0],[0], xstyle=1, ystyle=1,yrange=[binnedLcYMin,binnedLcYMax], xrange=xrange, xtickformat='(A1)',$   ;;/ystyle,
            ytitle='Normalized Transit Flux', position=position1BinLC, charsize = charsize
      oplot, binjd*24, binflux, psym=8, symsize=symsize

      ;; binned model
      bintime, alltime, allmodel, allerr, 5d0, binjdm, binmodel, binerr
      binjd = [-9d9,binjd,9d9]
      binmodel = [1d0, binmodel, 1d0]
      binflux = [1d0, binflux, 1d0]
      oplot, binjd*24d0, binmodel, thick=2, color=red, linestyle=0

      ;; residuals below
      plot, [0],[0], xstyle=1, ystyle=1, yrange=[-binnedLcResyRange,binnedLcResyRange], xrange=xrange, ytitle='O-C',$
       xtitle=xtitle, position=position2BinLC, /noerase, yticks=2,ytickv=[-binnedLcResyRange+0.001,0,binnedLcResyRange-0.001], yminor=binnedLcResyRange*1000, charsize=charsize 
      oplot, binjd*24, binflux-binmodel,psym=8,symsize=symsize
      oplot, [-9d9,9d9],[0,0], color=red, linestyle=2

   endif
endif  
   
   
   
   
;************* PLOT Secondary Transit Light Curves******************
;
;
if (tranfit[0] ne -1) and (nsecondaries gt 0) then begin 
   if keyword_set(psname) or keyword_set(debug) then begin
      if keyword_set(psname) then begin
         xsize = 20;12.5
         ysize=xsize/aspect_ratio + (nsecondaries)*2
         if ysize ge 25 then ysize = 25
         !p.font=0
         !p.multi=0
         device, xsize=xsize,ysize=ysize, yoffset=2.0,  /PORTRAIT
         loadct, 39, /silent
         red = 254
         symsize = 0.33
      endif else if keyword_set(debug) then begin
         !p.multi=0
         xsize = 600
         ysize=(xsize/aspect_ratio + (nsecondaries)*150) < screen[1]
         if win_state(3) then wset, 3 $
         else window, 3, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
         red = '0000ff'x
         symsize = 1         
      endif
      spacing = secondarySpacing

      ;; position keyword required for proper bounding box
      plot, [0],[0],yrange=[1.0d0-(spacing*1.0d0),1.00d0+(spacing*0.5d0)+spacing*((nsecondaries-1d0)>0)],$
            xrange=[-leftSecondaryPlotWidthHrs,rightSecondaryPlotWidthHrs],/xstyle,/ystyle,position=position1AllLCs,$
            ytitle='Normalized Secondary Flux + Constant',xtitle = TextoIDL('Time - T_S (hrs)'),charsize=charsize
   endif

   n = -1;
   for i=0, ntransits-1 do begin
      trandata = (*(transitptrs[tranfit[i]])) 
      if (trandata.secondary eq 0) then continue;
      n += 1;
       ;; If long exposures, create several model points and average 
       npoints = n_elements(trandata.bjd)                                     ;;********
       if options.ninterp gt 1 then begin                                     ;;********
          transitbjd = trandata.bjd#(dblarr(options.ninterp)+1d0) + $          ;;********
                       ((dindgen(options.ninterp)/(options.ninterp-1d0)-0.5d0)/$   ;;********
                        1440d*options.exptime)##(dblarr(npoints)+1d)                ;;********
          modelflux = dblarr(npoints,options.ninterp) + 1d0                       ;;********
          planetvisible = dblarr(npoints,options.ninterp) + 1d0                 ;;********
;;          print, 'interpolating ', npoints, 'points => ', n_elements(transitbjd), ' (exptime=',options.exptime,')'
       endif else begin                                                         ;;********
          transitbjd = trandata.bjd                                             ;;********
          modelflux = dblarr(npoints) + 1d0                                     ;;********
          planetvisible = dblarr(npoints) + 1d0                                 ;;********
;;         print, 'NO interpolation of ', npoints, ' points', ' (exptime=',options.exptime,')'
       endelse                                                                  ;;********

      ;; time of periastron corresponding to transit
      tpn = tp + (trandata.epoch)*period + pars[trandata.ndx] 

      ;; correct all times to target barycenter (?)
      transitbjd = bjd2target(transitbjd, inclination=inc, $
                              a=a/AU,tp=tpn, period=period, e=e,omega=omega) 
      transitflux = trandata.flux
      transiterr = trandata.err
      
      ;; get the quadratic limb-darkening for each transit 
      ldcoeffs = quadld(logg, teff, feh, trandata.band)
      u1 = ldcoeffs[0]
      u2 = ldcoeffs[1]
      if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity
      
      ;; the impact parameter for each BJD
      z = exofast_getb(transitbjd, i=inc, a=ar, tperiastron=tpn, $
                       period=period, e=e,omega=omega,z=depth)
      
      ;; the flux during transit
;;      npoints = n_elements(transitbjd)                                      ;;********
;;      modelflux = dblarr(npoints) + 1d0                                     ;;********
      primary = where(depth gt 0,complement=secondary)     
      if primary[0] ne -1 then begin
         exofast_occultquad, z[primary], u1, u2, p, mu1
         modelflux[primary] = mu1
      endif
      
      ;; add the flux from the planet 
      ;; i.e., everywhere but during secondary eclipse
;;      planetvisible = dblarr(npoints) + 1d0                                 ;;********
      if secondary[0] ne - 1 then begin
         exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
         planetvisible[secondary] = mu1
      endif
      modelflux += depth2*planetvisible - depth2 
      
       ;; now integrate the data points (before detrending)                            ;;******
       if options.ninterp gt 1 then modelflux = total(modelflux,2)/options.ninterp     ;;******      
      
      ;; detrending; this could be done analytically 
      ;; but what about covariances with non-linear parameters?
      if trandata.ndetrend gt 0 then begin
         d0 = trandata.detrend - (replicate(1,npoints)##total(trandata.detrend,2))/npoints ;; zero average
         detrend = transpose(pars[trandata.ndx+2:trandata.ndx+2+trandata.ndetrend-1]#d0)
         if nSinusoids gt 0 then begin
            for f=0, nSinusoids-1 do begin
               detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                          pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
            endfor
         endif
         modelflux += detrend ;; add all detrending
      endif else if nSinusoids gt 0 then begin
         detrend = transitbjd*0d0
         for f=0, nSinusoids-1 do begin
            detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                      pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
         endfor   
         modelflux += detrend ;; add all detrending        
      endif else detrend = 0d0
      
      
      ;; include the blended flux in the filter, blend=(companion flux in filter)/(companion+target flux in filter)
      blend = 0d0
      if useBlend eq 1 then begin ;;set to 1 to include all specified filter blends, 0 to NOT include any blends even if specifed
        if trandata.band eq 'U' then blend = blendVals[0] $
        else if trandata.band eq 'B' then blend = blendVals[1] $
        else if trandata.band eq 'V' then blend = blendVals[2] $
        else if trandata.band eq 'R' then blend = blendVals[3] $
        else if trandata.band eq 'I' then blend = blendVals[4] $
        else if trandata.band eq 'J' then blend = blendVals[5] $
        else if trandata.band eq 'H' then blend = blendVals[6] $
        else if trandata.band eq 'K' then blend = blendVals[7] $
        else if trandata.band eq 'Sloanu' then blend = blendVals[8] $
        else if trandata.band eq 'Sloang' then blend = blendVals[9] $
        else if trandata.band eq 'Sloanr' then blend = blendVals[10] $
        else if trandata.band eq 'Sloani' then blend = blendVals[11] $
        else if trandata.band eq 'Sloanz' then blend = blendVals[12] $
        else if trandata.band eq 'b' then blend = blendVals[13] $
        else if trandata.band eq 'v' then blend = blendVals[14] $
        else if trandata.band eq 'y' then blend = blendVals[15] $
        else if trandata.band eq 'Kepler' then blend = blendVals[16] $
        else if trandata.band eq 'CoRoT' then blend = blendVals[17] $
        else if trandata.band eq 'Spit36' then blend = blendVals[18] $
        else if trandata.band eq 'Spit45' then blend = blendVals[19] $
        else if trandata.band eq 'Spit58' then blend = blendVals[20] $
        else if trandata.band eq 'Spit80' then blend = blendVals[21] $
        else begin
           printf, log, 'ERROR: bandname not recognized: ' + trandata.band
           message, 'ERROR: bandname not recognized: ' + trandata.band
        endelse
      endif
      
      
      ;; normalize the model
      modelflux = pars[trandata.ndx+1]*(modelflux*(1d0-blend)+blend)

      chi2 += total(((transitflux - modelflux)/transiterr)^2)
      if not finite(chi2) then stop

      if keyword_set(debug) or keyword_set(psname) then begin

         ;; *** the fully-sampled model for pretty plots ***
         minbjd = min(trandata.bjd,max=maxbjd)
         minbjd -= 1.0d0 & maxbjd += 1.0d0
         npretty = ceil((maxbjd-minbjd)*1440d0) ;; 1 per minute
         prettybjd = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
         z = exofast_getb(prettybjd, i=inc, a=ar, tperiastron=tpn, $
                          period=period, e=e,omega=omega,z=depth)
         prettyflux = dblarr(npretty) + 1d0
         primary = where(depth gt 0,complement=secondary)
         if primary[0] ne -1 then begin
            exofast_occultquad, z[primary], u1, u2, p, mu1
            prettyflux[primary] = mu1
         endif
      
         ;; add the flux from the planet 
         ;; i.e., everywhere but during secondary eclipse
         planetvisible = dblarr(npretty) + 1d0
         if secondary[0] ne - 1 then begin
            exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
            planetvisible[secondary] = mu1
         endif
         prettyflux += depth2*planetvisible - depth2 
         prettytime = (prettybjd - ts - (trandata.epoch+((tc gt ts)?1.0d0:0.0d0))*period - pars[trandata.ndx])*24.d0
         ;; detrend/f0 are removed from data before plotting

         ;; ************************************

         ;; x axis is time - ts in hours
         time = (trandata.bjd - ts - (trandata.epoch+((tc gt ts)?1.0d0:0.0d0))*period - pars[trandata.ndx])*24.d0

         ;; subtract blend and normalize the data
         d = (transitflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         modelflux0 = (modelflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         
         ;; subtract the detrending from data
         d -= detrend
         modelflux0 -= detrend
                  
;         forprint, trandata.bjd, d, trandata.err, textout='transit.' + strtrim(i,2) + '.txt',/nocomment, format='(f14.6,x,f7.5,x,f7.5)'

         xmin=min(time,max=xmax)
         plotsym, 0, /fill
         ;; transit + model
         oplot, time, d + spacing*(nsecondaries-n-1), psym=8, symsize=symsize
         oplot, prettytime, prettyflux + spacing*(nsecondaries-n-1), thick=2, color=red, $
                linestyle=0
;         bintime, time,  d + spacing*(nsecondaries-n-1), trandata.err, 5d0*24d0, tbinjd, tbinflux, tbinerr
;         oplot, tbinjd, tbinflux, psym=8, symsize=symsize;color=red, symsize=symsize*1.5


         ;; residuals
         if not keyword_set(nosecondaryresiduals) then begin
            oplot, time, trandata.flux-modelflux + 1.0d0+residualOffsetFromSecondaryBaseline + spacing*(nsecondaries-n-1), $
                   psym=4, symsize=symsize, color=red
            oplot, [-9d9,9d9], [0,0] + 1.0d0+residualOffsetFromSecondaryBaseline + spacing*(nsecondaries-n-1), thick=2, $
                   color=black, linestyle=1
         endif
        ;; centered above transit
         xyouts, 0, 1.0d0 + secondarylabelOffset + spacing*(nsecondaries-n-1), trandata.label,charsize=1.0,alignment=0.5, color = red  ;1.0025 *********************
         
         if n eq 0 then begin
            alltime = time/24d0 
            allflux = d
            allerr = trandata.err
            allmodel = modelflux0
         endif else begin
            alltime = [alltime,time/24d0]
            allflux = [allflux,d]
            allerr = [allerr,trandata.err]
            allmodel = [allmodel,modelflux0]
         endelse
      endif
   endfor
endif   
   
;;********************* Separate Secondary Residual Plot **********************

if (tranfit[0] ne -1) and (nsecondaries gt 0) and (keyword_set(psname) or keyword_set(debug)) then begin
   if keyword_set(psname) or keyword_set(debug) then begin
      if keyword_set(psname) then begin
         xsize = 20;12.5
         ysize=xsize/aspect_ratio + (nsecondaries)*2
         if ysize ge 25 then ysize = 25
         !p.font=0
         !p.multi=0
         device, xsize=xsize,ysize=ysize, yoffset=2.0,  /PORTRAIT
         loadct, 39, /silent
         red = 254
         symsize = 0.33
      endif else if keyword_set(debug) then begin
         !p.multi=0
         xsize = 600
         ysize=(xsize/aspect_ratio + (nsecondaries)*150) < screen[1]
         if win_state(3) then wset, 3 $
         else window, 3, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
         red = '0000ff'x
         symsize = 1         
      endif
      spacing = secondaryResidualSpacing

      ;; position keyword required for proper bounding box
      plot, [0],[0],yrange=[1.0d0-(spacing*1.0d0),1.00d0+(spacing*0.5d0)+spacing*((nsecondaries-1d0)>0)],$
            xrange=[-leftSecondaryPlotWidthHrs,rightSecondaryPlotWidthHrs],/xstyle,/ystyle,position=position1AllLCs,$
            ytitle='Residual Secondary Flux + Constant',xtitle = TextoIDL('Time - T_S (hrs)'),charsize=charsize
   endif
   n = -1;
   for i=0, ntransits-1 do begin
       trandata = (*(transitptrs[tranfit[i]])) 
       if (trandata.secondary eq 0) then continue;
       n += 1;
      
       ;; If long exposures, create several model points and average 
       npoints = n_elements(trandata.bjd)                                     ;;********
       if options.ninterp gt 1 then begin                                     ;;********
          transitbjd = trandata.bjd#(dblarr(options.ninterp)+1d0) + $          ;;********
                       ((dindgen(options.ninterp)/(options.ninterp-1d0)-0.5d0)/$   ;;********
                        1440d*options.exptime)##(dblarr(npoints)+1d)                ;;********
          modelflux = dblarr(npoints,options.ninterp) + 1d0                       ;;********
          planetvisible = dblarr(npoints,options.ninterp) + 1d0                 ;;********
;;          print, 'interpolating ', npoints, 'points => ', n_elements(transitbjd), ' (exptime=',options.exptime,')'
       endif else begin                                                         ;;********
          transitbjd = trandata.bjd                                             ;;********
          modelflux = dblarr(npoints) + 1d0                                     ;;********
          planetvisible = dblarr(npoints) + 1d0                                 ;;********
;;         print, 'NO interpolation of ', npoints, ' points', ' (exptime=',options.exptime,')'
       endelse                                                                  ;;********

      ;; time of periastron corresponding to transit
      tpn = tp + (trandata.epoch)*period + pars[trandata.ndx] 

      ;; correct all times to target barycenter (?)
      transitbjd = bjd2target(transitbjd, inclination=inc, $
                              a=a/AU,tp=tpn, period=period, e=e,omega=omega) 
      transitflux = trandata.flux
      transiterr = trandata.err
      
      ;; get the quadratic limb-darkening for each transit 
      ldcoeffs = quadld(logg, teff, feh, trandata.band)
      u1 = ldcoeffs[0]
      u2 = ldcoeffs[1]
      if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity
      
      ;; the impact parameter for each BJD
      z = exofast_getb(transitbjd, i=inc, a=ar, tperiastron=tpn, $
                       period=period, e=e,omega=omega,z=depth)
      
      ;; the flux during transit
;;      npoints = n_elements(transitbjd)                                      ;;********
;;      modelflux = dblarr(npoints) + 1d0                                     ;;********
      primary = where(depth gt 0,complement=secondary)     
      if primary[0] ne -1 then begin
         exofast_occultquad, z[primary], u1, u2, p, mu1
         modelflux[primary] = mu1
      endif
      
      ;; add the flux from the planet 
      ;; i.e., everywhere but during secondary eclipse
;;      planetvisible = dblarr(npoints) + 1d0                                 ;;********
      if secondary[0] ne - 1 then begin
         exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
         planetvisible[secondary] = mu1
      endif
      modelflux += depth2*planetvisible - depth2 
      
       ;; now integrate the data points (before detrending)                            ;;******
       if options.ninterp gt 1 then modelflux = total(modelflux,2)/options.ninterp     ;;******      
      
      ;; detrending; this could be done analytically 
      ;; but what about covariances with non-linear parameters?
      if trandata.ndetrend gt 0 then begin
         d0 = trandata.detrend - (replicate(1,npoints)##total(trandata.detrend,2))/npoints ;; zero average
         detrend = transpose(pars[trandata.ndx+2:trandata.ndx+2+trandata.ndetrend-1]#d0)
         if nSinusoids gt 0 then begin
            for f=0, nSinusoids-1 do begin
               detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                          pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
            endfor
         endif
         modelflux += detrend ;; add all detrending
      endif else if nSinusoids gt 0 then begin
         detrend = transitbjd*0d0
         for f=0, nSinusoids-1 do begin
            detrend += pars[trandata.ndx+2+trandata.ndetrend + f*2]*sin(2*!dpi*frequencies[f]*transitbjd + $
                      pars[trandata.ndx+2+trandata.ndetrend + 1 + f*2])
         endfor   
         modelflux += detrend ;; add all detrending        
      endif else detrend = 0d0
      
      ;; include the blended flux in the filter, blend=(companion flux in filter)/(companion+target flux in filter)
      blend = 0d0
      if useBlend eq 1 then begin ;;set to 1 to include all specified filter blends, 0 to NOT include any blends even if specifed
        if trandata.band eq 'U' then blend = blendVals[0] $
        else if trandata.band eq 'B' then blend = blendVals[1] $
        else if trandata.band eq 'V' then blend = blendVals[2] $
        else if trandata.band eq 'R' then blend = blendVals[3] $
        else if trandata.band eq 'I' then blend = blendVals[4] $
        else if trandata.band eq 'J' then blend = blendVals[5] $
        else if trandata.band eq 'H' then blend = blendVals[6] $
        else if trandata.band eq 'K' then blend = blendVals[7] $
        else if trandata.band eq 'Sloanu' then blend = blendVals[8] $
        else if trandata.band eq 'Sloang' then blend = blendVals[9] $
        else if trandata.band eq 'Sloanr' then blend = blendVals[10] $
        else if trandata.band eq 'Sloani' then blend = blendVals[11] $
        else if trandata.band eq 'Sloanz' then blend = blendVals[12] $
        else if trandata.band eq 'b' then blend = blendVals[13] $
        else if trandata.band eq 'v' then blend = blendVals[14] $
        else if trandata.band eq 'y' then blend = blendVals[15] $
        else if trandata.band eq 'Kepler' then blend = blendVals[16] $
        else if trandata.band eq 'CoRoT' then blend = blendVals[17] $
        else if trandata.band eq 'Spit36' then blend = blendVals[18] $
        else if trandata.band eq 'Spit45' then blend = blendVals[19] $
        else if trandata.band eq 'Spit58' then blend = blendVals[20] $
        else if trandata.band eq 'Spit80' then blend = blendVals[21] $
        else begin
           printf, log, 'ERROR: bandname not recognized: ' + trandata.band
           message, 'ERROR: bandname not recognized: ' + trandata.band
        endelse
      endif
      
      
      ;; normalize the model
      modelflux = pars[trandata.ndx+1]*(modelflux*(1d0-blend)+blend)

;;      chi2 += total(((transitflux - modelflux)/transiterr)^2)

      if keyword_set(debug) or keyword_set(psname) then begin

         ;; *** the fully-sampled model for pretty plots ***
         minbjd = min(trandata.bjd,max=maxbjd)
         minbjd -= 1.0d0 & maxbjd += 1.0d0
         npretty = ceil((maxbjd-minbjd)*1440d0) ;; 1 per minute
         prettybjd = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
         z = exofast_getb(prettybjd, i=inc, a=ar, tperiastron=tpn, $
                          period=period, e=e,omega=omega,z=depth)
         prettyflux = dblarr(npretty) + 1d0
         primary = where(depth gt 0,complement=secondary)
         if primary[0] ne -1 then begin
            exofast_occultquad, z[primary], u1, u2, p, mu1
            prettyflux[primary] = mu1
         endif
      
         ;; add the flux from the planet 
         ;; i.e., everywhere but during secondary eclipse
         planetvisible = dblarr(npretty) + 1d0
         if secondary[0] ne - 1 then begin
            exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
            planetvisible[secondary] = mu1
         endif
         prettyflux += depth2*planetvisible - depth2 
         prettytime = (prettybjd - ts - (trandata.epoch+((tc gt ts)?1.0d0:0.0d0))*period - pars[trandata.ndx])*24.d0
         ;; detrend/f0 are removed from data before plotting

         ;; ************************************

         ;; x axis is time - ts in hours
         time = (trandata.bjd - ts - (trandata.epoch+((tc gt ts)?1.0d0:0.0d0))*period - pars[trandata.ndx])*24.d0
         ;; subtract blend and normalize the data
         d = (transitflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         modelflux0 = (modelflux/pars[trandata.ndx+1] - blend)/(1d0-blend)
         
         ;; subtract the detrending from data
         d -= detrend
         modelflux0 -= detrend
                  
;         forprint, trandata.bjd, d, trandata.err, textout='transit.' + strtrim(i,2) + '.txt',/nocomment, format='(f14.6,x,f7.5,x,f7.5)'

         xmin=min(time,max=xmax)
         plotsym, 0, /fill
         
         ;; residuals
         oplot, time, trandata.flux-modelflux + 1.0d0 + spacing*(nsecondaries-n-1), $
                   psym=8, symsize=symsize, color=black
         oplot, [-9d9,9d9], [0,0] + 1.0d0 + spacing*(nsecondaries-n-1), thick=2, $
                   color=red, linestyle=1
         xyouts, 0, 1.0d0+secondaryResidualLabelOffset + spacing*(nsecondaries-n-1), trandata.label,charsize=1.0,alignment=0.5, color = red
      endif
   endfor 
   
   if keyword_set(debug) or keyword_set(psname) then begin
      if keyword_set(psname) then begin
         aspect_ratio=1.5
         xsize=20;;10.5
         ysize=xsize/aspect_ratio
         device, xsize=xsize,ysize=ysize
      endif
      
      ;;; Binned Secondary Light Curve
      
      bintime, alltime, allflux, allerr, 5d0, binjd, binflux, binerr
      
      xrange = [-leftSecondaryPlotWidthHrs,rightSecondaryPlotWidthHrs] ;;[min(binjd)-0.01d0,max(binjd)+0.01d0]*24d0
      xtitle = TextoIDL('Time - T_S (hrs)')

      if not keyword_set(psname) then begin
         if win_state(4) then wset, 4 $
         else window, 4, xpos=2d0*screen[0]/3d0, ypos=2d0*screen[1]/3d0, $
                      xsize=screen[0]/3d0, ysize=screen[1]/3d0
      endif

      ;; binned data
      plot, [0],[0], xstyle=1, ystyle=1,yrange=[binnedSecondaryLcYMin,binnedSecondaryLcYMax], xrange=xrange, xtickformat='(A1)',$   ;;/ystyle,
            ytitle='Normalized Secondary Flux', position=position1BinLC,charsize=charsize
      oplot, binjd*24, binflux, psym=8, symsize=symsize

      ;; binned model
      bintime, alltime, allmodel, allerr, 5d0, binjdm, binmodel, binerr
      binjd = [-9d9,binjd,9d9]
      binmodel = [1d0, binmodel, 1d0]
      binflux = [1d0, binflux, 1d0]
      oplot, binjd*24d0, binmodel, thick=2, color=red, linestyle=0

      ;; residuals below
      plot, [0],[0], xstyle=1, ystyle=1, yrange=[-binnedSecondaryLcResyRange,binnedSecondaryLcResyRange], xrange=xrange, ytitle='O-C',$
       xtitle=xtitle, position=position2BinLC, /noerase, yticks=2,ytickv=[-binnedSecondaryLcResyRange,0,binnedSecondaryLcResyRange], yminor=binnedSecondaryLcResyRange*1000,charsize=charsize
      oplot, binjd*24, binflux-binmodel,psym=8,symsize=symsize
      oplot, [-9d9,9d9],[0,0], color=red, linestyle=2
   endif
endif 
 
 ;;;*****************END Secondary Light Curve Plots*******************

if useDoppTom and keyword_set(psname) then begin
  pars1 = pars
  ;pars1[14] = -47.5*(!dpi/180d0)
  ;pars1[6] = cos(89d0*(!dpi/180d0))
  ;junk = multifast_dopptom(pars1)  
  incl = acos(pars1[6])
  inc = incl * 180d0 /!dpi
  e = sqrtecosw^2 + sqrtesinw^2
  omega = atan(sqrtesinw,sqrtecosw)*180.d0/!dpi
  ecosw = e*cos(omega*!dpi/180d0)
  esinw = e*sin(omega*!dpi/180d0)
  bp = ar*cos(incl)*(1d0-e^2)/(1d0+esinw)
  ;; Winn 2010
  t14 = Period/!dpi*asin(sqrt((1d0+p)^2 - bp^2)/(sin(incl)*ar))*sqrt(1d0-e^2)/(1d0+esinw) ;; eq 14, 16  
  t23 = Period/!dpi*asin(sqrt((1d0-p)^2 - bp^2)/(sin(incl)*ar))*sqrt(1d0-e^2)/(1d0+esinw) ;; eq 15, 16
  tau = (t14-t23)/2d0
  Tfwhm = t14-tau
  egressphase = Tfwhm*0.5d0 / 10^(pars[2])
  xDTtitle='Velocity (km/s)'
  yDTtitle='Orbital Phase'
  lambdachar= '!9' + String("154B) + '!X'
  darkblue=60
  for d = 0, n_elements(dopptomptrs) - 1 do begin
    DTccf2d = (*(dopptomptrs[d])).ccf2d
    DTmodel = (*(dopptomptrs[d])).model
    DTresid = (*(dopptomptrs[d])).resid
    DTbjd = (*(dopptomptrs[d])).bjd
    DTvel = (*(dopptomptrs[d])).vel
    DTlabel = (*(dopptomptrs[d])).label
    ;set_plot, 'PS'
    xsize=14
    ysize=24
    !p.font=0
    !p.multi=[0,1,1]
    device, Bits_per_Pixel=8, /PORTRAIT
    device, xsize=xsize,ysize=ysize,xoffset=5.0,yoffset=2.0
    loadct, 0, /Silent
    
    nper = round((mean(DTbjd)-pars[1]) / 10^(pars[2]))
    localtc = pars[1] + 10^(pars[2])*nper
    phase = (DTbjd-localtc) / 10^(pars[2])
    vsini = pars[13]/1000.0d0
    
    velplotmax = max([abs(min(DTvel)),abs(max(DTvel))])
    phaseplotmax = max([abs(min(phase)),abs(max(phase))])
    minvel = min(DTvel)
    maxvel = max(DTvel)
    minphase = min(phase)
    maxphase = max(phase)

    modeltext = 'Model ('+lambdachar+'='+STRING((pars1[14]*180d0/!dpi), FORMAT='(D0.1)') + ', i='+STRING(inc, FORMAT='(D0.1)') + ')'
    DTrange = [DTrangelow[d], DTrangehigh[d]]
    DTscale = dblarr(256,2)
    step=(DTrange[1]-DTrange[0])/255
    for i=0, 255 do begin
      DTscale[i,0]=DTrange[0]+i*step
      DTscale[i,1]=DTscale[i,0]
    endfor
    
    multifast_DT_plotimage, -DTccf2d, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                            IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, ytitle = yDTtitle, $
                            XTICKFORMAT='(A1)', title = DTlabel + ' Doppler Data', charsize=charsize, position = [0.05, 0.70, 0.90, 0.95]
    

    loadct, 39, /Silent
    oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
    oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
    oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=4, color=darkblue, linestyle=0
    oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=4, color=darkblue, linestyle=0
    loadct, 0, /Silent
    multifast_DT_plotimage, -DTmodel, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                            IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, $
                            XTICKFORMAT='(A1)', ytitle = modeltext, charsize=charsize, /noerase, position = [0.05, 0.44, 0.90, 0.69]

    loadct, 39, /Silent
    oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
    oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
    oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=4, color=darkblue, linestyle=0
    oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=4, color=darkblue, linestyle=0
    loadct, 0, /Silent                            
    multifast_DT_plotimage, -DTresid, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                            IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, $
                            ytitle='Residuals', xtitle=xDTtitle, charsize=charsize, /noerase, position = [0.05, 0.18, 0.90, 0.43]
    
    loadct, 39, /Silent
    oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
    oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
    oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=4, color=darkblue, linestyle=0
    oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=4, color=darkblue, linestyle=0
    loadct, 0, /Silent       
    
    multifast_DT_plotimage, DTscale, XRANGE = DTrange, YRANGE = [-1,1], $
                        IMGXRANGE = DTrange, IMGYRANGE = [-1,1], RANGE = DTrange, $
                        YTICKS=1, YTickformat='(A1)', xticklen=0.2, xtitle='Fractional Variation', charsize=charsize, /noerase, position = [0.05, 0.08, 0.90, 0.10]
   
 ;   multifast_colorbar, min=DTrangelow[d], max=DTrangehigh[d],DIVISIONS=ceil((DTrange[1]-DTrange[0])*50d0), charsize=charsize, position = [0.05, 0.05, 0.90, 0.10]                     
                            
                            
;    for l=-18, 18 do begin
;      pars1[14] = l * 10 * !dpi / 180d0
;      junk = multifast_dopptom(pars1)
;      titletext = 'Model (lambda = ' + STRING((l*10), FORMAT='(I4)') + ', incl = ' + STRING((inc), FORMAT='(F5.2)') + ')'
;      
;      multifast_DT_plotimage, DTmodel, IMGXRANGE = [min(DTvel),max(DTvel)], IMGYRANGE = [min(phase),max(phase)], RANGE = [-0.06, 0.06], $;RANGE = [min(DTccf2d), max(DTccf2d)], $
;                            title = titletext; position = [0.05, 0.05, 0.90, 0.95]
;    endfor
  endfor
endif

;; add priors (priors[1,*] should be infinity for no prior)
if n_elements(priors) ne 0 then begin
   if n_elements(priors) ne n_elements(pars)*2 then begin
      printf, log, "ERROR: PRIORS must be an 2xNPARS array"
      message, "ERROR: PRIORS must be an 2xNPARS array"
   endif
   chi2 += total(((pars-priors[0,*])/priors[1,*])^2)
endif else begin
   printf, log, 'no priors!'
   print, 'no priors!'
endelse

if keyword_set(psname) then begin
   device, /color, bits=24,  /PORTRAIT
   device, xoffset=2.0,yoffset=2.0
   device, /close
   !p.font=-1
   !p.multi=0
   set_plot, mydevice
endif

;; print out the parameters and chi2 at each step
if keyword_set(debug) then begin
   format='('
   for i=0, n_elements(pars)-1 do format = format + 'f0.6,x,'
   format = format + 'f0.6)'
   printf, log, pars, chi2, format=format
   print, pars, chi2, format=format
endif

if not finite(chi2) then stop

return, chi2

end

