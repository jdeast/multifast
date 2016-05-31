;+
; NAME:
;   EXOFAST_RVTRAN
;
; PURPOSE: 
;   Simultaneously fits RV and transit data for a single planet. Please
;   cite Eastman et al (in prep) if you make use of this routine in
;   your research. Please report errors or bugs to
;   jeastman@lcogt.net
;
; CALLING SEQUENCE:
;   exofast_rvtran, rvfile, transitfile, band
;
; INPUTS:
;
;   RVFILE      - The filename of the RV data. Must have 3 columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) RV (m/s)
;                   3) err (m/s)
;                 NOTE: The units must be as specified, or the fit
;                 will be wrong.
;   TRANSITFILE - The filename of the transit data. Must have at
;                 least 3 columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) flux
;                   3) err
;                   4) Detrend parameter 1
;                   ....
;                   N+3) Detrend parameter N
;                 NOTE: The units must be as specified, or the fit
;                 will be wrong.
;   BAND        - The bandpass of the observed transit (see quadld.pro)
;
; OPTIONAL INPUTS:
;   PRIORS      - Priors on each of the parameters. Must be an N x 2
;                 elements array. If set, a gaussian prior is applied
;                 with a central value of priors[0,*] and width of
;                 priors[1,*]. For no prior on some parameters, set
;                 the corresponding width to infinity. i.e., 
;                 priors[1,*] = !values.d_infinity
;   MINPERIOD   - The minimum period to consider. The default is 1 day.
;   MAXPERIOD   - The maximum period to consider. The default is the
;                 range of input times.
;   NMIN        - The number of minima in the Lomb-Scargle Periodogram
;                 to fit a full Keplerian orbit. Default is 5. If the
;                 eccentricity is large, this may need to be
;                 increased. The execution time of the initial global
;                 fit is directly proportional to this value (though
;                 is a small fraction of the total time of the MCMC fit).
;   PREFIX      - Each of the output files will have this string as a
;                 prefix. Default is RVFILE without the extension.
;   PNAME       - If set, starting values from exoplanets.org for the
;                 planet PNAME will be used to seed the fit.
;   SIGCLIP     - If set, an iterative fit will be performed excluding
;                 data points more than SIGCLIP*error from the best
;                 fit. Default is no clipping.
;   NTHIN       - If set, only every NTHINth element will be
;                 kept. While this typically doesn't affect the
;                 resultant fit because there is a high correlation
;                 between adjacent steps, the only advantage to this
;                 is improved memory management. Therefore, this
;                 should be used only when the unthinned chain is
;                 not well-mixed with the maximum number of steps.
;
; OPTIONAL KEYWORDS:
;   CIRCULAR  - If set, the fit will be forced to be circular (e=0,
;               omega_star - pi/2)
;   NOSLOPE   - If set, it will not fit a slope to the data.
;   REFIT     - By default, if the output file already exists (see
;               PREFIX), the MCMC fit will not be redone. Set this
;               keyword to refit and overwrite all output files.
;   UPDATE    - Update the local copy of the exoplanets.org file (only
;               applied when PNAME is specified).
;
; OUTPUTS:
;
;   Each of the output files will be preceeded by PREFIX (defined
;   above). The second "?" will be either "f" for flat (if /noslope is
;   set) or "m" for slope. The first "?" will be either "c" for
;   circular if /circular is set or "e" for eccentric.
;
;   pars.?.?.txt   - A text file with the value of each parameter at
;                    each link in the chain.
;   pdfs.?.?.ps    - A postscript plot of each posterior distribution
;                    function. The second "?" will be either "f" for
;                    flat (if /noslope is set) or "m" for slope. The
;                    first "?" will be either "c" for circular if
;                    /circular is set or "e" for eccentric.
;   covars.?.?.ps  - A postscript plot of the covariances between each
;                    parameter. The title of each plot is the
;                    correlation coefficient.
;   median.?.?.tex - The LaTeX source code for a table of the median
;                    values and 68% confidence interval, rounded appropriately.
;   best.?.?.tex   - The LaTeX source code for a table of the best fit
;                    values and 68% confidence interval, rounded appropriately.
;   model.?.?.ps   - An postscript plot of the best-fit models and residuals.
;
; COMMON BLOCKS:
;   BLOCK:
;     RV      - A structure with three tags, corresponding to each
;               column of RVFILE.
;     TRANSIT - A structure with 3 + N tags, corresponding to each
;               column of the TRANSITFILE.
;     PRIORS  - The priors on each of the parameters
;     BAND    - The band of the observed transit
;     DEBUG   - A flag to specify if the data/fit should be
;               plotted. This is intended for debugging only; it
;               dramatically increases the execution time.
;   LD_BLOCK:
;     LOGGS   - The log(g) for all stars
;     TEFFS   - The T_eff for all stars
;     FEHS    - The [Fe/H] for all stars
;     QUADLDA - The linear limb darkening coefficients for all stars
;     QUADLDB - The quadratic limb darkening coefficient for all stars
; EXAMPLES:  
;
; MODIFICATION HISTORY
; 
;  2011/09 -- Public Release, Jason Eastman (OSU)
;-
pro multifast, rvpath=rvpath, tranpath=tranpath, refit=refit, priors=priors, $
             update=update,pname=pname,circular=circular,$
             noslope=noslope, nthin=nthin, secondary=secondary, $
             sigclip=sigclip, prefix=prefix,maxsteps=maxsteps,$
             rossiter=rossiter, debug=debug, display=display, ttv=ttv,redo=redo,$
             minperiod=minperiod,maxperiod=maxperiod, $
             transitfiles=transitfiles, replotmodels=replotmodels, replotall=replotall

COMMON rv_block, rvdata, rvdebug
COMMON chi2_block, rvptrs, transitptrs, priors0, debug0, rvfit, tranfit, log
COMMON dopptom_block, useDoppTom, dopptomptrs
COMMON rmerrscaling, rmerrscale
COMMON sinusoid_block, nSinusoids, frequencies, isDEMC, nAmoebaIterations
isDEMC = 0
nAmoebaIterations = 0L
if dopptomptrs EQ !null then dopptomptrs = ptrarr(1,/allocate_heap)
if rvdebug EQ !null then rvdebug = 0
if useDoppTom EQ !null then useDoppTom = 0
if nSinusoids EQ !null then nSinusoids = 0
if nSinusoids EQ !null then nSinusoids = 0
if frequencies EQ !null then frequencies = [0,0]

if keyword_set(debug) then debug0 = 1 $
else debug0=0
priors0 = priors
chi2func = 'multifast_chi2'
resolve_all, resolve_function=[chi2func,'exofast_random'],/cont,/quiet
G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2009

if keyword_set(circular) then circtxt = 'c.' $
else circtxt = 'e.'

if keyword_set(noslope) then slopetxt = 'f.' $
else slopetxt = 'm.'

;;print to log
if keyword_set(replotmodels) or keyword_set(replotall) then begin
   logname = prefix + circtxt + slopetxt + 'replot.log'
endif else begin
   logname = prefix + circtxt + slopetxt + 'log'
endelse
openw, log, logname, /get_lun 

npars = 16    ;;###
;; read all the RV data (one per telescope) 
;; into structures, referenced through a pointer array
rvfiles = file_search(rvpath, count=nrvtel)
rvptrs = ptrarr(nrvtel,/allocate_heap)
for i=0, nrvtel-1 do *(rvptrs[i]) = readrv(rvfiles[i])
npars += nrvtel

;; read all the transit data (one per transit) 
;; into structures, referenced through a pointer array
;; to reference, e.g., do
;; plot, (*(transitptrs[0])).bjd, (*(transitptrs[0])).flux,/ps,/ys
;
;transitfiles = file_search(tranpath,count=ntransits)
ntransits = n_elements(transitfiles)
nsecondaries = 0;
transitptrs = ptrarr(ntransits,/allocate_heap)
for i=0, ntransits-1 do begin
  *(transitptrs[i]) = readtransit3(transitfiles[i])
  (*(transitptrs[i])).ndx = npars
  npars = npars + 2 + (*(transitptrs[i])).ndetrend + 2*nSinusoids ;; TTV, F0, C0... Cn
  if ((*(transitptrs[i])).secondary eq 1) then nsecondaries += 1;
endfor


;; find unique bands
ntransits = n_elements(transitptrs)
bands = strarr(ntransits)
for j=0, ntransits-1 do bands[j] = (*(transitptrs[j])).band
bands = bands[UNIQ(bands, SORT(bands))]
nbands = n_elements(bands)

;; default is no sigma clipping
if n_elements(sigclip) ne 1 then sigclip = !values.d_infinity
tranfit = -1
rvfit = 0
;; do a global fit on RV only (no RM)
;; baseline file must come before RM file
;; RM filename must contain "rm"
;; multiple telescopes must be named such that IDL sorts them in order: 
;; baseline Tel#1, RM1 tel#1, RM2 tel#1, baseline Tel#2, RM1 tel#2)
period = 10^priors0[0,2]
perr = priors0[1,2] * alog(10.0d) * period
rvzeros = dblarr(nrvtel)
rvzerosscale = dblarr(nrvtel)
i_rm = 0;
if not keyword_set(replotall) and not keyword_set(replotmodels) then begin

    printf, log, 'Performing individual RV best fits'
    print, 'Performing individual RV best fits'
    for i=0, nrvtel-1 do begin
       if strpos(rvfiles[i],'rm') eq -1 then begin
          printf, log, 'RV dataset: ' + rvfiles[i]
          print, 'RV dataset: ' + rvfiles[i]
          rvdata = *(rvptrs[i])
          bestrv = multifast_prefitrv(circular=1,noslope=noslope,$
                                    minperiod=minperiod,maxperiod=maxperiod,$
                                    nmin=nmin,chisq=chi2, period=period, log=log)            
          if (n_elements(bestrv) eq 1) then begin   
            printf, log, '**** RV fit did not converge for ', rvfiles[i], ' ****'                      
            print, '**** RV fit did not converge for ', rvfiles[i], ' ****'
            stop
          endif
          rvzeros[i] = bestrv[5]
          rvzerosscale[i] = 0.0d0;
          ;; the number of degrees of freedom
          dof = n_elements((*(rvptrs[i])).err) - 7
          if keyword_set(noslope) then dof += 1
          if keyword_set(circular) then dof += 2
          
          
    
    ;;      if ((i eq 0) and (dof gt 1)) then begin  ;;for special TRES plot
          if (dof gt 1) then begin
            errscale = sqrt(chi2/chisqr_cvf(0.5,dof))   
          endif else begin
            errscale = 1d0  ;;change to 10 for special TRES plot
          endelse
          (*(rvptrs[i])).err *= errscale
          printf, log, 'Errors scaled by '  + strtrim(errscale,2) + ', mean scaled error' + strtrim(mean((*(rvptrs[i])).err))
          print, 'Errors scaled by '  + strtrim(errscale,2) + ', mean scaled error' + strtrim(mean((*(rvptrs[i])).err))
    ;;      printf, log, 'Scaling RV errors by ' + strtrim(errscale,2)
          
    ;;      print, 'errscale=', errscale
    ;;      printf, log, 'Mean RV error (scaled)' + strtrim(mean((*(rvptrs[i])).err))
    
       endif else begin
          printf, log, 'RM dataset: ' + rvfiles[i]
          print, 'RM dataset: ' + rvfiles[i]
          rvzeros[i] = median((*(rvptrs[i])).rv) ;bestrv[5]  ;;rvzeros = [rvzeros,bestrv[5]] 
          (*(rvptrs[i])).err *= rmerrscale[i_rm];   ;;RM error scaling values must be entered in the calling procedure, one for each RM dataset, and listed in order of processing
          printf, log, 'Errors manually scaled by '  + strtrim(rmerrscale[i_rm],2) + ', mean scaled error' + strtrim(mean((*(rvptrs[i])).err))
          print, 'Errors manually scaled by '  + strtrim(rmerrscale[i_rm],2) + ', mean scaled error' + strtrim(mean((*(rvptrs[i])).err))
          i_rm += 1;
       endelse
    rvfit++
    endfor
    
    printf, log, 'Total RV datasets: ', rvfit - i_rm
    print, 'Total RV datasets: ', rvfit - i_rm
    printf, log, 'Total RM datasets: ', i_rm
    print, 'Total RM datasets: ', i_rm
    
    printf, log, 'Performing individual light curve best fits'
    print, 'Performing individual light curve best fits'
    
    ;; now use RV fit to get a reasonable guess on transit data
    ;; convert to combined parameterization
    tc = priors[0,1];
    if tc eq 0 then tc = bestrv[0];
    utc = priors[1,1];
    if utc eq !values.d_infinity then utc = 0.05d0
    logp = priors[0,2]
    ulogp = priors[1,2]
    if ulogp eq !values.d_infinity then ulogp = 0.0001d0
    period = 10^priors[0,2]; bestrv[1];
    ecosw = bestrv[2]
    esinw = bestrv[3]
    e = sqrt(ecosw^2 + esinw^2)
    if e eq 0 then omega = !dpi/2d0 $
    else omega = atan(esinw, ecosw)
    sqrtecosw = sqrt(e)*cos(omega)
    sqrtesinw = sqrt(e)*sin(omega)
    
    logk = alog10(bestrv[4]);priors[0,5];
    gamma = bestrv[5];0;
    slope = bestrv[6];0;
    
    ;; use Hot Jupiter-like guesses
    cosi = priors[0,6];
    if cosi eq 0.0d0 then cosi = 0.05d0
    ucosi = priors[1,6];
    if ucosi eq !values.d_infinity then ucosi = 0.05d0
    p = priors[0,7];
    if p eq 0d0 then p = 0.1d0
    up = priors[1,7];
    if up eq !values.d_infinity then up = 0.01d0
    loga = priors[0,8];
    if loga eq 0d0 then loga = 1d0
    uloga = priors[1,8];
    if uloga eq !values.d_infinity then uloga = 0.05d0
    logg = priors[0,9]
    ulogg = priors[1,9]
    if ulogg eq !values.d_infinity then ulogg = 0.5d0
    teff = priors[0,10]
    uteff = priors[1,10]
    if uteff eq !values.d_infinity then uteff = 500d0
    feh = priors[0,11]
    ufeh = priors[1,11]
    if ufeh eq !values.d_infinity then ufeh = 1d0
    depth2 = priors[0,12]
    udepth2 = priors[1,12]
    if udepth2 eq !values.d_infinity then udepth2 = 0.003d0
    vsini = priors[0,13]
    uvsini = priors[1,13]
    if uvsini eq !values.d_infinity then uvsini = vsini * 0.2d0
    lambda = priors[0,14]
    ulambda = priors[1,14]
    if ulambda eq !values.d_infinity then ulambda = !dpi
    macturb = priors[0,15]   ;;###
    umacturb = priors[1,15]
    if umacturb eq !values.d_infinity then umacturb = 3000d0   ;;###
    
    dindex = 16+nrvtel  ;;###
    dwindex = 16+nrvtel ;;###
    coeffs = [priors[0, dindex++], priors[0, dindex++]]
    widths = [(priors[1, dwindex] eq !values.d_infinity) ? (0d0*dwindex++)+0.1 : priors[1, dwindex++], (priors[1, dwindex] eq !values.d_infinity) ? (0d0*dwindex++)+0.1 : priors[1, dwindex++]]
    for i=0, ntransits-1 do begin
       if i ne 0 then begin
          coeffs = [coeffs, priors[0, dindex++], priors[0, dindex++]]
          widths = [widths, (priors[1, dwindex] eq !values.d_infinity) ?  (0d0*dwindex++)+0.1 : priors[1, dwindex++], (priors[1, dwindex] eq !values.d_infinity) ?  (0d0*dwindex++)+0.1 : priors[1, dwindex++]]
       endif
       ndetrends = (*(transitptrs[i])).ndetrend + nSinusoids*2
       for j=1, ndetrends do begin
          coeffs = [coeffs, priors[0, dindex++]]
          widths = [widths, (priors[1, dwindex] eq !values.d_infinity) ?  (0d0*dwindex++)+0.1 : priors[1, dwindex++]]
       endfor 
    endfor
    
    ;;          0     1    2       3         4        5      6    7       8      9      10      11    12      13     14       15   16:15+nrvzeros 16+nrvzeros:npars-1    ;;###
    
    allpars = [slope, tc, logp  ,sqrtecosw,sqrtesinw,logk,  cosi,  p,    loga,   logg,  teff,   feh,  depth2, vsini, lambda,  macturb, rvzeros,  coeffs]    ;;###
    
    masterscale = [0, utc, ulogp , 0d0,      0d0,    0d0,  ucosi, up,   uloga, ulogg,  uteff,  ufeh,   0d0,     0d0,   0d0,     0d0, rvzerosscale, widths]  ;;###
    
    ;; list parameters for confirmation of setup
    ;;printf, log, '        RV_Slope         Tc           log(Period)      sqrtecosw        sqrtesinw        log(K)          cos(i)         rp/rstar       log(a/rstar)     log(g)          Teff            [Fe/H]      secondary_depth     vsini         orbit_angle    RV_zero_point       TTV          Baseline'
    ;;print, 'Prior Centers'
    ;;printf, log, allpars
    ;;print, 'Prior Widths'
    ;;printf, log, masterscale 
    
    
    
    scale = dblarr(npars)
    tofitrv = [0,1,2,3,4,5,14,15,lindgen(nrvtel)+16] ;;###
    
    scale[tofitrv] = masterscale[tofitrv]
    
    rvfit = lindgen(nrvtel) ;; fit all the data sets
    
    if keyword_set(circular) then begin
       dof += 2
       scale[3:4] = 0d0
    endif
    
    amobtol = 1d-8
    titleformat='(A24,x,A14,x,A16,x,A11,x,A8,x,A9,x,A5,x,A7,x,A7,x,A6,x,A6,x,A9,x,A7,x,A8,x,A8)'
    format='(A26,x,f14.6,x,f14.6,x,f11.6,x,f8.1,x,f9.2,x,f5.2,x,f8.6,x,f8.6,x,f6.4,x,f6.4,x,f9.3,x,I5,x,f8.3,x,f8.3)'
    holdUseDoppTom = useDoppTom  ;; save dopptom request
    useDoppTom = 0  ;;disable dopptom chi2 contribution during transit only fit
;    for i=0, ntransits-1 do begin
;       tranfit = i
;       nAmoebaIterations = 0L
;       ;; calculate the appropriate epoch for this transit
;
;       if ((*(transitptrs[i])).secondary eq 0) then (*(transitptrs[i])).epoch = round((mean((*(transitptrs[i])).bjd) - tc)/10^priors[0,2]) $
;       else (*(transitptrs[i])).epoch = round(((mean((*(transitptrs[i])).bjd) - tc)/10^priors[0,2])-0.5d0)
;       
;       
;    
;       ;; fit the Transit data alone
;       trandata  = *(transitptrs[i])
;    
;       ;; quadratic limb darkening parameters
;       
;       ldcoeffs = quadld(logg, teff, feh, trandata.band)
;       u1 = ldcoeffs[0]
;       u2 = ldcoeffs[1]
;    
;       tranpars = allpars
;       if tranpars[trandata.ndx+1] eq 0 then tranpars[trandata.ndx+1] = 1d0 ;; f0
;    
;       ;; if no TTVs, just fit F0, CN for each transit; otherwise, fit it too.
;       if (trandata.secondary eq 0) then begin
;          masterscale[12] = 0.0d0;
;          if keyword_set(ttv) then tofit = [6,7,8,9,10,11,lindgen(trandata.ndetrend+2+nSinusoids*2)+trandata.ndx] $
;          else tofit = [1,6,7,8,9,10,11,lindgen(trandata.ndetrend+1+nSinusoids*2)+trandata.ndx+1]
;       endif else begin
;          masterscale[12] = udepth2
;          if keyword_set(ttv) then tofit = [6,7,8,9,10,11,12,lindgen(trandata.ndetrend+2+nSinusoids*2)+trandata.ndx] $
;          else tofit = [1,6,7,8,9,10,11,12,lindgen(trandata.ndetrend+1+nSinusoids*2)+trandata.ndx+1]   
;       endelse      
;       
;       transcale = dblarr(npars)
;       transcale[tofit] = masterscale[tofit]
;    
;       rvfit = -1
;       
;       repeat begin
;          
;          besttran = multifast_amoeba(1d-6,function_name=chi2func,$
;                                    p0=tranpars,scale=transcale,nmax=500000,log=log)
;                                    
;          if besttran[0] eq -1 then stop
;          chi2 = call_function(chi2func,besttran,modelflux=modelflux,psname=prefix + 'transit.' + strtrim(i,2) + '.' + circtxt + slopetxt + 'ps')
;          
;          outliers = where(abs(trandata.flux - modelflux)/trandata.err gt $
;                           sigclip,noutliers)
;          if noutliers gt 0 then begin
;             trandata.err[outliers] = !values.d_infinity
;             (*(transitptr[i])).err[outliers] = !values.d_infinity
;          endif
;       endrep until noutliers eq 0
       
    for i=0, ntransits-1 do begin
       tranfit = i
       nAmoebaIterations = 0L
       
       transitndx = (*(transitptrs[i])).ndx
       (*(transitptrs[i])).ndx = 16 + nrvtel
       transitendndx = transitndx + 2 + (*(transitptrs[i])).ndetrend + 2*nSinusoids - 1
       
       ;; calculate the appropriate epoch for this transit
       if ((*(transitptrs[i])).secondary eq 0) then (*(transitptrs[i])).epoch = round((mean((*(transitptrs[i])).bjd) - tc)/10^priors[0,2]) $
       else (*(transitptrs[i])).epoch = round(((mean((*(transitptrs[i])).bjd) - tc)/10^priors[0,2])-0.5d0)
     
       ;; fit the Transit data alone
       trandata  = *(transitptrs[i])
    
       ;; quadratic limb darkening parameters
       ldcoeffs = quadld(logg, teff, feh, trandata.band)
       u1 = ldcoeffs[0]
       u2 = ldcoeffs[1]
    
       tranpars = [allpars[0:15+nrvtel],allpars[transitndx:transitendndx]]
       priors0 = [[priors[*,0:15+nrvtel]],[priors[*,transitndx:transitendndx]]]
       if tranpars[trandata.ndx+1] eq 0 then tranpars[trandata.ndx+1] = 1d0 ;; f0
    
       ;; if no TTVs, just fit F0, CN for each transit; otherwise, fit it too.
       if (trandata.secondary eq 0) then begin
          masterscale[12] = 0.0d0;
          if keyword_set(ttv) then tofit = [6,7,8,9,10,11,lindgen(trandata.ndetrend+2+nSinusoids*2)+trandata.ndx] $
          else tofit = [1,6,7,8,9,10,11,lindgen(trandata.ndetrend+1+nSinusoids*2)+trandata.ndx+1]
       endif else begin
          masterscale[12] = udepth2
          if keyword_set(ttv) then tofit = [6,7,8,9,10,11,12,lindgen(trandata.ndetrend+2+nSinusoids*2)+trandata.ndx] $
          else tofit = [1,6,7,8,9,10,11,12,lindgen(trandata.ndetrend+1+nSinusoids*2)+trandata.ndx+1]   
       endelse      
       
       transcale = dblarr(n_elements(tranpars))
       tranmasterscale = [masterscale[0:15+nrvtel],masterscale[transitndx:transitendndx]]
       transcale[tofit] = tranmasterscale[tofit]
    
       rvfit = -1
       
       repeat begin
          
          besttran = multifast_amoeba(1d-6,function_name=chi2func,$
                                    p0=tranpars,scale=transcale,nmax=500000,log=log)
                                    
          if besttran[0] eq -1 then stop
          chi2 = call_function(chi2func,besttran,modelflux=modelflux,psname=prefix + 'transit.' + strtrim(i,2) + '.' + circtxt + slopetxt + 'ps')
          
          outliers = where(abs(trandata.flux - modelflux)/trandata.err gt $
                           sigclip,noutliers)
          if noutliers gt 0 then begin
             trandata.err[outliers] = !values.d_infinity
             (*(transitptr[i])).err[outliers] = !values.d_infinity
          endif
       endrep until noutliers eq 0       
       
       ;; scale the transit errors so P(chi^2) = 0.5
       dof = n_elements(where(finite(trandata.err))) - n_elements(tofit) + 4 ;; 4 priors
       errscale = sqrt(chi2/chisqr_cvf(0.5,dof))
       (*(transitptrs[i])).err *= errscale
       
       ;; print results
       
       ;;titleformat='(                   A24                 A14               A16          A11        A8      A9     A5      A7     A7     A6    A6    A9     A7        A8        A8)'
       ;;format='(                        A26                 f14.6             f14.6        f11.6      f8.1    f9.2   f5.2    f8.6   f8.6   f6.4  f6.4  f9.3   I5       f8.3       f8.3)'
       
       if i eq 0 then printf, log, 'transit file name', '      T_0    ', '      T_C      ','period  ','epoch', 'a/R*','i  ','depth', 'RMS', 'u1', 'u2','chi2', 'dof', 'chi2/dof', 'errscale', format=titleformat
       if i eq 0 then print,       'transit file name', '      T_0    ', '      T_C      ','period  ','epoch', 'a/R*','i  ','depth', 'RMS', 'u1', 'u2','chi2', 'dof', 'chi2/dof', 'errscale', format=titleformat

       printf, log, transitfiles[i],besttran[1],besttran[1]+(10^besttran[2])*(trandata.epoch+((trandata.secondary eq 0) ? 0.0d0 : 0.5d0)), 10^besttran[2],$
                    trandata.epoch+((trandata.secondary eq 0) ? 0.0d0 : 0.5d0), 10^(besttran[8]),acos(besttran[6])*180d0/!dpi,(trandata.secondary eq 0) ? besttran[7]^2:besttran[12], stddev(trandata.flux - modelflux), $
                    ((trandata.secondary eq 0) ? u1 : 0.0d0), ((trandata.secondary eq 0) ? u2 : 0.0d0), chi2, dof, chi2/dof, errscale, format=format  
       print,       transitfiles[i],besttran[1],besttran[1]+(10^besttran[2])*(trandata.epoch+((trandata.secondary eq 0) ? 0.0d0 : 0.5d0)), 10^besttran[2],$
                    trandata.epoch+((trandata.secondary eq 0) ? 0.0d0 : 0.5d0), 10^(besttran[8]),acos(besttran[6])*180d0/!dpi,(trandata.secondary eq 0) ? besttran[7]^2:besttran[12], stddev(trandata.flux - modelflux), $
                    ((trandata.secondary eq 0) ? u1 : 0.0d0), ((trandata.secondary eq 0) ? u2 : 0.0d0), chi2, dof, chi2/dof, errscale, format=format  

       
       (*(transitptrs[i])).ndx = transitndx
    ;   multifast_demc, besttran, 'chi2_all6', pars, chi2=chi2, tofit=tofit,$
    ;                 nthin=nthin,maxsteps=maxsteps, burnndx=burnndx
    ;   save, pars, burnndx, chi2, filename='transit.' + strtrim(i,2) + '.idl'
       
    ;if trandata.ndetrend ne 0 then begin
       ;; update the combined parameter array for the Tc, F0, and C_i
       allpars[transitndx:transitndx+trandata.ndetrend+1+nSinusoids*2] = $    ;;ndetrend-1
          tranpars[16+nrvtel:16+nrvtel+trandata.ndetrend+1+nSinusoids*2]    ;;ndetrend-1
    ;endif
    
    endfor
    
    priors0 = priors
    
    printf, log, 'Total primary transits: ', ntransits - nsecondaries
    print, 'Total primary transits: ', ntransits - nsecondaries
    printf, log, 'Total secondary transits: ', nsecondaries
    print, 'Total secondary transits: ', nsecondaries
    
    ;; find the MCMC scale of the RV parameters
    
    ;; exclude parameters as desired
    tofit = [indgen(6),indgen(nrvtel)+16]  ;;###
    
    if keyword_set(circular) then tofit = tofit[where(tofit ne 3 and tofit ne 4)]
    if keyword_set(noslope) then tofit = tofit[where(tofit ne 0)]
    
    rvfit = indgen(nrvtel)
    
    tranfit = -1
    ;debug0=1
    rvscale = multifast_getmcmcscale(allpars,chi2func,tofit=tofit,log=log)
    
    for i=0, nrvtel-1 do begin
      rvzerosscale[i] = rvscale[where(tofit eq 16+i)]  ;;###
    endfor

    ;print, 'rvzeroscale=', rvzerosscale
    
    useDoppTom = holdUseDoppTom  ;;restore dopptom request
    
    allpars[1] = tc
    ;ugamma = rvscale[where(tofit eq 16)] ;;###
    
    ;; Replicate coefficient names repeated for each transit
    for i=0, ntransits-1 do begin
       if i eq 0 then begin
          coeffnames = [['TTV_{' +  strtrim(i,2) + '}','Transit Timing Variation (days)'],$
                        ['F_{' + strtrim(i,2) + ',0}','Baseline Flux']] 
          tcnames = ['T_{C,' +  strtrim(i,2) + '}','Mid-transit time (\bjdtdb)']
       endif else begin
          coeffnames = [[coeffnames],$
                        ['TTV_{' +  strtrim(i,2) + '}','Transit Timing Variation (days)'],$
                        ['F_{' + strtrim(i,2) + ',0}','Baseline Flux']]
          tcnames = [[tcnames],['T_{C,' +  strtrim(i,2) + '}','Mid-transit time (\bjdtdb)']]
       endelse
       nlin = (*(transitptrs[i])).ndetrend + nSinusoids*2
    ;   if nlin ge 1 then coeffnames=[[coeffnames],['C_{'+strtrim(i,2)+','+strtrim(indgen(nlin),2)+'}','Detrending parameter']]
       if nlin ge 1 then for j=1, nlin do coeffnames=[[coeffnames],['C_{'+strtrim(i,2)+','+strtrim(j,2)+'}','Detrending parameter']]
    endfor
    
    ;; each of the RV zero points
    gammanames = ['\gamma_{' + (*(rvptrs[0])).label + '}','m/s']
    for i=1, nrvtel-1 do gammanames = [[gammanames],['\gamma_{' + (*(rvptrs[i])).label + '}','m/s']]
    
    ;; names for limb darkening coefficients in each band (Claret, 2011)
    ldnames = [['u_{1' + bands[0] + '}','Linear Limb-darkening'],['u_{2' + bands[0] + '}','Quadratic Limb-darkening']]
    for j=1, nbands-1 do $
       ldnames = [[ldnames],['u_{1' + bands[j] + '}','Linear Limb-darkening'],['u_{2' + bands[j] + '}','Quadratic Limb-darkening']]
    
    ;; output labels and units
    latexnames = [['\dot{\gamma}', 'RV slope (m/s/day)'],$    ;; 0 Step parameters
                  ['T_C',          'Time of inferior conjunction (\bjdtdb)'  ],$ ;; 1
                  ['\log{P}',      'Log of period'                           ],$ ;; 2
                  ['\ecosw',       ''                                        ],$ ;; 3
                  ['\esinw',       ''                                        ],$ ;; 4
                  ['\log{K}',      'Log of the RV semi-amplitude'            ],$ ;; 5
                  ['\cos{i}',      'Cos of inclination'                      ],$ ;; 6
                  ['R_{P}/R_{*}',  'Radius of the planet in stellar radii'   ],$ ;; 7
                  ['a/R_*',        'Semi-major axis in stellar radii'        ],$ ;; 8
                  ['\log{g_*}',    'Surface gravity (cgs)'                   ],$ ;; 9
                  ['\teff',        'Effective temperature (K)'               ],$ ;; 10
                  ['\feh',         'Metallicity'                             ],$ ;; 11
                  ['\delta_S',     'Secondary eclipse depth'                 ],$ ;; 12
                  ['v\sin{I_*}',   'Rotational velocity (m/s)'               ],$ ;; 13
                  ['\lambda',      'Spin-orbit alignment (degrees)'          ],$ ;; 14
                  ['NR Vel. W.',   'Non-rotating line width (m/s)'           ],$ ;; 15                 ;;###
                  [gammanames                                                ],$ ;; 16:15+nrvtel       ;;###
                  [coeffnames                                                ],$ ;; 16+nrvtel:npars-1  ;;###
                  ['e',            'Eccentricity'                            ],$ ;; npars
                  ['\omega_*',     'Argument of periastron (degrees)'        ],$ ;; npars+1
                  ['T_{P}',        'Time of periastron (\bjdtdb)'            ],$ ;; npars+2
                  ['b_{S}',        'Impact parameter'                        ],$ ;; npars+3
                  ['T_{S,FWHM}',   'FWHM duration (days)'                    ],$ ;; npars+4
                  ['\tau_S',       'Ingress/egress duration (days)'          ],$ ;; npars+5
                  ['T_{S,14}',     'Total duration (days)'                   ],$ ;; npars+6
                  ['P_{S}',        'A priori non-grazing eclipse probability'],$ ;; npars+7
                  ['P_{S,G}',      'A priori eclipse probability'            ],$ ;; npars+8
                  ['T_{S}',        'Time of eclipse (\bjdtdb)'               ],$ ;; npars+9
                  ['M_{*}',        'Mass (\msun)'                            ],$ ;; npars+10
                  ['R_{*}',        'Radius (\rsun)'                          ],$ ;; npars+11
                  ['L_{*}',        'Luminosity (\lsun)'                      ],$ ;; npars+12
                  ['\rho_*',       'Density (cgs)'                           ],$ ;; npars+13
                  ['P',            'Period (days)'                           ],$ ;; npars+14
                  ['a',            'Semi-major axis (AU)'                    ],$ ;; npars+15
                  ['M_{P}',        'Mass (\mj)'                              ],$ ;; npars+16
                  ['R_{P}',        'Radius (\rj)'                            ],$ ;; npars+17
                  ['\rho_{P}',     'Density (cgs)'                           ],$ ;; npars+18
                  ['\log{g_{P}}',  'Surface gravity'                         ],$ ;; npars+19
                  ['T_{eq}',       'Equilibrium temperature (K)'             ],$ ;; npars+20
                  ['\Theta',       'Safronov number'                         ],$ ;; npars+21
                  ['\fave',        'Incident flux (\fluxcgs)'                ],$ ;; npars+22
                  ['K',            'RV semi-amplitude (m/s)'                 ],$ ;; npars+23
                  ['K_R',          'RM amplitude (m/s)'                      ],$ ;; npars+24
                  ['M_P\sin{i}',   'Minimum mass (\mj)'                      ],$ ;; npars+25
                  ['M_{P}/M_{*}',  'Mass ratio'                              ],$ ;; npars+26
                  ['i',            'Inclination (degrees)'                   ],$ ;; npars+27
                  ['b',            'Impact parameter'                        ],$ ;; npars+28
                  ['\delta',       'Transit depth'                           ],$ ;; npars+29
                  ['T_{FWHM}',     'FWHM duration (days)'                    ],$ ;; npars+30
                  ['\tau',         'Ingress/egress duration (days)'          ],$ ;; npars+31
                  ['T_{14}',       'Total duration (days)'                   ],$ ;; npars+32
                  ['P_{T}',        'A priori non-grazing transit probability'],$ ;; npars+33
                  ['P_{T,G}',      'A priori transit probability'             ],$ ;; npars+34
                  ['u',            'RM linear limb darkening'                ],$ ;; npars+35
                  [ldnames                                                   ],$ ;; npars+36:npars+2*nbands+35
                  ['f(m1,m2)',     'Mass function (\mj)'],$
                  [Tcnames]] ;; T_Cs
    
    
    ;; re-organize into more intuitive output table
    sidelabels = ['Stellar Parameters:','Planetary Parameters:','RV Parameters:',$
                  'Primary Transit Parameters:','Secondary Eclipse Parameters:']
    order = [-1,[10,11,12,13]+npars,9,10,11,13,14,15,$             ;; stellar pars  ;;###
             -1,[0,1,14,15,16,17,18,19,20,21,22]+npars,$           ;; planetary pars
             -1,1,[2,23,24,25,26,35]+npars,indgen(nrvtel)+16,0,3,4,npars+2*nbands+36,$ ;; RV Pars  ;;###
             -1,7,8,[27,28,29,30,31,32,33,34]+npars,$              ;; Primary pars
             indgen(npars-16-nrvtel)+16+nrvtel,npars+2*nbands+37+lindgen(ntransits),$                   ;; TTV, F0 and C0, T_Cs   ;;###
             indgen(2*nbands)+npars+36,$                           ;; Limb-dark pars
             -1,12,[9,3,4,5,6,7,8]+npars]                          ;; Secondary pars
    
    ;; all parameters to fit (exclude them later, if desired)
    tofit = [indgen(npars)]
    namendx = indgen(n_elements(latexnames[0,*]))
    
    ;; if no slope is desired, exclude it
    if keyword_set(noslope) then begin
       ;; exclude slope
       tofit = tofit[where(tofit ne 0)]
       namendx = namendx[where(namendx ne 0)]
       order = order[where(order ne 0)]
    endif else slopetxt = 'm.'
    
    ;; if we assume the orbit is circular, exclude ecosw, esinw
    if keyword_set(circular) then begin
       ;; exclude 3,4 (ecosw,esinw)
       tofit = tofit[where(tofit ne 3 and tofit ne 4)]
    
       ;; also skip the derived quantities e,omega,Tp,secondary eclipse pars
       namendx = namendx[where(namendx ne 3 and namendx ne 4 and $
                               namendx ne npars and $
                               namendx ne npars+1 and $
                               namendx ne npars+2 and $
                               namendx ne npars+3 and $
                               namendx ne npars+4 and $
                               namendx ne npars+5 and $
                               namendx ne npars+6 and $
                               namendx ne npars+7 and $
                               namendx ne npars+8)]
       order = order[where(order ne 3 and order ne 4 and $
                           order ne npars and $
                           order ne npars+1 and $
                           order ne npars+2 and $
                           order ne npars+3 and $
                           order ne npars+4 and $
                           order ne npars+5 and $
                           order ne npars+6 and $
                           order ne npars+7 and $
                           order ne npars+8)]
    endif
    
    ;; if no secondary depth is desired, exclude it
    if not keyword_set(secondary) then begin
       ;; exclude 12 (secondary depth)
       tofit = tofit[where(tofit ne 12)]
       namendx = namendx[where(namendx ne 12)]
       order = order[where(order ne 12)]
    endif
    
    ;; exclude the rossiter and dopp. tom. parameters vsini and lambda
    if not keyword_set(rossiter) and not useDoppTom then begin
       tofit = tofit[where(tofit ne 13 and tofit ne 14)]
       namendx = namendx[where(namendx ne 13 and namendx ne 14)]
       order = order[where(order ne 13 and order ne 14)]
       angular = []
    endif else begin
       angular = where(tofit eq 14)
    endelse
    
    if not useDoppTom then begin ;; include the line width
       tofit = tofit[where(tofit ne 15)]
       namendx = namendx[where(namendx ne 15)]
       order = order[where(order ne 15)]
    endif
    
    if not keyword_set(rossiter) then begin 
       namendx = namendx[where(namendx ne npars+24)] ;; and RM amplitude
       order = order[where(order ne npars+24)] ;; and RM amplitude
    endif
   
    ;; exclude TTVs and set sinusoidal detrending phases as angular.
    if not keyword_set(ttv) then begin
       for i=0, ntransits-1 do begin
          tofit = tofit[where(tofit ne (*(transitptrs[i])).ndx)]
          namendx = namendx[where(namendx ne (*(transitptrs[i])).ndx)]
          order = order[where(order ne (*(transitptrs[i])).ndx)]
       endfor
    endif
    if nSinusoids gt 0 then begin
       for i=0, ntransits-1 do begin
         for f=0, nSinusoids-1 do begin
           angular = [angular,where(tofit eq ((*(transitptrs[i])).ndx+2+(*(transitptrs[i])).ndetrend + 1 + f*2))]
         endfor
       endfor  
    endif  
;    printf, log, 'angular = ', angular
;    print, 'angular = ', angular    
;    printf, log, 'angular parameters = ', tofit[angular]
;    print, 'angular parameters = ', tofit[angular]
    printf, log, 'fitted parameters for combined amoeba best fit = '
    print, 'fitted parameters for combined amoeba best fit = '
    printf, log, tofit
    print, tofit
    
    ;tofit = tofit[where(tofit ne 35 and tofit ne 37)]
    ;namendx = namendx[where(namendx ne 35 and namendx ne 37)]
    ;order = order[where(order ne 35 and order ne 37)]
    
    ;; input labels to plotting program
    latexnames = latexnames[*,namendx]
    ;; printf, log, latexnames
    ;; map the static indices onto variable indices
    norder = n_elements(order)
    vorder = intarr(norder)
    for i=0, norder-1 do vorder[i] = where(namendx eq order[i])
    
    ;; trick TeXtoIDL
    mydevice = !d.name
    set_plot, 'PS'
    parnames = TeXtoIDL(transpose(latexnames[0,*]),font=0)
    set_plot, mydevice
    
    uslope = priors[1,0]
    if uslope eq !values.d_infinity then uslope = 1.0d0
    
    utc = priors[1,1];
    if utc eq !values.d_infinity then utc = 0.01d0
    
    ulogp = priors[1,2]
    if ulogp eq !values.d_infinity then ulogp = 0.0029d0
    
    usecosw = priors[1,3];
    if usecosw eq !values.d_infinity then usecosw = 0.05
    
    usesinw = priors[1,4];
    if usesinw eq !values.d_infinity then usesinw = 0.05
    
    ulogk = priors[1,5];
    if ulogk eq !values.d_infinity then ulogk = 0.014
    
    ucosi = priors[1,6];
    if ucosi eq !values.d_infinity then ucosi = 0.05d0
    
    up = priors[1,7];
    if up eq !values.d_infinity then up = 0.01d0
    
    ulogar = priors[1,8];
    if ulogar eq !values.d_infinity then ulogar = 0.7d0
    
    ulogg = priors[1,9]
    if ulogg eq !values.d_infinity then ulogg = 1.0d0
    
    uteff = priors[1,10]
    if uteff eq !values.d_infinity then uteff = 2000d0
    
    ufeh = priors[1,11]
    if ufeh eq !values.d_infinity then ufeh = 1d0
    
    depth2 = priors[0,12]
    udepth2 = priors[1,12]
    if udepth2 eq !values.d_infinity then udepth2 = 0.002d0
    vsini = priors[0,13]
    uvsini = priors[1,13]
    if uvsini eq !values.d_infinity then uvsini = vsini * 0.1d0
    if useDoppTom then uvsini = 0d0
    uvsini = 0d0
    
    lambda = priors[0,14]
    ulambda = priors[1,14]
    if ulambda eq !values.d_infinity then ulambda = !dpi
    umacturb = priors[1,15]
    if umacturb eq !values.d_infinity then umacturb = 3000d0   ;;###
    
    printf, log, 'Finding best fit for the combination of all RV, RM, DT, and light curve datasets'
    print, 'Finding best fit for the combination of all RV, RM, DT, and light curve datasets'
    
    ;; find the best fit of the two combined
    pars = allpars
    masterscale =[uslope,utc,ulogp,usecosw,usesinw,ulogk,ucosi,up,ulogar,ulogg,uteff,ufeh,udepth2,uvsini,ulambda,umacturb,rvzerosscale,widths];,ugamma,dblarr(npars-nrvtel-15)+1d-2]*3d0   ;;###
    
    scale = dblarr(npars)

    scale[tofit] = masterscale[tofit]
    ;printf, log, 'scale[tofit]', scale
    ;print, 'scale[tofit]', scale
    ;print, 'Rp/R* = ',pars[7]
    ;print, 'u Rp/R* = ',scale[7]
    
    nmax=1d5
    tranfit = indgen(n_elements(transitptrs))
    rvfit = indgen(n_elements(rvptrs))
    
    ;debug0=1
    nAmoebaIterations = 0L
    best = multifast_amoeba(1d-6,function_name=chi2func,p0=pars,scale=scale,$
                          ncalls=ncalls,nmax=1d8,log=log)

    if best[0] eq -1 then printf, log, 'ERROR: Could not find best combined fit'                      
    if best[0] eq -1 then message, 'ERROR: Could not find best combined fit'
    titleformat='(A14,x,A11,x,A9,x,A5,x,A7)'
    format='(f14.6,x,f11.6,x,f9.2,x,f5.2,x,f8.6)'
    printf, log, ' '
    print, ' '    
    printf, log, 'Combined bestfit transit parameters'
    print, 'Combined bestfit transit parameters'
    printf, log, '      T_0    ', 'period  ', 'a/R*', 'i  ','depth', format=titleformat
    print, '      T_0    ', 'period  ', 'a/R*', 'i  ','depth', format=titleformat

    printf, log,best[1], 10^best[2], 10^(best[8]),acos(best[6])*180d0/!dpi,best[7]^2, format=format  
    print,      best[1], 10^best[2], 10^(best[8]),acos(best[6])*180d0/!dpi,best[7]^2, format=format 
    
    ;;; generate plot of data, best fit, and residuals
    modelfile = prefix + 'model.amoebabestfit.' + circtxt + slopetxt + 'ps'
    chi2 = call_function(chi2func,best,psname=modelfile)
    texbestfile = prefix + 'amoebabestfit.' + circtxt + slopetxt + 'tex'
    label = "tab:" + prefix
    caption = prefix + '  Amoeba Best Fit +/- Prior'
    bestlatexnames = [['\dot{\gamma}', 'RV slope (m/s/day)'],$    ;; 0 Step parameters
                  ['T_C',          'Time of inferior conjunction (\bjdtdb)'  ],$ ;; 1
                  ['\log{P}',      'Log of period'                           ],$ ;; 2
                  ['\ecosw',       ''                                        ],$ ;; 3
                  ['\esinw',       ''                                        ],$ ;; 4
                  ['\log{K}',      'Log of the RV semi-amplitude'            ],$ ;; 5
                  ['\cos{i}',      'Cos of inclination'                      ],$ ;; 6
                  ['R_{P}/R_{*}',  'Radius of the planet in stellar radii'   ],$ ;; 7
                  ['a/R_*',        'Semi-major axis in stellar radii'        ],$ ;; 8
                  ['\log{g_*}',    'Surface gravity (cgs)'                   ],$ ;; 9
                  ['\teff',        'Effective temperature (K)'               ],$ ;; 10
                  ['\feh',         'Metallicity'                             ],$ ;; 11
                  ['\delta_S',     'Secondary eclipse depth'                 ],$ ;; 12
                  ['v\sin{I_*}',   'Rotational velocity (m/s)'               ],$ ;; 13
                  ['\lambda',      'Spin-orbit alignment (degrees)'          ],$ ;; 14
                  ['MacTurb',      'Macro turbulance'                        ],$ ;; 15    ;;###
                  [gammanames                                                ],$ ;; 16:15+nrvtel  ;;###
                  [coeffnames                                                ]] ;; 16+nrvtel:npars-1   ;;###
    exofast_latextab, [transpose(best),transpose(masterscale)], texbestfile, parnames=bestlatexnames[0,0:npars-1], $  //latexnames do not match bestfit pars, need to make names list from priors
                      units=bestlatexnames[1,0:npars-1], $ ;order=vorder, sidelabels=sidelabels, 
                      caption=caption, label=label
    
    printf, log, 'amoeba best fit chi2 = ', chi2
    print, 'amoeba best fit chi2 = ', chi2
    printf, log, ' '
    print, ' '
    printf, log, 'Performing MCMC for the combination of all RV, RM, and light curve datasets'
    print, 'Performing MCMC for the combination of all RV, RM, and light curve datasets'
    
    ;debug0=1
    ;; do the MCMC fit

    multifast_demc, best, chi2func, pars, chi2=chi2, tofit=tofit,$
                  nthin=nthin, maxsteps=maxsteps, angular=angular, burnndx=burnndx, log=log

    printf, log, 'Finished MCMC'
    print, 'Finished MCMC'
    
    ;; chains didn't converge, probably a test run
    burnndx = burnndx < (maxsteps-4)
    
    sz = size(pars)
    nallpars = sz[1]
    nsteps = sz[2]
    nchains = sz[3]
    pars = reform(pars,sz[1],sz[2]*sz[3])
    chi2 = reform(chi2,sz[2]*sz[3])
    minchi2 = min(chi2,bestndx)
    
    ;; the best fit parameters (lowest chi^2)
    bestamoeba = best
    best[tofit] = pars[*,bestndx]
    ;; generate the model fit from the best MCMC values, not AMOEBA
    modelfile = prefix + 'model.mcmcbestfit.' + circtxt + slopetxt + 'ps'
    bestchi2 = call_function(chi2func,best,psname=modelfile)
    
    ;; get the named parameters out of the parameter array (for each combination)
    if not keyword_set(noslope) then begin
       slope = pars[0,*]
       offset = 0
    endif else offset = -1
    tc = pars[1+offset,*]
    logp = pars[2+offset,*]
    period = 10^logp
    
    ;; find the tc with the lowest covariance with period
    ;np = 0
    ;mincovar = 1d0
    ;done = 0
    ;repeat begin
    ;   covar = abs(correlate(period,tc + period*np))
    ;   if covar lt mincovar then begin
    ;      mincovar = covar
    ;      bestnp = np
    ;      if np ge 0 then np++
    ;      if np lt 0 then np--
    ;   endif else begin
    ;      if np gt 0 then np = -1 $
    ;      else done = 1
    ;   endelse
    ;endrep until done
    ;tc = tc + period*bestnp
    ;pars[1+offset,*] = tc
    
    if not keyword_set(circular) then begin
       sqrtecosw = pars[3+offset,*]
       sqrtesinw = pars[4+offset,*]   
       e = sqrtecosw^2 + sqrtesinw^2
       omega = atan(sqrtesinw,sqrtecosw)*180.d0/!dpi
       ecosw = e*cos(omega*!dpi/180d0)
       esinw = e*sin(omega*!dpi/180d0)
       pars[3+offset,*] = ecosw
       pars[4+offset,*] = esinw
    endif else begin
       e = 0d0
       omega = !dpi/2d0
       offset -= 2
       esinw = 0d0
       ecosw = 0d0
    endelse
    logk = pars[5+offset,*]
    k = 10^logk
    cosi = pars[6+offset,*]
    i = acos(cosi)*180.d0/!dpi
    p = pars[7+offset,*]
    ar = 10^(pars[8+offset,*])
    pars[8+offset,*] = ar
    logg = pars[9+offset,*]
    teff = pars[10+offset,*]
    feh  = pars[11+offset,*]
    if keyword_set(secondary) then depth2 = pars[12+offset,*] $
    else offset -= 1
    
    if keyword_set(rossiter) or useDoppTom then begin
       vsini = pars[13+offset,*]
       ;; amplitude of the RM effect 
       ;; Eq 6., Gaudi & Winn 2007
       KR = vsini*p^2/(1d0-p^2)
    endif else offset -= 1
    
    if keyword_set(rossiter) or useDoppTom then begin
       pars[14+offset,*] *= 180d0/!dpi ;; convert lambda to degrees
    endif else offset -= 1
    
    if not useDoppTom then begin
      offset -= 1;
    endif
    
    ;; the central transit times for each transit
    tcs = dblarr(ntransits,sz[2]*sz[3])
    for j=0, ntransits-1 do begin
    
       if keyword_set(ttv) then begin
          timeoff = pars[(*(transitptrs[j])).ndx+offset,*]
       endif else timeoff = 0d0     
    
       tcs[j,*] = tc + (*(transitptrs[j])).epoch*period + timeoff
    
    endfor
    
    sini = sin(i*!dpi/180d0)
    G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
    a = (10^logg*(period*86400d0)^2/(4d0*!dpi^2*ar^2*100d0) + $
         K*(Period*86400d0)*sqrt(1d0-e^2)/(2d0*!dpi*sini))/6.9566d8         ;; R_sun
    rstar = a/ar                                                            ;; R_sun
    mp = 2d0*!dpi*(K*86400d0/6.9566d8)*a^2*sqrt(1d0-e^2)/(G*period*sini)    ;; M_sun
    mstar = 4d0*!dpi^2*a^3/(G*period^2) - mp                                ;; M_sun
    msini = mp*sini
    q = mp/mstar
    massfunction = msini^3/(mp+mstar)^2
    
    sigmab = 5.670373d-5/3.839d33*6.9566d10^2 ;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
    lstar = 4d0*!dpi*rstar^2*teff^4*sigmaB ;; L_sun
    
    bp = ar*cosi*(1d0-e^2)/(1d0+esinw)
    bs = ar*cosi*(1d0-e^2)/(1d0-esinw)
    rp = p*rstar
    delta = p^2
    loggp = alog10(G*mp/rp^2*9.31686171d0)      ;; cgs
    
    ;; densities of star and planet
    rhostar = mstar/(rstar^3)*1.41135837d0 ;; cgs
    rhop = mp/(rp^3)*1.41135837d0 ;; cgs
    
    corrmp = correlate(transpose(mp),transpose(rp),/double)
    
    ;; limb darkening coefficients for each band (Claret, 2011)
    for j=0, nbands-1 do begin                        
       ldcoeffs = quadld(logg, teff, feh, bands[j])
       if j eq 0 then allldcoeffs = ldcoeffs $
       else allldcoeffs = [allldcoeffs,ldcoeffs]
    endfor
    
    ;; linear limb darkening law for RM effect
    ;; TRES is more like B+V... does it matter?
    u1lin = linld(logg,teff,feh,'V')
    
    ;; Winn 2010
    t14 = period/!dpi*asin(sqrt((1d0+p)^2 - bp^2)/(sini*ar))*$
          sqrt(1d0-e^2)/(1d0+esinw) ;; eq 14, 16
    notransit = where(bp gt 1d0+p) 
    if notransit[0] ne -1 then t14[notransit] = 0d0
    t23 = period/!dpi*asin(sqrt((1d0-p)^2 - bp^2)/(sini*ar))*$
          sqrt(1d0-e^2)/(1d0+esinw) ;; eq 15, 16
    grazing = where(bp gt 1d0-p) 
    if grazing[0] ne -1 then t23[grazing] = 0d0
    tau = (t14-t23)/2d0
    Tfwhm = t14-tau
    primaryprobgraz = (rstar+rp)/a*(1d0 + esinw)/(1d0-e^2) ;; equation 9
    primaryprob     = (rstar-rp)/a*(1d0 + esinw)/(1d0-e^2) ;; equation 9
    
    ;; durations of secondary
    t14s = period/!dpi*asin(sqrt((1d0+p)^2 - bs^2)/(sini*ar))*$
          sqrt(1d0-e^2)/(1d0-esinw) ;; eq 14, 16
          
    ;;modification from Jason on 4/26/2013
    notransit = where(bs gt 1d0+p) 
    if notransit[0] ne -1 then t14s[notransit] = 0d0   
    ;;modification end   
          
    t23s = period/!dpi*asin(sqrt((1d0-p)^2 - bs^2)/(sini*ar))*$
          sqrt(1d0-e^2)/(1d0-esinw) ;; eq 15, 16
    grazing = where(bs gt 1d0-p) 
    if grazing[0] ne -1 then t23s[grazing] = 0d0
    taus = (t14s-t23s)/2d0
    Tfwhms = t14s-taus
    secondaryprobgraz = (rstar+rp)/a*(1d0 - esinw)/(1d0-e^2) ;; equation 10
    secondaryprob     = (rstar-rp)/a*(1d0 - esinw)/(1d0-e^2) ;; equation 10
    
    ;; phases for primary and secondary transits
    phase = exofast_getphase(e,omega*!dpi/180d0, /primary)
    phase2 = exofast_getphase(e,omega*!dpi/180d0, /secondary)
    
    ;; periastron passage and secondary eclipse times
    tp = tc - period*phase
    ts = tc + period*(phase2-phase)
    
    ;; it's possible tp/ts could be split down the middle 
    ;; then the median would be meaningless -- correct that
    ;; tp
    hist = histogram(tp,nbins=100,locations=x)
    max = max(hist,modendx)
    mode = x[modendx]
    toohigh = where(tp gt (mode + period/2d0))
    if toohigh[0] ne -1 then tp[toohigh] -= period
    toolow = where(tp lt (mode - period/2d0))
    if toolow[0] ne -1 then tp[toolow] += period
    
    ;; ts
    hist = histogram(ts,nbins=100,locations=x)
    max = max(hist,modendx)
    mode = x[modendx]
    toohigh = where(ts gt (mode + period/2d0))
    if toohigh[0] ne -1 then ts[toohigh] -= period
    toolow = where(ts lt (mode - period/2d0))
    if toolow[0] ne -1 then ts[toolow] += period
    
    ;; the planet's equilibrium temperature 
    ;; no albedo, perfect redistribution
    ;; Eq 2, Hansen & Barman, 2007
    teq = teff*sqrt(1d0/(2d0*ar))
    
    ;; Safronov Number eq 2, Hansen & Barman, 2007
    safronov = ar*q/p
    
    ;; <F>, the time-averaged flux incident on the planet
    ;; (ar*(1d0+e^2/2d0)) = time averaged distance to the planet
    sigmasb = 5.6704d-5 ;; stefan boltzmann constant 
    incidentflux = sigmasb*teff^4/(ar*(1d0+e^2/2d0))^2/1d9 ;; erg/s/cm^2
    
    mjup = 0.000954638698d0 ;; m_sun
    rjup = 0.102792236d0    ;; r_sun
    AU = 215.094177d0 ;; r_sun
    
    ;; convert times as measured in the SSB
    tcbjdtdb = target2bjd(tc,inclination=i*!dpi/180d0, a=a/AU, tp=tp, $
                          period=period, e=e,omega=omega*!dpi/180d0,q=q)
    tpbjdtdb = target2bjd(tp,inclination=i*!dpi/180d0, a=a/AU, tp=tp, $
                          period=period, e=e,omega=omega*!dpi/180d0,q=q)
    tsbjdtdb = target2bjd(ts,inclination=i*!dpi/180d0, a=a/AU, tp=tp, $
                          period=period, e=e,omega=omega*!dpi/180d0,q=q)
    
    ;; parameters to plot PDFs, covariances, and quote 68% confidence interval
    ;; must match output labels
    if not keyword_set(circular) then begin
       angular = n_elements(pars[*,0])+1
       pars = [pars,e,omega,tpbjdtdb,bs,tfwhms,taus,$
               t14s,secondaryprob,secondaryprobgraz]
    endif else angular = []
    
    pars = [pars,tsbjdtdb,mstar,rstar,lstar,rhostar,period,a/au,mp/mjup,rp/rjup,$
            rhop,loggp,teq,safronov,incidentflux,k]
    
    if keyword_set(rossiter) or useDoppTom then begin
       angular = [angular,where(tofit eq 14)]
    endif
    
    if keyword_set(rossiter) then begin
       pars = [pars,KR]
    endif
    
    pars = [pars,msini/mjup,q,i,bp,delta,tfwhm,tau,t14,primaryprob,$
            primaryprobgraz,u1lin,allldcoeffs,massfunction/mjup,tcs]

    if nSinusoids gt 0 then begin
       for i=0, ntransits-1 do begin
         for f=0, nSinusoids-1 do begin
           angular = [angular,where(tofit eq ((*(transitptrs[i])).ndx+2+(*(transitptrs[i])).ndetrend + 1 + f*2))]
           pars[where(tofit eq ((*(transitptrs[i])).ndx+2+(*(transitptrs[i])).ndetrend + 1 + f*2)),*] *= 180d0/!dpi 
         endfor
       endfor  
    endif 
    
    ;; the best fit parameters
    bestpars = pars[*,bestndx]
    
    ;; save the chain and chi2 the MCMC chain.
    idlfile = prefix + 'mcmc.' + circtxt + slopetxt + 'idl'
    nallpars = n_elements(pars[*,0])
    pars = reform(pars,nallpars,nsteps,nchains)
    chi2 = reform(chi2,nsteps,nchains)
    if (n_elements(angular) eq 0) then begin
        save, pars, chi2, latexnames, burnndx, nsteps, nchains, nallpars, bestpars, parnames, tofit, $
             vorder, sidelabels, bestamoeba, best, period, $
             rvdata, rvdebug, $
             rvptrs, transitptrs, priors0, debug0, rvfit, tranfit, $
             useDoppTom, dopptomptrs, $
             nSinusoids, frequencies, $
             rmerrscale, $
             filename=idlfile      
    endif else begin
        save, pars, chi2, latexnames, burnndx, nsteps, nchains, nallpars, bestpars, parnames, tofit, $
             angular, $
             vorder, sidelabels, bestamoeba, best, period, $
             rvdata, rvdebug, $
             rvptrs, transitptrs, priors0, debug0, rvfit, tranfit, $
             useDoppTom, dopptomptrs, $
             nSinusoids, frequencies, $
             rmerrscale, $
             filename=idlfile
    endelse
endif else begin
    printf, log, 'Restoring data from IDL save file'
    print, 'Restoring data from IDL save file'
    restore, prefix + 'mcmc.' + circtxt + slopetxt + 'idl'
    
    ;; generate the model fit from the best MCMC values
    printf, log, 'Creating model plots'
    print, 'Creating model plots'
    !P.MULTI=0
    modelfile = prefix + 'model.mcmcbestfit.' + circtxt + slopetxt + 'ps'
    bestchi2 = call_function(chi2func,best,psname=modelfile)    
endelse
;print, 'angular = ', angular
if not keyword_set(replotmodels) then begin
    ;; now trim the burn-in before making plots
    pars = reform(pars[*,burnndx:nsteps-1,*],nallpars,(nsteps-burnndx)*nchains)
    chi2 = reform(chi2[burnndx:nsteps-1,*],(nsteps-burnndx)*nchains)
    
    label = "tab:" + prefix
    caption = " values and 68\% confidence interval for " + prefix
    parfile = prefix + 'pdf.' + circtxt + slopetxt + 'ps'
    covarfile = prefix + 'covar.' + circtxt + slopetxt + 'ps'
    texmedfile = prefix + 'median.' + circtxt + slopetxt + 'tex'
    printf, log, 'Creating pdf and covariance plots'
    print, 'Creating pdf and covariance plots'
    multifast_plotdist, pars, finalpars, bestpars=bestpars, $
                      parnames=parnames, angular=angular, log=log,$
                      pdfname=parfile, covarname=covarfile,/degrees
    printf, log, 'Creating parameter table'
    print, 'Creating parameter table'
    exofast_latextab, finalpars, texmedfile, parnames=latexnames[0,*], $
                      units=latexnames[1,*], order=vorder, sidelabels=sidelabels, $
                      caption='Median' + caption, label=label+'median'
                      
    if useDoppTom then begin
      printf, log, 'Doppler Tomography Amoeba best chi2 = ', multifast_dopptom(bestamoeba)   
      print, 'Doppler Tomography Amoeba best chi2 = ', multifast_dopptom(bestamoeba)                  
      printf, log, 'Doppler Tomography MCMC best chi2 = ', multifast_dopptom(best)
      print, 'Doppler Tomography MCMC best chi2 = ', multifast_dopptom(best)
    endif
    printf, log, 'Total Amoeba best fit chi2 =', call_function(chi2func,bestamoeba)
    print, 'Total Amoeba best fit chi2 =', call_function(chi2func,bestamoeba)
    printf, log, 'Total MCMC best fit chi2 =', call_function(chi2func,best)
    print, 'Total MCMC best fit chi2 =', call_function(chi2func,best)
    best[tofit] = finalpars[0,0:(n_elements(tofit)-1)]
    if not keyword_set(circular) then begin
       ecosw = best[3]
       esinw = best[4]
       e = sqrt(ecosw^2+esinw^2)
       omega = atan(esinw,ecosw)
       best[3] = sqrt(e)*cos(omega)
       best[4] = sqrt(e)*sin(omega)
    endif else begin
       best[3]=0d0
       best[4]=0d0
    endelse
    best[8] = alog10(best[8])
    best[14] *= (!dpi/180.d0)
    if nSinusoids gt 0 then begin
       for i=0, ntransits-1 do begin
         for f=0, nSinusoids-1 do begin
           best[((*(transitptrs[i])).ndx+2+(*(transitptrs[i])).ndetrend + 1 + f*2)] *= !dpi/180d0 
         endfor
       endfor  
    endif     
    medianchi2modelfile = prefix + 'model.mcmcmedianfit.' + circtxt + slopetxt + 'ps'
    medianchi2 = call_function(chi2func,best,psname=medianchi2modelfile)
    printf, log, 'Total MCMC median parameter chi2 =', medianchi2
    print, 'Total MCMC median parameter chi2 =', medianchi2 
    ;; TTVs
    if keyword_set(ttv) then begin
       printf, log, 'Creating TTV results and plot'
       print, 'Creating TTV results and plot'
       openw, lun, prefix+'ttv.' + circtxt + slopetxt + 'txt', /get_lun
       npars = n_elements(pars[*,0])
       printf, log, ' '
       print, ' '
       printf, log, 'T_C (BJD_TDB)   ERR (d)    Observatory'
       print, 'T_C (BJD_TDB)   ERR (d)    Observatory'
       for i=0, ntransits-1 do begin
          if (((*(transitptrs[i])).secondary) eq 0) then begin
             printf, lun, finalpars[0,npars-ntransits+i],$
                     mean(finalpars[1:2,npars-ntransits+i]), $
                     (strsplit((*(transitptrs[i])).label,/extract))(0),format='(f14.6,x,f8.6,x,a10)'
             printf, log, finalpars[0,npars-ntransits+i],$
                     mean(finalpars[1:2,npars-ntransits+i]), $
                     (strsplit((*(transitptrs[i])).label,/extract))(0),format='(f14.6,x,f8.6,x,a14)'                 
             print, finalpars[0,npars-ntransits+i],$
                     mean(finalpars[1:2,npars-ntransits+i]), $
                     (strsplit((*(transitptrs[i])).label,/extract))(0),format='(f14.6,x,f8.6,x,a14)'
          endif
       endfor
       printf, log, ' '
       print, ' '
       close, lun
       free_lun, lun
       multifast_omc, prefix + 'ttv.' + circtxt + slopetxt + 'txt',period,epsname=prefix+'ttv.ps',log=log
    endif    
endif

;; display all the plots, if desired
if keyword_set(display) then begin
   spawn, 'gv ' + parfile + ' &'
   spawn, 'gv ' + covarfile + ' &'
   spawn, 'gv ' + modelfile + ' &'
   spawn, 'gv ' + prefix+'ttv.eps &'
   spawn, 'latpdf ' + prefix + 'median.' + circtxt + (strsplit(slopetxt,'.',/extract))(0)
endif

printf, log, 'Multifast finished'
print, 'Multifast finished'

;if log ne -1 then begin
   close, log
   free_lun, log
;endif

end
