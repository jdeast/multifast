function multifast_dopptom, pars

;; This code was written by Thomas Beatty at PSU and is an adaptation of the 
;; Collier-Cameron work on HD 189733b (https://arxiv.org/abs/0911.5361), 
;; altered to account for stars that are rotating *really* fast. The best description
;; that parallels this implementation is in:
;; http://adsabs.harvard.edu/abs/2013A%26A...550A..53B,
;; particularly section 4.

;; The DT model herein treats the planet shadow as a convolution of the average non-rotating line width
;; in the stellar spectrum, the instrumental resolution, and the rotation kernel. The integrated area
;; of the shadow is made to be the fraction of the limb-darkened star that is currently obscured.

;; The "macturb" parameter is the non-rotating line width. For super-fast rotating stars (vsini > 100km/s ??)
;; setting a near-zero prior on the macturb parameter may be okay, because the stellar rotation will completely
;; dominate the inherent line width, but for stars with slower rotation, this is inappropriate and the macturb
;; parameter should be fitted to give an appropriate model fit.

COMMON dopptom_block, useDoppTom, dopptomptrs

tc        = pars[1]       ;; transit center time
Period    = 10^(pars[2])  ;; Period of orbit
sqrtecosw = pars[3]       ;; ecosw
sqrtesinw = pars[4]       ;; esinw
cosi      = pars[6]       ;; cosi
p         = pars[7]       ;; rp/rstar
ar        = 10^(pars[8])
logg      = pars[9]       ;; density of the Primary (log(rho_star)/(rho_sun))
teff      = pars[10]      ;; effective temperature of the Primary
feh       = pars[11]      ;; metallicity of the Primary
vsini     = pars[13]/1000.      ;; vsini of the star
lambda    = pars[14]      ;; sky-projected spin-orbit angle
macturb   = pars[15]/1000.

if macturb lt 0d0 then return, !values.d_infinity ; invalid line width
if vsini lt 0d0 then return, !values.d_infinity ; invalid vsini
chisqr = 0.0d0;
convolve_limit = 5. ;; The gaussian broadening has to be less than this factor larger than vsini for it
                    ;; to try and include the rotation kernel. E.g., gaussian broadening ('GaussTerm') = 10km/s,
                    ;; vsini needs to be at least 2km/s. Otherwise just use a gaussian for the shadow.

inc = acos(cosi)
b = ar * cosi
e = sqrtecosw^2 + sqrtesinw^2
if e eq 0d0 then omega = !dpi/2d0 else omega = atan(sqrtesinw, sqrtecosw)

for d = 0, n_elements(dopptomptrs) - 1 do begin
    DTccf2d = (*(dopptomptrs[d])).ccf2d
    DTbjd = (*(dopptomptrs[d])).bjd
    DTvel = (*(dopptomptrs[d])).vel
    DTstepsize = (*(dopptomptrs[d])).stepsize 
    DTrms = (*(dopptomptrs[d])).rms
    DTchisqr0 = (*(dopptomptrs[d])).chisqr0
    DTRspec = (*(dopptomptrs[d])).Rspec
    DTblendFactor = (*(dopptomptrs[d])).blendFactor
    velsini = DTvel/vsini
    stepsizesini = DTstepsize / vsini
    meanstepsize = mean(stepsizesini)
    
    ; relaventVels = where((velsini gt -1.5) and (velsini lt 1.5)) ; we only care about these, to speed things up
    
    GaussTerm = sqrt(macturb^2. + (300000./DTRspec)^2.) ; so we can add in both the broadening from the instrument's spectral resolution, and the macturb parameter
    
    GaussRel = GaussTerm/vsini
    
    relvelslimit = 1.2 ; 1.0 + (GaussRel*5.)
    relaventVels = where((velsini gt -relvelslimit) and (velsini lt relvelslimit)) ; we only care about these, to speed things up
    
    IndepVels = (300000./DTRspec) / (meanstepsize*vsini) ; accounts for the fact that we're supersampling within the actual spectral res. of the spectrograph, to adjust chi^2 later
;    if (vsini ne 110.d0) Then Stop
    nper = round((mean(DTbjd)-tc) / period)
    localtc = tc + Period*nper
    phase = (DTbjd-localtc) / Period
    
    ;; init for planet shadow model
    velwidth = p
    c1 = 2 / (3.14159*velwidth)
    
    xp = ar * sin(phase*2.*!dpi)
    zp = ar * cos(phase*2.*!dpi) * cosi
    up = xp*cos(lambda) + zp*sin(lambda)
    
    ;; Calculate Collier-Cameron's "beta" (modelflux) as the actual transit LC in V
    ;; This presumes we're using an ~optical spectrograph for the observations
    
    modelflux = phase*0.0d0
    
    ldcoeffs = quadld(logg, teff, feh, 'V')
    u1 = ldcoeffs[0]
    u2 = ldcoeffs[1]
    if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity
    
    ;; original - only works correctly for circular orbits; z = exofast_getb(DTbjd, i=inc, a=ar, tperiastron=localtc, period=Period, e=e,omega=omega,z=depth)
    phasep = exofast_getphase(e,omega,/primary)
    localtp = localtc - period*phasep
    z = exofast_getb(DTbjd, i=inc, a=ar, tperiastron=localtp, period=Period, e=e,omega=omega,z=depth)
    
    primary = where(depth gt 0,complement=secondary)     
    if primary[0] ne -1 then begin
       exofast_occultquad, z[primary], u1, u2, p, mu1
       modelflux[primary] = mu1
    endif
    
    modelflux = 1-modelflux ; to just get the amount of light blocked, rather than the lightcurve
    
    DTmodel = DTccf2d * 0. + median(DTccf2d)
    nphase = n_elements(phase)
    for i = 0, nphase-1 do begin
    	if modelflux[i] gt 0.0d0 then begin ; only do this if the planet is on the star
    		if up[i] lt 1.0d0 then begin ; also only if the DT thinks it's on the star
    			if GaussRel lt convolve_limit then begin ; rotation is significant
    				G = velsini[relaventVels]*0.0d0
    				c2 = ((velsini[relaventVels]-up[i])/velwidth)^2.
    				validlocs = where(c2 lt 1.0d0)
    				;; Need to check if we have any valid locations
    				if validlocs[0] ne -1 then G[validlocs] = c1*sqrt(1-c2[validlocs])
    				normG = 1 / total(G*stepsizesini)
    				if normG ne !values.d_infinity then begin ; so we don't waste time on zeros
    					rotprofile = G * normG
    					unnormalized = multifast_gaus_convol(velsini[relaventVels], rotprofile, GaussRel)
    					normalization = 1 / total(unnormalized*stepsizesini)
    					DTmodel[relaventVels,i] += (modelflux[i] * normalization * unnormalized)
    				endif
    			endif else begin ; gaussian line width dominates
    				expterm = -1.*(velsini[relaventVels]-up[i])^2./(2*GaussRel^2.)
    				unnormalized = exp(expterm)
    				normalization = 1 / total(unnormalized*stepsizesini)
    				DTmodel[relaventVels,i] += (modelflux[i] * normalization * unnormalized)
    			endelse
    		endif
    	endif
    endfor
    
    DTresid = DTccf2d - DTmodel
    ;chisqr0 = DTchisqr0 / IndepVels
    chisqr += total(DTresid^2./DTrms^2.)/IndepVels; - chisqr0
    (*(dopptomptrs[d])).model = DTmodel
    (*(dopptomptrs[d])).resid = DTresid
endfor

return, chisqr

end