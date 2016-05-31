pro bintime, jd, flux, err, binsize, binjd, binflux, binerr

good = where(finite(jd))
jdfit = jd[good]
fluxfit = flux[good]
errfit = err[good]

jdmax = max(jdfit)
jdmin = min(jdfit)

;convert minutes to days
binday = binsize/(24.0d0*60.0d0) 

npoints = n_elements(jdfit)
nbins = (jdmax - jdmin)/binday

binjd = dblarr(nbins) + !values.d_nan
binflux = dblarr(nbins) + !values.d_nan
binerr = dblarr(nbins) + !values.d_nan

for i=0L,nbins-1 do begin
    tobin = where(jdfit ge (jdmin + i*binday) and jdfit lt (jdmin + (i+1)*binday))
    if tobin[0] ne -1 then begin
       ;; inverse variance weighted means
       binflux[i] = total(fluxfit[tobin]/(errfit[tobin]^2))/total(1d0/errfit[tobin]^2)
       binjd[i] =   total(jdfit[tobin]/(errfit[tobin]^2))/total(1d0/errfit[tobin]^2)
       binerr[i] = 1d0/(sqrt(total(1d0/errfit[tobin]^2)))
    endif
endfor

good = where(finite(binjd))
binjd = binjd[good]
binflux = binflux[good]
binerr = binerr[good]

end
