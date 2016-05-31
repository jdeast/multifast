pro findMassRadius
teff=6105d0
feh=-0.281d0
logg=4.058d0

massradius_torres, logg, teff, feh, mstar, rstar
;ulogm = 0.027d0
;ulogr = 0.014d0
logm=alog10(mstar)
errmp=10d0^(logm+0.027d0)
errmm=10d0^(logm-0.027d0)

logr=alog10(rstar)
errrp=10d0^(logr+0.014d0)
errrm=10d0^(logr-0.014d0)

print, "mstar = ", mstar,"       +",errmp-mstar, "          -",mstar-errmm
print, "rstar = ", rstar,"       +",errrp-rstar, "          -",rstar-errrm

end