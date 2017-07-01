function multifast_readdopptom, filenamebase, Rspec, DTblendFactor, lambdaRange, vsini 

;; First index is velocity, second index is BJD
ccf2d = -readfits(filenamebase+'dopptom.fits', /silent)
ccf2d = ccf2d/(1d0-DTblendFactor)
model = ccf2d * 0.0 + median(ccf2d)
resid = ccf2d - model
rdfloat, filenamebase+'bjds.dat', bjd,/double, /silent
rdfloat, filenamebase+'vels.dat', vel,/double, /silent

telescope = (strsplit(filenamebase,'.',/extract))(2)
night = strmid(filenamebase,1,4)+'-'+strmid(filenamebase,5,2)+'-'+strmid(filenamebase,7,2)
label = 'UT ' + night + ' ' + telescope
              
stepsize = vel*0.
nvels = n_elements(vel)
for i=0,(nvels-2) do begin
	stepsize[i] = vel[i+1]-vel[i]
endfor
baselinevels = where(abs(vel) gt vsini/1000d0)
nbaselinevels = n_elements(baselinevels)
if (nbaselinevels eq -1) then begin
    print, 'There are  no baseline Doppler tomography velocities outside the +/- vsini region.'
    print, 'Check the nominal setting of vsini in the priors and compare to the velocities listed in *.vels.dat'
    stop
endif

meanstepsize = mean(stepsize)
;print, 'mean stepsize = ', meanstepsize
stepsize[-1] = meanstepsize
;print, ccf2d[baselinevels,*]
baselineccf2d = ccf2d[baselinevels,*]
rms0 = stddev(baselineccf2d)
print, 'Doppler Tomography filename base = ', filenamebase
print, 'Number of baseline Doppler velocities = ', nbaselinevels
print, 'Doppler Tomography Baseline RMS = ', rms0
print, 'Doppler Tomography All Velocities RMS = ', stddev(ccf2d)
chisqr0 = total((baselineccf2d-median(baselineccf2d))^2./rms0^2.) / n_elements(baselineccf2d)
print, 'Doppler Tomography Baseline Chi2/dof = ', chisqr0
rms = rms0 * chisqr0
chisqr0 = total((ccf2d-median(ccf2d))^2./rms^2.)
print, 'Doppler Tomography All  Velocities Chi2/dof = ', chisqr0/n_elements(ccf2d)

dopptom=create_struct('ccf2d',ccf2d,'bjd',bjd,'vel',vel,'stepsize',stepsize,'rms',rms, $
                      'chisqr0',chisqr0,'Rspec',Rspec,'blendFactor',DTblendFactor,'model',model,'resid',resid, $
                      'lambdaRange',lambdaRange,'label',label,'tel',telescope,'night',night)
return, dopptom
end