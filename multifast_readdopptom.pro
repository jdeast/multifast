function multifast_readdopptom, filenamebase, Rspec, lambdaRange

;; First index is velocity, second index is BJD
ccf2d = readfits(filenamebase+'dopptom.fits', /silent)
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
meanstepsize = mean(stepsize)
stepsize[-1] = meanstepsize

rms0 = stddev(ccf2d)
chisqr0 = total(sqrt((ccf2d-median(ccf2d))^2./rms0^2.)) / n_elements(ccf2d)
rms = rms0 * chisqr0
chisqr0 = total(sqrt((ccf2d-median(ccf2d))^2./rms^2.))

dopptom=create_struct('ccf2d',ccf2d,'bjd',bjd,'vel',vel,'stepsize',stepsize,'rms',rms, $
                      'chisqr0',chisqr0,'Rspec',Rspec,'model',model,'resid',resid, $
                      'lambdaRange',lambdaRange,'label',label,'tel',telescope,'night',night)
return, dopptom
end