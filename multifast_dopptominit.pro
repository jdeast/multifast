pro multifast_dopptominit

COMMON dopptom_block, useDoppTom, DTccf2d, DTbjd, DTvel, DTstepsize, DTrms, DTchisqr0, Rspec, DTmodel, DTresid, lambdaRange

DTstepsize = DTvel*0.
nvels = n_elements(DTvel)
for i=0,(nvels-2) do begin
	DTstepsize[i] = DTvel[i+1]-DTvel[i]
endfor
DTmeanstepsize = mean(DTstepsize)
DTstepsize[-1] = DTmeanstepsize

rms0 = stddev(DTccf2d)
chisqr0 = total(sqrt((DTccf2d-median(DTccf2d))^2./rms0^2.)) / n_elements(DTccf2d)
DTrms = rms0 * chisqr0
DTchisqr0 = total(sqrt((DTccf2d-median(DTccf2d))^2./DTrms^2.))

end