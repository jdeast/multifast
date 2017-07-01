function readtransit3, filename, nmaxdetrends

;; Read the transit data file into a structure
;; (with an arbitary number of detrending variables)
band = (strsplit(filename,'.',/extract))(1)

if band eq 'u' then begin
   band = 'Sloanu'
   bandname = "u"
endif else if band eq 'g' then begin
   band = 'Sloang'
   bandname = "g"
endif else if band eq 'r' then begin
   band = 'Sloanr'
   bandname = "r"
endif else if band eq 'i' then begin
   band = 'Sloani'
   bandname = "i"
endif else if band eq 'z' then begin
   band = 'Sloanz'
   bandname = "z"
endif else bandname = band

line = ""
openr, lun, filename, /get_lun
readf, lun, line
entries = double(strsplit(line,/extract))
ncol = n_elements(entries)
nrow = file_lines(filename)

;; rewind file to beginning
point_lun, lun, 0
array = dblarr(ncol,nrow)
readf, lun, array
free_lun, lun

bjd = transpose(array[0,*])
flux = transpose(array[1,*])
err = transpose(array[2,*])

if (ncol gt 3) && (nmaxdetrends gt 0) then begin
   ndetrend = ncol-3 < nmaxdetrends
   d = array[3:3+ndetrend-1,*]
endif else begin
   d = 0d0
   ndetrend=0
endelse

secondary = 0;
if (strmid(filename,0,1) eq 's') then secondary = 1;

night = strmid(filename,1,4)+'-'+strmid(filename,5,2)+'-'+strmid(filename,7,2)
label = (strsplit(filename,'.',/extract))(2) + ' UT ' + night + ' ('+ bandname + ')'
transit=create_struct('bjd',bjd,'flux',flux,'err',err,'band',band,'ndx',0,$
                      'epoch',0.0,'detrend',d,'label',label,'ndetrend',ndetrend,'secondary',secondary)
return, transit

end
