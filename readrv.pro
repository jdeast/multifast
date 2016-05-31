function readrv, filename

label = (strsplit(filename,'.',/extract))(1)

;; Read the RV data file into the RV structure
readcol, filename, bjd, rv, err,bi,bierr, format='d,d,d,d,d',/silent,comment='#'

return, create_struct('bjd',bjd,'rv',rv,'err',err, 'bi',bi,'bierr',bierr,'label',label,'residuals',dblarr(n_elements(bjd)), 'rm',dblarr(n_elements(bjd)))

end
