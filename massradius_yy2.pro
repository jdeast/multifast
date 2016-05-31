function massradius_yy2, mstar, feh, teff, rstar, uteff,afe=afe,age=age,debug=debug,yyteff=yyteff, quad=quad, lsquad=lsquad,yytrack=yytrack

;; this is the approximate error in the YY isochrone modeling
uteff = 50d0 
if n_elements(afe) eq 0 then afe = 0d0

;; make these global for speed (restore is expensive)
COMMON YY_BLOCK, tracks

path = getenv('EXOFAST_PATH') + 'yy/'

if file_test(path + 'tracks.idl') then begin
   if n_elements(tracks) eq 0 then restore, path + 'tracks.idl' 
   sz = size(tracks)
   na = sz[1]
   nz = sz[2]
   nx = sz[3]
endif else begin
   ;; generate the tracks array and save it
   aname = ['a0o2','a2o2','a4o2']
   zname = ['x76997z00001', 'x7697z0001', 'x7688z0004',$
            'x767z001', 'x758z004', 'x749z007', 'x74z01',$
            'x71z02', 'x65z04', 'x59z06', 'x53z08']
   xmname = ['m04','m05','m06','m07','m08','m09','m10',$
             'm11','m12','m13','m14','m15','m16','m17','m18','m19','m20',$
             'm21','m22','m23','m24','m25','m26','m27','m28','m29','m30',$
             'm32','m34','m36','m38','m40','m42','m45','m50']
   na = n_elements(aname)
   nz = n_elements(zname)
   nx = n_elements(xmname)

   ;; nafe x nz x nmass x nage x 3 (teff, rstar, age) array 
   track1 = dblarr(na,nz,nx,24,3) + !values.d_nan
   track2 = dblarr(na,nz,nx,150,3) + !values.d_nan
   tracks = dblarr(na,nz,nx,174,3) + !values.d_nan

   ;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
   sigmab = 5.670373d-5/3.839d33*6.9566d10^2 

   for i=0,na-1 do begin
      for j=0,nz-1 do begin
         for k=0,nx-1 do begin 
            filename = path + aname[i] + path_sep() + zname[j] + path_sep() + xmname[k] + zname[j]
            
            ;; read all the first tracks
            if file_test(filename + '.track1') then begin
               readcol, path + aname[i] + path_sep() + zname[j] + path_sep() + xmname[k] + zname[j] + '.track1',$
                        n, age1,logteff1,loglstar1,format='i,d,d,d',skip=1,/silent
               track1[i,j,k,*,0] = 10^logteff1
               track1[i,j,k,*,1] = sqrt(10d0^loglstar1/(4d0*!dpi*track1[i,j,k,*,0]^4d0*sigmaB))
               track1[i,j,k,*,2] = age1
;               plot, logteff, loglstar
            endif
            
            ;; read all the second tracks
            if file_test(filename + '.track2') then begin
               readcol, path + aname[i] + path_sep() + zname[j] + path_sep() + xmname[k] + zname[j] + '.track2',$
                        n, age2,logteff2,loglstar2,format='i,d,d,d',skip=1,/silent
               track2[i,j,k,*,0] = 10^logteff2
               track2[i,j,k,*,1] = sqrt(10d0^loglstar2/(4d0*!dpi*track2[i,j,k,*,0]^4d0*sigmaB))
               track2[i,j,k,*,2] = age2
            endif

            if file_test(filename + '.track1') and file_test(filename + '.track2') then begin
               logteff = [logteff1,logteff2]
               loglstar = [loglstar1,loglstar2]
               age = [age1,age2]
               sorted = sort(age)
               tracks[i,j,k,*,0] = 10^logteff[sorted]
               tracks[i,j,k,*,1] = sqrt(10d0^loglstar[sorted]/(4d0*!dpi*tracks[i,j,k,*,0]^4d0*sigmaB))
               tracks[i,j,k,*,2] = age[sorted]
;               plot, tracks[i,j,k,*,0], tracks[i,j,k,*,1]
;               stop
            endif


         endfor
      endfor
   endfor
   save, filename=path + 'tracks.idl',tracks
endelse

;; convert [Fe/H] to Z (Table 2, Yi et al, 2001)
z = [[0.00001d0,-3.29d0],$ 
     [0.00010d0,-2.29d0],$  
     [0.00040d0,-1.69d0],$  
     [0.00100d0,-1.29d0],$  
     [0.00400d0,-0.68d0],$ 
     [0.00700d0,-0.43d0],$  
     [0.01000d0,-0.27d0],$  
     [0.02000d0, 0.05d0],$  
     [0.04000d0, 0.39d0],$  
     [0.06000d0, 0.60d0],$  
     [0.08000d0, 0.78d0]]

yyz = interpol(z[0,*],z[1,*],feh,quad=quad, lsquad=lsquad)

avalue = [0d0,0.3d0,0.6d0]
zvalue = [0.00001d0,0.0001d0,0.0004d0,0.001d0,0.004d0,$
          0.007d0,0.01d0,0.02d0,0.04d0,0.06d0,0.08d0]
xmass = [0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,$
         1.1d0,1.2d0,1.3d0,1.4d0,1.5d0,1.6d0,1.7d0,1.8d0,1.9d0,2.0d0,$
         2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0,$
         3.2d0,3.4d0,3.6d0,3.8d0,4.0d0,4.2d0,4.5d0,5.0d0]

z_one = interpol(indgen(nz),zvalue,yyz,quad=quad, lsquad=lsquad)
a_one = interpol(indgen(na),avalue,afe,quad=quad, lsquad=lsquad)
m_one = interpol(indgen(nx),xmass,mstar,quad=quad, lsquad=lsquad)

;; if the track never crosses through this rstar, return infinity
yyteff = !values.d_infinity 
mindiff = !values.d_infinity
age = !values.d_infinity

if keyword_set(debug) then print, mstar, rstar

;print, 'yyz=', yyz  ;***************************************************************************
;print, 'z_one', '(nz-1)', 'a_one', '(na-1)', 'm_one', '(nx-1)', 'z_one', 'a_one', 'm_one'
;print, z_one, (nz-1), a_one, (na-1), m_one, (nx-1), z_one, a_one, m_one, format='(9(f0.6,x))'
;print, ' ' 
;stop 
;
;; out of range
if z_one gt (nz-1) or a_one gt (na-1) or m_one gt (nx-1) or $
   z_one lt 0 or a_one lt 0 or m_one lt 0 then return, !values.d_infinity

;print, 'made it past out of range check' ;***************************************************************
;
;; this saves a lot of time (and imposes a solid prior) by not
;; interpolating/allowing values older than the universe
yytrack = dblarr(3,sz[4])
nages = 0

if afe eq 0d0 then begin
   repeat begin
      ;; 2D cubic convolution interpolation (metalicity and mass)
      yytrack[0,nages] = interpolate(tracks[0,*,*,nages,0],z_one,m_one,cubic=-0.5d0) ;; teff
      yytrack[1,nages] = interpolate(tracks[0,*,*,nages,1],z_one,m_one,cubic=-0.5d0) ;; rstar
      yytrack[2,nages] = interpolate(tracks[0,*,*,nages,2],z_one,m_one,cubic=-0.5d0) ;; age

      if ~finite(yytrack[0,nages]) then return, !values.d_infinity  ;;update from Jason 2015.08.31
      if ~finite(yytrack[1,nages]) then return, !values.d_infinity  ;;update from Jason 2015.08.31
      if ~finite(yytrack[2,nages]) then return, !values.d_infinity  ;;update from Jason 2015.08.31

      nages++
   endrep until yytrack[2,nages-1] ge 15 or nages eq sz[4]
endif else begin
   repeat begin
      ;; 3D trilinear interpolation (afe, metalicity, and mass)
      yytrack[0,nages] = interpolate(tracks[*,*,*,nages,0],a_one,z_one,m_one) ;; teff
      yytrack[1,nages] = interpolate(tracks[*,*,*,nages,1],a_one,z_one,m_one) ;; rstar
      yytrack[2,nages] = interpolate(tracks[*,*,*,nages,2],a_one,z_one,m_one) ;; age
      
      if ~finite(yytrack[0,nages]) then return, !values.d_infinity  ;;update from Jason 2015.08.31
      if ~finite(yytrack[1,nages]) then return, !values.d_infinity  ;;update from Jason 2015.08.31
      if ~finite(yytrack[2,nages]) then return, !values.d_infinity  ;;update from Jason 2015.08.31
      
      nages++
   endrep until yytrack[2,nages-1] ge 15 or nages eq sz[4]
endelse

yytrack = yytrack[*,0:nages-1]

;; locate all monotonic segments (inspired by lclxtrem.pro)
deriv = yytrack[0,1:nages-1] - yytrack[0,0:nages-2]
pos = where(deriv gt 0d0,complement=neg)
if pos[0] ne -1 then deriv[pos] = 1
zero = where(deriv[neg] eq 0d0)
if neg[0] ne -1 then deriv[neg] = -1
if zero[0] ne -1 then deriv[neg[zero]] = 0
deriv2 = deriv[1:nages-2] - deriv[0:nages-3]
ndx = where(deriv2 ne 0) + 1
if n_elements(ndx) eq 1 then ndx = [0,nages-1] $
else ndx = [0,ndx,nages-1]

for i=0, n_elements(ndx)-2 do begin
   ;; do not extrapolate, only interpolate
   if ((yytrack[1,ndx[i]] ge rstar) and (yytrack[1,ndx[i+1]] le rstar)) or $
      ((yytrack[1,ndx[i]] le rstar) and (yytrack[1,ndx[i+1]] ge rstar)) then begin

      yyteffnew = interpol(yytrack[0,ndx[i]:ndx[i+1]], yytrack[1,ndx[i]:ndx[i+1]],rstar,quad=quad, lsquad=lsquad)
      ;; if this segment is provides a better match to the observed Teff, use it.
      if abs(yyteffnew - teff) lt mindiff then begin
         yyteff = yyteffnew
         mindiff = abs(yyteff - teff)
         age = interpol(yytrack[2,ndx[i]:ndx[i+1]],yytrack[1,ndx[i]:ndx[i+1]],rstar,quad=quad, lsquad=lsquad)
      endif
      
;      if keyword_set(debug) then begin
;         ;print, yyteffnew, teff, mindiff, abs(yyteff - teff)
;         oplot, yytrack[0,ndx[i]:ndx[i+1]],yytrack[1,ndx[i]:ndx[i+1]],color=colors[i mod ncolors]
;         oplot, [yyteff],[rstar],psym=2,color=colors[i mod ncolors]
;      endif
   endif

endfor

chi2 = ((yyteff-teff)/uteff)^2

if keyword_set(debug) then begin 
   if !d.name ne 'PS' then begin
      device,window_state=win_state
      if win_state[2] eq 1 then wset, 2 $
      else window, 2
   endif
   print, 'mstar    ', 'rstar    ', 'feh      ', 'yyz      ', 'teff     ', 'yyteff   ', 'age      ', 'uteff    ', 'chi2'
   print, mstar, rstar, feh, yyz, teff, yyteff, age, uteff, chi2, format='(9(f0.6,x))'

   yyloggall = alog10(27443.4141d0*mstar/yytrack[1,*]^2)
   logg = alog10(27443.4141d0*mstar/rstar^2)
   plot, yytrack[0,*], yyloggall, xrange=[max(yytrack[0,*]),min(yytrack[0,*])],xtitle=textoidl('T_{eff} (K)'),ytitle='log g',/ynozero
   oplot, [teff], [logg], psym=2
;   plot, yytrack[0,*],yytrack[1,*],/ylog,xrange=[2500,8000],yrange=[0.1,10],$
;         ytitle=textoidl('R_*'),xtitle=textoidl('T_{eff}')
endif

return, chi2

stop
end
