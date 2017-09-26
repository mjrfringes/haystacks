;New Zodipic v1.0 - 2014-08-05
;Used to run ZODIPIC
;Andrew Lincowski
;NASA/Goddard Space Flight Center
;Greenbelt, MD 20770, USA

PRO new_zodipic,wls,inc,zodi,g,proj_dir,folder,dist,res,dt

  restore, proj_dir+'../IDL/zodipic.2.1/QabsVars.dat' ;required for zodipic

; Fixed ZODIPIC input values for Solar System model

  posang = 90.			 ; position angle for image (E of N)
  ring = 1				 ; scaling factor for Earth’s resonant ring dust structure. 
  bands = 1				 ; scaling factor for the dust bands from major asteroid families.
  blob = 1               ; scaling factor for the Earth’s trailing dust blob.
  radin = 0.03	         ; inner radius of the total dust structure in AU. 
  radout = 10.			 ; outer radius of the total dust structure in AU.
  outer_edge = 3.0	     ; outer edge of dust structure generated from Kelsall model. 
  						 ;   NOT USED.
  alpha = 4.5            ; ZODIPIC extrapolates from outer_edge to radout as r^-alpha.
  						 ;   NOT USED.

; Image handling

  nlam = n_elements(wls)  			 ; number of bins of wavelengths
  pixsize = res * 1000. / dist       ; pixel size in milli-arcseconds
  npix = ceil((2 * radout / res)/16.) * 16 ; image size in pixels, must be multiple of 16
  if npix lt 112 then npix = 112           ; mininum size for ZODIPIC

  fnu = fltarr(npix,npix,nlam)
  tempfnu = fltarr(npix,npix)

; Longitude of Earth
;
; Co-rotating planet orbital parameters.
; ZODIPIC specifies the location of the blob of Earth-trailing dust relative to the Earth.
; This data is From Table A.2 of Murray & Dermott's Solar System Dynamics.
; All angles get explicitly converted to radians by multiplying by pio180.
; earthlong is angle in degrees measured in the plane of the disk to the Earth
; position from the longitude of the ascending node.

  pio180 = !pi/180.
  a = [1.00000011] 						; semi-major axis
  ecc = [0.01671022] 					; eccentricity
  incl = pio180*[0.00005]               ; inclination
  longperi = pio180*[102.94719] 		; longitude of perihelion
  longnode = pio180*[348.73936] 		; longitude of ascending node
  meanlong = pio180*[100.46435] 		; mean longitude
  argperi = longperi - longnode 		; calculate argument of perihelion
  meananom = meanlong - longperi 		; calculate mean anomaly
  periods = a^(1.5)

  print,'Original mean anomaly: '+strtrim(meananom,2)
  meananom += 2*!PI*dt/periods
  print,'Years since epoch J2000: ' + strtrim(dt,2)
  print,'Mean anomoly: '+strtrim(meananom,2)

; Convert to cartesian coordinates

  GM = 4 * !pi * !pi   ; only matters for the values of vx, vy, and vz, which we ignore
  cartesian,GM,a,ecc,incl,longnode,argperi,meananom,x,y,z,vx,vy,vz

  print,'Earth x,y (AU): '+strtrim(x,2)+', '+strtrim(y)
  long = atan(y[0],x[0])/pio180
  long = 180 - long
  print,'Earth longitude: '+strtrim(long,2)

; Call ZODIPIC program

  FOR i = 0, nlam - 1 DO BEGIN
         zodipic,tempfnu,pixsize,wls(i),dist=dist,pixnum=npix,inclination=inc,positionangle=posang, earthlong=long,zodis=zodi,ring=ring,bands=bands,blob=blob,radin=radin,radout=radout,/hg,g=g,/nodisplay

     fnu(*,*,i) = tempfnu
     
  ENDFOR

  print,'Sum of ZODIPIC image flux: '
  sum = total(fnu)
  print,sum

; Define and write file for ZODIPIC
 minwl = strcompress(string(min(wls),format='(F0.2)'),/remove_all)
 maxwl = strcompress(string(max(wls),format='(F0.2)'),/remove_all)

  file = folder + 'ZODIPIC_zodi' + strcompress(string(fix(zodi)),/remove_all) + $
    'inc' + strcompress(string(fix(inc)),/remove_all) + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all)+'_'+ $ 
    minwl + '-' + maxwl+'um'
  save,filename=file+'.sav',fnu
  fits_write,file+'.fits',fnu

  print, 'Saved file '+file

end
