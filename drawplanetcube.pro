; DrawPlanetCube - 2014-08-06
; Reads in formatted spectra to generate planet cube for Solar System at 10 pc.
; Andrew Lincowski
; NASA/Goddard Space Flight Center
; Greenbelt, MD 20770, USA
;
; Previous versions in Python by Maxime Rizzo, Erika Nesvold, and Ashlee Wilkins.
; Edited by Aki Roberge, 2014-10-06.

pro drawplanetcube,zodi,inc,dist,wls,imgsize,res,proj_dir,folder,epoch_dt,sun,sun_spec_file,pl_names

; Establish variables  
  ref_dist = 10.			    ; Distance from observer to system in pc
  au_km = 1.49598d8			    ; 1 AU in km
  c_img = imgsize / 2. - 0.5    ; Image center 

  nlam = n_elements(wls)        ; Number of wavelengths
  ;pl_names = ['Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune']
  numpl = n_elements(pl_names)  ; Number of planets
 
  cos_i = cos(inc * !PI / 180.)
  sin_i = sin(inc * !PI / 180.)

 ; Strings for file saving
  zodis = strcompress(fix(zodi),/remove_all)
  incs = strcompress(fix(inc),/remove_all)
  minwl = strcompress(string(min(wls),format='(F0.2)'),/remove_all)
  maxwl = strcompress(string(max(wls),format='(F0.2)'),/remove_all)


SPECTRA:

; Create empty planet datacube.
  plcube = fltarr(imgsize,imgsize,nlam)
  
; Required spectral files
  spec_dir = proj_dir + 'Spectra_Final/'

; Sun spectrum
  sun_spec_file = 'SunKuruczSpectrum.txt'
  readcol,spec_dir+sun_spec_file,swl,sfnu ; wl in um, sfnu in Jy
  sunspec = interpol(sfnu,swl,wls)        ; Interpolate between wavelengths
  nlams = n_elements(sunspec)             ; Number of wavelengths in Sun spectrum

; Create planet spectrum array.
  plspec = fltarr(nlams,numpl) 

; Open planet spectral files
  FOR i = 0, numpl - 1 DO BEGIN
     file = pl_names[i]+'_geo_albedo.txt'
	 ; Read in planet spectra into tmp array 
     readcol,spec_dir+file,pwl,palb,FORMAT='F',/SILENT  
     plalb_tmp = interpol(palb,pwl,wls) ; Interpolate between wavelengths
     plspec(*,i) = plalb_tmp   ; Copy individual planet albedos into combined array
  ENDFOR

LOCATIONS:

; Planet orbital parameters.
; This is set up where each vector is a list of information in planetary order, 
;  i.e. element 0 = Venus, 1 = Earth, etc.
; This data is from Table A.2 of Murray & Dermott's Solar System Dynamics.
; All angles get explicitly converted to radians by multiplying by pio180.
; WARNING: Earth and Neptune are manually entered to a specific location after 
;  cartesian coordinates are calculated.

  pio180 = !pi/180.
  ; Planet radii in km, converted to AU
  p_radii = [6052., 6371.009, 3390., 69911., 58232., 25362., 24622.]/au_km
  ; Semi-major axis in AU
  a = [0.72333199,1.00000011,1.52366231,5.20336301,9.53707032,19.19126393,30.06896348]
  ; Eccentricity 
  ecc = [0.00677323,0.01671022,0.09341233,0.04839266,0.05415060,0.04716771,0.00858587]
  ;ecc = [0., 0, 0, 0, 0, 0, 0]  ; COMMENT THIS OUT FOR FINAL CUBES
  ; Inclination
  incl = pio180*[3.39471,0.00005,1.85061,1.30530,2.48446,0.76986,1.76917]
  ;incl = pio180*[0., 0, 0, 0, 0, 0, 0]  ; COMMENT THIS OUT FOR FINAL CUBES
  ; Longitude of perihelion     
  longperi = pio180*[131.53298,102.94719,336.04084,14.75385,92.43194,170.96424,44.97135]
  ; Longitude of ascending node 
  longnode = pio180*[76.73936180,348.73936,49.57854,100.55615,113.71504,74.22988,131.72169]
  ; Mean longitude
  meanlong = pio180*[181.97973,100.46435,355.45332,34.40438,34.40438,313.23218,-90] 
  ; Calculate argument of perihelion
  argperi = longperi - longnode 
  ; Calculate mean anomaly
  meananom = meanlong - longperi 
  periods = a^(1.5)
  meananom += 2*!PI*epoch_dt/periods
; Convert to cartesian coordinates
  GM = 4 * !pi * !pi ; This only matters for the values of vx, vy, and vz, which we ignore
  cartesian,GM,a,ecc,incl,longnode,argperi,meananom,x,y,z,vx,vy,vz
  pos = fltarr(3,numpl)         ; Num planets, num dimensions
  pos[0,*] = x
  pos[1,*] = y
  pos[2,*] = z

  plpos = fltarr(3,numpl)

; Rotating Earth and Neptune to align with dust features
  rho_earth = sqrt(pos[0,1]^2 + pos[1,1]^2)
  rho_nept = sqrt(pos[0,6]^2 + pos[1,6]^2)
  pos[0,1] = rho_earth
  pos[1,1] = 0
  pos[0,6] = rho_nept
  pos[1,6] = 0

; Position planets consistent with disk inclination
  plpos(0,*) = pos(0,*)
  plpos(1,*) = pos(1,*) * cos_i - pos(2,*) * sin_i
  plpos(2,*) = pos(1,*) * sin_i + pos(2,*) * cos_i

; Actual distances from star to planet
  pl_r = fltarr(numpl)		; Establish array for distances
  pl_r = sqrt(plpos(0,*)^2 + plpos(1,*)^2 + plpos(2,*)^2) 

; Need to consider light scattering angle due to planet's location above or below the disk
; Calculate star-planet-observer angle.
; 0 deg = planet between star and observer
; 180 deg = planet behind star
; 90 or 270 are in quadrature

; Vector joining planet to observer pointing toward -z-dir (-z points out of image) 
  plnt_obs_vec = [0.,0.,-ref_dist*206265] 
  plnt_obs_vec = transpose(plnt_obs_vec)

; Calculating scattering angle: the angle between the star-planet vector and the 
;   planet-observer vector.
  cos_obs = fltarr(numpl)
  cos_obs = plnt_obs_vec # plpos / (pl_r * sqrt(total(plnt_obs_vec^2)))
; Beta is the planet phase, which is 180 degrees out of phase from the scattering angle.
  beta = !PI - acos(cos_obs)
; Lambert phase function
  lambert_pf = (sin(beta) + (!PI - beta) * cos(beta)) / (!PI)

; Locate pixel locations of planets.
; Locations have been verified as within aphelion/perihelion distances.
  plpix = fltarr(2,numpl)
  plpix(0,*) = round(plpos(0,*) / res) + c_img
  plpix(1,*) = round(plpos(1,*) / res) + c_img

orbpar = [[p_radii],[a],[ecc],[incl],[longperi],[longnode],[meanlong],[argperi],[periods],[meananom],[transpose(plpos(0,*))],[transpose(plpos(1,*))],[transpose(pl_r)],[transpose(plpix(0,*))],[transpose(plpix(1,*))]]



  planetcube = fltarr(imgsize,imgsize,nlam)

; Calculate planet flux
  FOR i = 0, numpl - 1 do begin
   IF plpix(0,i) lt imgsize THEN BEGIN
	IF plpix(1,i) lt imgsize THEN BEGIN
    		planetcube(plpix(0,i),plpix(1,i),*) += plspec(*,i) * lambert_pf(i) * $
        		(p_radii[i]/pl_r[i])^2 * sunspec
	ENDIF
   ENDIF
  ENDFOR

; Add Sun spectrum to center of planet cube
  IF sun eq 1 THEN BEGIN
     planetcube(c_img,c_img,*) += sunspec
  ENDIF
  
; Scale planet fluxes by distance (unresolved sources)
  
  IF (dist ne ref_dist) then planetcube *= (ref_dist/dist)^2

; Print final information
  w = where(wls ge 0.55)
  print,'Sun location: ',c_img,c_img
  sun_flux = planetcube(c_img,c_img,[w[0]])
  print,'Sun flux at '+string(wls[w[0]],format='(F5.3)')+' um: '+ $
       strtrim(sun_flux,2) + ' Jy'
  FOR i = 0, numpl - 1 do begin
	IF plpix(0,i) lt imgsize THEN BEGIN
	IF plpix(1,i) lt imgsize THEN BEGIN
     print,pl_names[i] + ' location: ',plpix[*,i]
     ref_flux = planetcube(plpix(0,i),plpix(1,i),[w[0]])
     print,pl_names[i] + ' flux at '+string(wls[w[0]],format='(F5.3)')+' um: '+$
       strtrim(ref_flux,2)+' Jy'
	ENDIF
	ENDIF
  ENDFOR

; Save planet cube
  plfile = folder + 'planets_zodi' + zodis + 'inc' + incs + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all) + '_' + $ 
    minwl + '-' + maxwl+'um.fits'
  plspecfile = folder + 'planets_zodi' + zodis + 'inc' + incs + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all) + '_' + $ 
    minwl + '-' + maxwl+'um_spec.fits'
  orbparfile = folder + 'planets_zodi' + zodis + 'inc' + incs + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all) + '_' + $ 
    minwl + '-' + maxwl+'um_orbpar.fits'
  ; create header file

  mkhdr,hdr,plspec
 sxaddpar,hdr,'NAME','Planetary albedos',' '
  sxaddpar,hdr,'BUNIT','Dimensionless',' Albedo'
fits_write,plspecfile,plspec,hdr
mkhdr,hdr,orbpar
 sxaddpar,hdr,'NAME','Orbital Parameters',' '
fits_write,orbparfile,orbpar,hdr
  fits_write,plfile,planetcube
  delvarx,planetcube

END
