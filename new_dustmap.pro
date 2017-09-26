;Dust Map Compiler v1.1mini - 2014-08-05
;Used to run and compile the dust maps using Chris Stark's dustmap program
;Andrew Lincowski
;NASA/Goddard Space Flight Center
;Greenbelt, MD 20770, USA

;Previous versions by Christopher Stark & Ashlee Wilkins

PRO new_dustmap,lambda,inc,zodi,g,proj_dir,folder,dist,imgsize,res,epoch_dt

; The following are the text formatted scales used in the collisional
; grooming algorithm.  The string to the right indicates the geometric
; optical depth corresponding to that scale in "zodis."
 
  CASE zodi OF 
     1: BEGIN
        scale='1.00e+25' 
        zodis='1'
     END
     10: BEGIN
        scale='2.20e+26' 
        zodis='10'
     END
     100: BEGIN
        scale='7.00e+27' 
        zodis='100'
     END
  ENDCASE

LOCATIONS:

; Planet orbital parameters.
; Determines the longitude of co-rotating planet (Neptune).
; This data is From Table A.2 of Murray & Dermott's Solar System Dynamics.
; All angles get explicitly converted to radians by multiplying by pio180.

  pio180 = !pi/180.
  a = [30.06896348]              ; semi-major axis
  ecc = [0.00858587]             ; eccentricity
  inc_orb = pio180*[1.76917]     ; inclination
  longperi = pio180*[44.97135]   ; longitude of perihelion
  longnode = pio180*[131.72169]  ; longitude of ascending node
  meanlong = pio180*[304.88003]  ; mean longitude
  argperi = longperi - longnode  ; calculate argument of perihelion
  meananom = meanlong - longperi ; calculate mean anomaly
  periods = a^(1.5)

  meananom += 2*!PI*epoch_dt/periods

; Convert to cartesian coordinates

  GM = 4 * !pi * !pi       ; This only matters for vx, vy, and vz, which we ignore
  cartesian,GM,a,ecc,inc_orb,longnode,argperi,meananom,x,y,z,vx,vy,vz

  long = atan(y[0],x[0])/pio180
  print,'Longitude: '+strtrim(long,2)

; The path to the collisional Kuiper belt data files

  data_path = proj_dir + '../Datafiles/kuiper_COLL/'

; Parameters that control image

  nlambda = n_elements(lambda)
  dust_imgsize = imgsize - 1		; DUSTMAP only does even image sizes.
  dust_res = res
  fov = dust_res * 1000. / dist * dust_imgsize ; mas
  pixsize = dust_res * 1000. / dist ; mas
  
; Establish image cubes for summing flux from dustmap
  
  totimage = fltarr(dust_imgsize,dust_imgsize,nlambda)
  tempimage = fltarr(dust_imgsize,dust_imgsize,nlambda)

;----- THE FOLLOWING SHOULD NOT NEED TO BE CHANGED -----
; Some physical model parameters that need to be defined
; Changing these isn't simple since the output of the
; collisional grooming algorithm depended on them.

  betalist = reverse([0.01882,0.02503,0.03330,0.04428,0.05890,0.07833,0.10418,0.13856,0.18428,0.24510,0.32598,0.43355])

  numbetas = n_elements(betalist)
  files_per_beta = 40
  betatextlist = string(betalist,format='(%"%0.5f")')
  particledensity = 1.0                                ; in g cm^-3
  particlesize_mic = 0.57 / (particledensity * betalist) ; radius of grains in microns
  composition = 'waterice'
;-------------------------------------------------------

COLD:

  FOR i = 0,numbetas-1 DO BEGIN 	; iterate for all listed betas
     beta = betatextlist[i]
     particlesize = particlesize_mic[i]
     cold_files = data_path+'cold/collision_results-beta_'+beta+'-scale_'+scale+'-'+strcompress(string(indgen(files_per_beta)+i*files_per_beta),/remove_all)+'.dat'
         dustmap,cold_files,tempimage,distance=dist,inclination=inc,longitude=long,fov=fov,pixelsize=pixsize,datatype=3,/scattered,lambda=lambda,Tstar=5777.,lstar=1.,rdust=particlesize,composition=composition,/hg,g=g,/nodisp,/reduceprint
     
     totimage += tempimage      ; add flux from temp image for each particle size (beta)
     
  ENDFOR

HOT:

  FOR i = 0,numbetas-1 DO BEGIN
     beta = betatextlist[i]
     particlesize = particlesize_mic[i]
     hot_files = data_path+'hot/collision_results-beta_'+beta+'-scale_'+scale+'-'+strcompress(string(indgen(files_per_beta)+i*files_per_beta+1000),/remove_all)+'.dat'
     dustmap,hot_files,tempimage,distance=dist,inclination=inc,longitude=long,fov=fov,pixelsize=pixsize,datatype=3,/scattered,lambda=lambda,Tstar=5777.,lstar=1.0,rdust=particlesize,composition=composition,g=g,/hg,/nodisp,/reduceprint

     totimage += tempimage

  ENDFOR

PLUTINO:

  FOR i=0,numbetas-1 DO BEGIN
     beta = betatextlist[i]
     particlesize = particlesize_mic[i]
     plut_files = data_path+'plut/collision_results-beta_'+beta+'-scale_'+scale+'-'+strcompress(string(indgen(files_per_beta)+i*files_per_beta+2000),/remove_all)+'.dat'
     dustmap,plut_files,tempimage,distance=dist,inclination=inc,longitude=long,fov=fov,pixelsize=pixsize,datatype=3,/scattered,lambda=lambda,Tstar=5777.,lstar=1.0,rdust=particlesize,composition=composition,g=g,/hg,/nodisp,/reduceprint
     
     totimage += tempimage
     
  ENDFOR

FINAL:
minwl = strcompress(string(min(lambda),format='(F0.2)'),/remove_all)
maxwl = strcompress(string(max(lambda),format='(F0.2)'),/remove_all)

file = folder+'DUSTMAP_inc' + strcompress(string(long(inc)),/remove_all) $
  + 'zodi' + zodis + 'dist' + strcompress(fix(dist),/remove_all) + '_' $ 
  +composition+'_'+minwl + '-' + maxwl+'um'

save,filename=file+'.sav',totimage
fits_write,file+'.fits',totimage

END
