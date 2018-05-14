; 
; Haystacks v3.0 - 2018-05
; 
; v3.0 is a version of Haystacks best used for making all the cubes at once 
; Contributors: Andrew Lincowski, Aki Roberge
;
; Pieces of this code are copied, borrowed, or are some distant
; ancestor of code written in IDL or python by: Maxime Rizzo, Erika Nesvold, 
;   Ashlee Wilkins, and Chris Stark.
; Dustmap is from Kuchner & Stark (2010), revised by Chris Stark in 2014.
; ZODIPIC is from Kuchner (2007), revised by Chris Stark in 2014.
;

FUNCTION Haystacks,inc,systype,minlam,maxlam
; Usage:
; outvar=haystacks(inc,systype,minlam,maxlam)
; outvar returns the name of the file with its destination folder
; This function exports FITS cubes with all the required information in the header
; Data is organized into multiple FITS extensions for convenience
; Please refer to the Haystacks website for an explanation of the data structure
; https://asd.gsfc.nasa.gov/projects/haystacks/haystacks.html
; inc: inclination in degrees
; systype: either 'archean' or 'modern' depending on the desired epoch
; minlam: minimum wavelength in microns
; maxlam: maximum wavelength in microns


  t_start = systime(/sec)

; Project directory. Write out full path name (DUSTMAP cannot parse ~).
proj_dir="/local/data/nicolaus2/mrizzo/Haystacks/Solar_System/"
; Establish variables
  
;  inc = 0					; Default to face-on system
  zodi = 1					; Default to 1 zodi of exozodi
  const = 1					; Default to no local zodiacal emission
;  band = 1					; Default to Waveband 1 (0.3-0.6 um)
  multi = 0					; Default to no multithreading
  smoothcloud = 1           ; Default 1 = run Zodipic (inner dust)
  dust = 1                  ; Default 1 = run dustmap (outer dust)
  sun = 1                   ; Default 0 = don't add Sun spectrum  

 
;  sres = 300.				; Spectral resolution
  sres = 1.				; Spectral resolution
  g = 0.17                   ; H-G scattering phase function asymmetry parameter
;  gauss = 24               	; Kernel size
  gauss = 4               	; Kernel size
  dist = 10.               	; Distance in pc
  ref_dist = 10.           	; Reference Distance in pc
  imgsize = 601           	; Size of image in pixels. MUST BE AN ODD NUMBER
;  imgsize = 3333           	; Size of image in pixels. MUST BE AN ODD NUMBER
  res = 0.03          	; AU per pixel
  epoch_dt = 0.            	; Time from J2000, in years
  epoch_dt = 25.73 			; Epoch for Earth and Neptune in rough alignment along x-axis
  
;  date = (systime(/JULIAN)-365.25*6712.)/365.25 	; No day like today!
;  epoch_dt = date

  c_img = round(imgsize / 2. - 0.5)    ; Image center 

; Read in options
  
;  read,"Inclination of the system (in degrees): ",inc
;  read,"Zodi level of the system (1 zodi = Solar System, enter 1, 10, or 100): ",zodi
;  read,"System viewing distance (in pc): ",dist
;  read,"Desired spectral resolution: ",sres
;  read,"Run ZODIPIC? 1 = yes, 0 = no ",smoothcloud
;  read,"Run DUSTMAP? 1 = yes, 0 = no ",dust
;  read,'Add Sun? (1 = yes, 0 = no) ',sun
;  read,"Add local zodiacal dust flux?  (1 = yes, 0 = no) ",const  
;  read,"Waveband (enter 1 for 0.3-0.6 um, 2 for 0.6-1.2 um, 3 for 1.2-2.5 um): ",band

; Establish wavelength array in microns

; IF band eq 1 THEN wrange = [0.3,0.6]
; IF band eq 2 THEN wrange = [0.6,1.2]
; IF band eq 3 THEN wrange = [1.2,2.5]
; IF band eq 4 THEN wrange = [0.6,0.9]

 width = maxlam - minlam
 centwave = width/2. + minlam
 wbin = centwave/sres
 steps = round(width/wbin) + 1
 wls = findgen(steps)*wbin + minlam    ; Array of wavelengths in microns
 nlam = n_elements(wls)
; nlam = ALOG(maxlam/minlam)*sres
; end_vals = findgen(nlam+1)*ALOG(maxlam/minlam)/nlam + ALOG(minlam)
; mid_vals = 0.5*(end_vals[1:*]+end_vals[0:-2])
; endwls = EXP(end_vals)
; wls = EXP(end_vals)

 print,'Min/Max wavelenth: '+strcompress(string(min(wls),format='(F6.3)'),/remove_all) $ 
  	+' - '+strcompress(string(max(wls),format='(F6.3)'),/remove_all)+' um'
 print,'Number of wavelengths: '+strtrim(nlam,2)
 print,wls
; Make folders for output datafiles

; Strings for file saving
  incs = strcompress(fix(inc),/remove_all)
  zodis = strcompress(fix(zodi),/remove_all)
  minwl = strcompress(string(min(wls),format='(F0.2)'),/remove_all)
  maxwl = strcompress(string(max(wls),format='(F0.2)'),/remove_all)



; Establish spectra for Archean vs Modern

  IF systype EQ "modern" THEN BEGIN
	sun_spec_file = 'SunKuruczSpectrum.txt'
  	pl_names = ['Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune']
 	zodi_scl = 1
	folder1 = proj_dir + "Modern_zodi" + zodis + "inc" + incs + "/"
    
  ENDIF ELSE IF systype EQ "archean" THEN BEGIN
	sun_spec_file = 'SunKuruczArcheanSpectrum.txt'
  	pl_names = ['EarlyVenus','Hazy_ArcheanEarth','EarlyMars','Jupiter','Saturn','Uranus','Neptune']
 	zodi_scl = 3
	folder1 = proj_dir + "Archean_zodi" + zodis + "inc" + incs + "/"
  ENDIF


;  folder1 = proj_dir + "zodi" + zodis + "inc" + incs + "/" 
  spawn,'mkdir ' + folder1
  IF (sun eq 0) THEN sunstring = 'NoSun' ELSE sunstring = 'Sun'
  IF (const eq 0) THEN localstring = 'NoLocal' ELSE localstring = 'Local' 
  folder2 = folder1 + sunstring + '_' + localstring + '/'
  spawn,'mkdir ' + folder2

; Multithreading

  threads = !CPU.HW_NCPU        ; Number of CPU threads
;  if threads gt 2 then multi = 1
  if smoothcloud eq 0 or dust eq 0 then multi = 0
  print,'CPU Threads: '+string(threads)
  if multi eq 1 then print,'Multithreading...'

ZODIPIC:

; This section calls the wrapper that initializes ZODIPIC.
; This provides for multithreading concurrently with dustmap.
; Note zodipic multithreads itself too.

  z_start = systime(/sec)
  IF smoothcloud eq 1 THEN BEGIN

     ; Execute zodipic program initiator if no multithreading
     IF multi eq 0 THEN new_zodipic,wls,inc,zodi,g,proj_dir,folder1,ref_dist,res,epoch_dt 
     
     IF multi eq 1 THEN begin             	; If multithreading, start a new CPU thread
        obridge=obj_new("IDL_IDLBridge", output=' ')  ; Initiate IDL bridge/child process

        obridge->setvar, 'var1', wls                         
        obridge->setvar, 'var2', inc
        obridge->setvar, 'var3', zodi
        obridge->setvar, 'var4', g
        obridge->setvar, 'var5', proj_dir
        obridge->setvar, 'var6', folder1
        obridge->setvar, 'var7', ref_dist
        obridge->setvar, 'var8', res
        obridge->setvar, 'var9', epoch_dt

        obridge->execute, 'new_zodipic,var1,var2,var3,var4,var5,var6,var7,var8,var9', /nowait ; Run new_zodipic and continue
     ENDIF
  END
  z_end = systime(/sec)
  
DUSTMAP:
  
  d_start = systime(/sec)
  IF dust eq 1 THEN begin 
     ; Execute dustmap program initiator if no multithreading
     IF multi eq 0 THEN new_dustmap,wls,inc,zodi,g,proj_dir,folder1,ref_dist,imgsize,res,epoch_dt

     ; If multithreading, begin IDL procedures for child process
     IF multi eq 1 THEN BEGIN
        print,'Initiating concurrent dustmap set ...'
        dm=obj_new("IDL_IDLBridge", output=' ') ;initiate IDL bridge/child process
        dm->setvar, 'var1', wls
        dm->setvar, 'var2', inc
        dm->setvar, 'var3', zodi
        dm->setvar, 'var4', g
        dm->setvar, 'var5', proj_dir
        dm->setvar, 'var6', folder1
        dm->setvar, 'var7', ref_dist
        dm->setvar, 'var8', imgsize
        dm->setvar, 'var9', res
        dm->setvar, 'var10', epoch_dt
  
        dm->execute, 'new_dustmap,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10', /nowait
     ENDIF

  ENDIF
  d_end = systime(/sec)

  IF multi eq 0 THEN GOTO,FINALE_1

RUNNING:

  z_status = obridge->status(ERROR=z_err)
  d_status = dm->status(ERROR=d_err)
  print,'zodipic status: '+strcompress(string(z_status),/remove_all)
  print,'dustmap status: '+strcompress(string(d_status),/remove_all)

  IF (z_status eq 3) || (d_status eq 3) THEN BEGIN
     print,'ZODIPIC error, if any: '+strcompress(string(z_err))
     print,'dustmap error, if any: '+strcompress(string(d_err))
     GOTO,FINALE_1
  ENDIF

  print,'Both ZODIPIC and dustmap currently running...'
  
  WHILE (obridge->status() eq 1) && (dm->status() eq 1) DO wait,10

  z_status = obridge->status(ERROR=z_err)
  d_status = dm->status(ERROR=d_err)
  print,'zodipic status: '+strcompress(string(z_status),/remove_all)
  print,'dustmapstatus: '+strcompress(string(d_status),/remove_all)

  IF (dm->status() eq 0) THEN BEGIN
     d_end = systime(/sec)
     dtd = (d_end - d_start)/60. ; in minutes
     print,'Total dustmap run time: '+strcompress(string(dtd),/remove_all)+' minutes.'
     print,'ZODIPIC still running...'
     while obridge->status() eq 1 do wait,10
     z_end = systime(/sec)
     dtz = (z_end - z_start)/60. ; in minutes   
     print,'Total ZODIPIC run time: '+strcompress(string(dtz),/remove_all)+' minutes.' 
  ENDIF ELSE BEGIN
     IF (obridge->status() eq 2) || (obridge->status() eq 0) THEN BEGIN
        z_end = systime(/sec)
        dtz = (z_end - z_start)/60. ; in minutes   
        print,'Total ZODIPIC run time: '+strcompress(string(dtz),/remove_all)+' minutes.'
        print,'dustmap still running...'
        while dm->status() eq 1 do wait,10
        d_end = systime(/sec)
        dtd = (d_end - d_start)/60. ; in minutes
        print,'Total dustmap run time: '+strcompress(string(dtd),/remove_all)+' minutes.'
     ENDIF
  ENDELSE

  t_end = systime(/sec)
  dt = (t_end - t_start)/60.    ; in minutes
  print,'Total Run time: '+strcompress(string(dt),/remove_all)+' minutes.'
  
  z_status = obridge->status()
  d_status = dm->status()
  print,'zodipic status: '+strcompress(string(z_status),/remove_all)
  print,'dustmap status: '+strcompress(string(d_status),/remove_all)

  obj_destroy, obridge          ; Terminate IDL bridge/child process
  obj_destroy, dm

FINALE_1:

  IF multi eq 0 THEN BEGIN  
     t_end = systime(/sec)
     dt = (t_end - t_start)/60.  ; in minutes
     dtz = (z_end - z_start)/60. ; in minutes
     dtd = (d_end - d_start)/60. ; in minutes
     
     print,'Total Run time: '+strcompress(string(dt),/remove_all)+' minutes.'
     IF smoothcloud eq 1 THEN print,'Total ZODIPIC run time: '+strcompress(string(dtz),/remove_all)+' minutes.'
     IF dust eq 1 THEN print,'Total dustmap run time: '+strcompress(string(dtd),/remove_all)+' minutes.'
  ENDIF

SMOOTH:

; For smoothing dustmap
  IF dust eq 1 THEN BEGIN 
    s_start = systime(/sec)
    print,'Gaussian smoothing, kernel size ' + strcompress(gauss)

    d_file = 'DUSTMAP_inc' + incs + 'zodi' + zodis + 'dist' + $
      strcompress(fix(ref_dist),/remove_all) + '_waterice_' + minwl + '-' + maxwl + 'um.sav'
    restore,folder1 + d_file

    dustcube = totimage

    FOR k = 0, nlam - 1 DO BEGIN             ; Iterate over wavelengths
      cimage = totimage(*,*,k)               ; Cimage is a 2D array for each wavelength
      newimage = gauss_smooth(cimage,gauss,/edge_truncate)     ; Gaussian smoothing
      dustcube(*,*,k) = newimage             ; Copy smoothed image to cube
    ENDFOR

    file = 'DUSTMAP_inc' + incs + 'zodi' + zodis + 'dist' + $
      strcompress(fix(ref_dist),/remove_all) + '_waterice_'+ minwl + '-' + maxwl + 'um_smooth'
    save,file=folder1+file+'.sav',dustcube
    file_delete,folder1 + d_file
;    fits_write,folder1+file+'.fits',dustcube        ; Write cube FITS file

    print,"Dust smoothing complete."
    s_end = systime(/sec)
    dt_s = s_end - s_start
    print,'Dust smoothing runtime: ' + strcompress(dt_s) + ' seconds'
    
  ENDIF

DRAWDUSTCUBE:

  ddc_start = systime(/sec)
  
  print,'Drawing dust cube...'

  filezodipic = "ZODIPIC_zodi" + zodis + "inc" + incs + "dist" + $
    strcompress(fix(ref_dist),/remove_all) + "_" + minwl + '-' + maxwl + "um.sav" 
  restore,folder1+filezodipic   ; variable = fnu (x,y,lambda)


  filedustmap = 'DUSTMAP_inc' + incs + 'zodi' + zodis + 'dist' + $
    strcompress(fix(ref_dist),/remove_all) + '_waterice_'+ minwl + '-' + maxwl + 'um_smooth.sav'
  restore,folder1+filedustmap   ; variable = dustcube (x,y,lambda)
; Determine dimensions of dustcube
  npix = n_elements(dustcube(*,0,0))
  npix_z = n_elements(fnu(*,0,0))
  cpix_z = npix_z / 2
  cpix = npix / 2

; clean up files
  file_delete,folder1+filezodipic,folder1+filedustmap


; Patch zodipic into center of dustcube
  FOR k = 0, nlam - 1 DO BEGIN
     FOR i = 0, npix_z - 1 DO BEGIN
        FOR j = 0, npix_z - 1 DO BEGIN
			IF (cpix - cpix_z + i ge 0) AND (cpix - cpix_z + j ge 0) AND (cpix - cpix_z + i le npix-1) AND (cpix - cpix_z + j le npix-1) THEN BEGIN
				print,cpix - cpix_z + i,cpix - cpix_z + j
           		IF fnu(i,j,k) ge dustcube(cpix - cpix_z + i,cpix - cpix_z + j, k) THEN dustcube(cpix - cpix_z + i,cpix - cpix_z + j, k) = fnu(i,j,k)
			ENDIF
        ENDFOR
     ENDFOR
  ENDFOR
  
; Scale dust flux per pixel for distance so that surface brightness (flux / arcsec^2) 
; remains constant.
  IF (dist ne ref_dist) THEN BEGIN
  omega = (res/dist)^2		  	; Solid angle of pixel in arcsec^2
  ref_omega = (res/ref_dist)^2	; Solid angle of pixel in arcsec^2 for ref_dist.
  flux_scl = omega/ref_omega
  dustcube *= flux_scl
  ENDIF

; Scale the zodi level 
  dustcube *= zodi_scl  ; Scale up by 3 for young Solar System dust.


; Save file
;  file_new = folder2+'dust_zodi'+ zodis + 'inc' + incs + 'dist' + $ 
;    strcompress(string(fix(dist)),/remove_all)+'_' + minwl + '-' + maxwl + 'um'
;  fits_write,file_new+'.fits',dustcube
;  save,file=file_new+'.sav',dustcube
  
;  print,"Dust cube file "+file_new+" saved."
  ddc_end = systime(/sec)
  dt_ddc = ddc_end - ddc_start
  print,'Draw dust cube time: ' + strcompress(string(dt_ddc),/remove_all) + ' seconds'

STARSCALE:

; To rescale dust images from blackbody Sun to Kurucz Sun

  ss_start = systime(/sec)

  print, 'Scaling dust images to Sun model...'
  bbstar, 'Sun', wls, ref_dist, fstar_bb
  spec_dir = proj_dir + 'Spectra_Final/'
  
  readcol,spec_dir+sun_spec_file,swl,sfnu ; wl in um, sfnu in Jy
  sunspec = interpol(sfnu,swl,wls)        ; Interpolate between wavelengths
  ratio_z = sunspec / fstar_bb
  FOR i = 0, nlam -1 DO BEGIN
    dustcube[*,*,i] *= ratio_z[i]
  ENDFOR

; Save file
;  file_new = folder2+'dust_zodi'+ zodis + 'inc' + incs + 'dist' + $ 
;    strcompress(string(fix(dist)),/remove_all)+'_' + minwl + '-' + maxwl + 'um_scl'
;  fits_write,file_new+'.fits',dustcube
;  save,file=file_new+'.sav',dustcube
;  print,"Scaled dust cube file "+file_new+" saved."
  ss_end = systime(/sec)
  dt_ss = ss_end - ss_start
  print,'Scale dust cube time: ' + strcompress(string(dt_ss),/remove_all) + ' seconds'

LOCALZODI:

  IF const eq 1 THEN BEGIN
     print,"Now adding local zodiacal dust flux ..."

; Add zodiacal dust, i.e. observing from Earth means looking
; thru our own zodiacal dust. 23mag/arcsec^2 at 0.55 um, Johnson V filter = 
; 2.421e-6 Jy/arcsec^2. Calculated with http://www.gemini.edu/?q=node/11119
     fscale = 2.421d-6
     obs_scale = (res/dist)^2 * fscale  ; Jy/pix
     print,'Local zodiacal flux = ' + strtrim(obs_scale,2) + ' Jy/pix'
     
; Now pick out the Earth location from dustcube at 0.55 um.
; We will use this to determine local zodiacal flux at the other wavelengths.
     obs_flux = dustcube(c_img + 1/res,c_img,*)
     d = where(wls ge 0.55)
     lam_ref = wls[d[0]]
     print,'Reference wavelength for local zodi:'+string(lam_ref,format='(F6.3)')+' um'
     scale_ratio = (obs_scale/obs_flux[d[0]])
     FOR k = 0, nlam - 1 DO BEGIN    
        zodiflux = dustcube(c_img + 1/res,c_img,k) * scale_ratio
        print,'Zodi flux for '+string(wls[k],format='(F6.3)')+' um, in Jy per pix: ' $ 
          +strtrim(zodiflux,2)
        dustcube(*,*,k) += zodiflux
     ENDFOR

; Save file
;     file_new = folder2+'dust_zodi'+ zodis +'inc' + incs + 'dist' + $ 
;       strcompress(fix(dist),/remove_all) + '_'+ minwl + 'um_local'
;     fits_write,file_new+'.fits',float(dustcube)
;     save,file=file_new+'.sav',dustcube
;     print,"Local dust cube file "+file_new+" saved."
  ENDIF

PLANETCUBE:
  drawplanetcube,zodi,inc,dist,wls,imgsize,res,proj_dir,folder2,epoch_dt,sun,sun_spec_file,pl_names

BUILDCUBE:

  b_start = systime(/sec)

; Read in planets cube
  plfile = folder2 + 'planets_zodi' + zodis + 'inc' + incs + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all) + '_' + minwl + '-' + maxwl + 'um.fits'
  fits_read,plfile,plcube
  file_delete,plfile
  plspecfile = folder2 + 'planets_zodi' + zodis + 'inc' + incs + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all) + '_' + minwl + '-' + maxwl + 'um_spec.fits'
  fits_read,plspecfile,plspeccube,plspechdr
  file_delete,plspecfile
  orbparfile = folder2 + 'planets_zodi' + zodis + 'inc' + incs + 'dist' + $ 
    strcompress(string(fix(dist)),/remove_all) + '_' + minwl + '-' + maxwl + 'um_orbpar.fits'
  fits_read,orbparfile,orbparcube,orbparhdr
  file_delete,orbparfile

; Recenter dust image (due to its even size) so that center is one pixel
  print,'Recentering dust image...'

  FOR k = 0, nlam-1 DO BEGIN
    cimage = dustcube(*,*,k)
    cimage = rebin(cimage,npix * 2, npix * 2, /sample) ; Magnify image by 2
    cimage = shift(cimage,1,1)     ; Shift by 1/2 pixel of original image up and right
    cimage = rebin(cimage,npix,npix)  ; Shrink back to original image size
    dustcube(*,*,k) = cimage
  ENDFOR
  
; Put shifted dust image into array of same size as the planet image.
  npix_p = n_elements(plcube[0,*,0])
  tmpdustcube = dblarr(npix_p, npix_p, n_elements(wls))
  tmpdustcube[0:npix_p-2,0:npix_p-2,*] = dustcube[*,*,*]
  tmpdustcube[npix_p-1,*,*] = tmpdustcube[npix_p-2,*,*]
  tmpdustcube[*,npix_p-1,*] = tmpdustcube[*,npix_p-2,*]
  
; Add cubes together

  finalcube = plcube + tmpdustcube

  b_end = systime(/sec)
  dt_b = b_end - b_start
  print,'Final build runtime: ' + strcompress(string(dt_b),/remove_all) + ' seconds'

; Write files

  file_new = folder2 + systype+'_cube_zodi' + zodis +'inc' + incs + 'dist' + $
	strcompress(string(fix(dist)),/remove_all) + '_' + minwl + '-' + maxwl + 'um'
;  save,file=file_new+'.sav',wls,finalcube

;  writefits,file_new+'.fits',finalcube,head
  writefits,file_new+'.fits',0,head
  IF (sun eq 0) THEN h1 = 'F' ELSE h1 = 'T'
  sxaddpar,head,'SUN',h1,' Sun in central pixel',BEFORE='COMMENT'
  IF (const eq 0) THEN h2 = 'F' ELSE h2 = 'T'
  sxaddpar,head,'LOCAL',h2,' Local zodiacal background added',AFTER='SUN'
  sxaddpar,head,'COMMENT','Spectral image cube of the Solar System.',AFTER='LOCAL'
  sxaddpar,head,'TYPE',systype,'Modern or Archean?'
  sxaddpar,head,'BUNIT','Jy','Map units per pixel'
  sxaddpar,head,'INC',inc,'System inclination in degrees'
  sxaddpar,head,'ZODI',zodi*zodi_scl,'Number of zodis'
  sxaddpar,head,'MINLAM',minlam,'Minimum wavelength (um)'
  sxaddpar,head,'MAXLAM',maxlam,'Maximum wavelength (um)'
  IF (const eq 0) THEN h3 = 'F' ELSE h3 = 'T'
  sxaddpar,head,'SMOOTH',h3,'Whether to smooth the dustmap or not'
  sxaddpar,head,'GAUSS',gauss,'Smoothing kernel size'
  sxaddpar,head,'SRES',sres,'Spectral resolution Lam/dlam'
  sxaddpar,head,'G',g,'H-G scattering phase function asymmetry parameter'
  sxaddpar,head,'DIST',dist,'Distance to the system (pc)'
  sxaddpar,head,'PIXSCALE',res,'Pixel scale (AU)'
  sxaddpar,head,'EPOCH',epoch_dt,'Epoch in years since J2000'
  sxaddpar,head,'N_EXT',n_elements(wls),'Number of FITS extensions with cube data'
  

  modfits,file_new+'.fits',0,head
  
  FOR k = 0, nlam -1 DO BEGIN
    mkhdr,hdr,float(finalcube[*,*,k]),/IMAGE,/EXTEND
    sxaddpar,hdr,'BUNIT','Jy','Map units per pixel'
    sxaddpar,hdr,'WAVEL',wls[k],'Wavelength at the slice in microns'
    writefits,file_new+'.fits',float(finalcube[*,*,k]),hdr,/append
  ENDFOR

  mkhdr,hdr,float(wls),/IMAGE,/EXTEND
 sxaddpar,hdr,'BUNIT','microns',' Wavelength units'
  writefits,file_new+'.fits',float(wls),hdr,/append
  ;mwrfits,wls,file_new+'.fits',
  
  mkhdr,hdr2,float(sunspec),/IMAGE,/EXTEND
 sxaddpar,hdr2,'BUNIT','Jy',' Spectrum units'
  writefits,file_new+'.fits',float(sunspec),hdr2,/append
  ;mwrfits,subspec,file_new+'.fits'

;  file_wave = folder1 + 'wavelengths_' + minwl + 'um.fits'
; make new header to add wavelengths to previous FITS files
;  mkhdr,hwls,wls
;  sxaddpar,hwls,'BUNIT','microns',' Wavelength units'
; add wavelengths as a separate extension
;  mwrfits,wls,file_new+'.fits'


;  writefits,file_wave,wls,headw
;  sxaddpar,headw,'WUNIT','microns',' Wavelength units',BEFORE='COMMENT'
;  sxaddpar,headw,'COMMENT', $ 
;    ' Array of wavelengths for images in spectral image  cube.',AFTER='WUNIT'
;  modfits,file_wave,0,headw

  ; write planet spectra cube as an extension
mwrfits,plspeccube,file_new+'.fits',plspechdr
mwrfits,orbparcube,file_new+'.fits',orbparhdr


  print,'Dust cube ' + file_new + ' saved. Haystacks cube complete!'

  t_end = systime(/sec)
  dt = (t_end - t_start)/60.
  print,'Run time: '+strcompress(string(dt),/remove_all)+' minutes'
print,wls

RETURN,file_new+'.fits'
  
END
