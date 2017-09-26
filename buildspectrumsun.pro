; Build Sun Spectra - 2014-07-30
; Andrew Lincowski
; NASA/Goddard Space Flight Center
; Greenbelt, MD 20770, USA
;
; Previous versions in Python by Maxime Rizzo, Erika Nesvold, and Ashlee 
; Wilkins. Small edits by Aki Roberge, 2014-10-03.

pro buildSpectrumSun

  sun_name = 'SunKurucz'

; File structure and load files
  proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'
  spec_input_dir = proj_dir + 'Spectra_Input/'
  spec_final_dir = proj_dir + 'Spectra_Final/'
  spawn,'mkdir ' + spec_final_dir
 
; Kurucz Sun spectrum. wl in nm, sflam is 1st flux moment (Eddington flux) 
;  at Sun's surface in ergs/cm**2/s/ster/nm, sflam_cont is Eddington flux 
;  for solar continuum only.

  sun_file = 'kurucz_sun.txt'    
  readcol,spec_input_dir+sun_file,swl,sflam,sflam_cont   

; Set required parameters
  dist = 10 * 3.08568d13        ; 10 pc in km
  Rstar = 6.955d5               ; Radius of Sun in km

; Convert to flux at observer 10 pc from Sun.
  sflam *= 4 * !PI * (Rstar/dist)^2  	; erg/s/cm**2/nm

; Convert to fnu from flam
  c = 2.9979e17                 ; Speed of light (nm/s)
  sfnu = sflam * swl^2 / c      ; Displacement law (erg/s/cm**2/Hz)
  sfnu *= 1e23                  ; Convert to Janskys

  swl /= 1d3                   ; Convert nm to um

; Now for fun lets also output the spectrum
  tmp1 = plot(swl,sfnu, xtitle='Wavelength ('+!mu+'m)', $
    ytitle='Flux (Jy)', $
    xrange=[0.1,3.0],font_name='Helvetica',$
    font_size=18,margin=[0.14,0.12,0.04,0.04])
  text1 = text(0.75,0.8,'Sun',/norm,/overplot, font_name='Helvetica',$
    font_size=20,font_style='Bold')
  datetag_new

; Save spectrum files
  save,filename=spec_final_dir+sun_name+'.sav',swl,sfnu
  writecol,spec_final_dir+sun_name+'Spectrum.txt',swl,sfnu
  tmp1.Save, spec_final_dir+'Plots/Sun.png',resolution=200

END


