; Build Neptune Spectrum - 2014-10-07
; Aki Roberge
; NASA/Goddard Space Flight Center
; Greenbelt, MD 20770, USA
;
; Previous versions in Python by Maxime Rizzo, Erika Nesvold, Ashlee Wilkins,
; and Andrew Lincowski.

pro buildSpectrumNeptune

  mycolors
  planet = 'Neptune'

; ##### Optical #####

; File structure and load files
  proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'
  spec_input_dir = proj_dir + 'Spectra_Input/'
  spec_final_dir = proj_dir + 'Spectra_Final/'
  pl_dir = 'Neptune/'

  spawn,'mkdir ' + spec_final_dir

  pl_file = 'jovian_titan_karkoschka_1998.txt'     

  ; Planet spectrum: wavelength in nm, poch4 is the methane absorption 
  ;   coefficient, and the poalbX columns are various albedos for
  ;   Jupiter, Saturn, Uranus, Neptune, and Titan.
  
  readcol,spec_input_dir+pl_dir+pl_file,powl,powla,poch4,poalbJ,poalbS,$
    poalbU,poalbN,poalbT,format='D' 		
  	
  powl /= 1d3                   ; Convert nm to um

  ; Convert to geometric albedo.
  
  poalb = poalbN
  
; Now for fun lets also output the spectrum
  tmp1 = plot(powl,poalb, xtitle='Wavelength ('+!mu+'m)', $
    ytitle='Geometric albedo', $
    yrange=[-0.1,1.0],xrange=[0.1,3.0],/nodata,font_name='Helvetica',$
    font_size=18,margin=[0.14,0.12,0.04,0.04])
  tmp2 = plot([0.01,1d2],[0.0,0.0],linestyle=1,/overplot)
  text1 = text(0.75,0.8,planet,/norm,/overplot, font_name='Helvetica',$
    font_size=20,font_style='Bold')
  datetag_new
  tmp3 = plot(powl,poalb,/overplot,color='orange')

; ##### Near-IR #####

  pl_file = 'plnt_Neptune.txt'     

  ; Planet spectrum: wavelength in um, pflux is the planet flux in W m-2 um-1,
  ;   and perr is the flux error in W m-2 um-1.

  readcol,spec_input_dir+pl_dir+pl_file,piwl,piflux,pierr,format='D' 
  
  z = where(piflux lt 0.0)       ; Replace negative numbers with 0
  piflux[z] = 0d0

  ; Sun spectrum 
  
  sun_name = 'SunKurucz'
  readcol,spec_final_dir+sun_name+'Spectrum.txt',swl,sfnu
  sfnu = interpol(sfnu,swl,piwl)

  rangeo = where((powl gt 0.960) and (powl lt 1.020))
  rangei = where((piwl gt 0.960) and (piwl lt 1.020))
  tmp_refl = (piflux/sfnu)
  scl = mean(poalb[rangeo]/tmp_refl[rangei])
  pialb = tmp_refl*scl
  tmp4 = plot(piwl,pialb,/overplot,color='r')
    
  ; ##### Total #####

  cut = 1.0
  a = where(powl le cut)
  b = where(piwl ge cut)
  tmp_pwl = [powl[a],piwl[b]]
  tmp_palb = [poalb[a],pialb[b]]
    
  pwl = findgen(5501)*4d-4 + 0.3
  palb = interpol(tmp_palb,tmp_pwl,pwl)
  tmp5 = plot(pwl,palb,/overplot)

  ; Patch over NIR gap with a Burrows 15 AU Jupiter model.
  
  readcol,'~/Other_Projects/Haystacks/Planet_Spectra/Burrows/jup_15AU.txt', $ 
	w_b,f_b,ref_b
	
  rangeb = [1.8178, 1.8809]
  sclb = 2.07d8
  a = where((w_b ge rangeb[0]) and (w_b le rangeb[1]))
  b = where((pwl ge rangeb[0]) and (pwl le rangeb[1]))
  patch = interpol(ref_b[a]*sclb,w_b[a],pwl[b])
  palb[b] = patch

  x = where(pwl ge 0.55)
  ref_alb = palb[x[0]]
  print,'Geometric albedo at '+string(pwl[x[0]],format='(F4.2)')+' um = '$
	+string(ref_alb,format='(F5.3)')

 ; Save spectrum files and plot
  save,filename=spec_final_dir+planet+'_geo_albedo.sav',palb
  writecol,spec_final_dir+planet+'_geo_albedo.txt',pwl,palb
  tmp1.Save, spec_final_dir+'Plots/'+planet+'.png',resolution=200
 
END
