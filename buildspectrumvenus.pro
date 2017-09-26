; Build Venus Spectrum - 2014-10-07
; Aki Roberge
; NASA/Goddard Space Flight Center
; Greenbelt, MD 20770, USA
;
; Previous versions in Python by Maxime Rizzo, Erika Nesvold, Ashlee Wilkins, 
; and Andrew Lincowski.

pro buildSpectrumVenus

  planet = 'Venus'

; File structure and load files
  proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'
  spec_input_dir = proj_dir + 'Spectra_Input/'
  spec_final_dir = proj_dir + 'Spectra_Final/'
  pl_dir = 'Venus/'

  spawn,'mkdir ' + spec_final_dir

  pl_file = 'venus_flx_refl.dat'     

  ; Planet spectrum: wavelength in um, pflux is the top-of-the-atmosphere flux from 
  ;   Venus, psol is the solar flux at Venus, and pref is the ``reflectance''.
  readcol,spec_input_dir+pl_dir+pl_file,pwl,pflux,psol,pref 
  	
  ; Convert to geometric albedo.
  
  palb = (2/3.) * pflux / (psol/2.)
  
; Now for fun lets also output the spectrum
  tmp1 = plot(pwl,palb, xtitle='Wavelength ('+!mu+'m)', $
    ytitle='Geometric albedo', $
    yrange=[-0.1,1.0],xrange=[0.1,3.0],font_name='Helvetica',$
    font_size=18,margin=[0.14,0.12,0.04,0.04])
  tmp2 = plot([0.01,1d2],[0.0,0.0],linestyle=1,/overplot)
  text1 = text(0.75,0.8,planet,/norm,/overplot, font_name='Helvetica',$
    font_size=20,font_style='Bold')
  datetag_new

x = where(pwl ge 0.55)
ref_alb = palb[x[0]]
print,'Geometric albedo at '+string(pwl[x[0]],format='(F4.2)')+' um = '$
	+string(ref_alb,format='(F5.3)')

; Save spectrum files
  save,filename=spec_final_dir+planet+'_geo_albedo.sav',palb
  writecol,spec_final_dir+planet+'_geo_albedo.txt',pwl,palb
  tmp1.Save, spec_final_dir+'Plots/'+planet+'.png',resolution=200

END
