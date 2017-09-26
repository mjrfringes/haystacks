; Solar System zoo plot
; Aki Roberge
; NASA/Goddard Space Flight Center
; Greenbelt, MD 20770, USA
;

au_km = 1.49598d8			    ; 1 AU in km
planets = ['Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune']
p_radii = [6052., 6371.009, 3390., 69911., 58232., 25362., 24622.]/au_km
a = [0.72333199,1.00000011,1.52366231,5.20336301,9.53707032,19.19126393,30.06896348]

proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'
spec_final_dir = proj_dir + 'Spectra_Final/'
paper_fig_dir = '~/Other_Projects/Haystacks/Paper_2015/Haystacks/'

readcol,spec_final_dir+planets[0]+'_geo_albedo.txt',pwl_v,palb_v
readcol,spec_final_dir+planets[1]+'_geo_albedo.txt',pwl_e,palb_e
readcol,spec_final_dir+planets[2]+'_geo_albedo.txt',pwl_m,palb_m
readcol,spec_final_dir+planets[3]+'_geo_albedo.txt',pwl_j,palb_j
readcol,spec_final_dir+planets[4]+'_geo_albedo.txt',pwl_s,palb_s
readcol,spec_final_dir+planets[5]+'_geo_albedo.txt',pwl_u,palb_u
readcol,spec_final_dir+planets[6]+'_geo_albedo.txt',pwl_n,palb_n

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
; Paper figures
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!p.multi = [0,1,2]

; Geometric albedo

cv = cgcolor('Orange')
ce = cgcolor('Turquoise')
cm = cgcolor('Orange Red')
cj = cgcolor('Purple')
cs = cgcolor('Dark Goldenrod')
cu = cgcolor('Teal')
cn = cgcolor('Blue')
csol = cgcolor('Gold')

set_plot,'ps'
device,filename=paper_fig_dir+'Solar_System_albedo.eps',/color,/inches,$
	xsize=8,ysize=8,/encapsulated,/portrait,/helvetica,/bold,/isolatin1
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=1.8
plot,pwl_e,palb_e,xtitle='Wavelength ('+!mu+'m)', /nodata, $
    ytitle='Geometric albedo',yrange=[-0.05,1.0],xrange=[0.3,2.5], $
    xmargin=[6.5,1.5], ymargin=[3.25,1.0]
oplot,[0.01,1d2],[0.0,0.0],linestyle=1
; Venus
oplot,pwl_v,palb_v,color=cv
; Earth
oplot,pwl_e,palb_e,color=ce
; Mars
oplot,pwl_m,palb_m,color=cm
al_legend,[planets[0],planets[1],planets[2]],/top,/left,thick=4,$
	line=[0,0,0],color=[cv,ce,cm],linsize=0.35,pos=[0.70,0.95],/norm,box=0,$
	charsize=1.5
;
plot,pwl_j,palb_j,xtitle='Wavelength ('+!mu+'m)', /nodata, $
    ytitle='Geometric albedo',yrange=[-0.05,1.0],xrange=[0.3,2.5], $
    xmargin=[6.5,1.5], ymargin=[3.25,1.0]
oplot,[0.01,1d2],[0.0,0.0],linestyle=1
; Jupiter
oplot,pwl_j,palb_j,color=cj
; Saturn
oplot,pwl_s,palb_s, color=cs
; Uranus
oplot,pwl_u,palb_u,color=cu
; Neptune
oplot,pwl_n,palb_n,color=cn
al_legend,[planets[3],planets[4],planets[5],planets[6]],/top,/left,thick=4,$
	line=[0,0,0,0],color=[cj,cs,cu,cn],linsize=0.35,pos=[0.70,0.45],/norm,box=0,$
	charsize=1.5
device,/close
set_plot,'x'

; Reflectance

!p.multi = [0,1,2]

width = 2.5 - 0.3
centwave = width/2. + 0.3
wbin = centwave/500.
steps = round(width/wbin) + 1
wave = findgen(steps)*wbin + 0.3    ; Array of wavelengths in microns
nlam = n_elements(wave)
print,'Min/Max wavelenth: '+strcompress(string(min(wave),format='(F6.3)'),/remove_all) $ 
  	+' - '+strcompress(string(max(wave),format='(F6.3)'),/remove_all)+' um'
print,'Number of wavelengths: '+strtrim(nlam,2)

beta = !pi/2.
lambert_pf = (sin(beta) + (!PI - beta) * cos(beta)) / (!PI)

venus_tmp = interpol(palb_v,pwl_v,wave)
venus = venus_tmp * lambert_pf * (p_radii[0]/a[0])^2 
earth_tmp = interpol(palb_e,pwl_e,wave)
earth = earth_tmp * lambert_pf * (p_radii[1]/a[1])^2 
mars_tmp = interpol(palb_m,pwl_m,wave)
mars = mars_tmp * lambert_pf * (p_radii[2]/a[2])^2 
jupiter_tmp = interpol(palb_j,pwl_j,wave)
jupiter = jupiter_tmp * lambert_pf * (p_radii[3]/a[3])^2 
saturn_tmp = interpol(palb_s,pwl_s,wave)
saturn = saturn_tmp * lambert_pf * (p_radii[4]/a[4])^2 
uranus_tmp = interpol(palb_u,pwl_u,wave)
uranus = uranus_tmp * lambert_pf * (p_radii[5]/a[5])^2 
neptune_tmp = interpol(palb_n,pwl_n,wave)
neptune = neptune_tmp * lambert_pf * (p_radii[6]/a[6])^2 

set_plot,'ps'
device,filename=paper_fig_dir+'Solar_System_reflectance.eps',/color,/inches,$
	xsize=8,ysize=8,/encapsulated,/portrait,/helvetica,/bold,/isolatin1
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=1.8
plot,wave,earth,xtitle='Wavelength ('+!mu+'m)', /nodata,/ylog,$
    ytitle='Reflectance (F!Dp!N / F!DSun!N)',xrange=[0.3,2.5], yrange=[1d-16,5d-9],$
    xmargin=[8.0,1.5], ymargin=[3.25,1.0]
; Venus
oplot,wave,venus,color=cv
; Earth
oplot,wave,earth,color=ce
; Mars
oplot,wave,mars,color=cm
al_legend,[planets[0],planets[1],planets[2]],/top,/left,thick=4,$
	line=[0,0,0],color=[cv,ce,cm],linsize=0.35,pos=[0.17,0.715],/norm,box=0,$
	charsize=1.5
;
plot,wave,jupiter,xtitle='Wavelength ('+!mu+'m)', /nodata, /ylog, $
    ytitle='Reflectance (F!Dp!N / F!DSun!N)',xrange=[0.3,2.5], yrange=[1d-16,5d-9], $
    xmargin=[8.0,1.5], ymargin=[3.25,1.0]
; Jupiter
oplot,wave,jupiter,color=cj
; Saturn
oplot,wave,saturn,color=cs
; Uranus
oplot,wave,uranus,color=cu
; Neptune
oplot,wave,neptune,color=cn
al_legend,[planets[3],planets[4],planets[5],planets[6]],/top,/left,thick=4,$
	line=[0,0,0,0],color=[cj,cs,cu,cn],linsize=0.35,pos=[0.17,0.25],/norm,box=0,$
	charsize=1.5
device,/close
set_plot,'x'

; Spectral flux density

!p.multi = 0

readcol,spec_final_dir+'SunKuruczSpectrum.txt',wl_sun,flux_sun
sun = interpol(flux_sun,wl_sun,wave)

set_plot,'ps'
device,filename=paper_fig_dir+'Solar_System_absolute.eps',/color,/inches,$
	xsize=8,ysize=9,/encapsulated,/portrait,/helvetica,/bold,/isolatin1
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=1.8
plot,wave,venus*sun, xtitle='Wavelength ('+!mu+'m)', $
	ytitle='Spectral Flux Density (Jy)', yrange=[1d-15,1d-5], $
	xrange=[0.3,2.5], xmargin=[8.0,1.5], ymargin=[3.25,1.0], /nodata, /ylog, $
	xticklen=0.02, yticklen=0.02, $
	ytickname=['10!U-16!N', '10!U-14!N', '10!U-12!N', '10!U-10!N', '10!U-8!N', '10!U1!N']
oplot,wave, venus*sun, color=cv
oplot,wave, mars*sun, color=cm
oplot,wave, jupiter*sun, color=cj
oplot,wave, saturn*sun, color=cs
oplot,wave, uranus*sun, color=cu
oplot,wave, neptune*sun, color=cn
oplot,wave, earth*sun, color=ce
oplot,wave, sun*1d-7, color=csol
axis, 0.3, 3d-7, xaxis=0, xtickname=replicate(' ',5)
;
al_legend,[planets[3],planets[0],planets[4],planets[1],planets[2],planets[5], $
	planets[6]],/top,/left,thick=4,$
	line=[0,0,0,0,0,0,0], color=[cj,cv,cs,ce,cm,cu,cn], linsize=0.35, $
	pos=[0.17,0.35],/norm,box=0,charsize=1.7
;
al_legend,['Sun'],/top,/left,thick=4,line=[0],color=[csol],linsize=0.35, $
	pos=[0.7,0.9],/norm,box=0,charsize=1.7
device,/close
set_plot,'x'

stop

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
; My figures
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

; Geometric albedo 

tmp0 = plot(pwl_e,palb_e, xtitle='Wavelength ('+!mu+'m)', $
    ytitle='Geometric albedo', dimensions=[800,900],$
    yrange=[-0.1,1.0],xrange=[0.1,3.0],font_name='Helvetica',$
    font_size=18,margin=[0.12,0.095,0.04,0.04], /nodata)
tmp1 = plot([0.01,1d2],[0.0,0.0],linestyle=1,/overplot)

; Venus
tmp2 = plot(pwl_v,palb_v,color='gold',/overplot,name=planets[0],thick=2)
; Earth
tmp3 = plot(pwl_e,palb_e,color='black',/overplot,name=planets[1],thick=2)
; Mars
tmp4 = plot(pwl_m,palb_m,color='orange red',/overplot,name=planets[2],thick=2)
; Jupiter
tmp5 = plot(pwl_j,palb_j,color='blue violet',/overplot,name=planets[3],thick=2)
; Saturn
tmp6 = plot(pwl_s,palb_s,color='peru',/overplot,name=planets[4],thick=2)
; Uranus
tmp7 = plot(pwl_u,palb_u,color='teal',/overplot,name=planets[5],thick=2)
; Neptune
tmp8 = plot(pwl_n,palb_n,color='blue',/overplot,name=planets[6],thick=2)

l = legend(TARGET=[tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8], $
  POSITION=[0.87,0.91],/norm, /AUTO_TEXT_COLOR, font_size=18)

datetag_new

tmp0.Save,spec_final_dir+'Plots/Solar_System_albedo.png',resolution=200

tmp0.Close

; Absolute flux

!psym = 10

tmp0 = plot(wave,venus*sun, xtitle='Wavelength ('+!mu+'m)', $
    ytitle='Spectral Flux Density (Jy)', dimensions=[800,900],$
    yrange=[1d-12,1d-5],xrange=[0.3,1.0],font_name='Helvetica',$
    font_size=18,margin=[0.14,0.095,0.04,0.04], /nodata, /ylog, $
    xticklen=0.02, yticklen=0.02, $
    ytickname=['10!U-12!N','10!U-11!N', '10!U-10!N', '10!U-9!N', '10!U-8!N', '10!U-7!N', '10!U1!N', '10!U2!N'])
ax = tmp0.AXES

; Venus
tmp2 = plot(wave, venus*sun,color='orange',/overplot,name=planets[0],thick=2,/histogram)
; Earth
tmp3 = plot(wave,earth*sun,color='black',/overplot,name=planets[1],thick=2,/histogram)
; Mars
tmp4 = plot(wave,mars*sun,color='orange red',/overplot,name=planets[2],thick=2,/histogram)
; Jupiter
tmp5 = plot(wave,jupiter*sun,color='blue violet',/overplot,name=planets[3],thick=2,/histogram)
; Saturn
tmp6 = plot(wave,saturn*sun,color='peru',/overplot,name=planets[4],thick=2,/histogram)
; Uranus
tmp7 = plot(wave,uranus*sun,color='teal',/overplot,name=planets[5],thick=2,/histogram)
; Neptune
tmp8 = plot(wave,neptune*sun,color='blue',/overplot,name=planets[6],thick=2,/histogram)
; Sun
xaxis = axis('X',location=3d-7,tickname=[' ',' ',' ',' ',' ',' ',' ',' '])
tmp9 = plot(wave,sun*1d-7,color='gold',/overplot,name='Sun',thick=2,/histogram)

l1 = legend(TARGET=[tmp9], POSITION=[0.90,0.85],/norm, /AUTO_TEXT_COLOR, font_size=14)
l1 = legend(TARGET=[tmp5,tmp2,tmp6,tmp3], $
  POSITION=[0.43,0.25],/norm, /AUTO_TEXT_COLOR, font_size=14)
l2 = legend(TARGET=[tmp4,tmp7,tmp8], $
  POSITION=[0.64,0.25],/norm, /AUTO_TEXT_COLOR, font_size=14)

tmp0.Save,spec_final_dir+'Plots/Solar_System_absolute.png',resolution=200

tmp0.Close

END