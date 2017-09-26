dist = 10.0		; System distance in pc
res = 0.1		; Pixel resolution in AU
sres = 150.		; Spectral resolution

; Project directory

proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'
folder = "zodi1inc0/"
dir1 = proj_dir + folder + 'Sun_NoLocal_e0i0/'

sunloc = [600,600]
venusloc = [598, 607] - 1
earthloc = [611, 601] - 1
marsloc = [593, 588] - 1
jupiterloc = [597, 653] - 1
saturnloc = [695, 583] - 1
uranusloc = [687, 773] - 1
neptuneloc = [902, 601] - 1

;
; 0.3 um to 0.6 um
;

print,'BLUE CUBE'

file_b = 'cube_zodi1inc0dist'+strtrim(fix(dist),2)+'_0.3um'
fits_read, dir1+file_b+'.fits',finalcube,head
fits_read, proj_dir + folder + 'wavelengths_0.3um.fits',wls_b,headw

sun_b = transpose(finalcube[sunloc[0],sunloc[1],*])

venus_b = transpose(finalcube[venusloc[0],venusloc[1],*])
dust_v_b = transpose(finalcube[venusloc[0],venusloc[1]+1,*])
earth_b = transpose(finalcube[earthloc[0],earthloc[1],*])
dust_e_b = transpose(finalcube[earthloc[0],earthloc[1]+1,*])
mars_b = transpose(finalcube[marsloc[0],marsloc[1],*])
dust_m_b = transpose(finalcube[marsloc[0]-1,marsloc[1],*])
jupiter_b = transpose(finalcube[jupiterloc[0],jupiterloc[1],*])
dust_j_b = transpose(finalcube[jupiterloc[0],jupiterloc[1]+1,*])
saturn_b = transpose(finalcube[saturnloc[0],saturnloc[1],*])
dust_s_b = transpose(finalcube[saturnloc[0],saturnloc[1]+1,*])
uranus_b = transpose(finalcube[uranusloc[0],uranusloc[1],*])
dust_u_b = transpose(finalcube[uranusloc[0],uranusloc[1]+1,*])
neptune_b = transpose(finalcube[neptuneloc[0],neptuneloc[1],*])
dust_n_b = transpose(finalcube[neptuneloc[0],neptuneloc[1]+1,*])

yn1 = strcompress([-11,-10,-9,-8,-7],/remove_all)
yn2 = ['10!U'+yn1,': !C!A. ','10','100']

mycolors
set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_300nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_b,venus_b,/nodata,yrange=[1d-11,1d-4],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75], $
	ytickname = yn2
oplot,wls_b,sun_b/1d6
oplot,wls_b,jupiter_b,color=7		; purple
oplot,wls_b,venus_b,color=5			; yellow
oplot,wls_b,earth_b,color=3			; green
oplot,wls_b,saturn_b,color=8		; orange
oplot,wls_b,uranus_b,color=11		; light blue
oplot,wls_b,mars_b,color=2			; red
oplot,wls_b,neptune_b,color=4		; blue
oplot,wls_b,dust_e_b,line=1
al_legend,['Sun'],line=[0],color=[0],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Jupiter','Venus','Saturn','Earth'], $ 
	line=[0,0,0,0],color=[7,5,8,3],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.78],/norm
al_legend,['Uranus','Dust per pix at Earth','Mars','Neptune'], $ 
	line=[0,1,0,0],color=[11,0,2,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.78],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/dust_spec_300nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_b,venus_b,/nodata,yrange=[3d-13,6d-9],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wls_b,dust_v_b,color=5		; yellow
oplot,wls_b,dust_e_b,color=3		; green
oplot,wls_b,dust_m_b,color=2		; red
oplot,wls_b,dust_j_b,color=7		; purple
oplot,wls_b,dust_s_b,color=8		; orange
oplot,wls_b,dust_u_b,color=11		; light blue
oplot,wls_b,dust_n_b,color=4		; blue
al_legend,['Dust per pix at Venus','Earth','Mars','Saturn'], $ 
	line=[0,0,0,0],color=[5,3,2,8],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Jupiter','Uranus','Neptune'], $ 
	line=[0,0,0],color=[7,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.61,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_nosun_300nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_b,jupiter_b/sun_b,/nodata,/ylog,yrange=[6d-13,1.7d-8], $ 
	ytitle='F!Dplanet!N / F!DSun!N', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wls_b,(jupiter_b-dust_j_b)/sun_b,color=7				; purple
oplot,wls_b,(venus_b-dust_v_b)/sun_b,color=5		; yellow
oplot,wls_b,(earth_b-dust_e_b)/sun_b,color=3		; green
oplot,wls_b,(saturn_b-dust_s_b)/sun_b,color=8		; orange
oplot,wls_b,(uranus_b-dust_u_b)/sun_b,color=11		; light blue
oplot,wls_b,(mars_b-dust_m_b)/sun_b,color=2			; red
oplot,wls_b,(neptune_b-dust_n_b)/sun_b,color=4		; blue
oplot,wls_b,dust_e_b/sun_b,line=1
al_legend,['Jupiter','Venus','Saturn','Earth'], $ 
	line=[0,0,0,0],color=[7,5,8,3],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Uranus','Dust per pix at Earth','Mars','Neptune'], $ 
	line=[0,1,0,0],color=[11,0,2,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

;
; 0.6 um - 1.2 um
;

print,'GREEN CUBE'

file_g = 'cube_zodi1inc0dist'+strtrim(fix(dist),2)+'_0.6um'

fits_read, dir1+file_g+'.fits',finalcube,head
fits_read, proj_dir + folder + 'wavelengths_0.6um.fits',wls_g,headw

sun_g = transpose(finalcube[sunloc[0],sunloc[1],*])

venus_g = transpose(finalcube[venusloc[0],venusloc[1],*])
dust_v_g = transpose(finalcube[venusloc[0],venusloc[1]+1,*])
earth_g = transpose(finalcube[earthloc[0],earthloc[1],*])
dust_e_g = transpose(finalcube[earthloc[0],earthloc[1]+1,*])
mars_g = transpose(finalcube[marsloc[0],marsloc[1],*])
dust_m_g = transpose(finalcube[marsloc[0]-1,marsloc[1],*])
jupiter_g = transpose(finalcube[jupiterloc[0],jupiterloc[1],*])
dust_j_g = transpose(finalcube[jupiterloc[0],jupiterloc[1]+1,*])
saturn_g = transpose(finalcube[saturnloc[0],saturnloc[1],*])
dust_s_g = transpose(finalcube[saturnloc[0],saturnloc[1]+1,*])
uranus_g = transpose(finalcube[uranusloc[0],uranusloc[1],*])
dust_u_g = transpose(finalcube[uranusloc[0],uranusloc[1]+1,*])
neptune_g = transpose(finalcube[neptuneloc[0],neptuneloc[1],*])
dust_n_g = transpose(finalcube[neptuneloc[0],neptuneloc[1]+1,*])

yn1 = strcompress([-12,-10,-8],/remove_all)
yn2 = ['10!U'+yn1,': !C!A. ','100']

mycolors
set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_600nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_g,venus_g,/nodata,yrange=[1d-12,1d-4],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75], ytickname = yn2
oplot,wls_g,sun_g/1d6
oplot,wls_g,jupiter_g,color=7		; purple
oplot,wls_g,venus_g,color=5			; yellow
oplot,wls_g,saturn_g,color=8		; orange
oplot,wls_g,earth_g,color=3			; green
oplot,wls_g,mars_g,color=2			; red
oplot,wls_g,uranus_g,color=11		; light blue
oplot,wls_g,neptune_g,color=4		; blue
oplot,wls_g,dust_e_g,line=1
al_legend,['Sun'],line=[0],color=[0],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Jupiter','Venus','Saturn','Earth'], $ 
	line=[0,0,0,0],color=[7,5,8,3],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.86],/norm
al_legend,['Mars','Dust per pix at Earth','Uranus','Neptune'], $ 
	line=[0,1,0,0],color=[2,0,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.86],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/dust_spec_600nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_g,dust_v_g,/nodata,yrange=[3d-12,6d-9],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wls_g,dust_v_g,color=5		; yellow
oplot,wls_g,dust_e_g,color=3		; green
oplot,wls_g,dust_m_g,color=2		; red
oplot,wls_g,dust_j_g,color=7		; purple
oplot,wls_g,dust_s_g,color=8		; orange
oplot,wls_g,dust_u_g,color=11		; light blue
oplot,wls_g,dust_n_g,color=4		; blue
al_legend,['Dust per pix at Venus','Earth','Mars','Saturn'], $ 
	line=[0,0,0,0],color=[5,3,2,8],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Jupiter','Uranus','Neptune'], $ 
	line=[0,0,0],color=[7,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.61,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_nosun_600nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_g,jupiter_g/sun_g,/nodata,/ylog,yrange=[1d-14,1d-7], $ 
	ytitle='F!Dplanet!N / F!DSun!N', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wls_g,(jupiter_g-dust_j_g)/sun_g,color=7				; purple
oplot,wls_g,(venus_g-dust_v_g)/sun_g,color=5		; yellow
oplot,wls_g,(saturn_g-dust_s_g)/sun_g,color=8		; orange
oplot,wls_g,(earth_g-dust_e_g)/sun_g,color=3		; green
oplot,wls_g,(mars_g-dust_m_g)/sun_g,color=2			; red
oplot,wls_g,(uranus_g-dust_u_g)/sun_g,color=11		; light blue
oplot,wls_g,(neptune_g-dust_n_g)/sun_g,color=4		; blue
oplot,wls_g,dust_e_g/sun_g,line=1
al_legend,['Jupiter','Venus','Saturn','Earth'], $ 
	line=[0,0,0,0],color=[7,5,8,3],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Mars','Dust per pix at Earth','Uranus','Neptune'], $ 
	line=[0,1,0,0],color=[2,0,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

; 1.2 um - 2.5 um

print,'RED CUBE'

file_r = 'cube_zodi1inc0dist'+strtrim(fix(dist),2)+'_1.2um'

fits_read, dir1+file_r+'.fits',finalcube,head
fits_read, proj_dir + folder + 'wavelengths_1.2um.fits',wls_r,headw

sun_r = transpose(finalcube[sunloc[0],sunloc[1],*])

venus_r = transpose(finalcube[venusloc[0],venusloc[1],*])
dust_v_r = transpose(finalcube[venusloc[0],venusloc[1]+1,*])
earth_r = transpose(finalcube[earthloc[0],earthloc[1],*])
dust_e_r = transpose(finalcube[earthloc[0],earthloc[1]+1,*])
mars_r = transpose(finalcube[marsloc[0],marsloc[1],*])
dust_m_r = transpose(finalcube[marsloc[0]-1,marsloc[1],*])
jupiter_r = transpose(finalcube[jupiterloc[0],jupiterloc[1],*])
dust_j_r = transpose(finalcube[jupiterloc[0],jupiterloc[1]+1,*])
saturn_r = transpose(finalcube[saturnloc[0],saturnloc[1],*])
dust_s_r = transpose(finalcube[saturnloc[0],saturnloc[1]+1,*])
uranus_r = transpose(finalcube[uranusloc[0],uranusloc[1],*])
dust_u_r = transpose(finalcube[uranusloc[0],uranusloc[1]+1,*])
neptune_r = transpose(finalcube[neptuneloc[0],neptuneloc[1],*])
dust_n_r = transpose(finalcube[neptuneloc[0],neptuneloc[1]+1,*])

yn1 = strcompress([-12,-10,-8],/remove_all)
yn2 = ['10!U'+yn1,': !C!A. ','100']

mycolors
set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_1200nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_r,venus_r,/nodata,yrange=[1d-12,1d-4],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75], ytickname = yn2
oplot,wls_r,sun_r/1d6
oplot,wls_r,venus_r,color=5			; yellow
oplot,wls_r,jupiter_r,color=7		; purple
oplot,wls_r,saturn_r,color=8		; orange
oplot,wls_r,earth_r,color=3			; green
oplot,wls_r,mars_r,color=2			; red
oplot,wls_r,uranus_r,color=11		; light blue
oplot,wls_r,neptune_r,color=4		; blue
oplot,wls_r,dust_e_r,line=1
al_legend,['Sun'],line=[0],color=[0],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Venus','Jupiter','Earth','Saturn'], $ 
	line=[0,0,0,0],color=[5,7,3,8],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.86],/norm
al_legend,['Mars','Dust per pix at Earth','Uranus','Neptune'], $ 
	line=[0,1,0,0],color=[2,0,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.86],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/dust_spec_1200nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_r,dust_v_r,/nodata,yrange=[1d-12,1d-8],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wls_r,dust_v_r,color=5		; yellow
oplot,wls_r,dust_e_r,color=3		; green
oplot,wls_r,dust_m_r,color=2		; red
oplot,wls_r,dust_j_r,color=7		; purple
oplot,wls_r,dust_s_r,color=8		; orange
oplot,wls_r,dust_u_r,color=11		; light blue
oplot,wls_r,dust_n_r,color=4		; blue
al_legend,['Dust per pix at Venus','Earth','Mars','Saturn'], $ 
	line=[0,0,0,0],color=[5,3,2,8],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Jupiter','Uranus','Neptune'], $ 
	line=[0,0,0],color=[7,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.61,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_nosun_1200nm.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_r,jupiter_r/sun_r,/nodata,/ylog,yrange=[3d-17,1.2d-7], $ 
	ytitle='F!Dplanet!N / F!DSun!N', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wls_r,(jupiter_r-dust_j_r)/sun_r,color=7				; purple
oplot,wls_r,(venus_r-dust_v_r)/sun_r,color=5		; yellow
oplot,wls_r,(saturn_r-dust_s_r)/sun_r,color=8		; orange
oplot,wls_r,(earth_r-dust_e_r)/sun_r,color=3		; green
oplot,wls_r,(mars_r-dust_m_r)/sun_r,color=2			; red
oplot,wls_r,(uranus_r-dust_u_r)/sun_r,color=11		; light blue
oplot,wls_r,(neptune_r-dust_n_r)/sun_r,color=4		; blue
oplot,wls_r,dust_e_r/sun_r,line=1
al_legend,['Venus','Jupiter','Earth','Saturn'], $ 
	line=[0,0,0,0],color=[5,7,3,8],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Mars','Dust per pix at Earth','Uranus','Neptune'], $ 
	line=[0,1,0,0],color=[2,0,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

;
; Whole bandpass
;

print,'WHOLE BANDPASS'

wave = [wls_b,wls_g,wls_r]
sun = [sun_b,sun_g,sun_r]
venus = [venus_b,venus_g,venus_r]
dust_v = [dust_v_b,dust_v_g,dust_v_r]
earth = [earth_b,earth_g,earth_r]
dust_e = [dust_e_b,dust_e_g,dust_e_r]
mars = [mars_b,mars_g,mars_r]
dust_m = [dust_m_b,dust_m_g,dust_m_r]
jupiter = [jupiter_b,jupiter_g,jupiter_r]
dust_j = [dust_j_b,dust_j_g,dust_j_r]
saturn = [saturn_b,saturn_g,saturn_r]
dust_s = [dust_s_b,dust_s_g,dust_s_r]
uranus = [uranus_b,uranus_g,uranus_r]
dust_u = [dust_u_b,dust_u_g,dust_u_r]
neptune = [neptune_b,neptune_g,neptune_r]
dust_n = [dust_n_b,dust_n_g,dust_n_r]

yn1 = strcompress([-12,-10,-8],/remove_all)
yn2 = ['10!U'+yn1,': !C!A. ','100']

mycolors
set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_all.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wave,venus,/nodata,yrange=[1d-12,1d-4],/ylog,ytitle='Flux density (Jy)', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75], ytickname = yn2
oplot,wave,sun/1d6
oplot,wave,venus,color=5		; yellow
oplot,wave,jupiter,color=7		; purple
oplot,wave,saturn,color=8		; orange
oplot,wave,earth,color=3		; green
oplot,wave,mars,color=2			; red
oplot,wave,uranus,color=11		; light blue
oplot,wave,neptune,color=4		; blue
oplot,wave,dust_e,line=1
al_legend,['Sun'],line=[0],color=[0],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.22,0.92],/norm
al_legend,['Jupiter','Venus','Saturn','Earth'], $ 
	line=[0,0,0,0],color=[7,5,8,3],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.22,0.86],/norm
al_legend,['Uranus','Dust per pix at Earth','Mars','Neptune'], $ 
	line=[0,1,0,0],color=[11,0,2,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.45,0.86],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/dust_spec_nosun_all.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wave,dust_v/sun,/nodata,yrange=[5d-14,2d-10],/ylog,ytitle='F!Ddust!N / F!DSun!N', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wave,dust_v/sun,color=5		; yellow
oplot,wave,dust_e/sun,color=3		; green
oplot,wave,dust_m/sun,color=2		; red
oplot,wave,dust_j/sun,color=7		; purple
oplot,wave,dust_s/sun,color=8		; orange
oplot,wave,dust_u/sun,color=11		; light blue
oplot,wave,dust_n/sun,color=4		; blue
al_legend,['Dust per pix at Venus','Earth','Mars','Saturn'], $ 
	line=[0,0,0,0],color=[5,3,2,8],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.22,0.94],/norm
al_legend,['Jupiter','Uranus','Neptune'], $ 
	line=[0,0,0],color=[7,11,4],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.65,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_spec_nosun_all.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wave,jupiter/sun,/nodata,/ylog,yrange=[3d-17,1d-6], $ 
	ytitle='F!Dplanet!N / F!DSun!N', $
	xtitle='Wavelength ('+!mu+'m)',xmargin=[8.0,2.2],ymargin=[3.25,0.75]
oplot,wave,(jupiter-dust_j)/sun,color=7		; purple
oplot,wave,(venus-dust_v)/sun,color=5		; yellow
oplot,wave,(saturn-dust_s)/sun,color=8		; orange
oplot,wave,(earth-dust_e)/sun,color=3		; green
oplot,wave,(mars-dust_m)/sun,color=2		; red
oplot,wave,(uranus-dust_u)/sun,color=11		; light blue
oplot,wave,(neptune-dust_n)/sun,color=4		; blue
oplot,wave,dust_e/sun,line=1
al_legend,['Jupiter','Venus','Saturn','Earth'], $ 
	line=[0,0,0,0],color=[7,5,8,3],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.18,0.94],/norm
al_legend,['Uranus','Dust per pix at Earth','Neptune','Mars'], $ 
	line=[0,1,0,0],color=[11,0,4,2],/top, $
	/left,linsiz=0.25,charsize=1.7,box=0,pos=[0.41,0.94],/norm
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

a = where(wave ge 0.55)

print,'Earth reflectance at '+strcompress(wave[a[0]],/remove_all)+' = ' + $
	strcompress((earth[a[0]]-dust_e[a[0]])/sun[a[0]],/remove_all)

END