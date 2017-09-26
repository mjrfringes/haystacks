; Project directory. Write out full path name. DUSTMAP cannot parse ~
  proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'

band = 1
sres = 150.

wrange = [0.3,2.5]

width = wrange[1] - wrange[0]
centwave = width/2. + wrange[0]
wbin = centwave/sres
steps = round(width/wbin) + 1
wls = findgen(steps)*wbin + wrange[0]    ; Array of wavelengths in microns
nlam = n_elements(wls)
print,'Min/Max wavelenth: '+strcompress(string(min(wls),format='(F6.3)'),/remove_all) $ 
 	+' - '+strcompress(string(max(wls),format='(F6.3)'),/remove_all)+' um'
print,'Number of wavelengths: '+strtrim(nlam,2)

bbstar, 'Sun', wls, ref_dist, fstar_bb

spec_dir = proj_dir + 'Spectra_Final/'
sun_spec_file = 'SunKuruczSpectrum.txt'
readcol,spec_dir+sun_spec_file,swl,sfnu ; wl in um, sfnu in Jy
sunspec = interpol(sfnu,swl,wls)        ; Interpolate between wavelengths

mycolors
set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/star_compare.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=9,ysize=7,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls,fstar_bb,ytitle='Flux density (Jy)',xtitle='Wavelength ('+!mu+'m)', $
	xmargin=[6.0,1.8],ymargin=[3.25,0.5]
oplot,wls,sunspec,color=2
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

a = where(wls ge 0.551)

print,'Flux at '+strcompress(wls[a[0]],/remove_all)+' = '+ $ 
	strcompress(sunspec(a[0]),/remove_all)+' Jy'

;
; Patching over the gaps in the observed NIR giant planet spectra.
;

spec_dir = proj_dir + 'Spectra_Final/'
readcol,spec_dir+'Jupiter_geo_albedo_old.txt',wls_j,alb_j
readcol,spec_dir+'Saturn_geo_albedo_old.txt',wls_s,alb_s
readcol,spec_dir+'Uranus_geo_albedo_old.txt',wls_u,alb_u
readcol,spec_dir+'Neptune_geo_albedo_old.txt',wls_n,alb_n

; Burrows jup_4.0AU.txt for Jupiter
; jup_10AU.txt for Saturn
; jup_15AU.txt for Uranus and Neptune

readcol,'~/Other_Projects/Haystacks/Planet_Spectra/Burrows/jup_15AU.txt',wave,f_j,ref_j

rangej = [1.8127,1.8732]
ranges = [1.8180,1.8805]
rangeu = [1.832, 1.9148]
rangen = [1.8183, 1.8809]

scl = 2.07d8

b = where((wave ge rangen[0]) and (wave le rangen[1]))
chop_w = wave[b] & chop = ref_j[b]

plot,wls_n,alb_n,xrange=[1.73,1.95],yrange=[-0.0005,0.004]
oplot,wave,ref_j*scl,color=2
oplot,wave[b],chop*scl,color=3

c = where((wls_n ge rangen[0]) and (wls_n le rangen[1]))
new_chop = interpol(chop*scl,chop_w,wls_n[c])

oplot,wls_n[c],new_chop,color=8

new = alb_n
new[c] = new_chop

plot,wls_n,new,xrange=[1.8,1.9],yrange=[-0.0005,0.004]

readcol,spec_dir+'Jupiter_geo_albedo.txt',newwls_j,newalb_j
readcol,spec_dir+'Saturn_geo_albedo.txt',newwls_s,newalb_s
readcol,spec_dir+'Uranus_geo_albedo.txt',newwls_u,newalb_u
readcol,spec_dir+'Neptune_geo_albedo.txt',newwls_n,newalb_n

oplot,newwls_n,newalb_n,color=3,line=2

!p.multi = [0,1,4]

mycolors
set_plot,'ps'
device,file='~/Other_Projects/Haystacks/Notes/Mine/Figures/planet_patch.eps', $
	/encapsulated,/portrait,/helvetica,/isolatin1,xsize=8.5,ysize=11,/inches,/color,/bold
!p.thick=4 & !x.thick=3 & !y.thick=3 & !p.font=0 & !p.charsize=2.2
plot,wls_j,alb_j,xrange=[1.75,1.95],yrange=[-0.004,0.04],ytitle=' ', $
	xmargin=[8,2.5], ymargin = [3.25,0.75], xtitle='Wavelength ('+!mu+'m)'
oplot,newwls_j,newalb_j,color=2
xyouts,0.13,0.95,'Jupiter',/norm,color=2,charsize=2.0
xyouts,rangej[0],0.02,string(rangej[0],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
xyouts,rangej[1],0.01,string(rangej[1],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
;
plot,wls_s,alb_s,xrange=[1.75,1.95],yrange=[-0.007,0.07],ytitle=' ', $
	xmargin=[8,2.5], ymargin = [3.25,0.75], xtitle='Wavelength ('+!mu+'m)'
oplot,newwls_s,newalb_s,color=8
xyouts,0.13,0.7,'Saturn',/norm,color=8,charsize=2.0
xyouts,ranges[0],0.04,string(ranges[0],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
xyouts,ranges[1],0.02,string(ranges[1],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
;
plot,wls_u,alb_u,xrange=[1.75,1.95],yrange=[-0.0002,0.0022],ytitle=' ', $
	xmargin=[8,2.5], ymargin = [3.25,0.75], xtitle='Wavelength ('+!mu+'m)'
oplot,newwls_u,newalb_u,color=3
xyouts,0.13,0.45,'Uranus',/norm,color=3,charsize=2.0
xyouts,rangeu[0],0.0015,string(rangeu[0],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
xyouts,rangeu[1],0.001,string(rangeu[1],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
;
plot,wls_n,alb_n,xrange=[1.75,1.95],yrange=[-0.0002,0.002],ytitle=' ', $
	xmargin=[8,2.5], ymargin = [3.25,0.75], xtitle='Wavelength ('+!mu+'m)'
oplot,newwls_n,newalb_n,color=11
xyouts,0.13,0.2,'Neptune',/norm,color=11,charsize=2.0
xyouts,rangen[0],0.00015,string(rangen[0],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
xyouts,rangen[1],0.00015,string(rangen[1],format='(F5.3)')+' '+!mu+'m', $ 
	align=0.5,charsize=1.5
;
!p.charsize=1.5
datetag
device,/close
set_plot,'x'

!p.multi=0

END