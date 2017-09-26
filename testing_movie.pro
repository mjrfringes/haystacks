;
; Write Haystacks cube to video file.
;

; Project directory

proj_dir = '~/Other_Projects/Haystacks/Current_Code/Solar_System/'
folder = "zodi1inc0/"
dir1 = proj_dir + folder + 'NoSun_NoLocal/'

GOTO, ALL_INNER

; BLUE CUBE

print,'BLUE CUBE...'
file_b = 'cube_zodi1inc0dist'+strtrim(fix(10),2)+'_0.3um'
fits_read, dir1+file_b+'.fits',finalcube_b,head
fits_read, proj_dir + folder + 'wavelengths_0.3um.fits',wave_b,headw

sz_b = size(finalcube_b)

loadct,3
;window,2,xsize=sz_b[1],ysize=sz_b[2]
tmp = finalcube_b[*,*,0]
index = round(sz_b[3] / 2.)
tmp2 = finalcube_b[*,*,index]

;cgImage, tmp2, stretch=1, maxvalue=max(tmp),minvalue=min(tmp)

minv = min(tmp2)
maxv = median(tmp2)*30.

!p.font=0 & !p.charsize=2.0
set_plot,'ps'
device,file='tmp.ps',/helvetica
LABELS_b = cgGreek('lambda')+' = '+string(wave_b,format='(F4.2)')+' '+!mu+'m'
device,/close
set_plot,'x'
File_Delete, 'tmp.ps'

cgPS_Open,'movie.ps',/helvetica,/quiet
tmp3 = cgImgScl(tmp2, stretch=1, minvalue=minv, maxvalue=maxv)
cgDisplay,sz_b[1]*0.5,sz_b[2]*0.5
cgImage,tmp3
cgText,0.08,0.07,labels[index],/Normal,color=255
cgPS_Close,/png,width=sz_b[1],/Delete_PS

make_video, 'movie_300nm.mp4', finalcube_b, LABELS=labels_b

; Inner region

cent = sz_b[1]/2. - 0.5
res = 0.1		; AU per pixel
x1 = cent - (10/res)
x2 = cent + (10/res)
tmp_inner_b = finalcube_b[x1:x2,x1:x2,*]
sm_sz_b = size(tmp_inner_b)
inner_b = rebin(tmp_inner_b,sm_sz_b[1] * 4, sm_sz_b[2] * 4, sm_sz_b[3], /sample)

make_video, 'movie_300nm_inner.mp4', inner_b, LABELS=labels_b

; GREEN CUBE

print,'GREEN CUBE...'
file_g = 'cube_zodi1inc0dist'+strtrim(fix(10),2)+'_0.6um'
fits_read, dir1+file_g+'.fits',finalcube_g,head
fits_read, proj_dir + folder + 'wavelengths_0.6um.fits',wave_g,headw

sz_g = size(finalcube_g)

!p.font=0 & !p.charsize=2.0
set_plot,'ps'
device,file='tmp.ps',/helvetica
LABELS_g = cgGreek('lambda')+' = '+string(wave_g,format='(F4.2)')+' '+!mu+'m'
device,/close
set_plot,'x'
File_Delete, 'tmp.ps'

make_video, 'movie_600nm.mp4', finalcube_g, labels=labels_g

GREEN:

; Inner region

cent = sz_g[1]/2. - 0.5
res = 0.1		; AU per pixel
x1 = cent - (10/res)
x2 = cent + (10/res)
tmp_inner_g = finalcube_g[x1:x2,x1:x2,*]
sm_sz_g = size(tmp_inner_g)
inner_g = rebin(tmp_inner_g,sm_sz_g[1] * 4, sm_sz_g[2] * 4, sm_sz_g[3], /sample)

make_video, 'movie_600nm_inner.mp4', inner_g, labels=labels_g

; RED CUBE

print,'RED CUBE...'
file_r = 'cube_zodi1inc0dist'+strtrim(fix(10),2)+'_1.2um'
fits_read, dir1+file_r+'.fits',finalcube_r,head
fits_read, proj_dir + folder + 'wavelengths_1.2um.fits',wave_r,headw

!p.font=0 & !p.charsize=2.0
set_plot,'ps'
device,file='tmp.ps',/helvetica
LABELS_r = cgGreek('lambda')+' = '+string(wave_r,format='(F4.2)')+' '+!mu+'m'
device,/close
set_plot,'x'
File_Delete, 'tmp.ps'

sz_r = size(finalcube_r)

make_video, 'movie_1200nm.mp4', finalcube_r, labels=labels_r

; Inner region

cent = sz_r[1]/2. - 0.5
res = 0.1		; AU per pixel
x1 = cent - (10/res)
x2 = cent + (10/res)
tmp_inner_r = finalcube_r[x1:x2,x1:x2,*]
sm_sz_r = size(tmp_inner_r)
inner_r = rebin(tmp_inner_r,sm_sz_r[1] * 4, sm_sz_r[2] * 4, sm_sz_r[3], /sample)

make_video, 'movie_1200nm_inner.mp4', inner_r, labels=labels_r

; TOTAL CUBE

ALL:

print,'TOTAL CUBE...'
slices = total( [sz_b[3],sz_g[3],sz_r[3]] )
total = fltarr(sz_b[1],sz_b[2],slices)
total[*,*,0:sz_b[3]-1] = finalcube_b
total[*,*,sz_b[3]:sz_b[3]+sz_g[3]-1] = finalcube_g
total[*,*,sz_b[3]+sz_g[3]:sz_b[3]+sz_g[3]+sz_r[3]-1] = finalcube_r

labels = [labels_b,labels_g,labels_r]

sz = size(total)

;window,2,xsize=sz[1],ysize=sz[2]
;tmp3 = cgImgScl(total[*,*,307], stretch=1, minvalue=minv, maxvalue=maxv)
;tv,tmp3

make_video, 'movie_all.mp4', total, labels=labels

; Inner region

ALL_INNER:

cent = sz[1]/2. - 0.5
res = 0.1		; AU per pixel
x1 = cent - (10/res)
x2 = cent + (10/res)
tmp_inner_t = total[x1:x2,x1:x2,*]
sm_sz_t = size(tmp_inner_t)
inner_t = rebin(tmp_inner_t,sm_sz_t[1] * 5, sm_sz_t[2] * 5, sm_sz_t[3], /sample)

sz_in = size(inner_t)

a = where(inner_t eq max(inner_t))
maxind = get_indices(a[0],[sz_in[1],sz_in[2], sz_in[3]])
minv = min(inner_t[*,*,maxind[2]])
maxv = median(inner_t[*,*,maxind[2]])*250.

window,2,xsize=sz_in[1],ysize=sz_in[2]
tmp3 = cgImgScl(inner_t[*,*,307], stretch=3, minv=minv, maxv=maxv, gamma=0.25)
tv,tmp3

make_video, 'movie_all_inner.mp4', inner_t, gamma=0.25, minv=minv, $ 
	maxv=maxv, labels=labels

END