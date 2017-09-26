PRO Haystacks_loop

incs = [0,60,90]
systype = ['modern','archean']

; name of the output file that will contain all of the cubes's names
; they can then be easily grabbed for post-processing
fname = 'haystacks_loop_names.txt'
openw,lun,fname,/get_lun
close,lun

; cut the full range of wavelengths into chunks of 20% bandwidth
minlam = 0.3
maxlam = 2.5
bw = 0.2
nlam = fix(ALOG(maxlam/minlam)/bw)+1
print,'Number of wavelength chunks:'+string(nlam)
end_vals = findgen(nlam)*ALOG(maxlam/minlam)/fix(nlam-1) + ALOG(minlam)
wls = EXP(end_vals)
print,'Fractional bandwidth for each cube'+string(100*2.*(wls[1]-wls[0])/(wls[1]+wls[0]))+'%'

print,'Number of total cubes'+string(n_elements(incs)*n_elements(systype)*(nlam-1))

; run the loop on all the cases
FOR k=0, nlam-2 DO BEGIN
	lowlam = wls[k]
	highlam = wls[k+1]
	FOR i=0,2 DO BEGIN
;		FOR j=0,1 DO BEGIN
;			print,'Working on '+systype(j)+ ' Solar System, band '+string(band)+' at inclination '+string(incs(i))
			;print,'Working on modern  Solar System, band '+string(band)+' at inclination '+string(incs(i))
			haystacks,incs(i),systype(j),lowlam,highlam,cubename
			openw,lun,fname,/append
			printf,lun,cubename
			close,lun
;		ENDFOR
	ENDFOR
ENDFOR

END
