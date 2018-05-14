PRO Haystacks_loop_split

; This function runs all the haystacks cubes in parallel
; no arguments are needed, just run the IDL program
; it will output a file (fname) with the list of all cubes created


incs = [0,60,90]
systype = ['modern','archean']

; name of the output file that will contain all of the cubes's names
; they can then be easily grabbed for post-processing
fname = 'haystacks_loop_names.txt'
openw,4,fname
close,4

; cut the full range of wavelengths into chunks of 20% bandwidth
minlam = 0.3
maxlam = 2.5
bw = 0.2
nlam = fix(ALOG(maxlam/minlam)/bw)+1
print,'Number of wavelength chunks:'+string(nlam)
end_vals = findgen(nlam)*ALOG(maxlam/minlam)/fix(nlam-1) + ALOG(minlam)
wls = EXP(end_vals)
print,'Fractional bandwidth for each cube'+string(100*2.*(wls[1]-wls[0])/(wls[1]+wls[0]))+'%'

print,'Total number of cubes'+string(n_elements(incs)*n_elements(systype)*(nlam-1))

; run the loop on all the cases
FOR k=0, nlam-2 DO BEGIN
	start = systime(/sec)
	lowlam = wls[k]
	highlam = wls[k+1]
;	FOR i=0,2 DO BEGIN
;		FOR j=0,1 DO BEGIN
;			print,'Working on '+systype(j)+ ' Solar System, band '+string(band)+' at inclination '+string(incs(i))
;			;print,'Working on modern  Solar System, band '+string(band)+' at inclination '+string(incs(i))
;			
;			haystacks,incs(i),systype(j),lowlam,highlam,cubename
;			openw,lun,fname,/append
;			printf,lun,cubename
;			close,lun
;		ENDFOR
;	ENDFOR
	print,string(lowlam)+', '+string(highlam)
	b1=obj_new("IDL_IDLBridge", output='')  ; Initiate IDL bridge/child process
    b1->setvar, 'var1', incs(0)
    b1->setvar, 'var2', systype(1)
    b1->setvar, 'var3', lowlam
    b1->setvar, 'var4', highlam
	b1->execute,'outvar=haystacks(var1,var2,var3,var4)',/nowait
	print,'Process launched'

	b2=obj_new("IDL_IDLBridge", output='')  ; Initiate IDL bridge/child process
    b2->setvar, 'var1', incs(1)
    b2->setvar, 'var2', systype(1)
    b2->setvar, 'var3', lowlam
    b2->setvar, 'var4', highlam

	b2->execute,'outvar=haystacks(var1,var2,var3,var4)',/nowait
	print,'Process launched'
	b3=obj_new("IDL_IDLBridge", output='')  ; Initiate IDL bridge/child process
    b3->setvar, 'var1', incs(2)
    b3->setvar, 'var2', systype(1)
    b3->setvar, 'var3', lowlam
    b3->setvar, 'var4', highlam
	b3->execute,'outvar=haystacks(var1,var2,var3,var4)',/nowait

	print,'Process launched'                  
	b4=obj_new("IDL_IDLBridge", output='')  ; Initiate IDL bridge/child process
    b4->setvar, 'var1', incs(0)
    b4->setvar, 'var2', systype(0)
    b4->setvar, 'var3', lowlam
    b4->setvar, 'var4', highlam
	b4->execute,'outvar=haystacks(var1,var2,var3,var4)',/nowait

	print,'Process launched'
	b5=obj_new("IDL_IDLBridge", output='')  ; Initiate IDL bridge/child process
    b5->setvar, 'var1', incs(1)
    b5->setvar, 'var2', systype(0)
    b5->setvar, 'var3', lowlam
    b5->setvar, 'var4', highlam
	b5->execute,'outvar=haystacks(var1,var2,var3,var4)',/nowait

	print,'Process launched'
	b6=obj_new("IDL_IDLBridge", output='')  ; Initiate IDL bridge/child process
    b6->setvar, 'var1', incs(2)
    b6->setvar, 'var2', systype(0)
    b6->setvar, 'var3', lowlam
    b6->setvar, 'var4', highlam
	b6->execute,'outvar=haystacks(var1,var2,var3,var4)',/nowait
	print,'Process launched'

	WHILE (b1->status() eq 1) || (b2->status() eq 1) || (b3->status() eq 1) || (b4->status() eq 1) || (b5->status() eq 1) || (b6->status() eq 1) DO wait,10

	openw,4,fname,/append
	name = b1->GetVar('outvar')
	printf,4,name
	name = b2->getvar('outvar')
	printf,4,name
	name = b3->getvar('outvar')
	printf,4,name
	name = b4->getvar('outvar')
	printf,4,name
	name = b5->getvar('outvar')
	printf,4,name
	name = b6->getvar('outvar')
	printf,4,name
	close,4
	tend = systime(/sec)
	dt = tend-start
	print,'Time elapsed: '+strcompress(dt)+' seconds'

	
ENDFOR

END
