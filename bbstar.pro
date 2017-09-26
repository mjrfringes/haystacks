;
; Code to generate the blackbody stellar model used by ZODIPIC.
; Written by Aki Roberge [2015-04-09].
;

PRO bbstar, starname, wls, dist, fstar_bb

stellarparam, starname, rstar, lstar, tstar, dist, gk, zk  ; Sun has T_eff = 5777 K

k = 1.38066d-16
c = 2.9979d10			; speed of light in cm/s
h = 6.62608d-27
radsolar = 6.96d10  	; radius of the Sun in cm
lsolar = 3.86d33		; solar luminosity in ergs/sec
rstarcm = radsolar * rstar

; Calculate flux from star from Planck spectrum in ZODIPIC
; Bnu (erg s^-1 cm^-2 ster^-1 Hz^-1)

nu = c / (wls * 1d-4)		; frequency in Hz (wavelengths in cm)
xb = h*nu / (k*tstar)
bnu = xb^3.0/(exp(xb)-1.0)
bnu = bnu*2.0*((k*tstar)^3.0)/((c*h)^2.0)

distcm = dist * 3.085678d18  ; distance to star in cm

fstar_bb = 1e23 * !pi * bnu * ((rstarcm/distcm)^2.0)  ; stellar flux density in Jy

END