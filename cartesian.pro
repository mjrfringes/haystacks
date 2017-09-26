;correctly solves for the cartesian coordinates in the limit of low
;eccentricity (< a few tenths)
;by Christopher Stark

pro cartesian, GM, a, ecc, inc_cart, longnode, argperi, meananom, x, y, z, vx, vy, vz
  machine_epsilon = 1e-5

  ;first, compute ecc. anom
  E0 = meananom
  E1 = meananom + ecc * sin(E0)
  E2 = meananom + ecc * sin(E1)
  while max(abs(e0-e2)) gt machine_epsilon do begin
     E1 = meananom + ecc * sin(E0)
     E2 = meananom + ecc * sin(E1)
     den = E2 - 2.0*E1 + E0     
     j = where(abs(den) gt machine_epsilon)
     if j[0] ne -1 then begin
        E0[j] = E0[j] - (E1[j]-E0[j])*(E1[j]-E0[j])/den[j]
     endif
     j = where(abs(den) le machine_epsilon)
     if j[0] ne -1 then begin
        E0[j] = E2[j]
        E2[j] = E1[j]
     endif
  endwhile
  cosE = cos(E0)
  sinE = sin(E0)


  ;compute unrotated positions and velocities
  foo = sqrt(1.0 - ecc*ecc)
  meanmotion = sqrt(GM/(a*a*a))
  x = a * (cosE - ecc)
  y = foo * (a * sinE)
  z = fltarr(n_elements(y))
  denom = 1. / (1.0 - ecc * cosE)
  xd = (-a * meanmotion * sinE) * denom
  yd = foo * (a * meanmotion * cosE * denom)
  zd = fltarr(n_elements(yd))

  ;rotate by argument of perihelion in orbit plane
  cosw = cos(argperi)
  sinw = sin(argperi)
  xp = x * cosw - y * sinw
  yp = x * sinw + y * cosw
  zp = z
  xdp = xd * cosw - yd * sinw
  ydp = xd * sinw + yd * cosw
  zdp = zd

  ;rotate by inclination about x axis
  cosi = cos(inc_cart)
  sini = sin(inc_cart)
  x = xp
  y = yp * cosi - zp * sini
  z = yp * sini + zp * cosi
  xd = xdp
  yd = ydp * cosi - zdp * sini
  zd = ydp * sini + zdp * cosi

  ;rotate by longitude of node about z axis
  cosnode = cos(longnode)
  sinnode = sin(longnode)
  xf = x * cosnode - y * sinnode
  yf = x * sinnode + y * cosnode
  zf = z
  vx = xd * cosnode - yd * sinnode
  vy = xd * sinnode + yd * cosnode
  vz = zd
  x = xf
  y = yf
  z = zf
end
