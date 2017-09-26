;Approximate value of eccentric anomaly knowing the mean anomaly

pro anomaly,e,mm

  error = 1e-8
  if mm lt !PI then E = mm + e/2 else E = mm - e/2
  ratio = 1
  while abs(ratio) gt error do begin
     ratio = (E - e *sin(E) - mm) / (1 - e * cos(E))
     E = E - ratio
  end


end
