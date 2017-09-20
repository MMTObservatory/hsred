PRO hs_manual_rough, spec, lam=lam, min=min, overwrite=overwrite
  
  if NOT keyword_set(lam) then lam = lindgen(n_elements(spec(*)))

  
  deriv = spec*0
  
  for i = 0, n_elements(spec) -1 do begin
     if i eq 0 then begin
        deriv(i) = 0
     endif else begin
        deriv(i) = (spec(i) - spec(i-1))/(lam(i) - lam(i-1))
     endelse
  endfor
  
  rderiv = deriv*0
  rderiv[0:n_elements(spec)-2] = deriv[1:n_elements(spec)-1]
  rderiv[n_elements(spec)-1] =0
  
  if NOT keyword_set(min) then min = 300
  
  k = where(deriv LT 0 and rderiv GT 0 and spec gt min)
  
  
  wavein = strarr(n_elements(k))
  inten = dblarr(n_elements(k))
  test = ''
  
  for j = 0, n_elements(k) -1 do begin
     i = (n_elements(k) - 1)-j
     djs_plot, lam, spec, xrange=[lam(k(i)) - 300, lam(k(i)) + 300]
     djs_oplot, [lam(k(i)), lam(k(i))], [-1e3,1e12], color='red', thick=2
     splog, 'Please Enter Wavelegnth (X for skip)'
     splog, 'Estimated at ' + string(lam(k(i)))
     read, test, prompt='Value:    '
     if test eq 'STOP' then stop
     wavein(i) = test
     inten(i) = spec(k(i))
     
  endfor
  
  openw, 1, 'wavcal.dat'
  u = where(wavein ne 'x' and wavein ne 'X')
  wavein = wavein(u)
  inten = inten(u)
  k = k(u)
  
  for i = 0, n_elements(k) -1 do begin & $
     
     printf, 1, wavein(i) + string(inten(i)) + string(4607-k(i)) & $
  endfor
  close, 1
  
END
