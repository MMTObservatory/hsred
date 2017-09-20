PRO hs_snr_ivar, lam, flux, ivar, plugmap
  
  
  temp = flux * sqrt(abs(ivar))
  
  snr = findgen(300)
  
  for i = 0, 299 do begin
     
     k = where(abs(lam[*,i] - 6500.00) eq min(abs(lam[*,i]-6500.00)))
     lmin = k-200
     lmax = k+200
     
     djs_iterstat, temp[lmin:lmax, i], median=med
     snr(i) = med
     
  endfor
  
  mag = plugmap.rapmag
  
  k = where(mag gt 0 and snr gt 0)
  plotsym, 0, /fill
  
  plot, mag(k), snr(k), /ylog, xrange=[25,15], ps=8
  
  fit = linfit(mag(k), alog10(snr(k)))

  x = -1 * findgen(1000)/1000*(26-14) + 26

  line = fit(0) + x * fit(1)

  djs_oplot, x, 10^line, color='red', thick=2

  crit =  fit(0) + 21.0 * fit(1)
  
  
  xyouts, 21, 5, 'S/N at r=21 : ' + string(10.^crit, format='(d8.3)'), $
     charsize=1.5

  
END


  
  
