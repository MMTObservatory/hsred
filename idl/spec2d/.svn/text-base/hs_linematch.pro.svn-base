lampdefault = '/d0/rcool/data/hectdata/hs_henear.dat'
      lampfilename = lampdefault
  splog, 'Reading lamp file ', lampfilename
   readcol, lampfilename, lampwave, name, format='D, A'
   lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0.0d0}
   lamps = replicate(lamps, N_elements(lampwave))
   lamps.lambda = lampwave
   
   inten = lamps.lambda*0
   pixel = lamps.lambda*0
   for i = 0, n_elements(lamps.lambda) -1 do begin & $
      
      k = where( abs(lamps[i].lambda - lam[*,30]) EQ min(abs(lamps[i].lambda - lam[*,30]))) &$
      if k(0) gt 6 and k(0) lt 4608-6 then inten(i) = max(flux[(k(0)-5):(k(0)+5),odd(30)])   & $
       if k(0) gt 6 and k(0) lt 4608-6 then subset = flux[(k(0)-5):(k(0)+5),odd(30)] & $
       if k(0) gt 6 and k(0) lt 4608-6 then pixsub = lam[(k(0)-5):(k(0)+5),(30)] & $
       if k(0) gt 6 and k(0) lt 4608-6 then j = where(subset eq inten(i)) & $
       if k(0) gt 6 and k(0) lt 4608-6 then pixel(i) = pixsub(j) & $
   endfor
   
   
