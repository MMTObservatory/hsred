PRO hs_skyscale, wave, flux, plugmap, scale=scale, supersky=supersky, $
                 newwave=newwave, debug=debug
  
  
  k = where(plugmap.objtype eq 'SKY')
  newwave = 3900 + findgen(4200)*1.2
  newflux = dblarr(n_elements(newwave), n_elements(flux(0,*)))
  nfib = n_elements(flux(0,*))
  for i = 0, nfib-1 do newflux(*,i) = interpol(flux(*,i), wave(*,i), newwave)
  
  norm = newflux*0.0
  for i = 0,nfib-1 do norm[*,i] = median(newflux(*,i), 100)
  testspec = newflux/norm
  
  sky = djs_median(newflux(*,k), 2)
  snorm = median(sky, 100)
  ss = sky/snorm
  
  scaleall = flux*0.0
  junk = min(abs(newwave-7200), j1)
  junk = min(abs(newwave-7300), j2)
  
  for i = 0, nfib-1 do begin
     scaleval = findgen(40)/40*2+0.1
     sigma = scaleval*0.0
     for j = 0, n_elements(scaleval)-1 do begin
        test = newflux(*,i)/scaleval(j)-sky
        test = test-djs_median(test(j1:j2))
        djs_iterstat, test(j1:j2), mask=mask, sigrej=5
        test = test(j1:j2)
        k = where(mask eq 1, ct)
        sigma(j) = sqrt(total(test(k)^2))
     endfor
     
     x =  findgen(1000)/1000.*2.+0.1  
     plotsym, 0, /fill
     
     sigma1 = interpol(sigma, scaleval, x)
     junk = min(sigma1, m)
     if m eq 0 or m eq 999 then junk = min(abs(x-1), m)
     
     
        
     
     scaleall(*,i) = replicate(x(m), n_elements(flux(*,i)))
     
     if keyword_set(debug) then begin
        !P.multi=[0,1,2]
        !P.position=[0.1, 0.1, 0.3, 0.9]
        plot, scaleval, sigma, ps=8
        djs_oplot, x(m)*[1,1], [-1e4, 1e4], color='red'
        !P.position=[0.35, 0.1, 0.9, 0.9]
        plot, newwave, smooth(newflux(*,i)/scaleall(*,i)-sky,11),$
           xrange=[4000,8500], /xstyle, /ystyle
        djs_oplot, newwave, $
           djs_median(newflux(*,i)/scaleall(*,i)-sky, width=50, $
                      boundary='reflect'), color='blue', thick=2
        wait, 0.5
     endif
     
     
  endfor
  
  scale=scaleall
  ;;Check for ungodly answers
  k = where(scale(0,*) GT 1.1, ct)
  FOR i = 0, ct-1 DO $
     scale(*,k(i)) = replicate(djs_median(scaleall), n_elements(flux(*,i)))
 
  return
  
  
END
