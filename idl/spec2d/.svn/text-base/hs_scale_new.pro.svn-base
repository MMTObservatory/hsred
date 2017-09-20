FUNCTION gfit, X, P
  RETURN, P(0) + GAUSS1(X, P(1:3))
END

Pro hs_scale_new,  lam, flux, plugmap, sigma=sigma, scale=scale, $
                   debug=debug
  
  if not KEYWORD_SET(sigma) THEN sigma=2.0
  lbegin = [5450, 7335.2, 7987.0, 8390.0]
  lend = [5470, 7346, 7999, 8407]
  lwave = [ 4358, 4981.8, 5085, 5460.73, 5577.345, 6300.32, $
            6864.0, 8885.83, 8919.61, 7369.2,$
            7402.12, 7713.29, $
            7750.88, 7794.48, 7822.06, 7853.2, $
            7914.59, 8023.79, 8344.30, 8398.57, $
            8430.60,  8465.29, 8885.83, 8919.61, $
            8777.4183,8793.2410,8887.1229, $
            8919.8233, 8943.0301,8957.7980, 9002.1018, $
            6363, 5889]
  lbegin=lwave
  lend = lwave
  nfib = n_elements(flux(0,*))
  norm = flux * 0.0
  skynorm = norm*0.0
  
  for i = 0, nfib -1 do begin
     norm[*,i] = median(flux[*,i],100)
  endfor
  
  width = fltarr(n_elements(lbegin), nfib)
  dist = sqrt(plugmap.xfocal^2+plugmap.yfocal^2)
  c = min(abs(dist-250), crit)
  for i = 0, nfib -1 do begin
   
     for j = 0,n_elements(lbegin) -1 do begin
        
        junk = min(abs(lam(*,i)-lwave(j)), t)
        t1 = t-10
        t2 = t+10
        
        x = lam(t1:t2, i)
        y = flux(t1:t2,i)
        
        cont = y*0.0 + djs_median([y[0:5],y[16:19]])
        width(j,i) = int_tabulated(x, y)/int_tabulated(x, cont)
        
     endfor
  endfor
  
  scale = lam*0.0
  for i = 0, n_elements(lam(0,*))-1 do begin
     coeff = $
        robust_poly_fit( lwave, width(*,i)/width(*,crit), 3, yfit)
     
     scale(*,i) = 0
     for j = 0,n_elements(coeff)-1 do $
        scale(*,i) = scale(*,i) + coeff(j)*lam(*,i)^(j)
  endfor
  
  newwave = 3900 + findgen(4200)*1.2
  newflux = dblarr(n_elements(newwave), nfib)
  for i = 0, nfib-1 do $
     newflux(*,i) = interpol(flux(*,i)/scale(*,i), lam(*,i), newwave)
  
  if keyword_set(debug) then begin
     for i = 0, 149 do begin
        djs_plot, newwave, smooth(newflux(*,i)-sky,11), $
           xrange=[4000,8500], /xstyle, /ystyle
        wait, 0.5
     endfor
     
  endif
  return
END


        
        
        
     
