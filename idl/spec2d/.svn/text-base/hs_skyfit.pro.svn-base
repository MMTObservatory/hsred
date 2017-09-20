;Fit for fiber-to-fiber transmission correction by scaling the sky
;line height in each fiber to match that in the average 'supersky'
; based on SKYFIT module of SPECROAD pipeline, by J. Mink
PRO hs_skyfit, wave, flux, scale=scale, skyset=skyset
  
             sky=bspline_valu(wave,skyset)
             tmpsmooth=flux*0.0
             skysmooth=flux*0.0
             nsmooth=n_elements(flux[*,0])/50
             for q=0, 149 do tmpsmooth[*,q]=djs_median(flux[*,q],width=nsmooth, boundary='reflect') 
             fluxscale=flux-tmpsmooth
             for q=0, 149 do skysmooth[*,q]=djs_median(sky[*,q],width=nsmooth, boundary='reflect') 
             skyscale=sky-skysmooth

  
  nfib = n_elements(flux[0,*])
   
  
  scaleall = flux*0.0
;  scaleall2 = scaleall
  
  for i = 0, nfib-1 do begin
     scaleval = findgen(160)/160*2+0.1
     sigma = scaleval*0.0
;     sigma2 = scaleval*0.0
     for j = 0, n_elements(scaleval)-1 do begin
        test = fluxscale[*,i]/scaleval[j]-skyscale[*,i]
;        rval=r_correlate(skyscale[*,i],test)
;        sigma[j] = rval[0]
        rval=correlate(skyscale[*,i],test)
        sigma[j] = rval
     endfor
     
     x =  findgen(4000)/4000.*2.+0.1  
     
     sigma1 = interpol(sigma, scaleval, x)
     junk = min(sigma1^2.0, m)
   ;  sigma1a = interpol(sigma2, scaleval, x)
   ;  junk2 = min(sigma1a^2.0, m2)
;     print, m, m2, x[m], x[m2], junk, junk2
;    if reform(x[m]) lt 0.9 then stop
;     if m eq 0 or m eq 999 then junk = min(abs(x-1), m)
        
     
     scaleall[*,i] = replicate(x[m], n_elements(flux[*,i]))    
   ;  scaleall2[*,i] = replicate(x[m2], n_elements(flux[*,i]))    
     ;stop    
  endfor
;stop  
  scale=scaleall
  ;;Check for ungodly answers
  k = where((scale[0,*] GT 1.7) or (scale[0,*] lt 0.3), ct)
  FOR i = 0, ct-1 DO $
     scale[*,k[i]] = replicate(djs_median(scaleall), n_elements(flux[*,i]))
 
  return
  
  
END
