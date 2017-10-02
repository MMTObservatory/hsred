;Fit for fiber-to-fiber transmission correction by scaling the sky
;line height in each fiber to match that in the average 'supersky'
; based on SKYFIT module of SPECROAD pipeline, by J. Mink
PRO hs_skyfit, wave, flux, scale=scale, skyset=skyset, skyflux=skyflux
  
             nfib=(size(flux, /dimen))[1]
             npix=(size(flux, /dimen))[0]
             if not keyword_set(skyflux) then sky=bspline_valu(wave,skyset) else sky=skyflux
             scaleall = flux*0.0
            
             tmpsmooth=flux*0.0
             skysmooth=flux*0.0
             nsmooth=n_elements(flux[*,0])/50
             for q=0, nfib-1 do tmpsmooth[*,q]=djs_median(flux[*,q],width=nsmooth, boundary='reflect') 
             fluxscale=flux-tmpsmooth
             for q=0, nfib-1 do skysmooth[*,q]=djs_median(sky[*,q],width=nsmooth, boundary='reflect') 
             skyscale=sky-skysmooth
   
  
;  nfib = n_elements(flux[0,*])
   
  
  
;  scaleall2 = scaleall
  bestvals=fltarr(nfib)

  for i = 0, nfib-1 do begin
     scaleval = findgen(240)/240*3+0.3
     sigma = scaleval*0.0
     ;lets ignore wavelengths longer than 9010 or shorter than 3900
     ww=where((wave[*,i] ge 3900.) and (wave[*,i] lt 8800.),/null, cnt) 
     if cnt eq 0 then begin
        print, 'no wavelengths in range for skyfit'
        scale=scaleall
        scale[*,*]=1.0
        return
     endif
     for j = 0, n_elements(scaleval)-1 do begin
        test = fluxscale[ww,i]/scaleval[j]-skyscale[ww,i]
;        rval=r_correlate(skyscale[*,i],test)
;        sigma[j] = rval[0]
        rval=correlate(skyscale[ww,i],test)
        sigma[j] = rval
     endfor
     ;passing in a sky spectrum as both 'flux' and 'sky' results in
     ;this bug, as the sky is precisely subtracted and we correlate against
     ;an array of zeros. Setting this nan to 0.0 is what we want
     bad=where(finite(sigma,/nan),cnt)
     if cnt ge 1 then sigma[bad]=0.0
     x =  findgen(6000)/6000.*3.+0.3  
     
     sigma1 = interpol(sigma, scaleval, x)
     bestvals[i] = min(sigma1^2.0, m)
 ;    if bestvals[i] gt 0.1 then stop     
   ;  sigma1a = interpol(sigma2, scaleval, x)
   ;  junk2 = min(sigma1a^2.0, m2)
;     print, m, m2, x[m], x[m2], junk, junk2
;    if reform(x[m]) lt 0.9 then stop
;     if m eq 0 or m eq 999 then junk = min(abs(x-1), m)
        
     
     ;if we didn't find a local minimum, just skip this fiber
     if (m ne 0) and (m ne n_elements(x)-1) then scaleall[*,i] = replicate(x[m], npix) $ 
        else   scaleall[*,i] = replicate(1.0, npix)
   ;  scaleall2[*,i] = replicate(x[m2], n_elements(flux[*,i]))    
     ;stop    
  endfor
;stop  
  scale=scaleall
    ;check for fibers where the best correlation coefs were high,
    ;indicating we never passed through a minimum at zero, and discard these corrections
  k = where(bestvals gt 0.05, ct)
  FOR i = 0, ct-1 DO $
     scale[*,k[i]] = replicate(djs_median(scaleall), npix)
 
  return
  
  
END
