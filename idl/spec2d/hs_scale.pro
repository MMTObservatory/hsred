;+
; NAME:
;   hs_scale
;
; PURPOSE:
;   Find the scales for the fibers based on sky lines
;
; CALLING SEQUENCE:
;  hs_scale, lam, flux, scale=, 
;
; INPUTS:
;    lam - Input wavelenth [npix, nfib]
;    flux - input flux [npix, nfib]
;
; OPTIONAL KEYWORDS:
;  
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;   scale = scale for each of the fibers
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA

;------------------------------------------------------------------------------
Pro hs_scale,  lam, flux, sigma=sigma, scale=scale, norm=norm
  
  if not KEYWORD_SET(sigma) THEN sigma=2.0
  
  ;;Lines to fit
  lwave = [4046.56, 5460.73, 5577.345, 6300.32, $
           6864.0, 8885.83, 8919.61,7247.88, 7369.2,$
           7402.12, 7713.29, $
           7750.88, 7794.48, 7822.06, 7853.2, $
           7914.59, 8023.79, 8344.30, 8398.57, $
           8430.60,  8465.29];, 8885.83, 8919.61, $
           ;8838.5997,8777.4183,8793.2410,8887.1229, $
          ; 8919.8233, 8943.0301,8957.7980, 9002.1018]

  
  boxarea = 5*sigma
  hpix= fix(boxarea-1)/2
  intarea= 2*hpix+1
  boxcar = replicate(1,intarea)
  if boxarea - intarea GT 0 then $
     boxcar= [0.5 *(boxarea-intarea), boxcar, 0.5*(boxarea-intarea)]
  
  nfib = n_elements(flux(0,*))
  norm = flux * 0.0
  
  for i = 0, nfib -1 do begin
     norm[*,i] = median(flux[*,i],100)
  endfor
  
  width = fltarr(n_elements(lwave), nfib)
  
  for i = 0, nfib -1 do begin
     for j = 0,n_elements(lwave) -1 do begin
        
        k = where(abs(lam[*,i]-lwave(j)) eq min(abs(lam[*,i]-lwave(j))))
        min = k(0) -10
        max = k(0) +10
        ew = convol(flux[min:max,i]-norm[min:max,i],boxcar, /edge_truncate)
        ew = ew(11)
        width(j,i) = max(flux[min:max,i]/norm[min:max,i])
     endfor
  endfor
  
  ew = fltarr(n_elements(lwave))
  for i = 0, n_elements(lwave)-1 do begin
     ew(i) = median(width(i,*))
  endfor
  
  scale = fltarr(n_elements(lwave), nfib)
  wave = scale
  for i = 0, nfib-1 do begin
     scale(*,i) = djs_median(width(*,i)/ew)
  endfor
  
  lwave1 = findgen(n_elements(lwave), nfib)
  for i = 0, nfib-1 do begin
     lwave1(*,i) = lwave
  endfor
  
  
  scale1 = lam*0.0
  for i = 0, nfib-1 do scale1(*,i) = replicate(djs_median(scale(*,i)), $
                                               n_elements(scale1(*,i)))
  scale=scale1
  
END


        
        
        
     
