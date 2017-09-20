;+
; NAME:
;   hs_snr
;
; PURPOSE:
;   Create the signal to noise plot for a file
;
; CALLING SEQUENCE:
;  hs_snr, lam, flux, plugmap, sn2=sn2, noplot=noplot, snr=snr
;
; INPUTS:
;    lam = lam [npix,nfib]
;    flux = flux [npix, nfib]
;    plugmap = plugmap [nfib]
;
; OPTIONAL KEYWORDS:
;   
;
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   sn2 = critical signal to noise for the observation
;   snr = output signal to noise for each of the fiber
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------
PRO hs_snr, lam, flux,  plugmap, sn2=sn2, noplot=noplot, snr=snr, $
            nomag=nomag, middle=middle

  
  flux1=flux*0.0
  npix = n_elements(flux(*,0))
  flux1(1:npix-1,*) = flux(0:npix-2,*)
  diff = flux-flux1
  sum = flux+flux1
  temp = sqrt(2) * diff/sum
  nfib = n_elements(lam(0,*))                                                 
  snr = findgen(nfib)
  
  IF NOT keyword_set(middle) THEN BEGIN
     start = 5700
     stop = 6700
  endIF ELSE BEGIN
     start = lam(n_elements(lam(*,0))/2.0-50, 0)
     stop = lam(n_elements(lam(*,0))/2.0+50, 0)
  ENDELSE
  

  for i = 0, nfib -1 do begin
     
     lmin = where(abs(lam[*,i] - start) eq min(abs(lam[*,i]-start)))
     lmax = where(abs(lam[*,i] - stop) eq min(abs(lam[*,i]-stop)))
     
     djs_iterstat, temp[lmin(0):lmax(0), i], sigma=sigma 
     snr(i) = sigma
     
  endfor
  
  snr = 1/snr
  
  IF NOT keyword_set(nomag) THEN begin
     mag = plugmap.rapmag
     
     k = where(mag gt 0 and snr gt 0 AND mag LT 30) 
     plotsym, 0, /fill
     
     
     if NOT keyword_set(noplot) then begin
        
  
        plot, mag(k), snr(k), /ylog, xrange=[25,15], ps=8, xtitle='Magnitude',$
           ytitle='log SN', yrange=[0.1,max(snr(k))]
  
        fit = linfit(mag(k), alog10(snr(k)))
        x = -1 * findgen(1000)/1000*(26-14) + 26
        line = fit(0) + x * fit(1)
        djs_oplot, x, 10^line, color='red', thick=2
        crit =  fit(0) + 21.0 * fit(1)
        xyouts, 20, 0.2, 'S/N at r=21 : ' + $
           string(10.^crit, format='(d8.3)'), $
           charsize=1.5
     endif else begin
        fit = linfit(mag(k), alog10(snr(k)))
        x = -1 * findgen(1000)/1000*(26-14) + 26
        line = fit(0) + x * fit(1)
        crit =  fit(0) + 21.0 * fit(1)
     endelse
     sn2 = crit
  endif else begin
     sn2 = djs_median(snr)
  endelse
  
  
   
   
   
END


  
  
