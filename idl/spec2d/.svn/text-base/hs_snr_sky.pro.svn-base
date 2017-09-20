;+
; NAME:
;   hs_snr_sky
;
; PURPOSE:
;   Create the signal to noise plot for a file
;
; CALLING SEQUENCE:
;  hs_snr_sky, lam, flux, plugmap
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
;  
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




PRO hs_snr_sky, lam, flux,  plugmap

  
  flux1=flux*0.0
  npix = n_elements(flux(*,0))
  flux1(1:npix-1,*) = flux(0:npix-2,*)
   
   
  diff = flux-flux1
  sum = flux+flux1
  
  temp = sqrt(2) * diff/sum
  
    
  snr = findgen(300)
  
  for i = 0, 299 do begin
     
     lmin = where(abs(lam[*,i] - 5700.00) eq min(abs(lam[*,i]-5700.00)))
     lmax = where(abs(lam[*,i] - 6700.00) eq min(abs(lam[*,i]-6700.00)))
     
     djs_iterstat, temp[lmin(0):lmax(0), i], sigma=sigma
     snr(i) = sigma
     
  endfor
  
  snr = 1/snr
  
  mag = plugmap.rapmag
  
  k = where(mag lt 20 and snr gt 0 and mag gt 0)
  plotsym, 0, /fill
  
  !P.multi=[0,2,1]
  j = where(snr(k) gt 1)
  
  plot, plugmap[k(j)].ra, plugmap[k(j)].dec, ps=8, title='Signal Noise Noise Greater than 1', /isotropic, yrange=[min(plugmap[k].dec-1/6.), max(plugmap[k].dec)], xtitle='Ra', ytitle='dec'
   j = where(snr(k) lt 1)
  
  plot, plugmap[k(j)].ra, plugmap[k(j)].dec, ps=8, title='Signal Noise Noise Less than 1', /isotropic, yrange=[min(plugmap[k].dec-1/6.), max(plugmap[k].dec)], xtitle='Ra', ytitle='dec'
  
  !P.multi=0
  
END


  
  
