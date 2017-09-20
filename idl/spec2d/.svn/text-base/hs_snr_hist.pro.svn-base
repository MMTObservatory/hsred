;+
; NAME:
;   hs_snr_hist
;
; PURPOSE:
;   Create the signal to noise histogram for the data
;
; CALLING SEQUENCE:
;  hs_snr_file, lam, flux, plugmap
;
; INPUTS:
;    lam - wavelength [npix, nfib]
;    flux - flux [npix,nfib]
;    plugmap - plugmap [ nfib]
;
; OPTIONAL KEYWORDS:
;
;
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
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

PRO hs_snr_hist, lam, flux,  plugmap, title=title
  
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
  
  k = where(strmatch(plugmap.objtype, 'target'))
  
  hist1= rc_hist(snr, 0, 100, 1)
  max = max(snr)
  plot, hist1.midpoint, hist1.number,ps=10, $
     xrange=[0,max], xtitle='Signal to Noise', ytitle='Number', title=title
  
  mean1 = djs_mean(snr)
  median1 = djs_median(snr)

  maxn= max(hist1.number)
  locate=max(snr)*0.5
  xyouts, locate, maxn-3, 'Mean = ' + string(mean1), charsize=1.5
  xyouts, locate, maxn-6, 'Median = ' + string(median1), charsize=1.5
  xyouts, locate, maxn-9,  'N : SNR>3 = ' +$
     string(n_elements(where(snr gt 3.0))), charsize=1.5
  xyouts, locate, maxn-12, 'Fraction of SNR > 3' + $
     string(float(n_elements(where(snr gt 3.0)))/$
            float(n_elements(snr))*100.), charsize=1.5
  
END


  
  
