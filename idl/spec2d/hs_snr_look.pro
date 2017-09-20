
PRO hs_snr_look, lam, flux, plugmap, nsmooth=nsmooth
  
 
  
  flux1=flux*0.0
  npix = n_elements(flux(*,0))
  flux1(1:npix-1,*) = flux(0:npix-2,*)
   
   
  diff = flux-flux1
  sum = flux+flux1
  
  temp = sqrt(2) * diff/sum
  
    
  snr = findgen(300)
  k=indgen(300)
  for i = 0, 299 do begin
     
     lmin = where(abs(lam[*,i] - 5700.00) eq min(abs(lam[*,i]-5700.00)))
     lmax = where(abs(lam[*,i] - 6700.00) eq min(abs(lam[*,i]-6700.00)))
     
     djs_iterstat, temp[lmin(0):lmax(0), i], sigma=sigma
     snr(i) = sigma
     
  endfor
  
  snr = 1/snr
  
  sort1 = sort(snr)
  hs_page,lam[*,sort1], flux[*,sort1], title=snr(sort1), nsmooth=nsmooth,plugmap=plugmap
  
  
END


