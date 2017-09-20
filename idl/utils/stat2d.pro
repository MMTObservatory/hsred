FUNCTION stat2d, value, limitarray, min, max, binsize, equal=equal, $
                 weights=weights
  
  ;;This program does statistics in two dimensions.  The value array is the
  ;;one for which all the stats are taken, and the limit array is used
  ;;to limit to statistocs
  
  if not keyword_set(weights) then $
     weights = replicate(1.0, n_elements(value))
  if not keyword_set(equal) then begin
     
     nbin = (max-min)/binsize
  
     ;;Make nbin the approporiate long
     if (nbin ne long(nbin)) then nbin = long(nbin) +1
  
  
     ;;Make bounds
     bound = dblarr(nbin)
     bound(0) = min
     for i = 1, nbin-1 do begin
        bound(i) = bound(i-1) + binsize
     endfor
  
     valmean = dblarr(nbin)
     valsigma = dblarr(nbin)
     valnumber = dblarr(nbin)
  
     ;;Now loop through the bins
     for i = 0, nbin - 1 do begin
        k = where(limitarray gt bound(i) and limitarray lt bound(i) + binsize)
        if k(0) ne -1 then begin
           djs_iterstat, value(k)*weights(k), mean=meantemp, sigma=sigmatemp
           valmean(i) = meantemp
           valsigma(i)=sigmatemp
           valnumber(i) = n_elements(k)
        endif else begin
           valmean(i) = 0
           valsigma(i) = 0
           valnumber(i) = 0
        endelse
     endfor
     midpoint = bound + binsize/2.0
     
  endif else begin
     
     k = where(limitarray ge min and limitarray le max)
     nbin = n_elements(k)/binsize
     k1 = sort(limitarray)
     
     
     ;;Make bounds
     bound = dblarr(nbin+1)
     bound(0) = min
     ndone = 0
     for i = 1, nbin-1 do begin
        max1 = max(limitarray(k1(indgen(binsize)+ndone)))
        bound(i) = max1
        ndone = ndone + binsize
     endfor
     bound(nbin) = max
     
     
            
     valmean = dblarr(nbin)
     valsigma = dblarr(nbin)
     valnumber = dblarr(nbin)
  
     ;;Now loop through the bins
     for i = 0, nbin - 1 do begin
        k = where(limitarray ge bound(i) and limitarray lt bound(i+1))
        if k(0) ne -1 then begin
           djs_iterstat, value(k)*weights(k), mean=meantemp, sigma=sigmatemp
           valmean(i) = meantemp
           valsigma(i)=sigmatemp
           valnumber(i) = n_elements(k)
        endif else begin
           valmean(i) = 0
           valsigma(i) = 0
           valnumber(i) = 0
        endelse
     endfor
     
     midpoint = fltarr(n_elements(bound)-1)
     for i = 0, n_elements(midpoint) -1 do begin
        midpoint(i) = (bound(i+1)+bound(i))/2.0
     endfor
     
     
     
  endelse
       
     
    output = create_struct('midpoint', midpoint, $
                         'mean', valmean, $
                         'sigma', valsigma, $
                         'number', valnumber)

  return, output
  
END
