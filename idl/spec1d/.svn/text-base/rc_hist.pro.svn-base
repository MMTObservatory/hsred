;+
; NAME:
;   rc_hist
;
; PURPOSE:
;   Creates a histogram for the data give
;
; CALLING SEQUENCE:
;  histstruct= rc_hist( value, min, max, binsize)
;
; INPUTS:
;  value  -  vector to be used to create the histogram
;  min    -  minimum value of the array to be used in the histogram
;  max    -  maximum value of the array to be used in the histrogram
;  binsize - width of the bin to use in the historgam
;

; OUTPUTS:
;  histstuct = stucture with the histogram information included
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
; PROCEDURES CALLED:
;   readfits()
;   sxpar()
;   sxaddpar
;   writefits
;   create_struct
;
; REVISION HISTORY:
;   May 2004 - Written ; Cool UofA
;-
;-----------------------------------------------------------------------------


  
  
  
  
FUNCTION rc_hist, value, min, max, binsize
  
  ;;This program takes the array 'value' and makes a histogram of the values 
  ;;and returns a structure with the bounds and the number
  
  
  nbin = (max-min)/binsize
  
  ;;make sure nbin is integer
  if (nbin ne long(nbin)) then nbin = long(nbin) + 1
  
  bounds = dblarr(nbin)
  bounds(0) = min
  for i = 1, nbin-1 do begin
     bounds(i) = bounds(i-1) + binsize
  endfor
  
  number =  lindgen(nbin)
  
  for i = 0 , nbin - 1 do begin
     ;;Find the elements that are in this bin
     
     k = where(value ge bounds(i) and value lt bounds(i) + binsize)
     if k(0) ne -1 then number(i) = n_elements(k)
     if k(0) eq -1 then number(i) = 0
     
  endfor
  
  
  output = create_struct('midpoint', bounds+binsize/2., $
                         'number', number)
  
  return, output
  
END


