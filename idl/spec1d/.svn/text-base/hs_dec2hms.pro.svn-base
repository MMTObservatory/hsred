
; Written by Douglas Finkbeiner in ancient times
; 29 Nov 2000 - double keyword added
; pass type double, or suffer 0.003 arcsec error
;   2005-Sep-16 - call self recursively for array - DPF
;----------------------------------------------------------------------
function hs_dec2hms,angle_in, double=double, sep=sep  

; -------- if angle_in is an array, call self recursively
  n_angle = n_elements(angle_in) 
  if n_angle GT 1 then begin 
     hmsarr = strarr(n_angle)
     for i=0L, n_angle-1 do begin 
        hmsarr[i] = dec2hms(angle_in[i], double=double, sep=sep)
     endfor 
     return, hmsarr
  endif

; -------- if there is only one angle, carry on with old routine
 
  eps = (machar(double=double)).eps ; machine precision

  angle = double(angle_in)
  neg   = angle LT 0.0d
  angle = abs(angle)
  h     = floor(angle+eps)
  angle = (angle-h)*60.0d
  m     = floor(angle+eps*60)
  s     = ((angle-m)*60.0d) > 0 ; can be slightly LT 0 due to roundoff
  
  if keyword_set(sep) then cc=sep else cc=':'
  if keyword_set(double) then begin 
     if ((s ge 59.999949d) and (s le 60.00d)) then s=59.999949d
     hms=strtrim(string(h, format='(I3.2)'),2)+cc+strtrim(string(m, format='(I3.2)'),2)+cc+strtrim(string(s, format='(F8.4)'),2)
  endif else begin 
     if ((s ge 59.949d) and (s le 60.00d)) then s=59.949d
     hms=strtrim(string(h, format='(I3.2)'),2)+cc+strtrim(string(m, format='(I3.2)'),2)+cc+strtrim(string(s, format='(F5.1)'),2)
;     hms=string(h,cc,m,cc,s,format='(I3.2,A1,I3.2,A1,F5.1)')
  endelse 

  if (strmid(hms,7,1) eq ' ') then strput,hms,'0',7
  if neg then hms='-'+hms

return,hms
end
