function dec2hms, angle_in    
; pass type double, or suffer 0.003 arcsec error

    nangle = n_elements(angle_in)
    if nangle gt 1L then begin
       out = strarr(nangle)
       for k = 0L, nangle-1L do out[k] = dec2hms(angle_in[k])
       return, out
    endif
   
    angle=double(angle_in)
    neg = angle lt 0.0d
    angle=abs(angle)
    h=floor(angle)
    angle=(angle-h)*60.0d
    m=floor(angle)
    s=(angle-m)*60.0d
    if ((s ge 59.9949d) and (s le 60.00d)) then s=59.9949d
    
    hms=string(h,m,s,format='(I3,I3,F5.1)')
    if (strmid(hms,1,1) eq ' ') then strput,hms,'0',1
    if (strmid(hms,4,1) eq ' ') then strput,hms,'0',4
    if (strmid(hms,7,1) eq ' ') then strput,hms,'0',7
    if neg then strput,hms,'-',0
    
    out = hms

return, out
end
