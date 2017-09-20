;;This program is designed to emulate the hselect program from IRAF
;;files = [N] vector of file names
;;values = [U] vector of header values to find
;;exten = extention number to use

PRO hselect, files, values, exten=exten
  if not keyword_set(exten) then exten = 0
  
;  output = strarr(n_elements(files), n_elements(values))
  output = strarr(n_elements(files))
  for ifile = 0, n_elements(files) -1 do begin
     
     hdr = headfits(files(ifile), exten=exten)
     
     
     for ival = 0, n_elements(values)-1 do begin
        
        IF values(ival) EQ '$I' OR strupcase(values(ival)) EQ $
           'FILENAME' THEN BEGIN
           output(ifile) = output(ifile) + ' ' + files(ifile)
        ENDIF ELSE BEGIN
           output(ifile) = output(ifile) + ' ' + string(sxpar(hdr, values(ival)))
        ENDELSE
        
        
         ;output(ifile, ival) = sxpar(hdr, values(ival))
         ;if  values(ival) eq '$I' or strupcase(values(ival)) eq $ 
         ;      'FILENAME' then $
                 ; output(ifile, ival) = files(ifile)
               
     endfor
     
  endfor

  
  niceprint_vect, output
  
END


     
