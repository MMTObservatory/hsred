Pro correctplugmap, dir=dir
  
  ;;This program will go to a given directory, look for the plugdir sub-d
  ;;It will then read all of the files in this directory and look for 
  ;;instances where the aperture is either 43 or 45.  It will then switch 
  ;;the aperture numbers of these entries to correct for the mixed fiber names
  
  if not keyword_set(dir) then dir = cwd()
  
  files = findfile(dir + '/plugdir/*_all')
  
  if files(0) ne '' then begin
  for ifile = 0, n_elements(files) -1 do begin
     
     junk = strsplit(files(ifile), ' ' , length=ct)
     root = strmid(files(ifile), 0, ct-4)
     
     extensions = ['_all', '_targets', '_unused', '_sky']
     
     for iexten = 0, 3 do begin
        
        file = rc_readfile(root+extensions(iexten), 0)
      if strmatch(file(0), '') eq 0L then begin 
      if strmid(file(0), 0, 5) ne '#NOTE' then begin
        
        ap = indgen(n_elements(file))
        for i = 0, n_elements(file) -1 do begin
           junk = strsplit(file(i), ' ' , /extract)
           ap(i) = long(junk(0))
          
        endfor
        
        k43 = where(ap eq 43, n43)
        k45 = where(ap eq 45, n45)
        
        if n43 gt 0 then begin
           junk = strsplit(file(k43), ' ' , /extract)
           junk(0)= '45'
           file(k43) = ''
           for i = 0, n_elements(junk) -1 do begin
              file(k43) = file(k43) + ' ' + junk(i)
           endfor
           
        endif
        if n45 gt 0 then begin
           junk = strsplit(file(k45), ' ' , /extract)
           junk(0) = '43'
           file(k45) = ''
           for i = 0, n_elements(junk) -1 do begin
              file(k45) = file(k45) + ' ' + junk(i)
           endfor
           
        endif
        
        output = strarr(n_elements(file) +1)
        output(0) = '#NOTE - this files was swap corrected ' + im_today()
        output(1:n_elements(output)-1) = file
        
        openw, 1, root+extensions(iexten)
        for i = 0, n_elements(output) -1 do begin
           printf, 1, output(i)
        endfor
        close, 1
     endif
     
  endif
endfor
endfor
endif

end


        
        
        
