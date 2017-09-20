PRO correctmap, dir=dir
  
  ;;This program will find all of the various map files in a directory, and
  ;;correct the fact that apertrues 43/45 (1 indexed) are swapped.  
  
  if not keyword_set(dir) then dir = cwd()
  
  files1 = findfile(dir + '/' + '*map')
  files = files1
  
  if files(0) ne '' then begin
     
  for ifile = 0, n_elements(files) -1 do begin
     
     ;;Try to read the file as if it were a standard file
     readcol, files(ifile), ap, beam, target, ra, dec, tflag, $
        fiber, xfocal, yfocal, format='L,A,A,A,A,A,D,A,A'
     
     if total(long(ap)) eq 0 then type = 'cal'
     if total(long(ap)) ne 0 then type = 'obj'
     
     if type eq 'obj' then begin
        
        k47 = where(fiber eq 47, n47)
        k48 = where(fiber eq 48, n48)
        
        file = rc_readfile(files(ifile),0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
           if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47)
              file(k47) = file(k48)
              file(k48) = temp
              junk = strsplit(file(k47), ' ', /extract)
              junk1 = strsplit(file(k48), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48) = string
              
              
              
           ENDIF
        
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
        
        if file_test(files(ifile)+'1') eq 1l then begin
        file = rc_readfile(files(ifile)+'1',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
             
           if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47)
              file(k47) = file(k48)
              file(k48) = temp
              junk = strsplit(file(k47), ' ', /extract)
              junk1 = strsplit(file(k48), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48) = string
              
              
              
           ENDIF
        
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'1'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
        
        endif
        if file_test(files(ifile)+'2') eq 1l then begin

         file = rc_readfile(files(ifile)+'2',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
           
               if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47)
              file(k47) = file(k48)
              file(k48) = temp
              junk = strsplit(file(k47), ' ', /extract)
              junk1 = strsplit(file(k48), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48) = string
              
              
              
           ENDIF
        
           
           
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'2'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
        
        endif
        if file_test(files(ifile)+'1.tab') eq 1l then begin

          file = rc_readfile(files(ifile)+'1.tab',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
         
                      
           
               if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47+53)
              file(k47+53) = file(k48+53)
              file(k48+53) = temp
              junk = strsplit(file(k47+53), ' ', /extract)
              junk1 = strsplit(file(k48+53), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47+53) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48+53) = string
              
              
              
           ENDIF
        
             
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'1.tab'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
     endif
     
       if file_test(files(ifile)+'2.tab') eq 1l then begin
 
         file = rc_readfile(files(ifile)+'2.tab',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
                      
           
               if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47+53)
              file(k47+53) = file(k48+53)
              file(k48+53) = temp
              junk = strsplit(file(k47+53), ' ', /extract)
              junk1 = strsplit(file(k48+53), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47+53) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48+53) = string
              
              
              
           ENDIF
        
           
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'2.tab'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
     endif
     
        
     endif else if type eq 'cal' then begin
        readcol, files(ifile), ap, beam, target, $
        fiber, xfocal, yfocal, format='A,A,A,D,A,A' 
     
                k47 = where(fiber eq 47, n47)
        k48 = where(fiber eq 48, n48)
        
        file = rc_readfile(files(ifile),0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
           if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47)
              file(k47) = file(k48)
              file(k48) = temp
              junk = strsplit(file(k47), ' ', /extract)
              junk1 = strsplit(file(k48), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48) = string
              
              
              
           ENDIF
        
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
        
        if file_test(files(ifile)+'1') eq 1l then begin

        file = rc_readfile(files(ifile)+'1',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
             
           if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47)
              file(k47) = file(k48)
              file(k48) = temp
              junk = strsplit(file(k47), ' ', /extract)
              junk1 = strsplit(file(k48), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48) = string
              
              
              
           ENDIF
        
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'1'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
        endif
        if file_test(files(ifile)+'2') eq 1l then begin

         file = rc_readfile(files(ifile)+'2',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
           
               if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47)
              file(k47) = file(k48)
              file(k48) = temp
              junk = strsplit(file(k47), ' ', /extract)
              junk1 = strsplit(file(k48), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48) = string
              
              
              
           ENDIF
        
           
           
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'2'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
        endif
        if file_test(files(ifile)+'1.tab') eq 1l then begin

          file = rc_readfile(files(ifile)+'1.tab',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
         
                      
           
               if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47+53)
              file(k47+53) = file(k48+53)
              file(k48+53) = temp
              junk = strsplit(file(k47+53), ' ', /extract)
              junk1 = strsplit(file(k48+53), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47+53) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48+53) = string
              
              
              
           ENDIF
        
             
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'1.tab'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
     endif
     
        
        if file_test(files(ifile)+'2.tab') eq 1l then begin

         file = rc_readfile(files(ifile)+'2.tab',0)
        if strmid(file(0), 0, 5) ne '#NOTE'then begin
           
                      
           
               if n47 gt 0 and n48 gt 0 then begin
              temp = file(k47+53)
              file(k47+53) = file(k48+53)
              file(k48+53) = temp
              junk = strsplit(file(k47+53), ' ', /extract)
              junk1 = strsplit(file(k48+53), ' ', /extract)
              temp = junk(0)
              junk(0) = junk1(0)
              junk1(0) = temp
              string = '' 
              for i = 0, n_elements(junk) -1 do begin
                 string = string +  '  ' + junk(i)
              endfor
              file(k47+53) = string
               string = '' 
              for i = 0, n_elements(junk1) -1 do begin
                 string = string +  '  ' + junk1(i)
              endfor
              file(k48+53) = string
              
              
              
           ENDIF
        
           
        
           output = strarr(n_elements(file)+1)
           output(0) =$
              '#NOTE -This file was corrected for the swap of fiber 47/48'
           output(1:n_elements(output)-1) = file
        
           openw, 1, files(ifile)+'2.tab'
           for i = 0, n_elements(output) -1 do begin
              printf, 1, output(i)
           endfor
           close, 1
        endif
     endif
     
        
        
        
        
     endif
     
  endfor
endif

end

