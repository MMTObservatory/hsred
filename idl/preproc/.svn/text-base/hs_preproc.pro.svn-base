;;This program looks at the *fits files in the provided directory, and then makes the lists needed for future reductions.
; edited 01/29/13 by smoran to fix bugs in string processing for
; making bias list, and for making cal.list from existing plugmap files
PRO hs_preproc, inpath=inpath, overwrite=overwrite, plugdir=plugdir
  
  
  if NOT keyword_set(inpath) then inpath = cwd()
  if NOT keyword_set(plugdir) then plugdir = inpath + 'plugdir/'
  
  

  ;;Cehck to see if lists is already present if not make it
  if direxist(inpath + 'lists') EQ 0L then  spawn, 'mkdir ' + inpath + 'lists'
  
  files = findfile(inpath + '*.????.fits', count=nfiles)
  files2 = findfile(inpath + '*.????.fits.gz', count=nfiles2)
  if nfiles eq 0 then begin
     nfiles = nfiles2
     files = files2
  endif else if nfiles2 ne 0 then begin
     files = [files, files2]
     nfiles = nfiles+nfiles2
  end



  object = strarr(nfiles)
  imagetyp = strarr(nfiles)
  mjd = strarr(nfiles)
  exptime = fltarr(nfiles)
  
  for i = 0, nfiles -1 do begin
     hdr = headfits(files(i),exten=1)
     IF n_elements(hdr) EQ 1 THEN $
        message, 'It appears as though ' + files(i) + ' has a faulty header'
     object(i) = sxpar(hdr, 'OBJECT')
     imagetyp(i) = sxpar(hdr, 'IMAGETYP')
     mjd(i) = long(sxpar(hdr, 'MJD'))
     exptime(i) = sxpar(hdr, 'EXPTIME')
  endfor
  
  ;;Get the biases and make bias_mjd.list if it does not exhist
  
  if file_test('lists/bias.list') NE 0 AND $
     NOT keyword_set(overwrite) then begin
     splog, 'Bias file lists/bias.list already exists.  Will not '
     splog, 'overwrite'
  endif else if (file_test('lists/bias.list') NE 0 AND $
                 NOT keyword_set(overwrite)) OR file_test('lists/bias.fits') $
     EQ 0 Then begin
     
     k = where(strmatch(strupcase(imagetyp), 'ZERO*')and not strmatch(strupcase(object), 'JUNK*'), ct)
     openw, 1, 'lists/bias.list'
     printf, 1, '#FILE      OBJECT     IMAGETYPE'
     
     FOR I = 0, CT-1 DO BEGIN
        printf, 1, files[k[i]] + '   ' + object[k[i]] + '   ' +  imagetyp[k[i]]
     endfor
     
     close, 1
     
  endif
  
  
  ;;Get the dlfats and make dflat_mjd.list if it does not exhist
  
  if file_test('lists/dflat.list') NE 0 AND $
     NOT keyword_set(overwrite) then begin
     
     splog, 'Dflat File lists/dflat.list already exists.  Will not '
     splog, 'overwrite'
  endif else if (file_test('lists/dflat.list') NE 0 AND $
                 NOT keyword_set(overwrite)) OR $
     file_test('lists/dflat.fits') $
     EQ 0 Then begin
     
     k = where(strmatch(strupcase(imagetyp), 'DOMEFLAT*'), ct)
     if ct eq 0 then begin
        k = where(strmatch(strupcase(object), 'DOMEFLAT*'), ct)
        if ct gt 0 then begin 
           splog, 'No images with IMAGETYP=DOMEFLAT FOUND. ' + $
              'USING OBJECT=DOMEFLAT'
           splog, 'You may want to make sure that dflat.list is correct'
        endif
     endif
     
     if ct eq 0 then begin
        splog, 'NO DOMEFLAT FILES FOUND'
     endif else if ct gt 0 then begin
        
        openw, 1, 'lists/dflat.list'
        printf, 1, '#FILE      OBJECT     IMAGETYPE'
        FOR I = 0, CT-1 DO BEGIN
           printf, 1, files(k(i)) + '   ' + object(k(i)) $
              + '   ' +  imagetyp(k(i))
        endfor
     
        close, 1
     endif
  endif
  
  ;;Get the darks if it does not exhist
  
  if file_test('lists/dark.list') NE 0 AND $
     NOT keyword_set(overwrite) then begin
     
     splog, 'Sflat File lists/dark.list already exists.  Will not '
     splog, 'overwrite'
  endif else if (file_test('lists/dark.list') NE 0 AND $
                 NOT keyword_set(overwrite)) OR file_test('lists/dark.fits') $
     EQ 0 Then begin
     
     k = where(strmatch(strupcase(imagetyp), 'DARK*'), ct)
     if ct eq 0 then begin
        k = where(strmatch(strupcase(object), 'DARK'), ct)
        if ct gt 0 then begin 
           splog, 'No images with IMAGETYP=SKYFLAT FOUND. USING OBJECT=SKYFLAT'
           splog, 'You may want to make sure that lists/sflat.list is correct'
        endif
        
     endif
     
     if ct eq 0 then begin
        splog, 'NO DARK FILES FOUND'
     endif else if ct gt 0 then begin
        openw, 1, 'lists/dark.list'
        printf, 1, '#FILE      OBJECT     IMAGETYPE'
     
        FOR I = 0, CT-1 DO BEGIN
           printf, 1, files(k(i)) + '   ' + object(k(i)) $
              + '   ' +  imagetyp(k(i))
        endfor
     
        close, 1
     endif
  endif
  
  ;;Get the sflats and make sflat.list if it does not exhist
  
  if file_test('lists/sflat.list') NE 0 AND $
     NOT keyword_set(overwrite) then begin
     
     splog, 'Sflat File lists/sflat.list already exists.  Will not '
     splog, 'overwrite'
  endif else if (file_test('lists/sflat.list') NE 0 AND $
     NOT keyword_set(overwrite)) OR file_test('lists/sflat.fits') $
     EQ 0 Then begin
     
     k = where(strmatch(strupcase(imagetyp), 'SKYFLAT*'), ct)
     if ct eq 0 then begin
        k = where(strmatch(strupcase(object), 'SKYFLAT*'), ct)
        if ct gt 0 then begin 
           splog, 'No images with IMAGETYP=SKYFLAT FOUND. USING OBJECT=SKYFLAT'
           splog, 'You may want to make sure that lists/sflat.list is correct'
        endif
        
     endif
     
     if ct eq 0 then begin
        splog, 'NO SKYFLAT FILES FOUND'
     endif else if ct gt 0 then begin
        openw, 1, 'lists/sflat.list'
        
        printf, 1, '#FILE      OBJECT     IMAGETYPE'
        FOR I = 0, CT-1 DO BEGIN
           printf, 1, files(k(i)) + '   ' + object(k(i)) + '   ' $
              +  imagetyp(k(i))
        endfor
        
        close, 1
     endif
  endif
  
  
  ;;Get the arcs and make arc_Mjd.list if it does not exhist
  
  if file_test('lists/arc.list') NE 0 AND $
     NOT keyword_set(overwrite) then begin
     
     splog, 'Arc File lists/arc.list already exists.  Will not '
     splog, 'overwrite'
  endif else if (file_test('lists/arc.list') NE 0 AND $
                 NOT keyword_set(overwrite)) OR file_test('lists/arc.fits') $
     EQ 0 Then begin
     
     k = where(strmatch(strupcase(imagetyp), 'COMP*'), ct)
     if ct eq 0 then begin
        splog, 'NO ARCLINE FILES FOUND'
     endif else if ct gt 0 then begin
          
        openw, 1, 'lists/arc.list'
        printf, 1, '#FILE      OBJECT     IMAGETYPE'
     
        FOR I = 0, CT-1 DO BEGIN
           printf, 1, files(k(i)) + '   ' + object(k(i)) + '   '$
              +  imagetyp(k(i))
        endfor
        
        close, 1
     endif
  endif
  
  
  plugmaps = findfile(plugdir+'*all', count=nmaps)
  used = indgen(n_elements(files))*0.0
  
  if file_test(inpath + 'lists/cal.list') NE 0L and $
     NOT keyword_set(overwrite) then begin
     splog, 'File lists/cal.list already exists'
     splog, 'Not overwriting'
  endif else  begin
     
     openw, 1, inpath + 'lists/cal.list'
     printf, 1, '#FILES                              PLUGMAP'
     junk = strsplit(inpath, ' ', length=incount)
     if nmaps ge 1 then begin
     for i = 0, nmaps -1 do begin
        junk = strsplit(plugdir, ' ', length=dirct)
        junk = strsplit(plugmaps(i), ' ', length=count)
        junk = strsplit(inpath, ' ', length=incount)
        filetrim = files
        
        for j =0, n_elements(files) -1 do begin
           junk = strsplit(files(j), ' ', length=N)
           filetrim(j) = strmid(files(j), incount, n-incount-10)
        endfor
        plug = strmid(plugmaps(i),dirct,count-dirct-4)
        ;plug = strsplit(plug, '.', /extract)
        ;plug = plug(0)
        k = where(strmatch(filetrim, plug), ct)
        
                
        if ct ne 0 then begin    
           used(k) = 1
           output = ''
           
           for j = 0, n_elements(k) -1 do begin
              if j gt 0 then output = output + ','
              if exptime(k(j)) gt 1 then output = output + files(k(j))
              
           endfor
           
           printf, 1, output + '     ' + $
              strmid(plugmaps(i), dirct, count-dirct-4)
        endif
        
     endfor
  endif else begin
        k = where(strmatch(strupcase(imagetyp), 'OBJECT*') and (exptime gt 1.), ct)
        if ct eq 0 then begin
           splog, 'NO OBJECT FILES FOUND'
       endif else if ct gt 0 then begin
          filetrim=files[k]          
          for j =0, n_elements(filetrim) -1 do begin
             junk = strsplit(files[k[j]], ' ', length=N)
             filetrim[j] = strmid(files[k[j]], incount, n-incount-10)
          endfor
          filetrim=filetrim[sort(filetrim)]
          filetrim=filetrim[uniq(filetrim)]
          for j =0, n_elements(filetrim) -1 do begin
             match = strmatch(files[k],'*'+filetrim[j]+'*')
             used[k[where(match)]]=1
             outstring=strjoin(files[k[where(match)]],',')
             printf, 1, outstring + '  a ' 
          endfor
       endif
  endelse

     k = where(used eq 0, ct)
     
     if ct gt 0 then begin
        for j = 0, ct -1 do begin
           
           if strmatch(strupcase(imagetyp(k(j))), $
                       'OBJECT*') then begin
              if NOT strmatch(strupcase(object(k(J))), $
                              '*FLAT*') then begin
                 printf, 1, files(k(j))
              endif
           endif
           
           
        endfor
     endif
     close, 1
  endelse
  
  END

  
  
