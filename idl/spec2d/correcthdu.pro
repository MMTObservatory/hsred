PRO correcthdu, dir=dir
  
  ;;This program will go to a given directory, find all of the spObs,
  ;;spHect, spZbest, spZline, and spZall files, and correct the HDUs in each.  
  ;;This correction swaps only the *target* information. That means that 
  ;;any information from the spectra must stay the same
  
  if not keyword_set(dir) then dir = cwd()
  
  ;;Do the spObs first
  
  files = findfile(dir+'/spObs*fits')
  if files(0) ne '' then begin
     
  for ifile = 0, n_elements(files) -1 do begin
     
     ;;spobs - 6 HDUs - #5 is the plugmap
     plugmap = mrdfits(files(ifile), 5)
     k = [42, 44]
     tempplug = plugmap[k]
     
     ;;Now we need to correct the fiberid entry (the only spectra dependant)
     tempplug[0].fiberid = 44
     tempplug[1].fiberid = 42
     
     ;;Now feed this back into the plugmap if it will let me
     plugmap[42] = tempplug[1]
     plugmap[44] = tempplug[0]
     
     spawn, 'cp ' + files(ifile) + '  /tmp/spec.fits'

     for i = 0, 4 do begin
    
        temp = mrdfits('/tmp/spec.fits',i)
        
        if i eq 0 then mwrfits, temp, files(ifile), /create
        if i ne 0 then mwrfits, temp, files(ifile)
     endfor
     
     mwrfits, plugmap, files(ifile)
     
  endfor
  endif
    
  files = findfile(dir+'/spHect*fits')
    if files(0) ne '' then begin

  for ifile = 0, n_elements(files) -1 do begin
     
     ;;spobs - 9 HDUs - #5 is the plugmap
     plugmap = mrdfits(files(ifile), 5)
     k = [42, 44]
     tempplug = plugmap[k]
     
     ;;Now we need to correct the fiberid entry (the only spectra dependant)
     tempplug[0].fiberid = 45
     tempplug[1].fiberid = 43
     
     ;;Now feed this back into the plugmap if it will let me
     plugmap[42] = tempplug[1]
     plugmap[44] = tempplug[0]
     
     spawn, 'cp ' + files(ifile) + '  /tmp/spec.fits'

     for i = 0, 4 do begin
    
        temp = mrdfits('/tmp/spec.fits',i)
        if i eq 0 then mwrfits, temp, files(ifile), /create
        if i ne 0 then mwrfits, temp, files(ifile)
     endfor
     
     mwrfits, plugmap, files(ifile)
     
     for i =6, 9 do begin
        temp = mrdfits('/tmp/spec.fits',i)
        
        if size(temp, /tname) ne 'INT' then $
           mwrfits, temp, files(ifile)
     endfor
     
  
     
  endfor
  endif
  
  files = findfile(dir+'/spZbest*fits')
    if files(0) ne '' then begin

  for ifile = 0, n_elements(files) -1 do begin
     
     ;;spobs - 3 HDUs - #1 is the structure - only need to change
     ;;objtype, plug_ra, plug_dec
     
     plugmap = mrdfits(files(ifile), 1)
     k = [42, 44]
     tempplug = plugmap[k]
     
     ;;Now we need to correct the fiberid entry (the only spectra dependant)
     plugmap[42].objtype = tempplug[1].objtype
     plugmap[42].plug_ra = tempplug[1].plug_ra
     plugmap[42].plug_dec = tempplug[1].plug_dec
     plugmap[44].objtype = tempplug[0].objtype
     plugmap[44].plug_ra = tempplug[0].plug_ra
     plugmap[44].plug_dec = tempplug[0].plug_dec
     
     
     spawn, 'cp ' + files(ifile) + '  /tmp/spec.fits'

          
     mwrfits, plugmap, files(ifile), /create
     
     for i =1,2 do begin
        temp = mrdfits('/tmp/spec.fits',i)
        
        if size(temp, /tname) ne 'INT' then $
           mwrfits, temp, files(ifile)
     endfor
  endfor
  endif
   files = findfile(dir+'/spZall*fits')
   if files(0) ne '' then begin

  for ifile = 0, n_elements(files) -1 do begin
     
     ;;spobs - 1 HDUs - #1 is the structure - only need to change
     ;;objtype, plug_ra, plug_dec
     
     plugmap = mrdfits(files(ifile), 1)
     k43 = where(plugmap.fiberid eq 43)
     k45 = where(plugmap.fiberid eq 45)
     tempplug = plugmap[k43]
     
     ;;Now we need to correct the fiberid entry (the only spectra dependant)
     plugmap[k43].objtype = plugmap[k45].objtype
     plugmap[k43].plug_ra = plugmap[k45].plug_ra
     plugmap[k43].plug_ra = plugmap[k45].plug_dec
     plugmap[k45].objtype = tempplug.objtype
     plugmap[k45].plug_ra = tempplug.plug_ra
     plugmap[k45].plug_ra = tempplug.plug_dec
     
     
     spawn, 'cp ' + files(ifile) + '  /tmp/spec.fits'

          
     mwrfits, plugmap, files(ifile), /create
     
    endfor
 endif
 
    
    
    
end

        
