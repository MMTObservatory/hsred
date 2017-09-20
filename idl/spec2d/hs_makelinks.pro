;+
; NAME:
;   hs_makelinks
;
; PURPOSE:
;   Link the plugplan to the reduced data
;
; 
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------



PRO hs_makelinks, basedir=basedir, plugplan=plugplan, redplan=redplan, $
                  rerun=rerun
  
  ;;Check for base directory
  
  
  if not keyword_set(rerun) then rerun = 0
  if size(rerun, /tname) NE 'STRING' then $
     rerun = string(rerun, format='(i4.4)')  
  basedir = basedir + string(rerun, format='(i4.4)') +'/'
  
  if direxist(basedir) eq 0l then spawn, 'mkdir ' + basedir
  
  readcol, plugplan, dir, plugmap, plugfield, plugpass, $
     format='A,A,A,A', comment='#'
  readcol, redplan, redir, format='A', comment='#'
  
  ;;Now loop through the reduction files and look for files
  
  redir = redir + '/'
  for i = 0, n_elements(redir) -1 do begin
     files = findfile(redir(i) + 'reduction/' + rerun + '/spHect-???-' +$
                      rerun + '.fits', count=nfiles)
     if nfiles gt 0 then begin
        junk = strsplit(redir(i) +'reduction/' + rerun +'/', ' ',$
                        length=pathcount)
        
        pass = strmid(files, pathcount+7,  1)
        field = strmid(files, pathcount+1+7, 2)
        
        for j = 0, n_elements(pass) -1 do begin
           
           files1 =  findfile(redir(i) + 'reduction/' + rerun + '/spZ???-' +$
                              pass(j)+field(j)+'-'+rerun+     '.fits')
           files2 =  findfile(redir(i) + 'reduction/' + rerun + '/spZ????-' +$
                              pass(j)+field(j)+'-'+rerun+     '.fits')
           files3 =  findfile(redir(i) + 'reduction/' + rerun + '/qaplot_' +$
                              pass(j)+field(j)+'-'+rerun+     '.ps')
           files4 = findfile(redir(i) + 'reduction/' + rerun + '/spObs-' + $
                             pass(j)+field(j)+'-'+rerun+'*.fits')
           expid = strarr(n_elements(files4))
           
           for ifile = 0, n_elements(files4) - 1 do begin
              junk = strsplit(files4(ifile), '-.', /extract)
              expid(ifile)=junk(n_elements(junk)-2)
           endfor
        
           for ifile = 0, n_elements(expid) -1 do begin
              files5temp = findfile(redir(i) +'reduction/' + $
                                    rerun + '/spFluxcorr*' + expid(ifile)+$
                                    '*.fits')
              if n_elements(files5) ne 0 then files5 = [files5, files5temp]
              if n_elements(files5) eq 0 then files5 = files5temp
           endfor
           
           for ifile = 0, n_elements(expid) -1 do begin
              files6temp = findfile(redir(i) +'reduction/' + $
                                    rerun + '/hsFluxcalib*' + $
                                    expid(ifile)+'*fits')
              
              if n_elements(files6) ne 0 then files6 = [files6, files6temp]
              if n_elements(files6) eq 0 then files6 = files6temp
           endfor 
           
           k = where(plugfield eq long(field(j)) and plugpass eq long(pass(j)))
           
           if k(0) ne -1 then begin
              
              if direxist(basedir + dir(k)) EQ 0L then begin
                 split = strsplit(dir(k), '/', /extract)
                 
                 spawn,  'mkdir ' + basedir + split(0)
                 spawn, 'mkdir ' + basedir + split(0) + '/' + split(1)
              endif
              
              todo = [files(j), files1, files2, files3, files4, files5, files6]
              delvarx, files5
              delvarx, files6
              
              kl = where(strmatch(todo, '/*'),  nempty)
              todo=todo[kl]
              for l = 0, n_elements(todo) -1 do begin
                 spawn, 'ln -s ' + todo(l) +' ' +  basedir + dir(k)
                 
                 
              endfor
           endif
        endfor
     endif
  endfor
end


