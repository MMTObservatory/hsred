;+
; NAME:
;   hs_batch1d
;
; PURPOSE:
;  Given a plugplan and reduction plans this will reduce 1d data
;
; CALLING SEQUENCE:
;  hs_batch, redplan=redplan, rerun=rerun, plugplan=plugplan,link=link, $
;              basedir=basedir, checkaverage=checkaverage, docosmic=docosmic, $
;              nofluxplan=nofluxplan, superfile=superfile
;
;
; INPUTS:
; ;
; OPTIONAL KEYWORDS:
;   redplan      - List of directories that need to be reduced (with fluxing)
;   nofluxplan   - List of directories that need to be reduced without fluxing
;   plugplan     - A translation file between config and field pass
;   checkaverage - Same as in hs_extract
;   superfile    - same as in hs_extract
;   docosmic     - same as in hs_extract
;   link         - if set, run hs_makelinks at the end
;
; OUTPUTS:
;;
; OPTIONAL OUTPUTS:
;;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;   
;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA         
;                
;-
;------------------------------------------------------------------------------


PRO hs_batch1d, redplan=redplan, rerun=rerun, plugplan=plugplan,link=link, $
              basedir=basedir, checkaverage=checkaverage, docosmic=docosmic, $
              nofluxplan=nofluxplan, superfile=superfile
  
  ;;nofluxplan is for files that are to be reduced without fluxing
  
  
  
  if not keyword_set(basedir) then basedir = cwd()
  
  if not keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if size(rerun, /tname) NE 'STRING' then $
     rerun = string(rerun, format='(i4.4)')   
  
  
  if keyword_set(redplan) then begin
     
     readcol, redplan, redir, format='A', comment='#'
     origdir = redir
  endif
  if keyword_set(nofluxplan) then begin
     readcol, nofluxplan, redir, format='A', comment='#'
     orignofluxdir = redir
  endif
  
  stop = 'go'
   if strmid(plugplan,0,1) ne '/' and strmid(plugplan,0,1) ne '.' $
        and strmid(plugplan,0,1) ne '~' then plugplan = cwd() + plugplan
   if keyword_set(redplan) then begin
      
      if strmid(redplan,0,1) ne '/' and strmid(redplan,0,1) ne '.' $
        and strmid(redplan,0,1) ne '~'  then redplan = cwd() + redplan
   endif else if keyword_set(nofluxplan) then begin
      if strmid(nofluxplan,0,1) ne '/' and strmid(nofluxplan,0,1) ne '.' $
        and strmid(nofluxplan,0,1) ne '~' then nofluxplan = cwd() + nofluxplan
   endif
   
   
 while stop ne 'stop' do begin
  
    
  ;;Now do the 1-D reductions
  
   if keyword_set(redplan) then begin
     
     readcol, redplan, redir, format='A', comment='#'
     print, redir
     originaldir = cwd()
    
   
     for idir = 0, n_elements(redir) -1 do begin
        
        print, redir(idir)
        cd, redir(idir)+'reduction/' + rerun+'/'
        redfiles = findfile('spHect*fits', count=nfile)
        root = redfiles
        for i = 0, nfile -1 do begin
           junk = strsplit(redfiles(i), ' ',  length=nchar)
           root(i) = strmid(redfiles(i), 7, nchar-7)
        endfor
       
        
        if nfile gt 0 then begin
           
           for ifile = 0, nfile-1 do begin
              badfile = findfile('spZbest-'+root(ifile), count=nfile)
              if nfile gt 0 then print, badfile
              
              if nfile eq 0 then $
                           hs_reduce1d, redfiles(ifile)
           endfor
        endif
        
        
           cd, originaldir
  
        
     endfor
     
   
     
  endif
  
  
  
   if keyword_set(nofluxplan) then begin
     
     readcol, nofluxplan, redir, format='A', comment='#'
  
     originaldir = cwd()
  
  
     for idir = 0, n_elements(redir) -1 do begin
   
        cd, redir(idir)+'reduction/' + rerun+'/'
        redfiles = findfile('spHect*fits', count=nfile)
        root = redfiles

        for i = 0, nfile -1 do begin
           junk = strsplit(redfiles(i), ' ',  length=nchar)
           root(i) = strmid(redfiles(i), 7, nchar-7)
        endfor
        
        
        if nfile gt 0 then begin
           for ifile = 0, nfile-1 do begin
              badfile = findfile('spZbest-'+root(ifile), count=nfile)           
           if nfile gt 0 then print, badfile                                 
             
           if nfile eq 0 then $                                              
                       hs_reduce1d, redfiles(ifile) ,/pseudoflux 
           endfor
        endif
        
        
        cd, originaldir

        
     endfor
     
      
  endif
  
  
  
  
   
  if keyword_set(redplan) then begin
     
     readcol, redplan, redir, format='A', comment='#'
     newdir = redir
     
     
  endif
  if keyword_set(nofluxplan) then begin
     readcol, nofluxplan, redir, format='A', comment='#'
     newnofluxdir = redir
  endif
  if keyword_set(redplan) and keyword_set(nofluxplan) then begin
     test = origdir eq newdir and orignofluxdir eq newnofluxdir
     k = where(test eq 0, ct)
     if ct eq 0 then stop='stop'
  endif
  if keyword_set(redplan) and not keyword_set(nofluxplan) then begin
     test = origdir eq newdir
     k = where(test eq 0, ct)
     if ct eq 0  then stop='stop'
  endif
  if keyword_set(nofluxplan) and not keyword_set(redplan) then begin
     test = orignofluxdir eq newnofluxdir
     k = where(test eq 0, ct)
     if ct eq 0 then stop='stop'
  endif
  
  print, stop
  
endwhile

  
  
  
  
  
  
end

  
