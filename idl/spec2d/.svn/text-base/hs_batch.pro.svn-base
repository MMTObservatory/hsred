;+
; NAME:
;   hs_batch
;
; PURPOSE:
;  Given a plugplan and reduction plans this will reduce data
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


PRO hs_batch, redplan=redplan, rerun=rerun, plugplan=plugplan,link=link, $
              basedir=basedir, checkaverage=checkaverage, docosmic=docosmic, $
              nofluxplan=nofluxplan, superfile=superfile, oned=oned
  
  ;;nofluxplan is for files that are to be reduced without fluxing
  
  
  
  if not keyword_set(basedir) then basedir = cwd()
  
  if not keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if size(rerun, /tname) NE 'STRING' then $
     rerun = string(rerun, format='(i4.4)')   
    
     if strmid(plugplan,0,1) ne '/' and strmid(plugplan,0,1) ne '.' $
        and strmid(plugplan,0,1) ne '~' then plugplan = cwd() + plugplan
     if keyword_set(redplan) then begin
         if strmid(redplan,0,1) ne '/' and strmid(redplan,0,1) ne '.' $
                and strmid(redplan,0,1) ne '~' then redplan = cwd() + redplan
     endif else if keyword_set(nofluxplan) then begin
        if strmid(nofluxplan,0,1) ne '/' and $
           strmid(nofluxplan,0,1) ne '.' $
        and strmid(nofluxplan,0,1) ne '~' then nofluxplan = cwd() + nofluxplan
     endif
  
  
  if keyword_set(redplan) then begin
     
     readcol, redplan, redir, format='A', comment='#'
     originaldir = cwd()
     for i = 0, n_elements(redir) -1 do begin
         cd, redir(i)
         hs_reduce2d, rerun=rerun, plugplan=plugplan, $
           checkaverage=checkaverage, docosmic=docosmic,$
           superfile=superfule
         cd, originaldir
         
     endfor
     
 endif
 
 if keyword_set(nofluxplan) then begin
     readcol, nofluxplan, redir, format='A', comment='#'
     originaldir = cwd()
     for i = 0, n_elements(redir) -1 do begin
         cd, redir(i)
         hs_reduce2d_noflux, rerun=rerun, plugplan=plugplan, $
           checkaverage=checkaverage, docosmic=docosmic
         cd, originaldir
     endfor
  endif
  
 
  ;;Now do the 1-D reductions
 
  if keyword_set(oned) then $
     hs_batch1d,   redplan=redplan, rerun=rerun, plugplan=plugplan,link=link, $
              basedir=basedir, checkaverage=checkaverage, docosmic=docosmic, $
              nofluxplan=nofluxplan, superfile=superfile
   
  
end

  
