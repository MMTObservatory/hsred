;+
; NAME:
;   hs_calibproc
;
; PURPOSE:
;   Do the intial processing of calibration frames
;
; CALLING SEQUENCE:
;  hs_calibproc, listpath=listpath, calpath=calpath,$
;                  mjd=mjd, rerun=rerun, wavetweak=wavetweak, path=path, $
;                  doall=doall, dobias=dobias, dodome=dodome, dosky=dosky,$
;                  doarc=doarc, dowave=dowave, doflatten=doflatten, qaplot=qaplot
;
;
; INPUTS:
;   
;
; OPTIONAL KEYWORDS:
;   listpath      - Directory containing lists (default to 'lists/')
;   calpath       - Where calibration will be written (default to 'lists/')
;   mjd           - integer modified julian date for the observations. If 
;                   not set, it will be found in the header
;   dobias        - DO the bias combination
;   dodome        - combine the domeflats
;   dosky         - combine the skyflats
;   doarc         - combine the arcfiles
;   dowave        - find the wavelength solution
;   doall         - sets all of the above processes
;   qaplot       -turn on/off plotting of wavelength solution diagnostics
;
; OUTPUTS:
;     dobias - same as hs_biasgen
;     doarc  - same as hs_arcgen
;     dodome - same as hs_dflatgen
;     dosky  - same as hs_sflatgen
;     dowave - same as hs_maketraceset, hs_findwave
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------
PRO hs_calibproc, listpath=listpath, calpath=calpath,$
                  mjd=mjd, rerun=rerun, wavetweak=wavetweak, path=path, $
                  doall=doall, dobias=dobias, dodome=dodome, dosky=dosky,$
                  doarc=doarc, dowave=dowave, $
                  dodark=dodark, chelle=chelle, qaplot=qaplot
  
  if keyword_set(doall) then begin
     hs_calibproc, listpath=listpath, calpath=calpath, mjd=mjd,$
        rerun=rerun, wavetweak=wavetweak, path=path, $
        /dobias, /dodome, /dosky, /doarc, /dowave, dodark=dodark, qaplot=qaplot
       return
    endif
        
  
    ;;This module does all the combining/processing of calibration frames
    
    if NOT keyword_set(path) then path=cwd()
    if NOT keyword_set(listpath) then listpath = path+'lists/'
    if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)') else $
       rerun = string(rerun, format='(i4.4)')
    
    if NOT keyword_set(calpath) then calpath = path $
       + 'calibration/' + rerun + '/'
    
    if size(rerun, /tname) ne 'STRING' then $
       rerun = string(rerun, format='(i4.4)')
    
    
    ;;Check to see if the calibration directory exists
    
    if direxist(path+'calibration') eq 0L then $
       spawn, 'mkdir ' + path+'calibration'
    if direxist(path+'calibration/' + rerun) EQ 0L then $
       spawn, 'mkdir ' + path+'calibration/' + rerun
    
    
    
;;Create the biasfiles
    
    biasfile = findfile(listpath+'bias*')
    junk = strsplit(listpath, ' ', length = pathcount)
     
    if keyword_set(dobias) then begin
       if file_test(calpath + 'bias.fits') EQ 0L then begin
          readcol, listpath + 'bias.list', file, $
             format='A', comment='#'
          hs_biasgen, file, outdir=calpath, mjd=mjd
       endif
    endif
    
    if keyword_set(dodark) then begin
       if file_test(calpath + 'dark.fits') EQ 0L then begin
          readcol, listpath + 'dark.list', file, $
             format='A', comment='#'
          hs_darkgen, file, outdir=calpath, mjd=mjd
       endif
    endif
    
    
    if keyword_set(dodome) then begin
       ;;Create the dflatfiles
       dflatfile = findfile(listpath+'dflat*')
       junk = strsplit(listpath, ' ', length = pathcount)
       if file_test(calpath + 'dflat.fits') EQ 0L then begin
          readcol, listpath + 'dflat.list', $
             file, junk, format='A,A', comment='#'
          hs_dflatgen, file,  outdir=calpath, mjd=mjd
       endif
    endif
    
    
    if keyword_set(dosky) AND $
       file_test('lists/sflat.list') EQ 1l then begin
       ;;Create the sflatfiles
       
       sflatfile = findfile(listpath+'sflat*')
       junk = strsplit(listpath, ' ', length = pathcount)
       
       if file_test(calpath + 'sflat.fits') EQ 0L then begin
          readcol,  listpath + 'sflat.list', $
             file, junk, format='A,A', comment='#'
          hs_sflatgen, file,  outdir=calpath, mjd=mjd
       endif
    endif
    
    if keyword_set(doarc) then begin
       
       ;;Create the arcfiles
       arcfile = findfile(listpath+'arc*')
       junk = strsplit(listpath, ' ', length = pathcount)
       if file_test(calpath + 'arc.fits') EQ 0L then begin
          readcol,  listpath + 'arc.list', file, $
             junk, format='A,A', comment='#'
          hs_arcgen, file,  outdir=calpath, mjd=mjd
       endif
       
    endif
    
    if keyword_set(dowave) then begin
       
       if file_test('calibration/'+rerun+'/traceset.fits') eq 0L then begin
          ;;Make the traceset
          hs_maketraceset, /write, rerun=rerun, mjd=root
       endif
       
       
       if file_test('calibration/'+rerun+'/'+'wset.fits') EQ 0L then begin
          if keyword_set(wavetweak) then hs_tweakwave, arcguess, $
             rerun=rerun, root='' 
          hs_findwave, rerun=rerun, root=mjd, arcguess=arcguess, doplot=qaplot
          
       endif
    endif
    
 
 
  
  
  
  return
  
END
