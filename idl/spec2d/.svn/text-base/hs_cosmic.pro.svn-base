;+
; NAME:
;   hs_cosmic
;
; PURPOSE:
;   Calls the cosmic ray rejection routines
;
; CALLING SEQUENCE:
;  hs_cosmic, flux, newflux, sigclip=sigclip, objlim=objlim
;       sigfrac=sigfrac, blocksize=blocksize, rerun=rerun
;       nodelete=nodelete
;
; INPUTS:
;   flux       - input image
;
; OPTIONAL KEYWORDS:
;   sigclip - sigclip setting for LAcosmic (default to 5.0)
;   objlim  - objlim setting for LAcosmic (default to 2.0)
;   sigfrac - sigfrac setting for LAcosmic (default to 0.5)
;   blocksize - subimage size to use in LAcosmic (default to 1024)
;   rerun - rerun for the reduction (default to 0)
;
; OUTPUTS:
;      newflux  - CR rejected image
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;    This is quite slow.  
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
PRO hs_cosmic, flux, newflux,  sigclip=sigclip, $
               objlim=objlim, sigfrac=sigfrac, $
               blocksize=blocksize, rerun=rerun, $
               nodelete=nodelete
  
  
  if NOT keyword_set(sigclip) then sigclip =6.0
  if NOT keyword_set(objlim) then objlim = 2.0
  if NOT keyword_set(sigfrac) then sigfrac = 0.5
  if NOT keyword_set(blocksize) then blocksize= 1024
  if NOT keyword_set(rerun) then rerun=00
  if size(rerun, /tname) NE 'STRING' then $
     rerun = string(rerun, format='(i4.4)')
seed=systime(/sec)
filenum=string(fix(randomu(seed)*10000),format='(I4)')
  ;;In order to get the interpolations correct, 
  ;;you need to have the spectral direction in the X
  ;; Flip the image to do this
  
  tmp = transpose(flux*0.0)
  
  for i = 0, n_elements(tmp(*,0)) -1 do begin
     
     tmp(i,*) = flux(*,i) 
     
  endfor
  
  ;;Write this to a temp file
  mwrfits, tmp, 'calibration/' + rerun + '/tmp'+filenum+'-cr.fits', /create
  
  
  ;;Call the cosmic ray rejection
  
  hs_lacosmic,  'calibration/' + rerun + '/tmp'+filenum+'-cr.fits', $
     sigclip=sigclip, objlim=objlim, sigfrac=sigfrac, $
     blocksize=blocksize, /isbig, /verbose, gain=1, readn=2.8
  
  
  ;;Now read in the file
  
  newtmp = mrdfits('calibration/' + rerun + '/tmp'+filenum+'-cr-out.fits',0)
  
  newflux = transpose(newtmp*0.0)
  
   for i = 0, n_elements(newtmp(*,0)) -1 do begin
     
     newflux(*,i) = newtmp(i,*) 
     
  endfor
  
  if NOT keyword_set(nodelete) then begin
     
     ;;Now delete the temp files
     
     rmfile, 'calibration/' + rerun + '/tmp'+filenum+'-cr.fits'
     rmfile, 'calibration/' + rerun + '/tmp'+filenum+'-cr-out.fits'
     rmfile, 'calibration/' + rerun + '/tmp'+filenum+'-cr-mask.fits'

  endif
  
END
