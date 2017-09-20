;+
; NAME:
;   hs_readtraceset
;
; PURPOSE:
;  Read in the traceset
;
; CALLING SEQUENCE:
;  hs_readtraceset, ccdnum, xsol, rerun, 
;
; INPUTS:
;  ccdnum - CCD number
;
; OPTIONAL KEYWORDS:
;   rerun - reduction rerun
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   xsol - traceset xsol
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

;------------------------------------------------------------------------------
;;This will read a hs_traceset file

PRO hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun, filename=filename
  
  if not keyword_Set(rerun) then rerun = '0100'
  
  if not keyword_set(filename) then filename='calibration/' + rerun + '/traceset.fits'
  if file_test(filename) eq 0l then begin
     splog, 'Traceset has not been create yet - ENDING'
     return
  endif
  
  xsol = mrdfits(filename, ccdnum-1l, /silent)
    
  return
end

  
  
