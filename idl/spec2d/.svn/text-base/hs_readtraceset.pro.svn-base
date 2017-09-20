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

PRO hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun
  
  if not keyword_Set(rerun) then rerun = string(0, format='(i4.4)')
  
  if file_test('calibration/' + rerun + '/traceset.fits') eq 0l then begin
     splog, 'Traceset has not been create yet - ENDING'
     return
  endif
  
  xsol = mrdfits('calibration/' + rerun + '/traceset.fits', ccdnum-1l, /silent)
    
  return
end

  
  
