;+
; NAME:
;   hs_readwset
;
; PURPOSE:
;   Read in a wset file
;
; CALLING SEQUENCE:
;  hs_readwset,ccdnum, wset=, rerun=, 

; INPUTS:
;  ccdnum - side of chip to read
;
; OPTIONAL KEYWORDS:
;  wset = wavelngth solution
;  rerun  = reduction rerun
;
; OUTPUTS:
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

;------------------------------------------------------------------------------

PRO hs_readwset, ccdnum, wset=wset,rerun=rerun
  
  
  wset = mrdfits('calibration/' + rerun + '/wset.fits', ccdnum)
  
  return
END
