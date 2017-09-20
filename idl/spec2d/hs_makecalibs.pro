;+
; NAME:
;   hs_makecalibs
;
; PURPOSE:
;   Complete more calibration steps
;
; CALLING SEQUENCE:
;  hs_biasgen, biasfiles, outdir=outdir, mjd=mjd
;
; INPUTS:
;   biasfiles   - a linst of bias files that are to be comined
;
; OPTIONAL KEYWORDS:
;   outdir     - Output directory for the final combined arc
;   mjd        - Integer part of the MJD for the observations
;                this will create the outname for the file if
;                 in the headers of the file
;
; OUTPUTS:
;      creates the /calibration/XXXX/bias-YYYY.fits file where
;      XXXX is the rerun number and YYYY is the MJD
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
;    Not sure the scaling for the images is appropriate
;    is the entrie chip being used or just the illuminated
;    region
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------
;;This program does all the extra little steps that are needed to calibrate the
;;data.  Note: hs_tweakwave may be needed if you want to provide an initial 
;;guess for the wavelength solution.

PRO hs_makecalibs, arcguess=arcguess, rerun=rerun, root=root, $
                   dowave=dowave, tweakwave=tweakwave
  
  if NOT keyword_set(path) then path=cwd()
  if NOT keyword_set(listpath) then listpath = path+'lists/'
  if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if size(rerun, /tname) ne 'STRING' then rerun = $
     string(rerun, format='(i4.4)')
  if NOT keyword_set(calpath) then calpath = path + $
     'calibration/' + rerun + '/'
  
  ;;Check to see if the calibration directory exists
  
  if direxist(path+'calibration') eq 0L then spawn,$
     'mkdir ' + path+'calibration'
  if direxist(path+'calibration/' + rerun) EQ 0L then $
     spawn, 'mkdir ' + path+'calibration/' + rerun
  
  ;;Create the biasfiles
  
  biasfile = findfile(listpath+'bias*')
  junk = strsplit(listpath, ' ', length = pathcount)
    
  ;;Create the traceset 
  if file_test('calibration/'+rerun+'/traceset.fits') eq 0L then begin
     
     ;;Make the traceset
     hs_maketraceset, /write, rerun=rerun, mjd=root
  endif
  
  
  ;;Create the wavelength solution
  
  if keyword_set(dowave) then begin
     
     if keyword_set(tweakwave) then begin
        hs_tweakwave, arcguess, root=root, rerun=rerun
     endif
     
     if file_test('calibration/'+rerun+'/'+'wset.fits') EQ 0L then begin
        hs_findwave, rerun=rerun, root=root, arcguess=arcguess
        
     endif
     
  endif
  
END
 
 
 
