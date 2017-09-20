;+
; NAME:
;   hs_findwave
;
; PURPOSE:
;   Wrapper for the wavelength solution 
;
; CALLING SEQUENCE:
;  hs_findwave, arcname=arcname, rerun=rerun, root=root, arcguess=arcguess
;
; INPUTS:
;   
;
; OPTIONAL KEYWORDS:
;   arcname    - The arc file to find the wavelenth solution from
;   rerun      - Rerun number for the reduction - in order to find
;                output directory
;   root       - Modified Julian date for the arc
;   arcguess   - input solution guess - I wouldn't use this
;
; OUTPUTS:
;     Creates a wset for the data
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
;  
; REVISION HISTORY:
;   Apr 2014 -Enabled and tested 600line grating wavelength calibrations
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------
PRO hs_findwave, arcname=arcname, rerun=rerun, root=root, arcguess=arcguess, $
                 guessfile=guessfile, chelle=chelle, doplot=doplot
  
  if not keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if NOT keyword_set(arcname) then begin
     file = findfile('*fit*')
     N = n_elements(file)/2
     header = headfits(file(N), exten=1,  /silent)
     if not keyword_set(root) then root = sxpar(header, 'ROOT')
       arcname = 'calibration/' + rerun + '/arc.fits'
  endif
  
  if file_test(arcname)  eq 0 then begin
     splog, 'This arcname does not exist'
     return
  endif
  
  if file_test('calibration/' + rerun + '/traceset.fits') eq 0L then begin
     splog, 'You have not yet created a trace set - EXITTING'
     return
  endif
  
;  IF keyword_set(guessfile) THEN BEGIN
;     arcguess = mrdfits(guessfile, 0)
;  ENDIF
  
  for ccdnum = 1, 2 do begin
     IF NOT keyword_set(chelle) THEN begin
        hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun
        sigma= 2.0
        proftype = 2.0
        highrej = 15
        lowrej = 15
        npoly = 1
        wfixed = [1,1,1]
        
        hs_proc, arcname, ccdnum, arcimg, arcivar, rerun=rerun, root=root, header=header
        hs_extract_image, arcimg, arcivar, xsol, sigma, flux, fluxivar, $
           proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
           npoly=npoly, relative=1
        if strtrim(sxpar(header, 'DISPERSE'),2) eq '600_gpm' then begin
           tiltpos=float(sxpar(header, 'TILTPOS'))
           centers=[4800,5300,5800,6300,6600,6800,7300,7800,7800]
           tilt_positions=[-5.36,-6.28,-7.2,-8.1,-8.676,-9.03,-9.95,-10.89,-11.83]
           match=where(abs(tiltpos-tilt_positions) eq min(abs(tiltpos-tilt_positions)))
           center600=(centers[match])[0]
           lamptype=strtrim(sxpar(header, 'COMPLAMP'),2)
        endif else lamptype='henear'       
     ENDIF ELSE BEGIN
        hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun                       
        hc_getspec, arcname, ccdnum, specflux=flux, specivar=fluxivar, $
           rerun=rerun
        lamptype = 'thar'
     ENDELSE
     
     hs_wavecal, flux, fluxivar, xpeak, $
        ypeak, wset, ncoeff=arccoeff, $
        ccdnum=ccdnum, arcguess=arcguess, $
        doplot=doplot, center600=center600, $
        chelle=chelle, guessfile=guessfile, lamptype=lamptype
     
     IF ccdnum EQ 1 then begin
        mwrfits, wset, 'calibration/'+rerun+'/wset.fits', /create
     endif
     if ccdnum EQ 2 then begin
        mwrfits, wset, 'calibration/'+rerun+'/wset.fits'
     endif
     
     ;;Also, make a map of location of lines
     xloc = trace_crude(flux, yset=ycen, nave=2, nmed=2, xstart=0, $
                        ystart=0, maxshifte=0.1d, maxshift0=0.1d)
     IF ccdnum EQ 1 then begin
        mwrfits, xloc, 'calibration/'+rerun+'/lineloc.fits', /create
     endif
     if ccdnum EQ 2 then begin
        mwrfits, xloc, 'calibration/'+rerun+'/lineloc.fits'
     endif
     
  endfor
 
  return
end

