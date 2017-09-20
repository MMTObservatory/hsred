;+
; NAME:
;   hs_tweaksky
;
; PURPOSE:
;   Tweak the wavelength solution based on sky lines
;
; CALLING SEQUENCE:
;  hs_tweaksky, flux, fluxivar, wset, ccdnum=, ariset=, xsky=, rerun=, 
;               filename=, skywaves=, arcshift=, 
; INPUTS:
;   flux = flux [npix, nfib]
;   fluxivar = inverse varriance [npix, nfib]
;   wset  = input wavelenth solution
;
; OPTIONAL KEYWORDS:
;  xsky  - location of skylines
;  ccdnum - side of ccd working with
;  filename  - output file name for tweaks
;  skywaves - wavelength of skylines
;
;OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;  airset = output solution set
;  skyshift - shifts needed to be applied to the skylines
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
;- April 2014 - Modified to use ncoeff=6 in reconstructing the
;  wavelength solution
;------------------------------------------------------------------------------
PRO hs_tweaksky, flux, fluxivar, wset, ccdnum=ccdnum, $
                 airset=airset, xsky=xsky, rerun=rerun, $
                 filename=filename, skywaves=skywaves, arcshift=arcshift
  
  if keyword_set(filename) then begin
     filename1 = strsplit(filename, '/', /extract)
     fullroot = filename1[n_elements(filename1)-1]
     junk = strsplit(fullroot, ' ', length=ct)
     junk = strsplit(fullroot, '.', /extract)
     if junk[n_elements(junk)-1] eq 'gz' then begin
        root = strmid(fullroot, 0, ct-8)
     endif else begin
        root = strmid(fullroot, 0, ct-5)
     endelse
     filename1 = root
     
  endif
  
  if NOT keyword_set(filename) then filename1 = 'test'
  
 
  
  if NOT keyword_set(ccdnum) then begin
     splog, 'You must set the ccdnum keyword, please'
     return
  endif
  
  if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  
  filename1 = 'reduction/' + rerun + '/' + filename1 + '-' +$
     rerun + '-tweaks.fits'
  
  if file_test('calibration/'+rerun+'/lineloc.fits') eq 0L then begin
     splog, 'The wavelength calibrator has not been run as'
     splog, 'lineloc.fits is not present'
     return
  endif
  
  
   xcen = mrdfits('calibration/'+rerun+'/lineloc.fits', ccdnum - 1)
   hs_readwset, ccdnum, wset=wset, rerun=rerun
     
  
   hs_locateskylines, $
       getenv("HSRED_DIR") + '/etc/skylines.dat',$
      flux, fluxivar, wset, xcen, arcshift=arcshift, $
      xsky=xsky, skywaves=skywaves
   
   mwrfits, arcshift, filename1, /create 
   
   

   traceset2xy, wset,  transpose(xcen), lambda

    
   xy2traceset, transpose(double(xcen+arcshift)), $
          (lambda) ,  airset, ncoeff= 6    , xmin= 0  , xmax= 4607
     
   
END
