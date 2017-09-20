;+
; NAME:
;      hs_remove_telluric
;
; PURPOSE:
;      Correct hectospec spectra for telluric absorption in the A- and B-bands
;
; CALLING SEQUENCE:
;      hs_remove_telluric, lambda,flux=flux,ivar=ivar, airmass=airmass
; 
; INPUTS:
;      lambda - wavelength array corresponding either to a single
;                     spectrum or and array of spectra. For arrays of
;                     spectra, both dimensions of lambda must match
;                     those of flux ([NLAMBDA,NSPEC])
;
; REQUIRED KEYWORDS:
;      flux - spectrum, or 2D array of spectra to be corrected
;               ([NLAMBDA,NSPEC]). Must have same length as
;               lambda array. Corrected spectrum is returned back by
;               modifying this keyword
; OPTIONAL KEYWORDS:
;      ivar - inverse variance array associated with the spectrum
;               if passed, the adjusted ivar array is returned through
;               this keyword
;      airmass - a scalar or [NSPEC] array giving the airmass of 
;                      observation for each spectrum to be corrected
;                      If omitted, airmass is assumed to be 1.0 and no
;                      telluric correction is performed
;
; OUTPUTS:
;      No return value. Modified flux and ivar are returned by 
;      modifying keywords
;
;
; MODIFICATION HISTORY:
;   smoran 06/2014 - full rewrite implementing Beer's law algorithm to 
;                                avoid using propietary code


pro hs_remove_telluric, lambda,flux=flux,ivar=ivar, airmass=airmass


if not keyword_set(flux) then begin
  splog, 'No flux passed to HS_REMOVE_TELLURIC: returning'
  return
endif   
if n_elements(flux) ne n_elements(lambda) then begin
  splog, 'Dimensions do not match: flux, lambda'
  return
endif
if not keyword_set(ivar) then ivar=flux*0.0+1.0
   nspec=(size(flux))[2]
   npixsp =(size(flux))[1]
if not keyword_set(airmass) then airmass=fltarr(nspec)+1.0 else $
   if n_elements(airmass) ne nspec then begin
     print, 'Number of airmass values does not match number of spectra'
     return
   endif

;;LOAD IN TEMPLATE SPECTRUM AND CONSTRUCT CORRECTIONS
;Read in template spectrum for A-band, and make a copy for B-band
a_lambda=mrdfits('$HSRED_DIR/etc/HD-192281.fits')
a_flux=mrdfits('$HSRED_DIR/etc/HD-192281.fits',1)
a_ivar=mrdfits('$HSRED_DIR/etc/HD-192281.fits',2)
b_lambda=a_lambda
b_flux=a_flux
b_ivar=a_ivar
ab_airmass = 1.0 ;appropriate for this template

;set some pixel ranges around the A-band
   aband_low = min(where(a_lambda gt 7575.))
   aband_hi = min(where(a_lambda gt 7710.))
   amin_sideband = min(where(a_lambda gt 7505.))
   amax_sideband = min(where(a_lambda gt 7780.))

;and similar ranges for the B-band
   bband_low = min(where(b_lambda gt 6862.))
   bband_hi = min(where(b_lambda gt 6920.))
   bmin_sideband = min(where(b_lambda gt 6762.))
   bmax_sideband = min(where(b_lambda gt 7020.))

;zero out bad areas as a way to mask for the smooth continuum fit
   bad=where((a_ivar eq 0.0) or finite(a_ivar,/nan),bad_cnt)
  if bad_cnt ge 1 then a_ivar[bad]=0.0
  a_ivar[aband_low:aband_hi]=0.0
  a_ivar[0:(amin_sideband-1)]=0.0
  a_ivar[(amax_sideband+1):*]=0.0

   bad=where((b_ivar eq 0.0) or finite(b_ivar,/nan),bad_cnt)
  if bad_cnt ge 1 then b_ivar[bad]=0.0
  b_ivar[bband_low:bband_hi]=0.0
  b_ivar[0:(bmin_sideband-1)]=0.0
  b_ivar[(bmax_sideband+1):*]=0.0


;;Finish constructing the A-band correction

;normalize the template spectrum, tweaked with a careful choice of 
; pixels 
  a_fluxn=float(a_flux)/(total(a_flux[aband_low-10:aband_low])/11.)
  a_ivarn=a_ivar*(total(a_flux[aband_low-10:aband_low])/11.)^2.
; fit the continuum shape outside of the A-band
; svdfit is used because it will ignore points where we have set ivar=0
   coeff = svdfit(a_lambda, a_fluxn, weights=a_ivarn,2, yfit=afit, /double) 

;now divide the shape out, and set everything outside the
;region of interest to 1.0 
   a_flux = a_fluxn/ afit 
   a_flux[where(a_flux gt 1.0)]=1.0
   a_flux[0:aband_low-1]=1.0
   a_flux[aband_hi+1:*]=1.0

   aband = {a_lambda:float(a_lambda[aband_low:aband_hi]), $
         a_flux:float(a_flux[aband_low:aband_hi]) }

;; Now finish constructing the B-band correction

    ; normalize around the B-band with a careful choice of pixels
   b_fluxn=float(b_flux)/((total(b_flux[bband_low-2:bband_low])/3.))
   b_ivarn=b_ivar*(total(b_flux[bband_low-2:bband_low])/3.)^2.

  ; fit the continuum shape outside of the B-band
  ; svdfit is used because it will ignore points where we have set ivar=0
   coeff = svdfit(b_lambda, b_fluxn,weights = b_ivarn, 2, yfit=bfit, /double) 

  ;now divide the shape out, and set everything outside the
  ;region of interest to 1.0 
   b_flux = b_fluxn/ bfit 
   b_flux[where(b_flux gt 1.0)]=1.0
   b_flux[0:bband_low-1]=1.0
   b_flux[bband_hi+1:*]=1.0

   bband = {b_lambda:float([b_lambda[bband_low:bband_hi]]), $
         b_flux:float([b_flux[bband_low:bband_hi]]) }

;; Now we make the corrections by interpolating the template
;; onto the real spectral wavlength grid, one spectrum at a time

; initialize spline interpolation
   a_spline = spl_init(aband.a_lambda, aband.a_flux)
   b_spline = spl_init(bband.b_lambda, bband.b_flux)

;loop through the individual fibers
for i=0, nspec-1 do begin
  lambda1=lambda[*,i]
  flux1=flux[*,i]
  ivar1=ivar[*,i]

   ;Do A-band
   arange = where(lambda1 ge aband.a_lambda[0] and $
                   lambda1 le max(aband.a_lambda), acnt)
   if acnt eq 0 then begin
       splog, 'A-band does not fall in spectral range' 
   endif else begin 
       ; interpolate A band template to same wavelength array as input
       ; spectrum and scale according to airmass (Beer's law)
           aband_flux = spl_interp(aband.a_lambda, aband.a_flux, a_spline, $
                          lambda1[arange])^((airmass[i]/ab_airmass)^0.55)

        ; correct flux and ivar
           flux[arange,i] = flux1[arange]/aband_flux 
           ivar[arange,i] = ivar1[arange]*aband_flux^2 
           
   endelse

   ;Now B-band
   brange = where(lambda1 ge bband.b_lambda[0] and $
                   lambda1 le max(bband.b_lambda), bcnt)
   if bcnt eq 0 then begin
       splog, 'B-band does not fall in spectral range' 
   endif else begin 
       ; interpolate A band template to same wavelength array as input
       ; spectrum and scale according to airmass (Beer's law)
           bband_flux = spl_interp(bband.b_lambda, bband.b_flux, b_spline, $
                          lambda1[brange])^((airmass[i]/ab_airmass)^0.55)

        ; correct flux and ivar
           flux[brange,i] = flux1[brange]/bband_flux 
           ivar[brange,i] = ivar1[brange]*bband_flux^2 
           
   endelse

  endfor
  return
end
