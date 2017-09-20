;+
; NAME:
;   hs_makesuperflat
;
; PURPOSE:
;   Create a flat field vector for a primary correction 
;   preferably from the dome flats
;
; CALLING SEQUENCE:
;  hs_makesuperflat, usesky=usesky, flatname=flatname, fflat=fflat,
;                    rerun=rerun, root=root
;
; INPUTS:
;   
;
; OPTIONAL KEYWORDS:
;   usesky - use the sky flat rather than the dome flats for the 
;            correction
;   flatfile - file to make the flat vectors from
;   rerun - reduction rerun
;   root  - integer modified julian date for the observation
;
; OUTPUTS:
;    traceflat.fits in the calibration/XXXX directory
;
; OPTIONAL OUTPUTS:
;    fflat  - flat fielding vector for the data
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
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------
PRO hs_makesuperflat,  sflat=sflat, usesky=usesky, $
                       flatname=flatname, rerun=rerun, root=root, chelle=chelle
  
  if not keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if size(rerun, /tname) ne 'STRING' then $
     rerun = string(rerun, format='(i4.4)')
  
  
  if NOT keyword_set(flatfile) then begin
     file = findfile('*fit*')
     N = n_elements(file)/2
     header = headfits(file(N), exten=1,  /silent)
     if NOT keyword_set(root) then root = sxpar(header, 'mjd')
     flatfile = 'calibration/' + rerun + $
        '/dflat.fits'
     if keyword_set(usesky) then flatfile = 'calibration/' + $
        rerun + '/sflat.fits'
  endif
  
  
  for ccdnum = 1, 2 do begin
     
     IF NOT keyword_set(chelle) THEN begin
        hs_readtraceset, ccdnum, xsol=xsol,  rerun=rerun
        sigma= 2.0
        proftype = 2.0
        highrej = 15
        lowrej = 15
        npoly = 1
        wfixed = [1,1,1]
        
        hs_proc, flatfile, ccdnum, flat, flativar, rerun=rerun, root=mj
        
        hs_extract_image, flat, flativar, xsol, sigma, domeflux, $
           domefluxivar, $
           proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
           npoly=npoly, relative=1
     ENDIF ELSE BEGIN
        hc_getspec, flatfile, ccdnum, rerun=rerun, $
           specflux=domeflux, specivar=domefluxivar
     endELSE
     
     hs_readwset, ccdnum, wset=wset, rerun=rerun
     
     sset = superflat(domeflux, domefluxivar, wset)
     smooth_fit = hs_smooth_superflat(sset, wset)
          
     sflat = domeflux/smooth_fit
     for i = 0, n_elements(sflat(0,*)) -1 do begin 
        sflat(*,i) = sflat(*,i) / djs_median(sflat(*,i), width=400, boundary='reflect') 
     endfor
     
     fflat = fiberflat(domeflux/sflat,domefluxivar*sflat^2, wset , ncoeff=6)
     sflat = sflat*fflat
     
     if ccdnum eq 1 then begin $
        mwrfits, sflat, 'calibration/' + rerun + '/traceflat.fits', /create
     endif else if ccdnum eq 2 then begin $
        mwrfits, sflat, 'calibration/' + rerun + '/traceflat.fits'
     endif
     
  endfor
  
  
END
