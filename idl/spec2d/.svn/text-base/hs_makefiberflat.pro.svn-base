;+
; NAME:
;   hs_makefiberflat
;
; PURPOSE:
;   Create a flat field vector for a secondary correction 
;   preferably from the sky flats
;
; CALLING SEQUENCE:
;  hs_makefiberflat, usesky=usesky, flatname=flatname, fflat=fflat,
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
;    fiberflat.fits in the calibration/XXXX directory
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
PRO hs_makefiberflat, usesky=usesky, flatname=flatfile, dofringes=dofringes, $
                      fflat=fflat, rerun=rerun, root=root, chelle=chelle
  
  
  if not keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if size(rerun, /tname) ne 'STRING' then rerun =$
     string(rerun, format='(i4.4)')
  
  
 if NOT keyword_set(flatfile) then begin
    file = findfile('*fit*')
    N = n_elements(file)/2
    header = headfits(file(N), exten=1,  /silent)
    if not keyword_set(root) then root = sxpar(header, 'mjd')
    flatfile = 'calibration/' +rerun + '/dflat.fits'
    if keyword_set(usesky) then flatfile = 'calibration/' + $
       rerun + '/sflat.fits'
 endif
 
 for ccdnum = 1, 2 do begin
    
    IF NOT keyword_set(chelle) THEN begin
       hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun
       sigma= 2.0
       proftype = 2.0
       highrej = 15
       lowrej = 15
       npoly = 1
       wfixed = [1,1,1]
       
       hs_proc, flatfile, ccdnum, flat, flativar,rerun=rerun, root=root
       
       hs_extract_image, flat, flativar, xsol, sigma, domeflux, domefluxivar, $
          proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
          npoly=npoly, relative=1
    ENDIF ELSE BEGIN
       hc_getspec, flatfile, ccdnum, rerun=rerun, $
          specflux=domeflux, specivar=domefluxivar
    endelse
    
    
    hs_readwset, ccdnum, wset=wset,rerun=rerun
    
    sflat = mrdfits('calibration/' + rerun+'/traceflat.fits',ccdnum-1)
    ;if we made a fringing correction, must apply it to the domeflat here
    ;before dividing it from the sflat
    if file_test('calibration/'+rerun+'/'+'fringes_pattern.fits') and keyword_set(dofringes) then begin 
             fringes = mrdfits('calibration/'+rerun+'/'+ $
                           'fringes_pattern.fits',ccdnum-1, /silent)            
             tmpsmooth=sflat*0.0
             nsmooth=n_elements(sflat[*,0])/50
             for q=0, 149 do tmpsmooth[*,q]=djs_median(domeflux[*,q],width=nsmooth, boundary='reflect') 
             domeflux=domeflux-fringes*tmpsmooth   
             for q=0, 149 do tmpsmooth[*,q]=djs_median(sflat[*,q],width=nsmooth, boundary='reflect') 
             sflat=sflat-fringes*tmpsmooth   
    endif
 divideflat, domeflux, sflat, invvar=domefluxivar
    
    
    ;;I think 6 coeffs may be too many - too many wiggles in sflat
    ;;Trying 4 
    ;;fflat = fiberflat(domeflux,domefluxivar, wset , ncoeff=6)
    fflat = fiberflat(domeflux, domefluxivar, wset, ncoeff=4)
    if ccdnum eq 1 then begin $
       mwrfits, fflat, 'calibration/'+rerun+'/fiberflat.fits', /create
    endif else if ccdnum eq 2 then begin $
       mwrfits, fflat, 'calibration/'+rerun+'/fiberflat.fits'
    endif
    
 endfor
 
 return
 
END
