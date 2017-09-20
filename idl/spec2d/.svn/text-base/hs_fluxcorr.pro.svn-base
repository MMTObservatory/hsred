;+
; NAME:
;   hs_fluxcorr
;
; PURPOSE:
;   Find the fluxing vector for data with F stars
;
; CALLING SEQUENCE:
;   flux_out = hs_fluxcorr, lam1, flux1, ivar1, plugmap1, mask1,
;           outname=outname,$
;            plateid=plateid, rerun=rerun, psplot=psplot,$
;            checkaverage=checkaverage, stand=stand, $
;            header=header, cut=cut, superfile=superfile
;
; INPUTS:
;   lam1 = [Nim, Npix, Nfib] wavelength
;   flux1 = flux [Nim, Npix, Nfib]
;   ivar1 = ivar [Nim, Npix, Nfib]
;   plugmap1 = plugmap structure for the data (Nfiber)
;
;
; OPTIONAL KEYWORDS:
;   rerun   - rerun reductions for the data
;   checkaverage - if set, then truncate the spectrum where the 
;                  standard stars have fewer than 100 counts.
;   superfile    - if set use provided spectral types for the stars
;   header       - object hearder
;
; OUTPUTS:
;     fluxout - Fluxing vectors for the data
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------


function hs_fluxcorr, lam1, flux1, ivar1, plugmap1, mask1, outname=outname,$
                      plateid=plateid, rerun=rerun, psplot=psplot,$
                      checkaverage=checkaverage, stand=stand, $
                      header=header, cut=cut, superfile=superfile, $
                      useaverage=useaverage, fcalfile=fcalfile, $
                      scalfile=scalfile
  lam = lam1
  flux = flux1
  ivar = ivar1
  plugmap = plugmap1
  mask = mask1
  
  if keyword_set(useaverage) then $
     input_calset = mrdfits( getenv('HSRED_DIR')+ $
                             '/etc/average_flux.fits',1)
  
  if keyword_set(superfile) then superfile = getenv('HSRED_DIR') + $
     '/etc/superstand.fits'
    
  if not keyword_set(rerun) then rerun=0
  if not keyword_set(plateid) then plateid = 0
  
  if size(rerun, /tname) NE 'STRING' then $
     rerun = string(rerun, format='(i4.4)')     
  if size(plateid, /tname) NE 'STRING' then $
     plateid = string(plateid, format='(i5.5)')   
  if keyword_set(outname) AND NOT keyword_set(fcalfile) then begin
     
     fcalfile = strarr(n_elements(outname))
     scalfile = fcalfile(0)
     
     for i = 0, n_elements(outname) -1 do begin
        test = strsplit(outname(i), '/,.', /extract)
        ct = n_elements(test)
        base = test(ct-2)
        outname1 = strsplit(outname(i), '/,.fits', /extract)
        ct = n_elements(outname1)
        outname1= outname1(ct-1)
        
        fcalfile(i) = 'reduction/' + rerun +'/hsFluxcalib-' + $
           outname1 + '-' + plateid + '-' + $
           rerun +'.fits'
        if i eq 0 then scalfile = 'reduction/' + rerun + '/hsStd-' + outname1 $
           + '-' + plateid + '-' + $
           rerun+'.fits'
     endfor
     
  endif
  
  if not keyword_set(outname) then outname = 'test-0-0-0' 
  k = where(strmatch(plugmap.objtype, 'SPE*'))
  
  ;;Get the quantities in the format you want
  
  ndim = size(lam, /n_dimensions)
  lamorig = lam
  airtovac, lamorig
  
  if ndim eq 3 then begin
     
     size1 = size(lam)
     nim = size1(1)
     npix = size1(2)
     nfib = size1(3)
     newlam = dblarr(npix, n_elements(k)*nim)
     newflux = newlam
     newivar = newlam
     newmask = newlam
     count = 0
     nk = n_elements(k)
     for i = 0, nim -1 do begin
        
        newlam(*,count:count+(nk-1L))  =  transpose(transpose(lam(i,  *, k)))
        newflux(*,count:count+(nk-1L)) =  transpose(transpose(flux(i, *, k)))
        newivar(*,count:count+(nk-1L)) =  transpose(transpose(ivar(i, *, k)))
        newmask(*,count:count+(nk-1L)) =  transpose(transpose(mask(i, *, k)))
        plugmap.frames = i
        
        if i eq 0 then newplug = plugmap[k]
        if i ne 0 then newplug = [newplug, plugmap[k]]
        count = nk*(i+1)
        
     endfor
     
     lam = newlam
     flux = newflux
     ivar = newivar
     mask = newmask
     plugmap = newplug
     
  endif
  
  
  vaclam = lam
  airtovac, vaclam
  calset = mrdfits(getenv('HSRED_DIR') +$
                   '/etc/average_flux.fits',1, /silent)
  
  hs_sphoto_calib, alog10(vaclam), flux, ivar, $
     mask, plugmap, fcalfile, scalfile,  superfile=superfile, /stype,$
     input_calibset=input_calset
  repeatgo = 'false'
  
  for i = 0, n_elements(outname) -1 do begin
     
     fcal = mrdfits(fcalfile(i),1)
     flux_factor = bspline_valu(alog10(lamorig[i,*,*]), fcal)
     
     if keyword_set(checkaverage)  and NOT keyword_set(cut) then begin
        
        average_fluxfile = getenv('HSRED_DIR') + '/etc/average_flux.fits'
        acal = mrdfits(average_fluxfile,1, /silent)
        average_factor = bspline_valu(alog10(lamorig), acal)
        
        ratio = flux_factor/average_factor
        
        k = where(ratio LE 0.05, count)
        
        if count gt 0 then begin
           fluxtemp=flux*0.0
           npix = n_elements(flux(*,0))
           fluxtemp(1:npix-7,*) = flux(0:npix-8,*)
           
           k = where(strmatch(plugmap.objtype, 'SPE*'))
           val = fltarr(n_elements(k))
           
           for j =0, n_elements(k) -1 do begin 
              t = where(djs_median(flux[*,k(j)],width=200, $
                                   boundary='nearest') gt 100 and $
                        djs_median(fluxtemp[*,k(j)],width=200,$
                                   boundary='nearest') gt 100)
              val(i) = t(0)
           endfor
           
           djs_iterstat, val, median=cut
           
           splog, 'Bad flux solution. Truncating the spectrum to begin at '
           splog, 'pixel ' + string(cut)
           splog, 'Corresponding to ~ ' +string( lam[cut,k(0)])
           
           repeatgo='true'
        endif
     endif
     
     if i eq 0 then flux_out = dblarr(nim, npix, nfib)
     flux_out[i,*,*] = flux_factor
     
  endfor
  
  if repeatgo eq 'true' then begin
     
     lam2 = fltarr(nim, npix-cut, nfib)
     flux2 = lam2
     ivar2 = lam2
     mask2 = lam2 
     
     for im = 0, nim -1 do begin
        lam2[im,*,*] = lam1[im,cut:npix-1,*]
        flux2[im,*,*] = flux1[im,cut:npix-1,*]
        ivar2[im,*,*] = ivar1[im,cut:npix-1,*]
        mask2[im,*,*] = mask1[im,cut:npix-1,*]
     endfor
     flux_out = hs_fluxcorr(lam2, flux2, ivar2, plugmap1, mask2, $
                            outname=outname, plateid=plateid, $
                            rerun=rerun, /psplot, header=objhdr, $
                            cut=cut, superfile=superfile,$
                            useaverage=useaverage)
     
  endif
  
  return, flux_out
  
END
