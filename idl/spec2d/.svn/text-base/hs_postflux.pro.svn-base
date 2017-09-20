Pro hs_postflux, image, outimage=outimage, useaverage=useaverage, $
                 superfile=superfile
  
  if not keyword_set(outimage) then outimage = 'F'+image
  outname = outimage
  
  lam = mrdfits(image, 0)
  flux = mrdfits(image, 1)
  ivar = mrdfits(image, 2)
  plugmap = mrdfits(image, 5)
  mask = mrdfits(image, 3)
  
  npix = n_elements(lam[*,0])
  nfib = n_elements(lam[0,*])
  
  xxtemp = indgen(4608,300)
  for ifib = 0, 299 do begin 
     xxtemp(*,ifib) = xxtemp(*,0) 
  endfor
  xy2traceset, xxtemp, lam[*,*], $ 
     tempset, ncoeff=6
  
  telluric_factor = hs_telluric_corr(flux, ivar, tempset,$
                                     plugmap)
  divideflat, flux, telluric_factor, invvar = ivar
  
  if keyword_set(useaverage) then input_calset =$
     mrdfits( getenv('HSRED_DIR')+ '/etc/average_flux.fits',1)
  
  if keyword_set(superfile) then superfile = getenv('HSRED_DIR') + $
     '/etc/superstand.fits'
  
  if keyword_set(outname) then begin
     fcalfile = strarr(n_elements(outname))
     scalfile = fcalfile(0)
     
     for i = 0, n_elements(outname) -1 do begin
        test = strsplit(outname(i), '/,.', /extract)
        ct = n_elements(test)
        base = test(ct-2)
        outname1 = strsplit(outname(i), '/,.fits', /extract)
        ct = n_elements(outname1)
        outname1= outname1(ct-1)
        
        fcalfile(i) = 'hsFluxcalib-' + outname
        if i eq 0 then scalfile = 'hsStd-' + outname
        
     endfor
  endif
  
  k = where(strmatch(plugmap.objtype, 'SPE*'))
  vaclam = lam
  airtovac, vaclam
  
  
  calset = mrdfits(getenv('HSRED_DIR') +$
                   '/etc/average_flux.fits',1, /silent)
  
  
  hs_sphoto_calib, alog10(vaclam(*,k)), flux(*,k), ivar(*,k), $
     mask(*,k), plugmap[k], fcalfile, scalfile,  superfile=superfile, /stype,$
     input_calibset=input_calset
  
  repeatgo = 'false'
  fcal = mrdfits(fcalfile(0),1)
  flux_factor = bspline_valu(alog10(lam[*,*]), fcal)
  divideflat, flux, flux_factor,$
     invvar=ivar, minval=0.05*median(flux_factor[*,0])
  
  
  mwrfits, lam, outimage, /create
  mwrfits, flux,outimage
  mwrfits, ivar, outimage
  mwrfits, mask, outimage
  mwrfits, mask, outimage
  mwrfits, plugmap, outimage
  
  
end


