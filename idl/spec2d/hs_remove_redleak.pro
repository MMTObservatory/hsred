;Routine to use pre-computed model to subtract red leak contamination,
;based on position of this fiber in the focal plane. Keyword to also
;remove red emission at 8400A from MMTCAM does not follow the same
;pattern, and so this keyword implementation is broken at this time. 
;
pro hs_remove_redleak, airset, plugmap, exptime, specflux=specflux, invvar=specfluxivar, mmtcam=mmtcam, redflux=redflux
;get the lambda array
traceset2xy, airset, xx, lambda
;get the fiber positions in r,phi coords
fiber_x=plugmap.xfocal
fiber_y=plugmap.yfocal
fiber_num=plugmap.fiberid
xyz_to_angles,fiber_x,fiber_y,fltarr(n_elements(fiber_x)),fiber_radius,fiber_phi,theta2
;read in model normalized profile for red leak
if mmtcam eq 1 then norm_model=mrdfits(getenv('HSRED_DIR')+'/etc/redleak_normalized_shape.mmtcamOFF.fits',1) else $
if mmtcam eq 2 then norm_model=mrdfits(getenv('HSRED_DIR')+'/etc/redleak_normalized_shape.mmtcamON.fits',1) else $
norm_model=mrdfits(getenv('HSRED_DIR')+'/etc/redleak_normalized_shape.2014.fits',1)
scales=fltarr(n_elements(fiber_x))
red_flux_all=specflux*0.0
specflux_orig=specflux
bestvals=scales*0.0
for i=0, n_elements(scales)-1 do begin
    scale_model=mrdfits(getenv('HSRED_DIR')+'/etc/redleak_perfiber_model.fits',fix(fiber_num[i]), /silent)
   scales[i]=10.^(bspline_valu(fiber_radius[i], scale_model))
   yy=where(lambda[*,i] lt 8335. and lambda[*,i] ge 8315.)
   objcounts=median(specflux[yy,i]) ;subtract off a crude guess at the object continuum
   ww=where(lambda[*,i] gt 8322.)
   zz=where((lambda[*,i] gt 9010.) and (lambda[*,i] lt 9150.))

   ;try the same technique as hs_skyfit
     scaleval = findgen(240)/240.*33.+0.01
     sigma = scaleval*0.0
     for j = 0, n_elements(scaleval)-1 do begin
        redflux=bspline_valu(lambda[zz,i], norm_model)*scales[i]*exptime*scaleval[j]
        test = specflux[zz,i]-objcounts-redflux
        rval=correlate(test,redflux)
        sigma[j] = rval
     endfor
     ;passing in a sky spectrum as both 'flux' and 'sky' results in
     ;this bug, as the sky is precisely subtracted and we correlate against
     ;an array of zeros. Setting this nan to 0.0 is what we want
     bad=where(finite(sigma,/nan),cnt)
     if cnt ge 1 then sigma[bad]=0.0
     x =  findgen(6000)/6000.*33.+0.01  
     sigma1 = interpol(sigma, scaleval, x)
     bestvals[i] = min(sigma1^2.0, m)
     ;print, bestvals[i],m
     ;if we didn't find a local minimum, just skip this fiber
     if (m ne n_elements(x)-1) then scales[i]=scales[i]*x[m]

   red_flux=bspline_valu(lambda[ww,i], norm_model)*scales[i]*exptime
   nn=where(red_flux gt max(specflux[ww,i]),/null) ;cap at max level in spectrum
;   print, max(red_flux)
   red_flux[nn]=max(specflux[ww,i])
;   if (i eq 4) then begin
;      set_plot, 'x'
;      !P.multi=0
;      splot, lambda[ww,i], specflux[ww,i]
;      soplot, lambda[ww,i],red_flux+objcounts, color='red'
;      stop
;   endif


   red_flux_short=bspline_valu(lambda[zz,i], norm_model)*scales[i]*exptime
   ;kludge to make sure that oversubtracted fibers at least do not
   ;go too far below zero counts in case of low object s/n
   if median(specflux[zz,i]-red_flux_short) lt 0.0 then begin
      red_flux=red_flux*(median(specflux[zz,i]/red_flux_short))
   endif
   specflux[ww,i]=specflux[ww,i]-red_flux
   ;correction is only good to about 7%, so include this systematic
   ;uncertainty on top of the poission noise that is already accounted for
   specfluxivar[ww,i]=specfluxivar[ww,i]/(1.+specfluxivar[ww,i]*(0.07*red_flux)^2.)
   red_flux_all[ww,i]=red_flux
endfor
redflux=red_flux_all
return
end
