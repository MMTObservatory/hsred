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
for i=0, n_elements(scales)-1 do begin
    scale_model=mrdfits(getenv('HSRED_DIR')+'/etc/redleak_perfiber_model.fits',fix(fiber_num[i]), /silent)
   scales[i]=10.^(bspline_valu(fiber_radius[i], scale_model))
   ww=where(lambda[*,i] gt 8322.)
   zz=where((lambda[*,i] gt 9010.) and (lambda[*,i] lt 9150.))
   red_flux=bspline_valu(lambda[ww,i], norm_model)*scales[i]*exptime
   red_flux_short=bspline_valu(lambda[zz,i], norm_model)*scales[i]*exptime
   ;kludge to make sure that oversubtracted fibers at least do not
   ;go too far below zero counts in case of low object s/n
   if median(specflux[zz,i]-red_flux_short) lt 0.0 then red_flux=red_flux*(median(specflux[zz,i]/red_flux_short))
   specflux[ww,i]=specflux[ww,i]-red_flux
   ;correction is only good to about 7%, so include this systematic
   ;uncertainty on top of the poission noise that is already accounted for
   specfluxivar[ww,i]=specfluxivar[ww,i]/(1.+specfluxivar[ww,i]*(0.07*red_flux)^2.)
   red_flux_all[ww,i]=red_flux
endfor
redflux=red_flux_all
return
end
