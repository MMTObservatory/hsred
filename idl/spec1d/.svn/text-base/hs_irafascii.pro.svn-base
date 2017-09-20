PRO hs_irafascii, iraffile, outname
   
   hdr = headfits(iraffile)
   flux = mrdfits(iraffile, 0)
   
   wave = sxpar(hdr,'crval1') + $
      findgen(n_elements(flux[*,0]))*sxpar(hdr,'cdelt1')
   
   output = string(wave, format='(F9.4)')
   
   FOR i = 0, 299 DO output = output + ' ' + string(flux[*,i], format='(F7.3)')
   
   openw, 1, outname
   FOR i = 0, 299 DO printf, 1, output(i)
   close, 1
   
   
END
