FUNCTION rjc_centline, xval, yval
   
   npix = n_elements(xval)
   I_bar = 1/npix * total(yval)
   
   k = where(yval GT I_bar,ct)
   IF ct GT 0 THEN $
      cent = total( (yval(k)-I_bar)*xval(k)) / total( (yval(k)-I_bar))
   IF ct EQ 0 THEN cent = -1
   
   return, cent
END

