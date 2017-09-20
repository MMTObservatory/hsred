FUNCTION hc_construct_wset, specflux, oddcent=oddcent, oddwave=oddwave, $
                            evencent=evencent, evenwave=evenwave
   
;   hc_getspec, file, ccdnum, specflux=specflux
   
   odd = findgen(60)*2+1
   even = findgen(60)*2
   npix = n_elements(specflux(*,0))
   sort1 = sort(evenwave)
   evenwave = evenwave(sort1)
   evencent = evencent(sort1)
   
   evenwave = evenwave(1:n_elements(evenwave)-1)
   evencent = evencent(1:n_elements(evenwave)-1)
   excen = hc_tracearc(specflux(*,even), xstart=evencent, ystart=0, yset=eycen)
   
   sort1 = sort(oddwave)                                                    
   oddwave = oddwave(sort1)                                                
   oddcent = oddcent(sort1)                                                
   
   oddwave = oddwave(1:n_elements(evenwave)-1)                              
   oddcent = oddcent(1:n_elements(evenwave)-1)                              
   
   
   oxcen = hc_tracearc(specflux(*,odd), xstart=oddcent, ystart=0, yset=oycen)
   
   xx = findgen(npix)
   
   FOR ii = 0, 59 DO BEGIN
      k = where(excen[ii,*] GT 0)
      xy2traceset, reform(excen[ii,k]), evenwave[k],  etset, ncoeff=6
      k = where(oxcen[ii,*] GT 0)
      xy2traceset, reform(oxcen[ii,k]), oddwave[k], otset, ncoeff=6
      
         
      IF ii EQ 0 THEN BEGIN
         tset = [etset, otset]
         traceset2xy, etset, xx, ewave
         traceset2xy, otset, xx, owave
         wave = [[ewave], [owave]]
      ENDIF ELSE BEGIN
         traceset2xy, etset, xx, ewave                                        
         traceset2xy, otset, xx, owave                                         
         wave = [[wave], [ewave], [owave]]     
         tset = [tset, etset, otset]
      ENDELSE
   ENDFOR
   
   mwrfits, wave, 'calibration/0000/wave.fits'
   delvarx, xx
   xx = findgen(npix, 120)*0.0
   FOR i = 0, 119 DO xx(*,i) = findgen(npix)
   
   xy2traceset, xx, wave, wset, ncoeff=6
   
   
   return, wset
   
   
   
END


   
