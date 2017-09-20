PRO hc_getarc
   
   file = 'calibration/0000/arc.fits'
   
   hc_getspec, file, 1, specflux=arc
   hc_getspec, file, 2, specflux=arc1
   
   spawn, 'mkdir reduction'
   spawn, 'mkdir reduction/0000'
   spawn, 'mkdir reduction/0000/arc'
;   mwrfits, arc(*,0), 'reduction/0000/arc/000.fits', /create
;   mwrfits, arc1(*,2), 'reduction/0000/arc/120.fits', /create
   
   files =findfile('hect*fits')
   p = mrdfits(files(0), 5)
   k = where(p.target GT 0, ct)
   
   FOR i = 0, ct -1 DO BEGIN
      
      file = 'reduction/0000/arc/' + string(k(i), format='(i3.3)') +'.fits'
      mwrfits, arc(*,k(i)), file, /create
      
   endfor
   
   
   p = mrdfits(files(0), 6)
   k = where(p.target GT 0, ct)
   
   FOR i = 0, ct -1 DO BEGIN
      
      file = 'reduction/0000/arc/' + string(k(i)+120, format='(i3.3)') +'.fits'
      mwrfits, arc1(*,k(i)), file, /create
      
   endfor
   
   
   
   
END
