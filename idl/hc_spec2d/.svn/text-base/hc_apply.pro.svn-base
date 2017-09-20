PRO hc_apply
   
   files = findfile('spObs*fits')
   
   FOR i = 0, n_elements(files) -1 DO BEGIN
      
      f = mrdfits(files(i), 1)
      in = mrdfits(files(i), 2)
      plug = mrdfits(files(i), 5)
      
      k = where(plug.target GT 0, ct)
      hdr = headfits(files(i))
      helio = sxpar(hdr, 'HELIO_RV')
      
      FOR ii = 0, ct -1 DO BEGIN
         
         file = 'arc/database/id' + string(k(ii), format='(i3.3)')
         readcol, file, p, w
         xy2traceset, p, w, tset, ncoeff=6
         x = indgen(n_elements(f(*,ii)))
         traceset2xy, tset, x, lam1
         lam1 = lam1 / (1+helio/299792.458)

         IF ii EQ 0 THEN lam = lam1
         IF ii NE 0 THEN lam = [[lam], [lam1]]
         
      ENDFOR
      
      outfile = 'obj-' + files(i)
      mwrfits, lam, outfile, /create
      mwrfits, f(*,k), outfile
      mwrfits, in(*,k), outfile
      mwrfits, in(*,k), outfile
      mwrfits, in(*,k), outfile
      mwrfits, plug(k), outfile
      
      
      
   ENDFOR
   
   
END

