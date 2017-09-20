PRO hs_easycoadd, filelist, outname=outname
   
   color = ['red','yellow', 'orange', 'green', 'blue', 'cyan', $
            'magneta']
   color = [color, color, color, color, color]
   
   IF NOT keyword_set(outname) THEN outname = 'coadd.fits'
   
   l = mrdfits(filelist(0), 0)
   f = mrdfits(filelist(0), 1)
   i = mrdfits(filelist(0), 2)
   dx = l(1000,0)-l(999,0)
;   npix = ((l(4100,0)) - (l(300,0))) / dx
;   LAMOUT = min(l(300,0)) + findgen(npix)*dx
   lamout = l(300:n_elements(l(*,0))-100, 0)


   lam = dblarr(n_elements(filelist), n_elements(lamout), n_elements(l(0,*)))
   flux = lam * 0.0
   invvar = lam * 0.0
   
   plugmap = mrdfits(filelist(0), 5)
   med = dblarr(n_elements(filelist), n_elements(plugmap.object))
   scale = med * 0.0
   
   
   
   FOR i = 0, n_elements(filelist) - 1 DO BEGIN
      
      l = mrdfits(filelist(i), 0)
      f = mrdfits(filelist(i), 1)
      in = mrdfits(filelist(i), 2)
      
      FOR ii = 0, n_elements(plugmap.object) -1 DO BEGIN
            
         flux[i,*,ii] = interpol(f(*,ii), l(*,ii), lamout)
         invvar[i,*,ii] = interpol(in(*,ii), l(*,ii), lamout)
            
         IF plugmap(ii).target GT 0 THEN begin   
            ratio = flux[i,*,ii] / flux[0,*,ii]
            test = bspline_iterfit(indgen(n_elements(ratio)), ratio, $
                                   bkspace=500, $
                                   yfit=yfit)
            flux[i,*,ii] = flux[i,*,ii] / yfit
            invvar[i,*,ii] = invvar[i,*,ii] * yfit ^2
         endIF
         
         
      ENDFOR
      
      
         
            
   endfor
   
   lam_done = reform(flux[0,*,*])*0.0
   flux_done = lam_done * 0.0
   inv_done = lam_done * 0.0
   
   
   
   FOR i = 0, n_elements(filelist) -1 DO BEGIN
      IF i EQ 0 THEN lin = [lamout]
      IF i NE 0 THEN lin = [[lin], [lamout]]
   ENDFOR
   
   
   
   FOR ii = 0, n_elements(plugmap.object) -1 DO BEGIN
      print, format='("Co-adding ",i4," of ",i4,a1,$)', $
         ii, n_elements(plugmap.object), string(13b)
      
      f = transpose(flux[*,*,ii])
      i = transpose(invvar[*,*,ii])
      
      IF plugmap(ii).target GT 0 THEN begin
         
         djs_plot, lamout, flux(0, *, ii), /xs, /ys, $
            xtitle='Wavelength [Angstroms]', $
            ytitle='Counts', charthick=2, charsize=1.5, $
            ps=10
         
         FOR iii = 1, n_elements(filelist) -1 DO BEGIN
            djs_oplot, lamout, flux(iii, *, ii), /xs, /ys, $
               xtitle='Wavelength [Angstroms]', $
               ytitle='Counts', charthick=2, charsize=1.5, $
               ps=10, color=color(iii-1)
         ENDFOR
      ENDIF
      
      
      
      nflux = djs_avsigclip(f) ;;total(f*i,2)/total(i,2)
      nivar = 1/total(i,2)
    
      
      
      ;hs_combine1fiber, lin, f, i, newlam=lamout, $
      ;newivar=nivar, newflux=nflux
      
      flux_done[*,ii] = nflux
      inv_done[*,ii] = nivar
      lam_done[*,ii] = lamout
      
   ENDFOR
   
   
   mwrfits, lam_done, outname, /create
   mwrfits, flux_done, outname
   mwrfits, inv_done, outname
   mwrfits, inv_done, outname
   mwrfits, inv_done, outname
   mwrfits, plugmap, outname
   



   
   
END


        
   
   
