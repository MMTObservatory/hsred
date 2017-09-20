PRO hs_sky_pca, wset, flux, ivar, mask, plug, sset, $
                newwset=newwset, newflux=newflux, $
                newivar=newivar, newmask=newmask
   
   traceset2xy, wset, xx, lam
   ksky = where(strmatch(plug.objtype, 'SKY*'), ct)
   npix = 3834
   outlam = 3900 + 1.2*findgen(npix)
   
   nspec = n_elements(lam(0,*))
   outflux = fltarr(npix, nspec)	
   outivar = fltarr(npix, nspec)
   outmask = fltarr(npix, nspec)
   outflux1 = outflux*0.0
   sky = bspline_valu(outlam, sset)
   
   for ispec = 0, nspec-1 do begin
      
      hs_combine1fiber, lam[*,ispec], flux[*,ispec], ivar[*,ispec], $
         newlam=outlam, newflux=newflux1, newivar = newivar1, $
         finalmask=mask[*,ispec], andmask=mask1
      
      outflux[*,ispec] = newflux1
      outivar[*,ispec] = newivar1
      outmask[*,ispec] = mask1
                  
      test = bspline_iterfit(outlam, outflux(*,ispec), bkspace=200, $
                             yfit=yfit, niter=5)
      
      outflux1(*,ispec) = outflux(*,ispec) - yfit
      
   endfor
   
   k = where(strmatch(plug.objtype, 'SKY*'))
   em_pca, outflux(*,k), 5, vect
   
   junk = bspline_iterfit(outlam, vect(*,0), bkspace=200, yfit=yfit, niter=5)
   vect1 = vect(*,0)-yfit
   junk = bspline_iterfit(outlam, vect(*,1), bkspace=200, yfit=yfit, niter=5)
   vect2 = vect(*,1)-yfit
   
   
;   junk = bspline_iterfit(outlam, sky, bkspace=200, yfit=yfit, niter=5) ;
;   mask = outlam*0.0
;   mask(where(outlam GT 6140 AND outlam LT 6660)) = 1
;   mask(where(outlam GT 6790 AND outlam LT 7089)) = 1
;   mask(where(outlam GT 7230 AND outlam LT 8130)) = 1
;   mask(where(outlam GT 8250)) = 1
;   mask(where(sky GT 350))= 1.0
;   m = where(mask EQ 0)
;   test = interpol(sky[m], outlam[m], outlam)
;   junk = bspline_iterfit(outlam, test, bkspace=100, yfit=yfit, niter=5) ;
;;   yfit = bspline_valu(outlam, junk)
;   vect2 = sky-yfit
   
  
   l = where(outlam gt 5700 and outlam lt 6100)
   l1 = where(outlam GT 7650 AND outlam LT 8150)
   l = [l,l1]
   for i = 0, nspec-1 do begin
      
      npnt = 20
      scale = findgen(npnt)/float(npnt)*5000-2500
;      scale1 = findgen(npnt)/float(npnt)*0.2-0.1
      scale1 = scale
      chi = fltarr(npnt,npnt)
 ;     chi = fltarr(npnt)
      jmin=0
      mmin=0
      for j = 0, n_elements(scale) -1 do begin
         for m = 0, n_elements(scale) -1 do begin
            diff = outflux1(*,i)-vect1*scale(j)-vect2*scale1(m)
            djs_iterstat, diff, mask=mask
            ll = where(mask EQ 1)
     ;       chi(j) = total(diff(l(ll))^2)
            chi(j,m) = total(diff(l(ll))^2)
            IF chi(j,m) LT (chi(jmin,mmin)) THEN BEGIN
               jmin = j
               mmin = m
            ENDIF
         endfor
      endfor
      
      fit=poly_fit(scale1, chi(jmin,*), 2)
      besty= -1*(fit(1))/(2*fit(2))
      fit=poly_fit(scale, chi(*,mmin), 2)    
      bestx= -1*(fit(1))/(2*fit(2))  
;      fit=poly_fit(scale, chi, 2)    
;      bestx= -1*(fit(1))/(2*fit(2))  
;      outflux(*,i) = outflux(*,i) - vect(*,0)*bestx - sky*besty
;      outivar(*,i) = 1d/(1/outivar[*,i]+(vect(*,0)*bestx - sky*besty))
      outflux[*,i] = outflux[*,i] - vect(*,0)*bestx - vect(*,1)*besty
      
      ;IF i EQ 74 THEN stop
   ENDFOR
   
   x = findgen(npix)
   xx = rebin(x, npix, nspec)
   lamout = rebin(outlam, npix, nspec)
   xy2traceset, xx, lamout, wset
   
   newflux = outflux
   newwset = wset
   newivar = outivar
   newmask = outmask
   
   

   return
 
END

		

		
