PRO scale_sky, lam, flux, ivar, scale=scale, doplot=doplot
  
;  lam = mrdfits('25*fits', 0)
;  flux = mrdfits('25*fits', 1)
;  plug = mrdfits('25*fits', 3)
;  ivar = 1/(flux)
;  k = where(finite(ivar) eq 0, ct)
;  if ct gt 0 then ivar(k) = 0.0
  
  linecent = [5460.5,7340.7,7993.2,8398.8]
  leftbox = [5450,7335.2,7987.0,8390]
  rightbox = [5470,7346.0,7999.0,8407]
  
  leftmed1 = [5437,7333.6,7985,8388]
  leftmed2 = [5448,7335.2,7988.2,8391]
  
  rightmed1 = [5476,7348,8000.0,8406]
  rightmed2 = [5500.0, 7352, 8005.0, 8408]
  
  nspec = n_elements(lam(0,*))
  lineflux = fltarr(nspec, 4)
  
  if keyword_set(doplot) then begin
     !P.position=0
     !P.multi=0
     !P.multi=[0,4,2]
  endif
  
  
  for ispec = 0, nspec -1 do begin
     for iline = 0, 3 do begin
        
        lam1 = lam(*,ispec)
        flux1 = flux(*,ispec)
        ivar1 = ivar(*,ispec) 
        
        ;;Get leftval
        k = where(lam1 gt leftmed1(iline) and $
                  lam1 lt leftmed2(iline))
        leftval = djs_median(flux1(k))
        leftlam = djs_median([leftmed1(iline), leftmed2(iline)])
        
        
        ;;Get rightval
        k = where(lam1 gt rightmed1(iline) and $
                  lam1 lt rightmed2(iline))
        rightval = djs_median(flux1(k))
        rightlam = djs_median([rightmed1(iline), rightmed2(iline)])
        
        ;;Fit range
        k = where(lam1 gt leftbox(iline) and $
                  lam1 lt rightbox(iline))
        lam1 = lam1(k)
        flux1 = flux1(k)
        ivar1 = ivar1(k)
        
        ;;Subtract the background
        fit = linfit([leftlam,rightlam]-linecent(iline), [leftval,rightval])
        flux1 = flux1 - fit(0) - fit(1)*(lam1-linecent(iline))
                
        err = 1/sqrt(ivar1)
        k = where(ivar1 le 0, ct)
        if ct gt 0 then err(k) = 1e9
        
        result = gaussfit(lam1, flux1, a, measure_error=err, nterms=3)
        
        xx = leftbox(iline) + $
           findgen(1000)/1000.*(rightbox(iline)-leftbox(iline))
        zz = (xx-a(1))/a(2)
        pp = a(0)*exp(-zz^2/2.)
        
        if keyword_set(doplot) then begin
           plot, lam(*,ispec), flux(*,ispec), /xs, /ys, ps=10, $
              thick=2, xthick=2, ythick=2, xrange=linecent(iline) +[-20,20]
           
           djs_oplot, xx, pp+fit(0)+fit(1)*(xx-linecent(iline)), $
              color='green', thick=2
        
        endif
        
        lineflux(ispec, iline) = int_tabulated(xx, pp)

     endfor
  endfor
  
  if keyword_set(doplot) then begin
      !P.multi=0
      !P.position=0 
  endif


  for iline =0, 3 do begin
     lineflux(*,iline) = lineflux(*,iline)/djs_median(lineflux(*,iline))
  endfor

  npix = n_elements(lam(*,0))
  scale = fltarr(nspec) 
  scale = djs_median(lineflux, 2)
  scale = transpose(rebin(scale, 150,npix))
  

  
  return
  
end
