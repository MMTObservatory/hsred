
;;:This program will take a hectospec output, find the zbest file, and plot 
;;everything


PRO hs_specplot, nsmooth=nsmooth, rerun=rerun, plateid=plateid, autoline=autoline, pseudoflux=pseudoflux, nstart=nstart, root=root, wavemin=wavemin, xrange=xrange
  
  
  
    
  if NOT keyword_set(wavemin) then wavemin=4100
  if NOT keyword_set(xrange) then xrange=[4000, 8500]
  
  
  if NOT keyword_set(root) then   begin
     
     specfile ='spHect-'+ string(plateid, format='(i3.3)') + '-' + $
        string(rerun, format='(i4.4)') + '.fits'
     zfile = 'spZbest-' +  string(plateid, format='(i3.3)') + '-' +$
          string(rerun, format='(i4.4)') + '.fits'
  endif else begin
     specfile = 'spHect-' + root + '.fits'
     zfile = 'spZbest-' + root + '.fits'
   
  endelse
  
  
  
  if not keyword_set(nsmooth) then nsmooth =1
  
  objlam = mrdfits(specfile,0)
  objflux = mrdfits(specfile,1)
  objivar = mrdfits(specfile,2)
  mask = mrdfits(specfile,3)
  plugmap = mrdfits(specfile, 5)
  zmap = mrdfits(zfile,1)
  synth = mrdfits(zfile,2)
  
  
  
  wavemin = max(objlam[2,*]) > wavemin
   
  
  
   if Keyword_set(pseudoflux) then begin
      
      fcalfile = getenv('HSRED_DIR') + '/etc/average_flux.fits'
      tcalfile = getenv('HSRED_DIR') + '/etc/average_tell.fits'
      
      fcal = mrdfits(fcalfile, 1)
      tcal = mrdfits(tcalfile, 1)
      
      flux_factor = bspline_valu(alog10(objlam), fcal)
      tell_factor = bspline_valu((objlam), tcal)
      
      
      divideflat, objflux, flux_factor, invvar=objivar
      divideflat, objflux, tell_factor, invvar=objivar
   endif
  
  
  
   synthlam = 10.0^(alog10(wavemin) + 1e-4*findgen(n_elements(synth(*,0))))
   
   stop = 0
   if not keyword_set(nstart) then nstart = 0
   ispec  = nstart
   
   while stop eq 0 do begin
      
      
      smoothivar = smooth(objivar(*,ispec),21)
      k = where( (objivar[*,ispec]/smoothivar)^(-1) gt 3)
      
      
      ivarout = objivar(*,ispec)
      if k(0) ne -1 then ivarout(k) = 0 
      
      lam = smooth(objlam(*,ispec),nsmooth)
      flux = smooth(objflux(*,ispec), nsmooth)
      ivar = ivarout
      plug = plugmap[ispec]
      eigen = synth[*,ispec]
      zbest = zmap[ispec]
      k = where(ivar gt 0)
      lam = lam(k)
      flux = flux(k)
      
      
      target=  strmatch(plug.objtype, '*SKY*') or strmatch(plug.objtype,$
               '*REJECT*') or strmatch(plug.objtype, '*UNUSED*')
      
      if target eq 0 then begin 
         
         min = min(abs(lam[*,0]-xrange(0)), nmin)
         max = min(abs(lam[*,0]-xrange(1)), nmax)
         
         ymax = max(flux[nmin:nmax])*1.5
         ymin = min(flux[nmin:nmax])
         
         ymin = ymin - (ymax-ymin)*0.1
         diff = ymax-ymin    
         
         
         djs_plot, lam, flux, ytitle = 'Flux [10^{-17} erg/s/cm^{2}]', $
            xtitle='Observed Wavelength (Ang)', /xstyle, $
            yrange=[ymin, ymax], /ystyle, xrange=xrange
         djs_oplot, lam, sqrt(1/ivar), color='red'
         djs_oplot, synthlam, eigen, color='blue', thick=2
         
         
         djs_xyouts, xrange(0)+0.05*(xrange(1)-xrange(0)), ymax-0.1*diff, $
            zbest.class, color='orange', $
            charsize=2, thick=2
         djs_xyouts, xrange(0)+0.05*(xrange(1)-xrange(0)), ymax-0.15*diff, $
            zbest.subclass, color='orange', $
            charsize=2, thick=2
         djs_xyouts, xrange(0)+0.05*(xrange(1)-xrange(0)), ymax-0.2*diff , $
            'z = ' + string(zbest.z, format='(f5.3)'), color='orange', $
            charsize=2, thick=2
         djs_xyouts, xrange(0)+0.05*(xrange(1)-xrange(0)), ymax -0.25*diff, $
            plug.objtype, color='orange', $
            charsize=2, thick=2
         djs_xyouts, xrange(0)+0.05*(xrange(1)-xrange(0)), $
            ymax-0.3*diff, 'Fiber #' + $
            string(ispec,format='(i3.3)') , color='orange', $
            charsize=2, thick=2
         
         
     
     
      if keyword_set(autoline) and keyword_set(zmap) then begin
         
         if strmid(zbest.class, 0, 6) eq 'GALAXY' then begin
            
            z = zbest.z
            
            
            
               
        lines= [3797.90,3834.00,3889.00,3933.70,3968.50,4101.70,$
                4304.40,4340.48,4861.34,5175.36,5268.98,5892.50,$
                4172.00,4226.74]
        lineid = ['H8', 'H7' , 'H6' ,  'K',  'H',  'H \delta',  $
                  'G',  'H \gamma',  'H \beta',  'Mg',  'CaFe', $
                  'Na',  'CaI',  'CaI']
         for j = 0, n_elements(lines) -1 do begin
               djs_oplot, (1+z) * [lines(j), lines(j)], [-1e5,1e5], $
                  color='red', linestyle=3
               if (j/2 eq j/2.) then $
                  djs_xyouts, (1+z)* lines(j), ymin+0.05*diff, $
                  lineid(j), color='red', charsize=1.5
               if (j/2 ne j/2.) then $
                  djs_xyouts, (1+z) * lines(j),ymin+0.1*diff, $
                  lineid(j), color='red', charsize=1.5
            endfor
            
            
             
        
        lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                      1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                      4862.68, 4960.3,5008.24, 6549.86,$
                      6564.61, 6585.27, 6718.29 ,6732.67]
        lineid = ['Ly \beta', 'Ly \alpha', 'NV', 'Si V', 'C IV', $
                  'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                  'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                  'N II', 'S II', 'S II']
        
 
        
        for j = 0, n_elements(lines) -1 do begin
               djs_oplot, (1+z) * [lines(j), lines(j)], [-1e5,1e5], $
                  color='yellow', linestyle=3
               if (j/2 eq j/2.) then $
                  djs_xyouts, (1+z)* lines(j), ymin+0.05*diff, $
                  lineid(j), color='yellow', charsize=1.5
               if (j/2 ne j/2.) then $
                  djs_xyouts, (1+z) * lines(j),ymin+0.1*diff, $
                  lineid(j), color='yellow', charsize=1.5
            endfor
             
                        
            
         endif else if strmid(zbest.class, 0, 3) EQ 'QSO' then begin
           z = zbest.z
            
             
        lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                      1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                      4862.68, 4960.3,5008.24, 6549.86,$
                      6564.61, 6585.27, 6718.29 ,6732.67]
        lineid = ['Ly \beta', 'Ly \alpha', 'NV', 'Si V', 'C IV', $
                  'He II', 'C III', 'C II', 'Mg II', 'O II', 'H \gamma', $
                  'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                  'N II', 'S II', 'S II']
                    
         
            for j = 0, n_elements(lines) -1 do begin
               djs_oplot, (1+z) * [lines(j), lines(j)], [-1e5,1e5], $
                  color='green', linestyle=3
               if (j/2 eq j/2.) then $
                  djs_xyouts, (1+z)* lines(j), ymin+0.05*diff, $
                  lineid(j), color='green', charsize=1.5
               if (j/2 ne j/2.) then $
                  djs_xyouts, (1+z) * lines(j),ymin+0.1*diff, $
                  lineid(j), color='green', charsize=1.5
            endfor
            
         endif
         
      endif
      
     
     
     
     
     test=''
     splog, 'Commands'
     splog, 'N - next'
     splog, 'P - previous'
     splog, 'Stop - stop'
     read, prompt='command: ', test
     
     if strupcase(test) eq 'N' then begin
        ispec = ispec+ 1
        
     endif
     
     if strupcase(test) eq 'P' then begin
        ispec = ispec - 1
      
     endif
     
     if strupcase(test) eq 'STOP' then begin
        
        stop = 1
    endif    
    
     if strupcase(test) eq 'N' then ispec = ispec + 1 
     if strupcase(test) eq 'P' then ispec = ispec - 1
 
  endif else ispec = ispec + 1
  
  
 endwhile
 
end

