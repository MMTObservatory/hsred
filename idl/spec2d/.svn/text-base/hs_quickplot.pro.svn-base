PRO hs_quickplot, flux=flux, lam=lam, plugmap=plugmap, rerun=rerun
  
  if not keyword_set(rerun) then rerun = 0
  if not keyword_set(plateid) then plateid=0
  ;;This is the prototype for a signal to noise calculator for the qaplots,
   if size(rerun, /tname) NE 'STRING' then $
     rerun = string(rerun, format='(i4.4)')  
   if size(plateid, /tname) NE 'STRING' then $
     plateid = string(plateid, format='(i3.3)')  
   
  if not keyword_set(flux) then  flux = mrdfits(plateid +'-'+ rerun +'.fits',1)
  if not keyword_set(plugmap) then  plugmap = mrdfits(plateid +'-'+ rerun +'.fits',5)
  if not keyword_set(lam) then  lam = mrdfits(plateid + '-' + rerun + '.fits', 0)
   
 
  flux1 = flux
  flux1(0,*) = 0
  
  npix = n_elements(flux(*,0))
  flux1(1:npix-1,*) = flux(0:npix-2,*)
  
  diff = flux-flux1
  sum = flux+flux1
  
  snr = abs(sum/diff)
  
  nfibers = n_elements(flux(0,*))
  
  SN = fltarr(nfibers)
  for ifib = 0, nfibers - 1 do begin
     k = where(abs(lam(*,ifib) - 6000.) eq min(abs(lam(*,ifib)-6000.)))
     
     lmin = k(0)-400
     lmax = k(0)+400
     
     SN(ifib) = djs_median(snr[lmin:lmax,ifib])
  endfor
  
  mag = plugmap.rapmag
  
  
  plotsym, 0, /fill
  
  title = 'Rerun ' + rerun + ' plateid ' + plateid
  
;  window,1
  
  plot, mag, sn, ps =8, xrange=[25,17], /ylog, $
     xtitle='R Aperture Magnitude', ytitle='Signal to Noise Ratio', $
     charsize=1.5, title= 'Rerun ' + rerun + ' plateid ' + plateid, /xstyle, $
     /ystyle
  
  k = where(mag gt 0)
  
  fit = linfit(mag(k), alog10(sn(k)))
  
  x = -1 * findgen(1000)/1000*(26-14) + 26
  
  line = fit(0) + x * fit(1)
  
  djs_oplot, x, 10^line, color='red', thick=2
  
  crit =  fit(0) + 21.0 * fit(1)
  
  xyouts, 21, 5, 'S/N at r=21 : ' + string(10.^crit, format='(d8.3)'), charsize=1.5
  
  
  ;;Now do the residuals
  
  expect = fit(0) + mag(k)* fit(1)
  
  resid = alog10(sn(k)) - expect
  
  minr = -1
  maxr = 0.4
  
  dc = (maxr - minr)/300.  
  
  xfocal  =(plugmap.xfocal)(k)
  yfocal = (plugmap.yfocal)(k)
  resid(10) = -0.4
  
  setup = ''
  loadct, 39
  for i = 0, 299 do begin
     
     
     bmin = minr + float(i*dc)
     bmax = bmin + dc 
     
     k = where(resid ge bmin and resid lt bmax)
      
     
     
     if setup eq '' then begin
        if k(0) ne -1 then begin
           
           plot, [xfocal[k]], [yfocal[k]], ps=8, xrange=[-400,400], yrange=[-400,400], /nodata
           oplot,  [xfocal[k]], [yfocal[k]], ps=8, color=360-i
           colorbar, datarange=[minr, maxr], format='(f4.1)',   /vertical, position=[0.90, 0.1, 0.95, 0.90]
         setup = 'dd'
      endif
      
      endif else begin
        
        if k(0) ne -1 then $
           oplot, [xfocal[k]], [yfocal[k]], ps=8, color=360-i
      
     endelse
     
     

     
     
  endfor
  
  
  
     
     
     
     
  
  
  
  
END
