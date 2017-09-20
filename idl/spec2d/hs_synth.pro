PRO hs_synth, hectfile
   
   
   lamout = mrdfits(hectfile, 0)
   fluxout = mrdfits(hectfile, 1)
   plugmapout = mrdfits(hectfile, 5)
   
   k = where(plugmapout.objtype EQ  'SPECTROPHOTO_STD')
   
   
   
   
   FOR ii = 0, n_elements(k) -1 DO BEGIN
      
      standfile = getenv('HSRED_DIR')+'/etc/standstar.dat'
      readcol, standfile, sra, sdec, u, g, r, i, z, $
         reddening, format='D,D,D,D,D,D,D,D', $
         comment='#', /silent
      
      dist = sqrt( (plugmapout(k(ii)).ra-sra)^2 + $
                   (plugmapout(k(ii)).dec-sdec)^2)*3600.0
      junk = min(dist, j)
  
      maggies = k_project_filters(lamout(*,k(ii)), $
                                  1d-17*fluxout(*,k(ii)))
      mag = -2.5*alog10(maggies)
      smag = [u(j), g(j), r(j), i(j), z(j)]
      leff = k_lambda_eff()
      
      diff = smag - mag
      IF ii EQ 0 THEN out = dblarr(3, n_elements(k))
      out(*,ii) = diff(1:3)
      plotsym, 0, /fill
      color = ['red', 'orange', 'green', 'blue', 'magenta']
      color = [color, color, color, color]
      IF ii EQ 0 THEN $
         djs_plot, leff(1:3), diff(1:3), ps=8, yrange=[-0.5, 0.5], $
         xrange=[4000, 8000], color=color(ii), $
         xtitle='Effective Wavelength', $
         ytitle='De-reddened SDSS Mag - AGES Synth Mag', $
         charthick=2, charsize=1.2
      
      IF ii NE 0 THEN djs_oplot, leff(1:3), diff(1:3), ps=8, color=color(ii)
      
   ENDFOR
   djs_iterstat, out, mean=median, sigma=sigma
   djs_oplot, [0, 1e9], median*[1,1], linestyle=3, thick=3, color='red'
   djs_oplot, [0, 1e9], median + sigma * [-1, -1], color='blue', $
      linestyle=2, thick=3
   
   djs_oplot, [0, 1e9], median + sigma * [1, 1], color='blue', $
      linestyle=2, thick=3
   
   djs_xyouts, 4500, 0.45, hectfile, charthick=2, charsize=1.2
   djs_xyouts, 4500, 0.4, 'Mean = ' + strn(median), charthick=2, charsize=1.2
   djs_xyouts, 4500, 0.35,'Sigma = ' + strn(sigma), charthick=2, charsize=1.2
   
   

end
