FUNCTION make_model_spec, wave, lineintensity, linewave, nsmooth=nsmooth
   
   IF NOT keyword_set(nsmooth) THEN nsmooth=2
   
   npix = n_elements(wave)
   nline = n_elements(linewave)
   
   nsz = 4 * fix(nsmooth) + 1   ; kernal size
   gausskern = exp( -( ((nsz-1)/2 - findgen(nsz))^2 ) / nsmooth^2)
   gausskern = gausskern / total(gausskern) ; normalized
   
   npad = fix((npix+1/2))
   model = fltarr(npix)
   
   FOR iline = 0, nline -1 DO BEGIN
      qless = ( linewave[iline] LT wave)
      iloc = (where(qless))[0] -1
      IF (iloc GE 1 AND iloc LE npix -2) THEN BEGIN
         dx = linewave[iline] - wave[iloc]
         dpix = wave[iloc+1] - wave[iloc]
         model[iloc:iloc+1] = model[iloc:iloc+1] + $
            lineintensity[iline]*[1-dx/dpix, dx/dpix]
      endIF
   ENDFOR
   
   model = convol( model, gausskern, /center, /edge_truncate)
   
   return, model / max(model)
END


   
   
   

PRO makewset, rerun=rerun
   
   outfile = 'calibration/' + rerun + '/wset.fits'
   test = mrdfits('/tmp/wset_test.fits', 1)
   coeff = test.coeff
   
   wset = {func : test[0].func, xmin: test[0].xmin, xmax: test[0].xmax, $
           coeff : test.coeff}
   mwrfits, wset, outfile, /create
   test = mrdfits('/tmp/wset_test.fits', 2)
   coeff = test.coeff
   
   wset = {func : test[0].func, xmin: test[0].xmin, xmax: test[0].xmax, $
           coeff : test.coeff}
   mwrfits, wset, outfile
END


   
   
PRO hc_initfit, specflux, lampfile=lampfile, rerun=rerun, $
                debug=debug, psplot=psplot, xrange=xrange, crit=crit, $
                oddcent=oddcent, oddwave=oddwave, evencent=evencent, $
                evenwave=evenwave
   
   guessfile = '/tmp/hsred_arcguess.fits'
   
   IF NOT keyword_set(rerun) THEN rerun = '0000'
   
   IF NOT keyword_set(lampfile) THEN $
      lampfile = filepath('lampthar.dat', $
                          root_dir = getenv('HSRED_DIR'), $
                          subdirectory='etc')
   IF NOT keyword_set(crit) THEN crit = 500
   
;   hc_getspec, file, ccdnum, specflux=specflux, specivar=specivar
      
   
   guess = mrdfits(guessfile, 0)
   fit = poly_fit(guess(*,0), guess(*,1), 1)
   pix = findgen(n_elements(specflux(*,0)))
   
   
   
   
   wave = 0
   FOR i = 0, n_elements(fit)-1 DO $
      wave = wave + fit(i)*pix^i
   
   
   readcol, '/home/rcool/idl/hsred/etc/lampthar.dat', $
      wave1orig, intenorig
   npix = n_elements(wave)
   print, min(wave), max(wave)
   k = where(wave1orig GT wave(0) AND wave1orig LT wave(npix-1))
   wave1 = wave1orig(k)
   inten = intenorig(k)
   
   nkeep = 50
   sort1 = sort(-1*inten)
   cent = findgen(nkeep)*0.0
   wave1temp = wave1(sort1(0:nkeep-1))
   index = indgen(nkeep)
   window = 15
   radius = 3
   IF max(pix) GT 4000 THEN window = 30
   IF max(pix) GT 4000 THEN radius = 6
   
   FOR i = 0, n_elements(cent)-1 DO BEGIN
      junk = min(abs(wave-wave1temp(i)),j)
      IF j LT window+radius OR j GT n_elements(wave) -window-radius THEN BEGIN
         cent(i) = -1
      endIF ELSE BEGIN
         x0 = 0
         junk = max(specflux(J-window:j+window, 0), x0)
         x0 = x0 + j-window
         xval = pix[x0-radius:x0+radius]
         yval = specflux(x0-radius:x0+radius, 0)
         cent(i) = rjc_centline(xval, yval)
      ENDELSE
   ENDFOR
   
   k = where(cent GT 0)
   cent = cent(k)
   wave1temp = wave1temp(k)
   index = index(k)
   
   fit = poly_fit(cent, wave1temp, 3, yfit=yfit)
         
   IF NOT keyword_set(order) THEN order = 6
         
   wave = 0
   FOR i = 0, n_elements(fit)-1 DO $
      wave = wave + fit(i)*pix^i
   
   plotsym, 0, /fill
   djs_plot, cent, wave1temp-yfit, ps=8, /xs, /ys, $
      color='blue', symsize=2
   stop = ''
   in = ''
   x1 = cent
   y1 = wave1temp-yfit
   stop = ''
   superstop = ''
   WHILE superstop EQ '' DO begin
      WHILE stop EQ '' DO BEGIN
         splog, 'type stop to stop rejection - '
         read, in, prompt='Stop -- stop reject : '
         IF strupcase(in) EQ "STOP" THEN stop = '1'
         IF stop NE '1' THEN BEGIN
            
            splog, 'Click on data point to reject it'
            cursor,x, y, 4
            xr = max(x1)-min(x1)
            yr = max(y1)-min(y1)
            dist = sqrt( (x1-x)^2/xr^2 + (y1-y)^2/yr^2)
            
            min1 = min(dist,j)
            djs_oplot, [x1(j)], [y1(j)], ps=8, symsize=1.5
            k = where(dist NE min(dist))
            x1 = x1(k)
            y1 = y1(k)
            cent = cent(k)
            wave1temp = wave1temp(k)
            index = index(k)
         ENDIF
      ENDWHILE
      stop =''
      
      fit = poly_fit(cent, wave1temp, order, yfit=yfit)
      
      wave = 0
      FOR i = 0, n_elements(fit)-1 DO BEGIN
         wave = wave + fit(i)*pix^i
      ENDFOR
      plotsym, 0, /fill
      djs_plot, cent, wave1temp-yfit, ps=8, /xs, $
         /ys, symsize=2, color='blue'
      in = ''
      x1 = cent
      y1 = wave1temp-yfit
      splog, 'Type DONE to end fitting'
      read, in, prompt='Done? - '
      IF strupcase(in) EQ 'DONE' THEN superstop = 'l'
      IF strupcase(in) EQ 'ORDER' THEN BEGIN
         read, order, prompt = 'NEW order - - '
      ENDIF
      
      
   ENDWHILE
   wave = 0
   FOR i = 0, n_elements(fit)-1 DO BEGIN
      wave = wave + fit(i)*pix^i
   ENDFOR
   
   lags = findgen(2000)/10.-100
   waveorig = wave
   fit1_range = fit(1) + (findgen(20)/20.*4.0-2.0)*fit(1)
   fit2_range = fit(2) + (findgen(10)/10.*1.2-0.6)*fit(2)
   fit3_range = fit(3) + (findgen(6 )/6.*1.2-0.6)*fit(3)
   fit4_range = fit(4) ; (findgen(6)/6.*1.0-0.5)*fit(4)                       
  
   lags = findgen(320)/5. - 32
   
   bestcor = -1
   fitattempt = fit
   FOR iloop1 = 0, n_elements(fit1_range) -1 DO BEGIN
      FOR iloop2 = 0, n_elements(fit2_range) -1 DO BEGIN
         FOR iloop3 = 0, n_elements(fit3_range) -1 DO begin
            FOR iloop4 = 0, n_elements(fit4_range) -1 DO begin
               fitattempt(1) = fit1_range(iloop1)
               fitattempt(2) = fit2_range(iloop2)
               fitattempt(3) = fit3_range(iloop3)
               fitattempt(4) = fit4_range(iloop4)
               wavetest = 0
               FOR i = 0, n_elements(fitattempt)-1 DO BEGIN
                  wavetest = wavetest + fitattempt(i)*pix^i
               ENDFOR         
               k = where(inten GT crit)
               model = make_model_spec(wavetest, inten(k), wave1(k))
               
               ccor = c_correlate((specflux(*,0)/max(specflux(*,0))), $
                                  (model), lags)
               mm = max(ccor, j)
               IF mm GT bestcor THEN BEGIN
                  bestcor = mm
                  bestfit = fitattempt
                  print, iloop1, iloop2, iloop3
               ;;correct for the best fit lag
                  bestfit(0) = bestfit(0) - lags(j)*bestfit(1)
               ENDIF
               
            ENDFOR
         ENDFOR
      ENDFOR
   ENDFOR
   
   fit = bestfit
   
   wave = 0    
   FOR i = 0, n_elements(fit)-1 DO BEGIN              
      wave = wave + fit(i)*pix^i                         
   ENDFOR      
   
   xy2traceset, pix, wave, tset, ncoeff=6
   tset1 = tset
   
   IF keyword_set(debug) THEN begin
      djs_plot, wave, specflux(*,0)/max(specflux(*,0)), $
         /xs, /ys, xrange=xrange
      djs_oplot, wave, make_model_spec(wave, inten, wave1), $
         color='red', thick=2
      wait, 1
   ENDIF
   
   
   
   nkeep = 90
   sort1 = sort(-1*inten)
   cent = findgen(nkeep)*0.0
   wave1temp = wave1(sort1(0:nkeep-1))
   index = indgen(nkeep)
   
   FOR i = 0, n_elements(cent)-1 DO BEGIN
      junk = min(abs(wave-wave1temp(i)),j)
      IF j LT window+radius OR j GT n_elements(wave) - window-radius THEN BEGIN
         cent(i) = -1
      endIF ELSE BEGIN
         x0 = 0
         junk = max(specflux(J-window:j+window, 0), x0)
         x0 = x0 + j-window
         xval = pix[x0-radius:x0+radius]
         yval = specflux(x0-radius:x0+radius, 0)
         cent(i) = rjc_centline(xval, yval)
      ENDELSE
      
   ENDFOR
   
   
   k = where(cent GT 0)
   cent = cent(k)
   wave1temp = wave1temp(k)
         
   evencent = cent
   evenwave = wave1temp
   
   ispec = 1
   
   lags = findgen(2000)/10.-100

   
   junk = max(c_correlate(specflux(*,ispec-1), $
                          specflux(*,ispec), lags), ilag)
   fit(0) = fit(0) - fit(1)*lags(ilag)
   
   ;;new
   oddcent = evencent + lags(ilag)
   oddwave = evenwave
   
   IF keyword_set(debug) THEN $
      plot, lags,c_correlate(specflux(*,ispec-1), $     
                              specflux(*,ispec), lags), ps=10
   
   wave = 0
   FOR i = 0, n_elements(fit)-1 DO $
      wave = wave + fit(i)*pix^i
   
   
;   k = where(wave1orig GT wave(0) AND wave1orig LT wave(npix-1))
;   wave1 = wave1orig(k)
;   inten = intenorig(k)      
 ;  
 ;  wave = 0
 ;  FOR i = 0, n_elements(fit)-1 DO $
 ;     wave = wave + fit(i)*pix^i
 ;  
 ;  fit1_range = fit(1) + (findgen(10)/10.*1.0-0.5)*fit(1)                   
 ;  fit2_range = fit(2) + (findgen(10)/10.*1.2-0.6)*fit(2)                
 ;  fit3_range = fit(3) + (findgen(6 )/6.*1.2-0.6)*fit(3)                   
 ;  fit4_range = fit(4) ; (findgen(6)/6.*1.0-0.5)*fit(4)
 ;  lags = findgen(320)/5. - 32          
 ;  bestcor = -1
 ;  fitattempt = fit
 ;  
 ;  FOR iloop1 = 0, n_elements(fit1_range) -1 DO BEGIN
 ;     FOR iloop2 = 0, n_elements(fit2_range) -1 DO BEGIN
 ;        FOR iloop3 = 0, n_elements(fit3_range) -1 DO begin
 ;           FOR iloop4 = 0, n_elements(fit4_range)-1 DO BEGIN
 ;              fitattempt(1) = fit1_range(iloop1)
 ;              fitattempt(2) = fit2_range(iloop2)
 ;              fitattempt(3) = fit3_range(iloop3)
 ;              fitattempt(4) = fit4_range(iloop4)
 ;              wavetest = 0
 ;              FOR i = 0, n_elements(fitattempt)-1 DO BEGIN
 ;                 wavetest = wavetest + fitattempt(i)*pix^i
 ;              ENDFOR        
 ;;              
 ;              k = where(inten GT crit)
 ;              model = make_model_spec(wavetest, inten(k), wave1(k), $
 ;                                      nsmooth=2)
 ;              ccor =$
 ;                 c_correlate(sqrt(specflux(*,0)/max(specflux(*,0))>0), $
 ;                             sqrt(model), lags)
 ;              
 ;              mm = max(ccor, j)
 ;              IF mm GT bestcor THEN BEGIN
 ;                 bestcor = mm
 ;                 bestfit = fitattempt
 ;                 wavetemp = 0 
 ;                 
 ;                 print, iloop1, iloop2, iloop3, iloop4
 ;                 FOR i = 0, n_elements(fitattempt)-1 DO $
 ;                    wavetemp = wavetemp + bestfit(i)*pix^i  
 ;              ENDIF
 ;           ENDFOR
 ;        ENDFOR
 ;     ENDFOR
 ;  ENDFOR
   
 ;  fit = bestfit
 ;  
 ;  lags = findgen(1000)/100.-5
 ;  model = make_model_spec(wavetemp, inten, wave1)
 ;  junk = max(c_correlate(specflux(*,ispec), $
 ;                         model, lags), ilag)
 ;  fit(0) = fit(0) + fit(1)*lags(ilag)
 ;  
 ;  wave = 0
 ;  FOR i = 0, n_elements(fit)-1 DO $
 ;     wave = wave + fit(i)*pix^i
 ;  
 ;  nkeep = 90
 ;  sort1 = sort(-1*inten)
 ;  cent = findgen(nkeep)*0.0
 ;  wave1temp = wave1(sort1(0:nkeep-1))
 ;  index = indgen(nkeep)
 ;  
 ;  FOR i = 0, n_elements(cent)-1 DO BEGIN
 ;     junk = min(abs(wave-wave1temp(i)),j)
 ;     IF j LT window+radius OR j GT n_elements(wave)-window-radius THEN BEGIN
 ;        cent(i) = -1
 ;     endIF ELSE BEGIN
 ;        x0 = 0
 ;        junk = max(specflux(J-window:j+window, 0), x0)
 ;        x0 = x0 + j-window
 ;        xval = pix[x0-radius:x0+radius]
 ;        yval = specflux(x0-radius:x0+radius, 0)
 ;        cent(i) = rjc_centline(xval, yval)
 ;     ENDELSE
 ;     
 ;  ENDFOR
   
   
   
;   k = where(cent GT 0)
;   cent = cent(k)
;   wave1temp = wave1temp(k)
         
;   oddcent = cent
;   oddwave = wave1temp
   
   
   ;;new
   k = where(oddcent GT 0 AND oddcent LT n_elements(specflux(*,1))-1)
   oddcent = oddcent(k)
   oddwave = oddwave(k)
   
END

