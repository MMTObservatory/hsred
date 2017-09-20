PRO hc_initialwave, specflux, $
                    rerun=rerun, debug=debug
   
   
   IF NOT keyword_set(rerun) THEN rerun = 0
   IF size(rerun, /tname) NE 'STRING' THEN $
      rerun=string(rerun, format='(i4.4)')
   
   tmpfile = '/tmp/hsred_arcguess.fits'
         
;   hc_getspec, file, ccdnum, specflux=specflux, specivar=specivar, $
;      rerun=rerun
   splog, 'DOING INITIAL GUESS FOR APERTURE 0'
   tracenum=0
   
   ;;This is an interactive program to get an initial guess $
   ;;for the wavelength
   ;;solution of the image
   npix = n_elements(specflux(*,0))
   arc = specflux(*,0)
   pix = findgen(n_elements(arc))
   
   djs_plot, pix, arc, /xs, /ys, xtitle='Pixel Number', $
      ytitle='Counts', $
      charsize=2
   
   splog, 'Please choose a line to zoom in my clicking:'
   cursor, x, y , 4
   x = long(x)
   djs_plot, pix(x-40:x+40), arc(x-40:x+40), /xs, /ys, $
      xtitle='Pixel Number', $
      ytitle='Counts', charsize=2                               
   splog, 'Please click near the middle of the spectral line:'
   cursor, x, y, 4
   x = long(x)
   I_bar = 1/(11.0) * total(arc(x-5:x+5))
   temp = arc(x-5:x+5)
   xpix = pix(x-5:x+5)
   k = where( (temp - I_bar) GT 0)
   cent = total( (temp(k) - I_bar)*xpix(k) )/ total( (temp(k)-I_bar))
   
   djs_oplot, cent*[1,1], [-1e5, 1e5], linestyle=2, thick=2
   splog, 'Please enter the wavelength (in angstroms) of this line:'
   read, inwave, prompt='Wavelength : ' 
   
   IF n_elements(wave1) EQ 0 THEN BEGIN
      wave1 = inwave
      x1 = cent
   ENDIF ELSE BEGIN
      wave1 = [wave1, inwave]
      x1 = [x1, cent]
   ENDELSE
   
   stop = ''
   WHILE STOP EQ '' DO BEGIN
      DJS_PLOT, PIX, ARC, /XS, /YS, XTITLE='PIXEL NUMBER', $
         YTITLE='COUNTS', $
         CHARSIZE=2
      
      SPLOG, 'PLEASE CHOOSE A LINE TO ZOOM IN MY CLICKING:'
      CURSOR, X, Y ,4
      X = LONG(X)
      DJS_PLOT, PIX(X-40:X+40), ARC(X-40:X+40), /XS, /YS, $
         XTITLE='PIXEL NUMBER', $
         YTITLE='COUNTS', CHARSIZE=2                               
      SPLOG, 'PLEASE CLICK NEAR THE MIDDLE OF THE SPECTRAL LINE:'
      CURSOR, X, Y, 4
      X = LONG(X)
      I_BAR = 1/(11.0) * TOTAL(ARC(X-5:X+5))
      TEMP = ARC(X-5:X+5)
      XPIX = PIX(X-5:X+5)
      K = WHERE( (TEMP - I_BAR) GT 0)
      CENT = TOTAL( (TEMP(K) - I_BAR)*XPIX(K) )/ TOTAL( (TEMP(K)-I_BAR))
      DJS_OPLOT, CENT*[1,1], [-1E5, 1E5], LINESTYLE=2, THICK=2
      SPLOG, 'PLEASE ENTER THE WAVELENGTH (IN ANGSTROMS) OF THIS LINE:'
      READ, INWAVE, PROMPT='WAVELENGTH : ' 
      
      IF N_ELEMENTS(WAVE1) EQ 0 THEN BEGIN
         WAVE1 = INWAVE
         X1 = CENT
      ENDIF ELSE BEGIN
         WAVE1 = [WAVE1, INWAVE]
         X1 = [X1, CENT]
      ENDELSE
      
      SPLOG, 'PLEASE TYPE STOP IF YOU WOULD LIKE TO STOP FITTING LINES'
      dostop = ''
      READ, DOSTOP, PROMPT=' STOP - TO QUIT -- Any key to continue : '
      IF STRUPCASE(DOSTOP) EQ "STOP" THEN STOP = 'STOP'
   ENDWHILE
   
   fit = linfit(x1, wave1)
   tset = [fit(0)+fit(1)*(n_elements(arc))/2.0, $
           fit(1)*n_elements(arc)/2.0, 0, 0]
   y = [[x1],[wave1]]
   
   mwrfits, y, tmpfile, /create
   
   delvarx, wave1
   delvarx, x1
         
   return
   
   
   
END



   
   
   
