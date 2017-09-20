FUNCTION hc_tracearc, fimage, fivar, xstart=xstart, ystart=ystart,$
                      window=window, yset=ycen
   
   IF NOT keyword_set(window) THEN window = 10
   IF NOT keyword_set(ystart) THEN ystart = 0L
   
   nrows = n_elements(fimage(0,*))
   npoints = n_elements(xstart)
   xcen = fltarr(nrows,npoints)*0.0
   ycen = xcen*0.0
   npix = n_elements(fimage(*,0))
   
   ;;Do the ystart row
   row = fimage(*,ystart)
   pix = findgen(n_elements(fimage(*,0)))
   
   FOR i = 0, npoints -1 DO BEGIN
      x0 = xstart(i)
      IF x0 LT 15 or x0 GT max(pix)-15 THEN BEGIN
         xcen(ystart, i) = -1
      endIF ELSE $
         xcen(ystart, i) = xstart(i)
         ;cen(ystart, i) = rjc_centline(pix(x0-window/2:x0+window/2), $
         ;                             row(x0-window/2:x0+window/2))
      
   ENDFOR
   
   lags = findgen(2000)/10 -100
   FOR irow = ystart+1, nrows -1 DO BEGIN
      row = fimage(*,irow)
      
      ccor = c_correlate(fimage(*,irow-1), row, lags)
      junk = max(ccor, ilag)
      shift = lags(ilag)
            
      FOR i =0, npoints -1 DO BEGIN
         x0 = xcen(irow-1, i)+shift
        
         
         IF x0 LT window OR x0 GT npix-window THEN BEGIN
            xcen(irow, i) = -1
         ENDIF ELSE BEGIN
            xcen(irow, i) = rjc_centline(pix(x0-window/2:x0+window/2), $
                                         row(x0-window/2:x0+window/2))
         ENDELSE
         
      
      ENDFOR
;      plot, row, /xs
;      FOR ii = 0, npoints -1 DO $
;         djs_oplot, xcen(irow,ii)*[1,1], 1e4*[-1,1], color='red', thick=2
;      wait, 0.2
   
   
   ENDFOR
   
   
   
   
   FOR i = 0, nrows -1 DO ycen(i,*) = replicate(i, n_elements(ycen(i,*)))
   
   
   return, xcen
   
END

