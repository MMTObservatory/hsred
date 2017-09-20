PRO hs_comment, hectfile, nsmooth=nsmooth, $
                nstart=nstart, xrange=xrange, $
                nsynth=nsynth, wavemin=wavemin
   
   IF NOT keyword_set(psym) THEN psym=10
   IF NOT keyword_set(xrange) THEN xrange = [3800, 8500]
   IF NOT keyword_set(clipfact) THEN clipfact = 0.6
   
   IF NOT keyword_set(nsmooth) THEN nsmooth=11
   IF NOT keyword_set(nstart) THEN nstart = 0
   
   ;;Parse the hectfile to get the Zbest file and the comment file
   
   junk = strsplit(hectfile, ' ', length=ct)
   root = strmid(hectfile, 7, ct-7-5)
   commentfile = root + '_comment.txt'
   spzbest = 'spZbest-' + root + '.fits'
   
   IF file_test(commentfile) EQ 1 THEN BEGIN
      cdat = rc_readfile(commentfile, 0)
      cz = findgen(n_elements(cdat))*0.0
      comment = strarr(n_elements(cdat))
      
      FOR i = 0, 299 DO BEGIN
         junk = strsplit(cdat(i), '|', /extract)
         cz(i) = float(junk(0))
         comment(i) = junk(1)
      ENDFOR
      
   ENDIF ELSE BEGIN
      cz = replicate(-999, 300)
      comment = replicate(' ... ', 300)
   ENDELSE
   
   ;;Read the file
   filename = hectfile
   
   
   ;;Read the file
   lam = mrdfits(filename, 0)
   object = mrdfits(filename, 1)
   ivar = mrdfits(filename, 2)
   mask = mrdfits(filename, 4)*0.0
   plugmap = mrdfits(filename, 5)
   
   zdat = mrdfits(spzbest, 1)
   IF NOT keyword_set(nsynth) THEN nsynth = 2
   synth = mrdfits(spzbest, nsynth)
   
   IF NOT keyword_set(wavemin) THEN wavemin = 4100
   wavemin = max(lam[2,*]) > wavemin
   synthlam = 10.0^(alog10(wavemin) + 1e-4*findgen(n_elements(synth(*,0))))

   
   
   
   test = ''
   splog, 'Commands'
   splog, 'N - next'
   splog, 'P - previous'
   splog, 'Stop - stop'
   splog, 'A - accept redshift'
   splog, 'L - plot lines'
   splog, 'R - clear lines'
   splog, 'Z - Add user redshift'
   splog, 'COM - Add user comment'
   splog, 's - plot synthetic spectrum'
   splog, 'ul - plot lines with current user z'
   splog, 'sl - plot lines with spZbest lines'
   
   i = nstart
   stop = 'go'
  
   
      
      
    while stop ne 'stop'  AND i LT 300 do begin
     
     ;;Mask pixels were the ivar is lt 0.6 the smoothed value
     temp = ivar[*,i] / smooth(ivar[*,i], 21)
     newmask = temp lt clipfact
     junk = min(abs(lam[*,i]-5588),pix)
     mask[(pix-10):(pix+10),i]=mask[(pix-10):(pix+10),i]+$
        newmask[(pix-10):(pix+10)]
     junk = min(abs(lam[*,i]-5577), pix)
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     junk = min(abs(lam[*,i]-6300), pix)
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     junk = min(abs(lam[*,i]-6363), pix)
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     pixels = where(mask[*,i] eq 0)  
     
     l1 =  min(abs(lam[pixels,i]-xrange(0)),l)
     m1 =  min(abs(lam[pixels,i]-xrange(1)),m)
     
     l = l(0)
     m = m(0)
     
     ymax = max(smooth(object[pixels(l):pixels(m),i],nsmooth)) *1.1
     ymin = min(smooth(object[pixels(l):pixels(m),i],nsmooth))
     diff = ymax-ymin  
     
      djs_plot, lam[pixels,i], smooth(object[pixels,i],nsmooth),$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=2
     djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
        ymin + (ymax-ymin)*0.95, 'Aperture  ' + strn([i]), $
        color='yellow', charthick=2, charsize=2, /isolatin1
     
     
     djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
        ymin + (ymax-ymin)*0.90  ,  plugmap[i].objtype, $
        color='yellow', charthick=2, charsize=2, /isolatin1
     djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
        ymin + (ymax-ymin)*0.95, $
        'icode = ' + strn(plugmap[i].icode), $
        color='yellow', charthick=2, charsize=2, /isolatin1
     djs_xyouts, xrange(0)+(xrange(1)-xrange(0))*0.05, $
        ymin + (ymax-ymin)*0.85, $
        'spZbest = ' + strn(zdat[i].z), $
        color='green', charsize=2, thick=2, charthick=2, /isolatin1
      djs_xyouts, xrange(0)+(xrange(1)-xrange(0))*0.40, $
        ymin + (ymax-ymin)*0.85, $
        'User z = ' + strn(cz(i)), $
        color='green', charsize=2, thick=2, charthick=2, /isolatin1
     
	z = zdat(i).z
        
        lines= [3797.90,3834.00,3889.00,3933.70,3968.50,4101.70,$
                4304.40,4340.48,4861.34,5175.36,5268.98,5892.50,$
                4172.00,4226.74]
        labels = ['H8', 'H7' , 'H6' ,  'K',  'H',  'H \delta',  $
                  'G',  'H \gamma',  'H \beta',  'Mg',  'CaFe', $
                  'Na',  'CaI',  'CaI']
        for j = 0, n_elements(lines) -1 do begin
           djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
              color='magenta', linestyle=2
           djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
              labels(j), charsize=2, color='yellow'
        endfor
        
        lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                 1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                 4862.68, 4960.3,5008.24, 6549.86,$
                 6564.61, 6585.27, 6718.29 ,6732.67]
        labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si V', 'C IV', $
                  'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                  'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                  'N II', 'S II', 'S II']
        for j = 0, n_elements(lines) -1 do begin
           djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
              color='yellow', linestyle=2
           djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
              labels(j), charsize=2, color='yellow'
        endfor 	
        
        


     loop = 0 
     
     while loop eq 0 do begin
        read, prompt='Command: ', test
        
        if strupcase(test) eq 'N' then begin
           i = i+ 1
           loop = 1
        endif
        
        if strupcase(test) eq 'P' then begin
           i = i - 1
           loop = 1
        endif
        
        if strupcase(test) eq 'R' then begin
           i = i 
           loop = 1
        endIF
        
        if strupcase(test) eq 'STOP' then begin
           loop = 1
           stop = 'stop'
        ENDIF
        
        if strupcase(test) eq 'L'  $
           OR strupcase(test) EQ 'SL' $
           OR strupcase(test) EQ 'UL' then begin
           
           IF strupcase(test) EQ 'L' THEN BEGIN
              splog, 'Enter the Redshift you would like plotted'
              read, prompt='Redshift : ', z
           ENDIF
           IF strupcase(test) EQ 'SL' THEN z = zdat(i).z
           IF strupcase(test) EQ 'UL' THEN z = cz(i)
           
           lines= [3797.90,3834.00,3889.00,3933.70,3968.50,4101.70,$
                   4304.40,4340.48,4861.34,5175.36,5268.98,5892.50,$
                   4172.00,4226.74]
           labels = ['H8', 'H7' , 'H6' ,  'K',  'H',  'H \delta',  $
                     'G',  'H \gamma',  'H \beta',  'Mg',  'CaFe', $
                     'Na',  'CaI',  'CaI']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='magenta', linestyle=2
              djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
                 labels(j), charsize=2, color='yellow'
           endfor
           
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si V', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='yellow', linestyle=2
              djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
                 labels(j), charsize=2, color='yellow'
           endfor 
           
                      
        endif
        
        IF strupcase(test) EQ 'S' THEN BEGIN
           djs_oplot, synthlam, synth[*,i], color='blue', $
              thick=2, ps=10
        ENDIF
        
        IF strupcase(test) EQ 'Z' THEN BEGIN
           splog, 'Please enter the redshift for this spectrum'
           splog, 'Common notations:'
           splog, '---- z     - the redshift'
           splog, '---- -z    - and uncertain redshift'
           splog, '---- -999  - no attempt / no clue'
           splog, '---- 100   - spZbest looks right'
           
           test1 = ' ' 
           read, prompt='Command:', test1
           cz(i) = float(test1)
           
        ENDIF
        
   	IF strupcase(test) EQ 'A' THEN BEGIN
          cz(i) =  zdat(i).z
           
        ENDIF
        
        IF strupcase(test) EQ 'C' THEN BEGIN
           read, prompt='Comment:', test1
           comment(i) = test1
        ENDIF
        
        
        output = string(cz) + '   |   ' + comment
        
        openw, 1, commentfile
        FOR ii = 0, n_elements(output) -1 DO $
           printf, 1, output(ii)
        close, 1
        
        
        
     ENDWHILE
     
  ENDWHILE
  
  
  output = string(cz) + '   |   ' + comment
  
  openw, 1, commentfile
  FOR i = 0, n_elements(output) -1 DO $
     printf, 1, output(i)
  close, 1
        
             
        
  END

