PRO hs_manual, spec, lam
  
  
  lampdefault = filepath('lamphenear.dat', $
                         root_dir=getenv('IDLSPEC2D_DIR'), $
                         subdirectory='etc')
  lampfilename = lampdefault
  
  splog, 'Reading lamp file ', lampfilename
  readcol, lampfilename, lampwave, lampinten, lampquality, name, $
     format='D,F,A,A'
  lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0.0d0}
  lamps = replicate(lamps, N_elements(lampwave))
  lamps.lambda = lampwave
  
  nline = n_elements(lamps.lambda)
  lambda = lamps.lambda
  
  wave = lambda
  inten = wave
  quality = lampquality
  name = name
  keep = strarr(nline)
  test = ''
  count = 0
  go = 0
  col = wave*0
  
  for i = 0, nline -1 do begin
     
     if lambda(i) gt 3650 and lambda(i) lt 9150 then begin
        go = 0
        while go ne 1 do begin
           k = where( abs(lam-wave(i)) eq min(abs(lam-wave(i))))
           index = k(0)
           subset = spec(index-3:index+3)
           lset = lam(index-10:index+10)
           u = where(subset eq max(subset))
           inten(i) = subset(u)
           lj = where(spec eq inten(i) and abs(lam-wave(i)) lt 6)
           col(i) = lj(0)
           print, lj(0)
           
           djs_plot, lam, spec, xrange=[wave(i) - 500, wave(i) + 500], $
              yrange=[-10, inten(i) + 500]
           djs_oplot, [lam(lj(0)), lam(lj(0))], [-1e5,1e5], color='red'
           splog, 'Do you want to keep this line? (at wave = ' $
              + string(wave(i))
           splog, 'Yes (Y) or No (N)'
           read, test, prompt='Y or N :    '
           test = strupcase(test)
           if test ne 'Y' and test ne 'N' then begin
              splog, 'Y or N please'
              go = 0
           endif else begin
              go = 1
              keep(i) = test
               count = count + 1
            endelse
            
         endwhile
      endif else keep(i) = 'N'
      
   endfor
   
   
   k = where(keep eq 'Y')
   wave = wave(k)
   inten = inten(k)
   quality = quality(k)
   name = name(k)
   col = col(k)
   
   output = string(wave) + string(inten) + ' '+ quality + ' ' + $
      name + string(col)
   
   openw, 1, 'solution.dat'
   for i = 0, n_elements(output) -1 do begin
      printf, 1, output(i)
   endfor
   close, 1
   
   
END

   
      
     
            
