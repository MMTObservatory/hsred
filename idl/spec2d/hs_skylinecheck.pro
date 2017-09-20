PRO hs_skylinecheck, spec, lam
  
  readcol, getenv("HSRED_DIR") + '/etc/skylines.dat', wave, name, $
     format='D,A', comment='#'
  
  nline = n_elements(wave)
  keep = strarr(nline)
  test= ''
  
  for i = 0, nline -1 do begin
     
     if wave(i) gt 3650 and wave(i) lt 9150 then begin
        go = 0
        while go ne 1 do begin
           k = where( abs(lam-wave(i)) eq min(abs(lam-wave(i))))
           
           djs_plot, lam, spec, xrange=[wave(i) - 500, wave(i) + 500]
           djs_oplot, [lam(k(0)), lam(k(0))], [-1e5,1e5], color='red'
           splog, 'Do you want to keep this line? (at wave = ' + $
              string(wave(i))
           splog, 'Yes (Y) or No (N)'
           read, test, prompt='Y or N :    '
           test = strupcase(test)
           if test ne 'Y' and test ne 'N' then begin
              splog, 'Y or N please'
              go = 0
           endif else begin
              go = 1
              keep(i) = test
           endelse
           
        endwhile
     endif else keep(i) = 'N'
     
  endfor
  k = where(keep eq 'Y')
  wave = wave(k)
  name = name(k)
  output = string(wave) +  ' ' + name 
  
  openw, 1, 'new_skylines.dat'
  for i = 0, n_elements(output) -1 do begin
     printf, 1, output(i)
  endfor
  close, 1
  
  
END

   
      
     
            
