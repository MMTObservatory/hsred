PRO hs_skyplot, color=color
  
  if NOT keyword_set(color) then color='green'
  
  readcol, getenv("HSRED_DIR") + '/etc/skyline_plot.dat', wave, name,$
     format='D,A', comment='#'
  
  for i = 0, n_elements(wave) -1 do begin
     
     djs_oplot, [wave(i), wave(i)], [-1e5,1e6], color=color
     
  endfor
  
return  
END

  
  
  
 
