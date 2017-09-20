Pro hs_plugmap_to_ascii, file, outfile
  
  z = mrdfits(file, 5)
  
  output = replicate('', n_elements(z))
   
  for i = 0, n_elements(output) -1 do begin
     temp = string(z[i].fiberid-1, format='(i3)')
     output(i) = string(temp, format='(a6)') + ' ' + $
        string(z[i].ra, format='(d10.6)') + $
        string(z[i].dec, format='(d10.6)') + $
        string(z[i].fiberid, format='(a12)') +$
        string(z[i].icode, format='(a12)')  + $
        string(z[i].newicode, format='(a12)') + $
        string(z[i].xfocal, format='(a12)')  +  $
        string(z[i].yfocal, format='(a12)')
  endfor
  
  
     
  
  openw, 1, outfile
  
  printf, 1, '#1) Aperture Number'
  printf, 1, '#2-3) Ra/dec
  printf, 1, '#4) Fiberid
  printf, 1, '#5-6) Icode - newicode
  printf, 1, '#7-8) Xfocal/Yfocal'
  for i = 0, n_elements(output) -1 do begin
     printf, 1, output(i)
  endfor
  
  close, 1
  
  
END
