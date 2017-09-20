Pro hs_spZbest_to_ascii, file, outfile
  
  z = mrdfits(file, 1)
  
  output = z.class
  k = where(strmatch(z.subclass, ' *'))

  z(k).subclass='...'
  
  for i = 0, n_elements(output) -1 do begin
     temp = string(z[i].fiberid-1, format='(i3)')
     output(i) = string(temp, format='(a6)') + $
        string(z[i].objtype, format='(a19)') + $
        string(z[i].plug_ra, format='(a12)') + $
        string(z[i].plug_dec, format='(a12)') +$
        string(z[i].class, format='(a12)')  + $
        string(z[i].subclass, format='(a12)') + $
        string(z[i].z, format='(a12)')  +  $
        string(z[i].z_err, format='(a12)') + $
        string(z[i].rchi2, format='(a12)') + $
        string(z[i].zwarning, format='(a12)')
  endfor
  
  
     
  header = 'Fiberid       Objtype       Ra       Dec      Class     ' +$
     'Subclass     Z     Z_err       Chi^2        Zwarning'
  
  header1 ='-------------------------------------' + $
     '----------------------------------------------------------------------'
  
  
  openw, 1, outfile
  
  printf, 1, header
  printf, 1, header1
  for i = 0, n_elements(output) -1 do begin
     printf, 1, output(i)
  endfor
  
  close, 1
  
  
END
