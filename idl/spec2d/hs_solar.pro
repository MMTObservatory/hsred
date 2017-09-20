PRO hs_solar, z=z
 
  if NOT keyword_set(z) then z = 0
  lamH=[4340,     4101,     3970,    3889,    3835, 6563] * (1+z)
  txtH=['H\gg\n', 'H\gd\n', 'H\ge\n', 'H8',   'H9', textoidl('H\alpha')]
  for i = 0, 5 do begin
     
     djs_oplot, [lamH(i), lamH(i)], [-1e5,1e10], color='red'
  endfor
  
  ; set Ca line labels and text
  lamCa=[3934,  3968,   4227 ]*(1+z)
  
    for i = 0, 2 do begin
     
     djs_oplot, [lamCa(i), lamCa(i)], [-1e5,1e10], color='blue'
  endfor
  
return  
END

  
  
  
