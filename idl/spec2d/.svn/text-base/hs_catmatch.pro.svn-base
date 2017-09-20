PRO hs_catmatch, mapfile, catalog=catalog,outname=outname
  
  if NOT keyword_set(catalog) then catalog='/home/rcool/bootes.cat'
  
  ;;Read the targets
  readcol, catalog, ra, dec, type, rank, rmag, rapmag, bcode, $
     icode, rcode, junk, junk,junk, format='D,D,A,D,D,D,A,L,L,A,A,A'
  
  ;;Read the skies
  readcol, catalog+'_sky', ras, decs, junk, junk, junk, junk, $
     junk1, format='D,D,A,A,A,A,A'
      
   ;;Read the map
  readcol, mapfile, fiber1, junk, target1, ra1, dec1, junk, junk, $
     junk,junk, format='L,A,A,A,A,A,A,A,A'
  
  k = where(strmatch(target1,'target'))
  fiber1 = fiber1(k)
  target1=target1(k)
  ra1=ra1(k)
  dec1=dec1(k)
  ra1 = hms2dec(ra1)
  dec1 = hms2dec(dec1)
  
  
  ;;Now match the targets
  spherematch, ra, dec, ra1, dec1, 1.0/3600., match1, match2,dist
  ra2 = dindgen(300)
  dec2 = dindgen(300)
  target2 = strarr(300)
  dmin2= lindgen(300)*0.0
  rmag2 = dindgen(300)
  rapmag2 = dindgen(300)
  bcode2 = strarr(300)
  icode2 = lindgen(300)
  rcode2 = lindgen(300)
  target = strarr(300)
  fiber = lindgen(300)
  
  ct = n_elements(match1)
  
  ra2(0:ct-1) = ra(match1)
  dec2(0:ct-1) = dec(match1)
  fiber(0:ct-1) = fiber1(match2)
  rmag2(0:ct-1) = rmag(match1)
  rapmag2(0:ct-1) = rapmag(match1)
  bcode2(0:ct-1) = bcode(match1)
  icode2(0:ct-1) = icode(match1)
  rcode2(0:ct-1) = rcode(match1)
  dmin2(0:ct-1) = 0.0
  target(0:ct-1) = 'target'
  
  
  ;;Now match the skies
  spherematch, ra1, dec1, ras, decs, 1.0/3600., match1, match2, dist
  cts = n_elements(match1)
  
  ra2(ct:cts+ct-1) = ra(match2)
  dec2(ct:cts+ct-1) = dec(match2)
  fiber(ct:cts+ct-1) = fiber1(match1)
  rmag2(ct:cts+ct-1) = 0
  rapmag2(ct:cts+ct-1) = 0
  bcode2(ct:cts+ct-1) = 0
  icode2(ct:cts+ct-1) = 0
  rcode2(ct:cts+ct-1) = 0
  dmin2(ct:cts+ct-1) = 0
  target(ct:cts+ct-1) = 'sdsssky'

  
  ;;Now find the unused
    ;;Read the map
  
  readcol, mapfile, fiber1, junk, target1, ra1, dec1, junk, junk, $
     junk,junk, format='L,A,A,A,A,A,A,A,A'
  
  k = where(strmatch(target1,'unused'))
  fiber1 = fiber1(k)
  target1=target1(k)
  ra1=ra1(k)
  dec1=dec1(k)
  
  k = where(strmatch(target1, 'unused*'))
  
  if k(0) ne -1 then begin
     ctu = n_elements(k)
  
     ra2(cts+ct:cts+ct+ctu-1) = 0.0
     dec2(cts+ct:cts+ct+ctu-1) = 0.0
     fiber(cts+ct:cts+ct+ctu-1) = fiber1(k)
     rmag2(cts+ct:cts+ct+ctu-1) = 0
     rapmag2(cts+ct:cts+ct+ctu-1) = 0
     bcode2(cts+ct:cts+ct+ctu-1) = 0
     icode2(cts+ct:cts+ct+ctu-1) = 0
     rcode2(cts+ct:cts+ct+ctu-1) = 0
     dmin2(cts+ct:cts+ct+ctu-1) = 0
     target(cts+ct:cts+ct+ctu-1) = 'unused'
     
  endif else ctu = 0
                                ;Read the map
  readcol, mapfile, fiber1, junk, target1, ra1, dec1, junk, junk, junk,$
     junk, format='L,A,A,A,A,A,A,A,A'
  
  k = where(strmatch(target1,'sky*'))
  fiber1 = fiber1(k)
  target1=target1(k)
  ra1=ra1(k)
  dec1=dec1(k)
  
  k = where(strmatch(target1,'sky*'))
  

  
  
  if k(0) ne -1 then begin
     cts1 = n_elements(k)
     ra2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0.0
     dec2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0.0
     fiber(ctu+cts+ct:ctu+ct+cts-1+cts1) = fiber1(k)
     rmag2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
     rapmag2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
     bcode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
     icode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
     rcode2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
     dmin2(ctu+cts+ct:ctu+ct+cts-1+cts1) = 0
     target(ctu+cts+ct:ctu+ct+cts-1+cts1) = 'sky'
     
  endif else ctu = 0
  
  sort1 = sort(fiber)
  
  fiber = fiber(sort1)
  dmin = dmin2(sort1)
  ra = ra2(sort1)
  dec = dec2(sort1)
  rmag = rmag2(sort1)
  rapmag = rapmag2(sort1)
  bcode = bcode2(sort1)
  icode = icode2(sort1)
  rcode = rcode2(sort1)
  target = target(sort1)
  
  
  
  output = string(fiber) + ' ' + string(dmin) + ' '+ string(ra) + ' ' + $
     string(dec) + ' ' + string(rmag) + ' ' + string(rapmag) + ' ' +$
     string(bcode) + ' ' + string(icode) + ' ' +string(rcode)
  
  
  openw, 1, 'plugdir/' + outname+'_all'
  for i = 0, n_elements(output) -1 do begin
     printf, 1, output(i)
  endfor
  
  close, 1
  
  k = where(strmatch(target,'target'))
  if k(0) ne -1 then begin $
     openw, 1, 'plugdir/' + outname+'_targets'
     for i = 0, n_elements(k) -1 do begin
        printf, 1, output(k(i))
     endfor
  endif
  
  output = string(fiber) +'  '+ target
  
  
  close, 1
  k = where(strmatch(target,'*sky*'))
  if k(0) ne -1 then begin $
     openw, 1, 'plugdir/' + outname+'_sky'
     for i = 0, n_elements(k) -1 do begin
        printf, 1, output(k(i))
     endfor
  endif
  close, 1
  
  k = where(strmatch(target,'unused'))
  if k(0) ne -1 then begin $
     openw, 1, 'plugdir/' + outname+'_unused'
     for i = 0, n_elements(k) -1 do begin
        printf, 1, output(k(i))
     endfor
  endif
  close, 1


END


