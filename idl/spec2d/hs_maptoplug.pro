;+
; NAME:
;   hs_maptoplug
;
; PURPOSE:
;  Create a plugmap structure for the Hectospec generated map files
;
; CALLING SEQUENCE:
;  hs_maptoplug, mapfile, plugmap, ccdnum=ccdnum
;
; INPUTS:
;   mapfile - file to use to create plugmap
;   
;
; OPTIONAL KEYWORDS:
;   ccdnum  - ccdnumber to make plugmap for
;
; OUTPUTS:
;    plugmap - output stucture
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;   June 2014  - Revised to read everything needed, including
;                standard star info from standardstar.dat
;-
;------------------------------------------------------------------------------
PRO hs_maptoplug, mapfile, plugmap, ccdnum=ccdnum, dummy=dummy, chelle=chelle
  
  if keyword_set(chelle) then begin
   c1_max=119
   c2_min=120
   c2_max=239
  endif else begin
   c1_max=149
   c2_min=150
   c2_max=299
  endelse
  if not keyword_set(dummy) then $
    readcol, mapfile, aperture, id, target, rain, decin, $
     catid, fiber, xpos, ypos, $
     format='L,L,A,A,A,A,L,D,D', /silent else begin
;  readcol, mapfile, aperture, id, target, fiber, xpos, ypos, $
;     format='L,L,A,L,D,D', /silent
    readcol, mapfile, aperture, id, target, rain, decin, $
     catid, fiber, xpos, ypos, $
     format='L,L,A,A,A,A,L,D,D', /silent
    if n_elements(aperture) eq 0 then $
       readcol, mapfile, aperture, id, target, fiber, xpos, ypos, $
       format='L,L,A,L,D,D', /silent
    rain=strarr(c2_max+1)  
    decin=rain
    catid=intarr(c2_max+1)
    rain[*]='00:00:00.0'
    decin[*]='00:00:00.0'
  endelse
  
  if ccdnum eq 1 then begin
     aperture = aperture[0:c1_max]
     rain = rain[0:c1_max]
     id = id[0:c1_max]
     decin = decin[0:c1_max]
     fiber = fiber[0:c1_max]
     catid = catid[0:c1_max]
     xpos = xpos[0:c1_max]
     ypos = ypos[0:c1_max]
     target = target[0:c1_max]
  endif else if ccdnum eq 2 then begin
     aperture = aperture[c2_min:c2_max]
     rain = rain[c2_min:c2_max]
     id = id[c2_min:c2_max]
     decin = decin[c2_min:c2_max]
     fiber = fiber[c2_min:c2_max]
     catid = catid[c2_min:c2_max]
     xpos = xpos[c2_min:c2_max]
     ypos = ypos[c2_min:c2_max] 
     target=target[c2_min:c2_max]
  endif else begin
    print,  'CCDNUM is 1 or 2'
  endelse
  
  targetup=target
  ra = dindgen(n_elements(rain))
  dec = dindgen(n_elements(decin))
  n = n_elements(rain)
  
  for i = 0, n_elements(rain) -1 do begin
     targetup[i] = strupcase(target[i])
     split = strsplit(rain[i], ':', /extract)
     ra[i] = ten(split[0], split[1], split[2])*15
     split = strsplit(decin[i], ':', /extract)
     dec[i] = ten(split[0], split[1], split[2])
  endfor


;;Changed to look at standardstar file
readcol, getenv('HSRED_DIR') + '/etc/standstar.dat', sra, sdec, format='D,D'
     
spherematch, sra, sdec, ra, dec, 1.0/3600., m1, m2, dist
  
  IF m2[0] NE -1 THEN  targetup[m2] = 'SPECTROPHOTO_STD'
  k = where(strmatch(targetup, 'SPE*'), sp_cnt)
  mag = dblarr(150,5)
  BD17_color = 9.35 - [10.56, 9.64,9.35,9.25, 9.23]
  
  ;;For the standard stars, loop through them and
  ;;find the appropriate standstar 
  
  if sp_cnt ge 1 then begin
     standfile = getenv('HSRED_DIR')+'/etc/standstar.dat'
     readcol, standfile, sra, sdec, u, g, r, i, z, $
        reddening, format='D,D,D,D,D,D,D,D', $
        comment='#', /silent
     
     spherematch, ra[k], dec[k], sra, sdec, 1.0/3600., match1, match2, dist
          
     reddening2 = reddening
     reddening0 = reddening2*1.874
     reddening1 = reddening2*1.379
     reddening3 = reddening2*0.758
     reddening4 = reddening2*0.538
     
     mag[k[match1], 0] = u[match2] + reddening0[match2]
     mag[k[match1], 1] = g[match2] + reddening1[match2]
     mag[k[match1], 2] = r[match2] + reddening2[match2]
     mag[k[match1], 3] = i[match2] + reddening3[match2]
     mag[k[match1], 4] = z[match2] + reddening4[match2]     
  endif

  
  create_struct, plugmap, 'plugmap', $
     ['objtype', 'ra', 'dec', 'fiberid', $
      'rmag', 'rapmag', 'icode', 'rcode', 'bcode',$
      'mag', 'xfocal','yfocal', 'frames'], $
     'A,D,D,D,D,D,D,D,A,D(5),D,D,L', dimen=c2_min
  
  for i = 0, n  -1 do begin
     plugmap[i].objtype = targetup[i]
     plugmap[i].ra = ra[i]
     plugmap[i].dec = dec[i]
     plugmap[i].fiberid = aperture[i]
     plugmap[i].xfocal = xpos[i]
     plugmap[i].yfocal = ypos[i]
     plugmap[i].rmag = mag[i,2]
     plugmap[i].mag = mag[i,*]
     plugmap[i].rapmag = mag[i,2]
     plugmap[i].icode =catid[i]
     plugmap[i].rcode =fiber[i]
     plugmap[i].bcode =target[i]
     plugmap[i].frames= 0.0
  endfor
  
  
return
END
