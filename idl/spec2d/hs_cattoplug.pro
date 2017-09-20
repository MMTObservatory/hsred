;+
; NAME:
;   hs_cattoplug
;
; PURPOSE:
;   Make the plumap structure from the plugcat files
;
; CALLING SEQUENCE:
;  hs_cattoplug, root, plugmap, ccdnum=ccdnum
;
; INPUTS:
;   root = the root name for the files containing the plugcat 
;
;
; OPTIONAL KEYWORDS:
;   ccdnum    - set for which ccd number to use
;
; OUTPUTS:
;   plugmap    - output plugmap structure
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;    You also need the _map files from the hectospec output as well
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
;-
;------------------------------------------------------------------------------
PRO hs_cattoplug, root, plugmap, ccdnum=ccdnum, oldmap=oldmap
   
   IF keyword_set(oldmap) THEN BEGIN
      hs_cattoplug_orig, root, plugmap, ccdnum=ccdnum
      return
   ENDIF
   
   
  if NOT keyword_set(ccdnum) then begin
     splog, 'You must specify the ccdnumber'
     return
  endif
  
  
  ap = long(fltarr(150))
  dmin = dblarr(150)
  ra = dmin
  dec = dmin
  rmag = dmin
  rapmag = dmin
  bcode = (strarr(150))
  icode = long(dblarr(150))
  newicode = long(dblarr(150))
  rcode = long(dblarr(150))
  objtype = strarr(150)
  
  
  
  ;;get the target
  
  test = rc_readfile(root+'_targets', 0)
  testline = test(2)
  junk = strsplit(testline, ' ', /extract)
  
  IF n_elements(junk) EQ 10 THEN begin
     
     readcol, root + '_targets', ap1, dmin1, ra1, dec1, rmag1,$
        rapmag1, bcode1, icode1, rcode1,newicode1, $
        format='L,D,D,D,D,D,A,D,D,D', comment='#', /silent
  endIF ELSE BEGIN
     readcol, root + '_targets', ap1, dmin1, ra1, dec1, rmag1,$                                             
        rapmag1, bcode1, icode1, rcode1, $                                                         
        format='L,D,D,D,D,D,A,D,D', comment='#', /silent     
     newicode1 = icode1
  endelse
  
  
  
  if ccdnum eq 2 then k = where(ap1 gt 150)
  if ccdnum eq 1 then k = where(ap1 le 150)
  
  if k(0) ne -1 then begin
     ntarg = n_elements(ap1(k))
     ap(0:ntarg-1) = ap1(k)
     dmin(0:ntarg-1) = dmin1(k)
     ra(0:ntarg-1) = ra1(k)
     dec(0:ntarg-1) = dec1(k)
     rmag(0:ntarg-1) = rmag1(k)
     rapmag(0:ntarg-1) = rapmag1(k)
     bcode(0:ntarg-1) = bcode1(k)
     icode(0:ntarg-1) = icode1(k)
     newicode(0:ntarg-1) = newicode1(k)
     rcode(0:ntarg-1) = rcode1(k)
     objtype(0:ntarg-1)= 'TARGET'
  endif else ntarg = 0
  
  
  ;;get the sky fibers
  
  readcol, root+ '_sky', ap2, junk, format='L,A,D,D', comment='#', /silent
  
  
  if ccdnum eq 2 then k = where(ap2 gt 150)
  if ccdnum eq 1 then k = where(ap2 le 150)
  
  ap1 = ap2
  if k(0) ne -1 then begin
     nsky = n_elements(ap1(k))
  
     
     ap(ntarg:ntarg+nsky-1) = ap1(k)
     dmin(ntarg:ntarg+nsky-1) =0
     ra(ntarg:ntarg+nsky-1) = 0
     dec(ntarg:ntarg+nsky-1) = 0
     rmag(ntarg:ntarg+nsky-1) = 0
     rapmag(ntarg:ntarg+nsky-1) = 0
     bcode(ntarg:ntarg+nsky-1) = '000000000000'
     icode(ntarg:ntarg+nsky-1) = 0
     newicode(ntarg:ntarg+nsky-1) = 0
     rcode(ntarg:ntarg+nsky-1) = 0
     objtype(ntarg:ntarg+nsky-1) = 'SKY'
     
     
  endif else nsky = 0 
  
  ;;get the unused fibers
  
   readcol, root+ '_unused', ap3, junk, format='L,A', comment='#', /silent
  
  if n_elements(ap3) gt 0 then begin
  if ccdnum eq 2 then k = where(ap3 gt 150)
  if ccdnum eq 1 then k = where(ap3 le 150)
  ap1 = ap3
  
  if k(0) ne -1 then begin
     none = n_elements(ap1(k))
     
     ap(ntarg+nsky:ntarg+nsky+none-1) = ap1(k)
     
     dmin(ntarg+nsky:ntarg+nsky+none-1) =0
     ra(ntarg+nsky:ntarg+nsky+none-1) = 0
     dec(ntarg+nsky:ntarg+nsky+none-1) = 0
     rmag(ntarg+nsky:ntarg+nsky+none-1) = 0
     rapmag(ntarg+nsky:ntarg+nsky+none-1) =0
     bcode(ntarg+nsky:ntarg+nsky+none-1) = '00000000000000'
     icode(ntarg+nsky:ntarg+nsky+none-1) = 0
     newicode(ntarg+nsky:ntarg+nsky+none-1) = 0
     rcode(ntarg+nsky:ntarg+nsky+none-1) = 0
     objtype(ntarg+nsky:ntarg+nsky+none-1) = 'UNUSED'
     
     
  endif else none = 0 
endif

 
  index = sort(ap)
 
  ap = ap(index)
  
  dmin = dmin(index)
  ra = ra(index)
  dec = dec(index)
  rmag = rmag(index)
  rapmag = rapmag(index)
  bcode = bcode(index)
  icode = icode(index)
  newicode = newicode(index)
  rcode = rcode(index)
  objtype = objtype(index)    
  
  
  ;;Changed to look at standardstar file
  readcol, getenv('HSRED_DIR') + '/etc/standstar.dat', sra, sdec, format='D,D'
     
  spherematch, sra, sdec, ra*15.0, dec, 1.0/3600., m1, m2, dist
  
  IF m2(0) NE -1 THEN  objtype(m2) = 'SPECTROPHOTO_STD'
  k = where(strmatch(objtype, 'SPE*'))
  mag = dblarr(150,5)
  BD17_color = 9.35 - [10.56, 9.64,9.35,9.25, 9.23]
  
  
  ;;For the standard stars, loop through them and
  ;;find the appropriate standstar 
  
  if k(0) ne -1 then begin
     standfile = getenv('HSRED_DIR')+'/etc/standstar.dat'
     readcol, standfile, sra, sdec, u, g, r, i, z, $
        reddening, format='D,D,D,D,D,D,D,D', $
        comment='#', /silent
     
     spherematch, ra(k)*15, dec(k), sra, sdec, 1.0/3600., match1, match2, dist
     
     reddening2 = reddening
     reddening0 = reddening2*1.874
     reddening1 = reddening2*1.379
     reddening3 = reddening2*0.758
     reddening4 = reddening2*0.538
     
     
     mag(k(match1), 0) = u(match2) + reddening0(match2)
     mag(k(match1), 1) = g(match2) + reddening1(match2)
     mag(k(match1), 2) = r(match2) + reddening2(match2)
     mag(k(match1), 3) = i(match2) + reddening3(match2)
     mag(k(match1), 4) = z(match2) + reddening4(match2)
     
  endif
 
 
  split = strsplit(root, '/', /extract)
  mapname = findfile(split(n_elements(split) -1) + $
                     '*map*' +string(ccdnum, format='(i1.1)'))
  mapname = mapname(0)
  
  readcol, mapname, apjunk, beamjunk, objectjunk, rajunk, $
     decjunk, targetjunk, fiberjunk, platex, $
     platey, format='L,L,A,A,A,L,L,F,F', /silent
  
  
  
  delvarx, plugmap
  plugmap = {objtype : objtype(0), $
             ra : ra(0), $
             dec : dec(0), $
             fiberid : ap(0), $
             rmag : rmag(0), $
             mag : reform(mag(0,*)), $
             rapmag : rapmag(0), $
             icode : icode(0), $
             newicode : newicode(0), $
             rcode : rcode(0), $
             bcode : bcode(0), $
             xfocal : platex(0), $
             yfocal : platey(0), $
             frames : 0.0}
  
  plugmap = replicate(plugmap, 150)
  
  plugmap.objtype = objtype
  plugmap.ra = ra*15
  plugmap.dec = dec
  plugmap.fiberid = ap
  plugmap.rmag = rmag
  plugmap.mag = transpose(mag)
  plugmap.rapmag = rapmag
  plugmap.icode = icode
  plugmap.newicode = newicode
  plugmap.rcode = rcode
  plugmap.bcode = bcode
  plugmap.xfocal = platex
  plugmap.yfocal = platey
  plugmap.frames = replicate(0.0, 150)
    
  
  return
end

