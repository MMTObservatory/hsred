;+
; NAME:
;   hs_maketraceset
;
; PURPOSE:
;   Create a trace set for hectospec data set
;
; CALLING SEQUENCE:
;  hs_maketraceset, flatfile=flatfile, write=write, xcen=xcen, 
;                   tset=tset, xsol=xsol, rerun=rerun, mjd=mjd
;
; INPUTS:
;   
;
; OPTIONAL KEYWORDS:
;   flatfile - file to make the flat vectors from
;   write    - Turn on if you want to write the data set out
;
; OUTPUTS:
;    
;
; OPTIONAL OUTPUTS:
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
;-
;------------------------------------------------------------------------------
PRO hs_maketraceset, flatfile=flatfile, rootsearch=rootsearch, $
                     write=write, xcen=xcen, $
                     tset=tset, xsol=xsol,rerun=rerun, mjd=mjd, $
                     ystart=ystart,nfiber=nfiber, chelle=chelle
  
  if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if NOT keyword_set(flatfile) then begin
     file = findfile('*fit*')
     N = n_elements(file)/2
     header = headfits(file(N), exten=1,  /silent)
     if not keyword_set(mjd) then mjd = long(sxpar(header, 'MJD'))
     flatfile = 'calibration/'+rerun+'/dflat.fits'
     if file_test(flatfile) eq 0 then begin
        flatfile = findfile('calibration/' + rerun + '/dflat-*fits')
     endif   
  endif
  
  if file_test(flatfile) eq 0 then begin
     splog, 'This flatfile does not exist'
     return
  endif
  
  write = 1.0
  IF keyword_set(chelle) THEN nfiber = 125
  IF keyword_set(chelle) THEN ystart = 200
  for ccdnum = 1, 2 do begin
     
     IF ccdnum EQ 1 AND keyword_set(chelle) THEN BEGIN
        ystart = 2000
        nfiber = 124
       npbundle= 1
       deltax = 13
     endif
     
     IF ccdnum EQ 2 AND keyword_set(chelle) THEN BEGIN
        nfiber = 120
        ystart=2000
        delvarx, npbundle
        delvarx, deltax
     endif
     
     
        
     hs_proc, flatfile, ccdnum, flatim, flativar, rerun=rerun, root=mjd, $
        chelle=chelle
     xsol = trace150crude(flatim, flativar, yset=ycen, ystart=ystart, $
                          nfiber=nfiber, chelle=chelle, deltax=deltax, $
                         npbundle=npbundle)
   ;  xy2traceset, ycen, xsol, tset, ncoeff=5
   ;  traceset2xy, tset, ycen, xsol
     output = xsol
     
     if keyword_set(write) then begin
        if ccdnum eq 1 then begin
           mwrfits, output, 'calibration/'+rerun+'/traceset.fits', /create
           mwrfits, ycen, 'calibration/' + rerun + '/traceset_ycen.fits', /create
        endif
        
        if ccdnum eq 2 then begin
           mwrfits, output, 'calibration/'+rerun+'/traceset.fits'
           mwrfits, ycen, 'calibration/'+rerun+'/traceset_ycen.fits'
        endif
     endif
  endfor
  return
  
END

  
  
  
