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
                     ystart=ystart,nfiber=nfiber, chelle=chelle, outfile=outfile
  
  if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if NOT keyword_set(flatfile) then begin

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
     
     hs_proc, flatfile, ccdnum, flatim, flativar, rerun=rerun, $
        chelle=chelle, header=header

     IF ccdnum EQ 1 AND keyword_set(chelle) THEN BEGIN
        binning=strsplit(sxpar(header,'CCDSUM'), ' ')
        binx=fix(binning[0])
        biny=fix(binning[1])
        ystart = 2000/biny
        nfiber = 120
        ;transpose the flatim so that we start counting 120 good traces
        ;from the center, not the edge
        flatim=rotate(flatim,5)
        flativar=rotate(flativar,5)
       ;npbundle= 1
       ;deltax = 13/binx
     endif
     
     IF ccdnum EQ 2 AND keyword_set(chelle) THEN BEGIN
        binning=strsplit(sxpar(header,'CCDSUM'), ' ')
        binx=fix(binning[0])
        biny=fix(binning[1])
        nfiber = 120
        ystart=2000/biny
        delvarx, npbundle
        delvarx, deltax
     endif
     
     xsol = trace150crude(flatim, flativar, yset=ycen, ystart=ystart, $
                          nfiber=nfiber, chelle=chelle, deltax=deltax, $
                         npbundle=npbundle)
   ;  xy2traceset, ycen, xsol, tset, ncoeff=5
   ;  traceset2xy, tset, ycen, xsol
     IF ccdnum EQ 1 AND keyword_set(chelle) THEN BEGIN
        ;if we transposed chelle chip 1, need to flip the xvals back
        output=reverse((size(flatim))[1]-xsol-1,2)
     endif else output = xsol

     if not keyword_set(outfile) then begin
        outfile='calibration/'+rerun+'/traceset.fits'
        outfile_ycen='calibration/'+rerun+'/traceset_ycen.fits'
     endif else begin
        outfile=outfile
        outfile_ycen=(strsplit(outfile,'.',/extract))[0]+'_ycen.fits'
     endelse     
     if keyword_set(write) then begin
        if ccdnum eq 1 then begin
           mwrfits, output, outfile, /create
           mwrfits, ycen, outfile_ycen, /create
        endif
        
        if ccdnum eq 2 then begin
           mwrfits, output, outfile
           mwrfits, ycen, outfile_ycen
        endif
     endif
  endfor
  return
  
END

  
  
  
