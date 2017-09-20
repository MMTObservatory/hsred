;+
; NAME:
;   hs_arcgen
;
; PURPOSE:
;   Create a super arc image from all of the single arc images
;
; CALLING SEQUENCE:
;  hs_arcgen, arcfiles, outdir=outdir, mjd=mjd
;
; INPUTS:
;   arcfiles   - a linst of arc files that are to be comined
;
; OPTIONAL KEYWORDS:
;   outdir     - Output directory for the final combined arc
;   mjd        - Integer part of the MJD for the observations 
;                this will create the outname for the file if
;                 in the headers of the file
;
; OUTPUTS:
;      creates the /calibration/XXXX/arc-YYYY.fits file where
;      XXXX is the rerun number and YYYY is the MJD
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
;    Not sure the scaling for the images is appropriate
;    is the entrie chip being used or just the illuminated
;    region
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA         
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------

pro hs_arcgen1, files, outfile=outfile, outdir=outdir, $
                 sigrej=sigrej, maxiter=maxiter, gzip=gzip
  
  nfile = n_elements(files)
  if (NOT keyword_set(sigrej)) then begin
     if (nfile) LE 2 then sigrej = 1.0 $
     else if (nfile EQ 3) then sigrej=1.1 $
     else if (nfile EQ 4) then sigrej = 1.3 $
     else if (nfile EQ 5) then sigrej = 1.6 $
     else if (nfile EQ 6) then sigrej = 1.9 $
     else sigrej = 2.0
  endif
  
  if (NOT keyword_set(maxiter)) then maxiter = 3.0
 
  for iamp = 1, 4 do begin
  for ifile =0, nfile-1 do begin
     
        splog, 'Reading file #', ifile+1, ' of ', nfile 
        thisimage = readfits(files[ifile], exten_no=iamp, header, /silent)
        thisivar = 1/thisimage
             
        ;;For the first image, setup the arrays for the combination 
        if (ifile EQ 0) then begin
           hdr0 = header
           imgarr = make_array(dimension=[size(thisimage, /dimens), nfile], $
                               /float)
           inmask = make_array(dimension=[size(thisimage, /dimens), nfile], $
                               /byte)
           sxaddpar, hdr0, 'NEDP', nfile, 'Number of exposures in this file'
        endif
        
        if nfile gt 1 then begin
          imgarr[*,*,ifile]= thisimage
          inmask[*,*,ifile] = thisivar LT 0
        endif else begin
          imgarr[*,*]= thisimage
          inmask[*,*] = thisivar LT 0
        endelse
     endfor
     
     if nfile gt 1 then mnimg = djs_avsigclip(imgarr, sigrej=sigrej, $
                 maxiter=maxiter, inmask=inmask) else mnimg=imgarr
     splog, 'Writing combined arc ', $
       djs_filepath(outfile, root_dir = outdir), $
       ' amplifier', string(iamp)
     
     sxdelpar, hdr0, 'EXTNAME'
     
     
     if (iamp eq 1) then begin
        junk  = headfits(files(0), exten=0)
        sxaddpar, junk, 'bitpix', long(-32)
        writefits, djs_filepath(outfile, root_dir = outdir), 1, junk
        writefits, djs_filepath(outfile, root_dir = outdir), mnimg, $
          hdr0, /append
     endif else begin
        writefits, djs_filepath(outfile, root_dir = outdir), mnimg, $
          hdr0, /append
     endelse 
     
  endfor
  
  if keyword_set(gzip) then begin
     
     ;;You gotta GZIP it
     spawn, 'gzip ' + djs_filepath(outfile, root_dir = outdir)
  endif
  
  
  return
  
end

  

PRO hs_arcgen, arcfiles,  outdir=outdir, indir=indir, gzip=gzip, mjd=mjd
  
  
  if not keyword_set(outdir) then outdir = cwd()
  if not keyword_set(indir) then indir = cwd()
  
  ;;Grab what we need out of the FITS headers
  
  
  ;;Setup the subregions of the chips
  chip = findgen(4)+1
  stop = 'go'
  
  imagetyp = strarr(n_elements(arcfiles))

  
  
  ;;Check to make sure that the files are all arces
  for ifile = 0, n_elements(arcfiles) - 1l do begin
     head = headfits(arcfiles(ifile), exten=1)

     imtype = strmid((sxpar(head, 'imagetyp')),0, 4)
     if imtype ne 'arc' and imtype ne 'comp' then begin
        splog, 'File ' + (arcfiles(ifile)) + ' is not a arc file'
        stop = 'stop'
     endif
     
  endfor
  

  ;;Check to be sure that there are proper numbers of arces in the directory
  if n_elements(arcfiles) GT 50 then begin
     splog, 'Too many arces (' + strn(n_elements(arcfiles)) + $
        ') to load into memory'
     stop = 'stop'
  endif
  
  
    
  
  if stop ne 'stop' then begin
     
     arcname = 'arc.fits'
     splog, 'Generating Arc ' + arcname
     splog, 'Output directory ' + outdir
     hs_arcgen1, arcfiles, outfile=arcname, outdir=outdir, $
        sigrej = sigrej, maxiter=maxiter, gzip=gzip
     
  endif
  
  return
  
end
;--------------------------------------------------------------
     
     
