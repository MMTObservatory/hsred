;+
; NAME:
;   hs_dflatgen
;
; PURPOSE:
;   Create a super dflat image from all of the single dflat images
;
; CALLING SEQUENCE:
;  hs_flatgen, dflatfiles, outdir=outdir, mjd=mjd
;
; INPUTS:
;   dflatfiles   - a linst of dflat files that are to be comined
;
; OPTIONAL KEYWORDS:
;   outdir     - Output directory for the final combined arc
;
; OUTPUTS:
;      creates the /calibration/XXXX/dflat.fits file where
;      XXXX is the rerun number
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
pro hs_dflatgen1, files, outfile=outfile, outdir=outdir, $
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
           imgarr = make_array(dimension=[size(thisimage, $
                                               /dimens), nfile], /float)
           inmask = make_array(dimension=[size(thisimage, $
                                               /dimens), nfile], /byte)
           sxaddpar, hdr0, 'NEDP', nfile, 'Number of exposures in this file'
        endif
        
        imgarr[*,*,ifile]= thisimage
        inmask[*,*,ifile] = thisivar LT 0
        
     endfor
     
     mnimg = djs_avsigclip(imgarr, sigrej=sigrej, maxiter=maxiter, $
                           inmask=inmask)
     splog, 'Writing combined dflat ', $
        djs_filepath(outfile, root_dir = outdir), $
        ' amplifier', string(iamp)
     
     sxdelpar, hdr0, 'EXTNAME'
   
     
     if (iamp eq 1) then begin
        junk  = headfits(files(0), exten=0)
        sxaddpar, junk, 'bitpix', long(-32)
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), 1, junk
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), mnimg, hdr0, /append
     endif else begin
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), mnimg, hdr0, /append
     endelse 
     
  endfor
  
  if keyword_set(gzip) then begin
     
     ;;You gotta GZIP it
     spawn, 'gzip ' + djs_filepath(outfile, root_dir = outdir)
  endif
  
  
  return
  
end

  

PRO hs_dflatgen, dflatfiles,  outdir=outdir, indir=indir, gzip=gzip, mjd=mjd
  
  
  if not keyword_set(outdir) then outdir = cwd()
  if not keyword_set(indir) then indir = cwd()
  
  ;;Grab what we need out of the FITS headers
  
  
  ;;Setup the subregions of the chips
  chip = findgen(4)+1
  stop = 'go'
  
  imagetyp = strarr(n_elements(dflatfiles))

  
  
  ;;Check to make sure that the files are all dflates
  for ifile = 0, n_elements(dflatfiles) - 1l do begin
     head = headfits(dflatfiles(ifile), exten=1)

     imtype = strmid((sxpar(head, 'imagetyp')),0, 8)
     object = strmid((sxpar(head, 'object')), 0, 8)
     if imtype ne 'domeflat' and object ne 'domeflat' then begin
        splog, 'File ' + (dflatfiles(ifile)) + ' is not a domeflat file'
        stop = 'stop'
     endif
     
  endfor
  

  ;;Check to be sure that there are proper numbers of dflates in the directory
  if n_elements(dflatfiles) GT 50 then begin
     splog, 'Too many dflates (' + strn(n_elements(dflatfiles)) + $
        ') to load into memory'
     stop = 'stop'
  endif
  
  
    
  
  if stop ne 'stop' then begin
     dflatname = 'dflat.fits'
     splog, 'Generating Dflat ' + dflatname
     splog, 'Output directory ' + outdir
     hs_dflatgen1, dflatfiles, outfile=dflatname, outdir=outdir, $
        sigrej = sigrej, maxiter=maxiter, gzip=gzip
     
  endif
  
  return
  
end
;--------------------------------------------------------------
     
     
