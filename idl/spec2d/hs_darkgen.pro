;+
; NAME:
;   hs_darkgen
;
; PURPOSE:
;   Create a super dark image from all of the single dark images
;
; CALLING SEQUENCE:
;  hs_darkgen, darkfiles, outdir=outdir, mjd=mjd
;
; INPUTS:
;   darkfiles   - a linst of dark files that are to be comined
;
; OPTIONAL KEYWORDS:
;   outdir     - Output directory for the final combined dark
;   model      - set to fit a median-smoothed model to the dark, rather than
;                saving the actual dark frame (lowers noise in the
;                dark subtraction)
;
; OUTPUTS:
;      creates the /calibration/XXXX/dark.fits file where
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
pro hs_darkgen1, files, outfile=outfile, outdir=outdir, model=model
  
  nfile = n_elements(files)
  if (NOT keyword_set(sigrej)) then begin
     if (nfile) LE 2 then sigrej = 1.0 $
     else if (nfile EQ 3) then sigrej=1.1 $
     else if (nfile EQ 4) then sigrej = 1.3 $
     else if (nfile EQ 5) then sigrej = 1.6 $
     else if (nfile EQ 6) then sigrej = 1.9 $
     else sigrej = 2.0
  endif
  

  if (NOT keyword_set(maxiter)) then maxiter = 10.0

  for iamp = 1, 4 do begin
     
     darktime = findgen(nfile)
     for ifile =0, nfile-1 do begin
     
        splog, 'Reading file #', ifile+1, ' of ', nfile 
        thisimage = readfits(files[ifile], exten_no=iamp, header, /silent)
        
        ;;Overscan and trim correct the bias image
        
        thisivar = 1/thisimage
        darktime(ifile) = sxpar(header, 'DARKTIME')
        
        ;;For the first image, setup the arrays for the combination 
        if (ifile EQ 0) then begin
           hdr0 = header
           imgarr = make_array(dimension=[size(thisimage, $
                                               /dimens), nfile], /float)
           imgarr2=imgarr
           inmask = make_array(dimension=$
                               [size(thisimage, /dimens), nfile], /byte)
           sxaddpar, hdr0, 'NEDP', nfile, 'Number of exposures in this file'
        endif
        
        imgarr[*,*,ifile]= thisimage 
        imgarr2[*,*,ifile]= thisimage 
        inmask[*,*,ifile] = thisivar LT 0
        
     endfor
     
     ;;Mean dark time
     dtime_mean = djs_mean(darktime)
     scale = dtime_mean/darktime
     for ifile = 0, nfile -1 do begin
        if keyword_set(model) then imgarr[*,*,ifile] = median(reform(imgarr[*,*,ifile]*scale(ifile)),4) else imgarr[*,*,ifile] = imgarr[*,*,ifile]*scale(ifile)
        imgarr2[*,*,ifile] = imgarr2[*,*,ifile]*scale(ifile)
     endfor
     
     mnimg = djs_avsigclip(imgarr, sigrej=sigrej, maxiter=maxiter, $
                           inmask=inmask)
     mnimg2 = djs_avsigclip(imgarr2, sigrej=sigrej, maxiter=maxiter, $
                           inmask=inmask)
     ;chelle dark creation per NCaldwell proceedure
     if keyword_set(model) then begin
        mnimg=median(mnimg,4)
        dmodel=mnimg*0.0
        sz=size(mnimg)
        nrow=sz[1]
        ncol=sz[2]
        xvals=findgen(ncol)
       ; stop
        for zz=0,nrow-1 do begin
           yvals=reform(mnimg[zz,*])
           err=bspline_iterfit(xvals,yvals,yfit=yfit, everyn=20, /eachgroup)
           dmodel[zz,*]=yfit
        endfor
        mnimg=dmodel
     endif
     ;stop
     splog, 'Writing combined bias ', $
        djs_filepath(outfile, root_dir = outdir), $
        ' amplifier', string(iamp)
     
     sxdelpar, hdr0, 'EXTNAME'
     sxdelpar, hdr0, 'BSCALE'
     sxdelpar, hdr0, 'BZERO'
     
     if (iamp eq 1) then begin
        junk  = headfits(files(0), exten=0)
        sxaddpar, junk, 'bitpix', long(-32)
        sxaddpar, hdr0, 'darktime', dtime_mean
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), 1, junk
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), mnimg, hdr0, /append
     endif else begin
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), mnimg, hdr0, /append
     endelse 
     
  endfor
  
   
  return
  
end

  

PRO hs_darkgen, biasfiles,  outdir=outdir, indir=indir, mjd=mjd, model=model
  
  
  if not keyword_set(outdir) then outdir = cwd()
  if not keyword_set(indir) then indir = cwd()
  
  ;;Grab what we need out of the FITS headers
  
  
  ;;Setup the subregions of the chips
  chip = findgen(4)+1
  stop = 'go'
  imagetyp = strarr(n_elements(biasfiles))

  ;;Check to make sure that the files are all biases
  for ifile = 0, n_elements(biasfiles) - 1l do begin
     head = headfits(biasfiles(ifile), exten=1)

     imtype = strmid((sxpar(head, 'imagetyp')),0, 4)
  endfor
  
  
  ;;Check to be sure that there are proper numbers of biases in the directory
  if n_elements(biasfiles) GT 50 then begin
     splog, 'Too many darks (' + strn(n_elements(biasfiles)) + $
        ') to load into memory'
     stop = 'stop'
  endif
  
  
  if n_elements(biasfiles) LT 3 then begin
     splog, 'Expected at least 3 darks, got ' + $ 
        strn(n_elements(biasfiles))
     stop = 'stop'
  endif
  
  
  if stop ne 'stop' then begin
     biasname = 'dark.fits'
     splog, 'Generating Dark ' + biasname
     splog, 'Output directory ' + outdir
     hs_darkgen1, biasfiles, outfile=biasname, outdir=outdir, model=model
     
  endif
  
  return
  
end
;--------------------------------------------------------------
     
     
