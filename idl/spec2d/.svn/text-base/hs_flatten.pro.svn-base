;+
; NAME:
;   hs_flatten
;
; PURPOSE:
;   Create pixel-to-pixel flat-field from a stack of Hectospec spectral flats.
;
; CALLING SEQUENCE:
;   spflatten2, flatname, allflats, [ pixflat, sigrej=, maxiter=, $
;    oldflat=, outfile=, indir=, outdir=, tmpdir=, $
;    pixspace=, nord=, lower=, upper=, mincounts=, /nodelete ]
;
; INPUTS:
;   flatname   - Name of flat image for tracing arc
;  
;   allflats   - Name(s) of raw SDSS flat-field image(s).
;                Note that many flats from many nights can be combined.
;
; OPTIONAL INPUTS:
;   sigrej     - Sigma rejection level; default to 1, 1, 1.1, 1.3, 1.6 or 1.9
;                for 1,2,3,4,5 or 6 flats.  For more then 6 flats, default
;                to 2.0.
;   maxiter    - Number of rejection iterations; default to 2.
;   oldflat    - Name of old flat-field from which to select pixels to mask;
;                if not set, then read the masks in this source distribution.
;   outfile    - Write the image PIXFLAT to this file.
;   indir      - Input directory for FLATNAME; default to current directory
;   outdir     - Output directory for OUTFILE; default to current directory
;   tmpdir     - Directory for temporary files; default to current directory
;   pixspace   - Approximate spacing in pixels for break points in the
;                spline fits to individual fibers; default to 50 pixels.
;                Note that this spacing need be small enough to fit out
;                the flux variations of multiple fibers crossing a single
;                column, but it need not fit out the actual flat-field
;                variations.
;   mincounts  - The summed model image must have a minimum of this many counts
;                at each pixel, or the output flat is set to zero for that
;                pixel; default to 100 counts.
;   nodelete   - If set, then do not delete temporary files.
;
; PARAMETERS FOR SLATEC_SPLINEFIT:
;   nord       - Default to 4
;   lower      - Default to 2
;   upper      - Default to 2
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   pixflat    - Image containing all the information about pixel-to-pixel
;                variations.  Illumination variations are removed.
;
; COMMENTS:
;   This program writes 2*nflat temporary files to disk to save internal memory.
;   But it still creates an two arrays (FLATARR (float) and OUTMASK (byte))
;   that is as large as all of the input flats.
;   Thus, if you are using ten 16 MB flat images,
;   that array will be 10*(16+4)=200 MB (in addition to a few other images).
;
;   The header contains the number of files used in each camera (NEXP),
;   and an identifier for each of those files (EXPID*).  Those identifiers
;   contain the camera name, MJD, flat exposure number, and arc exposure
;   number, all dash-separated.
;
; EXAMPLES:
;
; BUGS:
;  Since the traces move around this might night work too well.  This is 
;  not a trivial problem as the wavelength solution, and thus the fringe 
;  patterns are not continuous fiber-to-fiber. 
;
;  Also, this is very traumatic on the CPU. It is very possible, if the wrong
;  settings are used for the memory settings, will crash the computer.
;
; PROCEDURES CALLED:
;   bspline_iterfit()
;   bspline_valu()
;   calcscatimage()
;   djs_avsigclip()
;   djs_filepath()
;   extract_image
;   fileandpath()
;   genflatmask()
;   readfits()
;   repstr()
;   rmfile
;   sdssproc
;   sxaddpar
;   sxpar()
;   superflat()
;   trace320crude()
;   traceset2xy
;   writefits
;   xy2traceset
;
; REVISION HISTORY:
;   13-Oct-1999  Written by D. Schlegel, APO
;   April 2004 - Modified hs_spflatten2 to work with HS  data.  
;                R Cool U of A
;-
;------------------------------------------------------------------------------
pro hs_flatten, flatname,  allflats, pixflat, $
                sigrej=sigrej, maxiter=maxiter, $
                oldflat=oldflat, outfile=outfile, indir=indir, $
                outdir=outdir, tmpdir=tmpdir, $
                pixspace=pixspace, nord=nord, lower=lower, $
                upper=upper, mincounts=mincounts, $
                nodelete=nodelete, rerun=rerun
  
  for ccdnum=1,2 do begin
  
   if (N_elements(pixspace) EQ 0) then pixspace = 50
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 2
   if (N_elements(upper) EQ 0) then upper = 2
   if (NOT keyword_set(mincounts)) then mincounts = 0.

   nflat = N_elements(allflats)
   ngrow = 2

   if (NOT keyword_set(sigrej)) then begin
      if (nflat LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 flats
       else if (nflat EQ 3) then sigrej = 1.1 $
       else if (nflat EQ 4) then sigrej = 1.3 $
       else if (nflat EQ 5) then sigrej = 1.6 $
       else if (nflat EQ 6) then sigrej = 1.9 $
       else sigrej = 2.0
   endif
   if (NOT keyword_set(maxiter)) then maxiter = 2

   tmpname1 = 'tmp-pix-' + string(findgen(n_elements(allflats)),$
                                  format='(i2.2)')+'.fits.gz'
   tmpname2 = 'tmp-ymodel-' + string(findgen(n_elements(allflats)),$
                                     format='(i2.2)')+'.fits.gz'
   

   ;---------------------------------------------------------------------------
   ; First find the wavelength solution
   ;---------------------------------------------------------------------------

   ;------
   ; Read flat-field image that corresponds to the arc

   splog, 'Reading flat ', flatname
   hs_proc, flatname[0], ccdnum, flatimg, flativar, indir=indir, $
      header=flathdr, $
      rerun=rerun
   
       

   dims = size(flatimg, /dimens)
   nx = dims[0]
   ny = dims[1]

   ;------
   ; Create spatial tracing from flat-field image

   splog, 'Tracing 150 fibers in ',  flatname
   xsol = trace150crude(flatimg, flativar, yset=ycen, maxdev=0.15)

   splog, 'Fitting traces in ',  flatname
   xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
   traceset2xy, tset, ycen, xsol


   ; Compute wavelength calibration for arc lamp

   arccoeff = 5

   splog, 'Searching for wavelength solution'
   if file_test('calibration/' + rerun + '/wset.fits') eq 0L then begin
        hs_findwave, rerun=rerun, root=root

   endif
   wset = mrdfits('calibration/'+rerun+'/wset.fits', ccdnum)

   ;---------------------------------------------------------------------------
   ; Construct wavelength image
   ;---------------------------------------------------------------------------

   ;------
   ; Compute wavelength at every pixel on the CCD
   ; WAVEIMG will be in units of log10(lambda)

   waveimg = fltarr(nx,ny)

   traceset2xy, wset, xx, loglam
   loglam = alog10(loglam)
   
   xx = loglam*0.0
   for i = 0, n_elements(xx(0,*)) - 1 do begin
      xx(*,i) = indgen(n_elements(xx(*,0)))
   endfor
   xy2traceset, xx, loglam, wset
   delvarx, xx
   
   xy2traceset, transpose(xsol), transpose(loglam), tmpset, $
    func='legendre', ncoeff=5, xmin=0, xmax=nx-1, upper=2.0, lower=2.0
   xtmp = 0
   traceset2xy, tmpset, xtmp, waveimg

   ;---------------------------------------------------------------------------
   ; Read or create a mask image
   ;---------------------------------------------------------------------------
   
   
   maskimg = flatimg * 0.0

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      hs_proc, allflats[iflat], ccdnum, flatimg, flativar, $
         indir=indir, header=flathdr, $
          rerun=rerun
      
      if (iflat EQ 0) then begin
         hdr0 = flathdr
         sxaddpar, hdr0, 'NEXP', nflat, 'Number of exposures in this file', $
          before='EXPTIME'
      endif

      sxaddpar, hdr0, string(iflat+1,format='("EXPID",i2.2)'), $
       string( sxpar(flathdr,'CAMERAS'), $
        sxpar(flathdr,'EXPOSURE'), sxpar(flathdr,'EXPOSURE'), $
        format='(a2,"-",i8.8,"-",i8.8)'), $
        'ID string for exposure ' + strtrim(iflat+1,2), before='EXPTIME'

      ;----------------------
      ; Create spatial tracing from flat-field image

      xsol = trace150crude(flatimg, flativar, yset=ycen, maxdev=0.15)

      xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
      traceset2xy, tset, ycen, xsol

      ;----------------------
      ; Extract the flat-field vectors
      ; Apply the bad pixel mask, MASKIMG, so that we avoid those bad regions
      ; when computing the supersky vector.

      sigma = 2.0
      proftype = 3
      highrej = 15
      lowrej = 15
      npoly = 9 ; Try fitting 9 polynomial background terms to each row
      wfixed = [1,1,1] ; Just fit the first gaussian term

      extract_image, flatimg, flativar * (1-maskimg), $
       xsol, sigma, flux, fluxivar, $
       proftype=proftype, wfixed=wfixed, $
       highrej=highrej, lowrej=lowrej, npoly=npoly, /relative, $
       ansimage=ansimage

      ;----------------------
      ; Construct the "superflat" vector for this particular frame

      x2 = xsol ; CCD-X position
      asset = superflat(flux, fluxivar, wset, x2=x2, $
       fibermask=fibermask, minval=0.0, lower=lower, upper=upper, $
       nord=4, npoly=3)
      
      flux = 0
      fluxivar = 0
      
      delvarx, maskimg
      
      x2 = float(djs_laxisgen([nx,ny],iaxis=0)) ; CCD-X position
      fitimg = bspline_valu(waveimg, asset, x2=x2)
      ;fitimg = fitimg > 0.02    ; ???

      ;----------------------
      ; Construct the scattered light image

      dims = size(xsol,/dimens)
      nrow = dims[0]
      ntrace = dims[1]
      nterms = n_elements(wfixed)
      yfnorm = (2.0*findgen(nrow)-nrow)/nrow
      ;scatfit = fchebyshev(yfnorm,npoly) # ansimage[ntrace*nterms:*,*]

      ;----------------------
      ; divide by the superflat image.
      ; Then we're left with an image that should show a flat response
      ; if every pixel were the equally sensitive, except we'll still
      ; see the 320 fiber profiles.

      flatimg = (flatimg) / fitimg
      flativar = flativar * fitimg^2
      fitimg = 0


      ;----------------------
      ; Create the array of preliminary flats

      if (iflat EQ 0) then begin
         pixflatarr = fltarr(nx,ny,nflat)
         inmask = bytarr(nx,ny,nflat)
      endif

      ;----------------------
      ; Determine YMODEL image

      splog, 'Solving for YMODEL'



      ymodel = 0.0 * flatimg
      for i=0, nx-1 do begin
         ; Burles counter for column number...
         print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
         
         
         indx = where(flativar[i,*] GT 0, ngind)

         ; If fewer than 5 good points, then fit all points in this column.
         ; Also, change the values of FLATIVAR to prevent SLATEC_SPLINEFIT()
         ; from crashing.
         if (ngind LT 5) then begin
            indx = lindgen(ny)
            flativar[i,*] = 1.0
         endif


         ;------
         ; The following spline fit chooses breaks points evenly separated
         ; in row number.

         yaxis = findgen(ny)

         sset = bspline_iterfit(yaxis[indx], transpose(flatimg[i,indx]), $
          invvar=transpose(flativar[i,indx]), nord=nord, bkspace=pixspace, $
          upper=4, lower=4, maxiter=3, maxrej=ceil(0.05*n_elements(indx)))
         
         ymodel[i,*] = bspline_valu(yaxis, sset)

      endfor

      ;----------------------
      ; Set a pixel mask that is bad for YMODEL counts <= 1, or for pixels
      ; more than 2 pixels to the left (right) of the first (last) fiber.

      inmask[*,*,iflat] = ymodel GT 1 ; =1 for good
      for iy=0, ny-1 do begin
         x1 = long(xsol[iy,0]) - 3
         x2 = long(xsol[iy,ntrace-1]) + 3
         if (x1 GE 0) then inmask[0:x1,iy,iflat] = 0
         if (x2 LE nx-1) then inmask[x2:nx-1,iy,iflat] = 0
      endfor

      pixflatarr[*,*,iflat] = (flatimg * inmask[*,*,iflat]) $
       / (ymodel * inmask[*,*,iflat] + (inmask[*,*,iflat] EQ 0))

      ;----------------------
      ; Write FLATIMG and YMODEL to disk
      if (nflat GT 1) then begin
         writefits, tmpname1[iflat], flatimg
         writefits, tmpname2[iflat], ymodel
      endif
flatimg = 0
ymodel = 0

   endfor

   if (nflat EQ 1) then begin

      pixflat = temporary(pixflatarr)

   endif else begin
      ; Find deviant pixels in each pixflat
      inmask = temporary(1b - inmask) ; Change to =0 for good
      meanimg = djs_avsigclip(pixflatarr, 3, sigrej=sigrej, maxiter=maxiter, $
       inmask=inmask, outmask=outmask)
meanimg = 0
pixflatarr = 0


      outmask = inmask
      outmask = temporary(1-outmask)  ; Change to 0=bad, 1=good
      
      
      flatimgsum = 0
      ymodelsum = 0
      for iflat=0, nflat-1 do begin

         junk = where(outmask[*,*,iflat] EQ 0, ct)
         splog, 'Number of pixels masked in ', allflats[iflat], ' = ', ct

         flatimgsum = flatimgsum + outmask[*,*,iflat] * $
            readfits(tmpname1[iflat])
         ymodelsum = ymodelsum + outmask[*,*,iflat] * readfits(tmpname2[iflat])

         if (NOT keyword_set(nodelete)) then begin
            rmfile, tmpname1[iflat]
            rmfile, tmpname2[iflat]
         endif

      endfor

      qgood = ymodelsum GE mincounts
      pixflat = qgood * (flatimgsum > 0) / (ymodelsum + (1-qgood))

      
   endelse

   if (keyword_set(outfile)) then begin
      sxaddpar, hdr0, 'EXTEND', 'T'
      if ccdnum eq 1 then $
         writefits, djs_filepath(outfile, root_dir=outdir), pixflat, hdr0
      if ccdnum eq 2 then $
         writefits,  djs_filepath(outfile, root_dir=outdir),$
         pixflat, hdr0,/append
   endif
   
endfor


   return
end
;------------------------------------------------------------------------------
