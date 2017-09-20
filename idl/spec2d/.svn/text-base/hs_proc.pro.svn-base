;+
; NAME:
;   hs_proc
;
; PURPOSE:
;  Read the Hectospec data
;
; CALLING SEQUENCE:
;  hs_proc, infile, ccdnum, image, invvar, header=haeder, indir, $
;           pixelmask=, rerun=, root=, docosmic=, 
;
; INPUTS:
;     infile  - file containing image
;     ccdnum  - ccdnum of the side of the image you would like
;     
; OPTIONAL KEYWORDS:
;   rerun = reduction rerun
;   root  = reduction MJD
;   docosmic - do CR rejection
;
; OUTPUTS:
;
;   image - output image
;   invvar - output inverse variance
;
; OPTIONAL OUTPUTS:
;
;   header - image header
;   pixelmask - badpixel mask
;
; COMMENTS:
;   This code will read the bad pixel mask and bias correct the data
;   automatically (as well as trim and overscan)
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA

;------------------------------------------------------------------------------;

;------------------------------------------------------------------------------
;  Create the bad column mask (1 for a masked pixel) with image size nc,nr
;  If the operation doesn't work just return 0 for no masked pixels

function make_badpix_mask, bcfile, camrow, camcol, nc=nc, nr=nr, $
                           silent=silent,rerun=rerun, root=root
  
  if NOT keyword_set(nc) then nc=4096L
  if NOT keyword_set(nr) then nr=4096L
  
  readcol, bcfile, xbegin, xend, ystart, yend,  format = 'L,L', $
     comment='#', /silent
  bcmask = long(fltarr(nc,nr))*0.0
  if n_elements(xbegin) gt 0 then  begin
     for i = 0, n_elements(xbegin) -1 do begin
        bcmask[xbegin(i):xend(i), ystart(i):yend(i)] = 1
        
     endfor
  endif
  return, bcmask
end


FUNCTION gfit, X, P
  RETURN, P(0) + GAUSS1(X, P(1:3))
END

PRO hs_proc, infile, ccdnum, image, invvar, header=header, $
             indir=indir, outfile=outfile, pixelmask=pixelmask, $
             rerun=rerun, root=root, applyflat=applyflat, docosmic=docosmic, $
             nobias=nobias, chelle=chelle
  
  
  if not keyword_set(indir) then indir = cwd()
  if not keyword_set(outdir) then outdir = cwd()
  if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  
  if ccdnum ne 1 and ccdnum ne 2 then begin
     splog, 'CCD Number 1 and 2 are the allowed values for ccdnum'
     return
  endif
  if ccdnum eq 1 then amp = [1,2] 
  if ccdnum eq 2 then amp = [3,4]
  
  
  ;;Construct IMAGE
  ;;read in the first image to get the right sizes
  
  
  for iamp = 0, 1 do begin
     ;;Read the first image
     
     splog, 'Reading in file ' + infile 
     splog, 'Amplifier ' + $
        string(amp[iamp], format='(i1.1)') + ', please be patient'
     
     thisimage = readfits(infile, exten_no=amp[iamp], header, /silent)
     
     ;;Find the bad columns
     bcmask = thisimage * 0.0
     bcfile = getenv("HSRED_DIR") + '/etc/badpix-' + $
        string(ccdnum, format='(i1.1)')
     IF keyword_set(chelle) THEN  bcfile = getenv("HSRED_DIR") + '/etc/hc_badpix-' + $        
        string(ccdnum, format='(i1.1)')     
     if ccdnum eq 1 and iamp eq 0 then bcfile = bcfile + '.0.dat'
     if ccdnum eq 1 and iamp eq 1 then bcfile = bcfile + '.1.dat'
     if ccdnum eq 2 and iamp eq 0 then bcfile = bcfile + '.0.dat'
     if ccdnum eq 2 and iamp eq 1 then bcfile = bcfile + '.1.dat'  
     
     dimen = size(thisimage, /dimen)
     nc = dimen[0]
     nr = dimen[1]
 
     bcmask = make_badpix_mask(bcfile, camrow, camcol, $
                               nc=nc, nr=nr, silent=silent)
     pixelmask =bcmask
     hs_overscan, thisimage, header
     hs_trim, thisimage, header=header, immask=bcmask
     
     if (iamp EQ 0) then begin
        imsize = size(thisimage, /dimension)
        ncols = imsize[0]
        nrows = imsize[1]
        output = dblarr(ncols, nrows, 2)
        outputivar = output
     endif
     
     if not keyword_set(root) then root = (long(sxpar(header, 'mjd')))
     
     ;;Apply the bias correction
     biasfile = 'calibration/'+rerun+'/bias.fits'
     
     ;;Check to see if the biasfile exists in this directory
          
     pixfile = lookforgzip(djs_filepath(biasfile, root_dir=indir))
     if file_test(pixfile) eq 0 then begin
        splog, 'Cannot find the bias'
        stop
     endif
     
     
     pixbias = readfits(pixfile(0), $
                        exten_no=amp[iamp], biashead, /silent)
     
     ;;Now interpolate over bad pixels
     
     ;;As of june 23, 04 I have changed this to loop through the rows 
     ;;and create an interpolation on the pixels based on the pixels 
     ;;that are not flagged in the bpm
     
     ;;Only do the interpolation if there are any badpixels
     bad = where(bcmask eq 1, badcount)
     
     if badcount gt 0 then begin
        ;;Loop the rows
        for i = 0, n_elements(bcmask(0,*))-1 do begin
           ;;choose the GOOD pixels
           goodpix = where(bcmask(*,i) eq 0)
           badpix = where(bcmask(*,i) eq 1, badcount)
           
           if badcount gt 0 then begin
              row = indgen(n_elements(bcmask(*,i)))
              temprow1 = interpol(thisimage(goodpix,i),goodpix, badpix)
              m1 = max(thisimage(badpix(0)-20:badpix(0)+20, i))
              
              temprow = interpol(thisimage(goodpix,i), goodpix, row)
              thisimage(*,i) = temprow 
              
           endif
        endfor   
     endif
     
     if keyword_set(nobias) then pixbias=pixbias*0.0
     
     
     thisimage = thisimage - pixbias
     delvarx, pixbias
                                ;        
     ;;Add pixbias to header
     sxaddpar, header, 'PIXBIAS', pixfile(0)
     
     ;;Now correct for the gain
     gain = float(sxpar(header, 'gain'))
     
     IF gain GT 1.5 THEN BEGIN 
        
        ;;This is the hardwired gain
     
        ;;From Nelson Cawldwell, the gain and readnoise in the headers
        ;;are incorrect.  I am supplying the absolute value here
        gain = [1.0, 0.99, 0.949, 0.937]
        gain = gain(amp(iamp)-1)
        sxaddpar, header, 'GAINCOR', 'Applied gain of ' + strn(gain)
     ENDIF
     
     
     ;;Now the images are in electrons, not DN
     thisimage = thisimage * gain
     
     
     ;; This is the second amplifer on the chip, then reverse the image
     if iamp eq 1 then begin
        columns = lindgen(ncols)
        thisimage[columns, *] = reverse(thisimage[columns, *])
        bcmask[columns,*]=reverse(bcmask[columns,*])
     endif
     ;;Construct the inverse variance map for this image
     
     thisinvvar = thisimage * 0.0
     ;;The flux term 
     expr1 = abs(thisimage)
     
     ;;The readnoise/gain term
     ;;rdnoise = double(sxpar(header, 'rdnoise'))
     ;;As above, the header readnoise is not correct
     rdnoise = 2.8
     expr3 = (rdnoise*gain)^2.0
     thisivar = 1.0/(expr1 + expr3) *(1-bcmask)
     delvarx, expr1
     delvarx, expr3
     output(*,*,iamp) = thisimage
     outputivar(*,*,iamp) = thisivar
     pixelmask1 = output*0.0
     
     pixelmask1(*,*,iamp)= bcmask
     delvarx,thisimage
     delvarx,thisivar
  endfor
  
  ;;Now mesh them together
  
  image = dblarr(2*ncols, nrows)
  invvar = image
  pixelmask2 = image
  
  image(0:(ncols-1), * ) = output(*,*,0)
  invvar(0:(ncols-1), * ) = outputivar(*,*,0)
  image(ncols:(2*ncols-1), * ) = output(*,*,1) 
  invvar(ncols:(2*ncols-1), * ) = outputivar(*,*,1)
  pixelmask2(ncols:(2*ncols-1), * ) = pixelmask1(*,*,1)
  pixelmask2(0:(ncols-1), * ) = pixelmask1(*,*,0)
  pixelmask=pixelmask2
  
  ;;Now do the cosmic ray rejection if so requested
  if keyword_set(docosmic) then begin
     
     hs_cosmic, image, newflux, rerun=rerun
     image = newflux
     delvarx, newflux
     ;;Check for infinite pixels from this
     mask = finite(image)
     invvar = invvar * mask
     k = where(mask eq 0)
     if k(0) ne -1 then image(k) = 0.0
     sxaddpar, objhdr, 'COSMIC', 'CR-rejected using LACosmic'
  endif
  if keyword_set(applyflat) then begin
     if file_test('calibration/' + rerun + '/pixflat.fits') NE 1 then begin
        splog, 'The pixel-to-pixel flat has not been generated!!!!'
     endif else begin
        
        
        pixflat = mrdfits('calibration/' + rerun + '/pixflat.fits', ccdnum-1)
        mask = pixflat LE 0.5
        
        pixflat = pixflat*(1-mask) + mask
        mask = abs(pixflat) gt 1
        pixflat= pixflat*(1-mask) + mask
        
        image = image / pixflat
        
        invvar = invvar*pixflat^2
        delvarx, mask
        delvarx, pixflat
     endelse
     
  endif
  
         
     
  sxaddpar, header, 'naxis1', imsize[0]
  sxaddpar, header, 'naxis2', imsize[1]
  sxaddpar, header, 'extname', strn(ccdnum)
  sxaddpar, header, 'namps', 2
  sxaddpar, header, 'datasec', '[1:' + strn(imsize[0]) + $
     ','+ '1:' + strn(imsize[1]) +']'
  sxaddpar, header, 'ccdsec', '[1:' + strn(imsize[0]) + $
     ','+ '1:' + strn(imsize[1]) +']'
  sxaddpar, header, 'detsize',   '[1:' + strn(imsize[0]) + $
     ','+ '1:' + strn(2*imsize[1]) +']'
  sxdelpar, header, 'biassec'
  sxdelpar, header, 'trimsec'
  sxdelpar, header, 'origsec'
  sxaddpar, header, 'EXTNAME', strn(ccdnum)
  sxdelpar, header, 'BSCALE'
  sxdelpar, header, 'BZERO'
  sxaddpar, header, 'imageid', strn(ccdnum)
  if ccdnum eq 2 then begin
     sxaddpar, header, 'detsec', '[1:' +strn(imsize[0])$
        + ',' + strn(imsize[1]+1) $
        + ':' + strn(2*imsize[1]) +']'
  endif else begin
     sxaddpar, header, 'detsec', '[1:' +strn(imsize[0]) + ',' + strn(1) $
        + ':' + strn(imsize[1])     +']'   
  endelse
  
  ;;reverse the image
  thisimage1 = image * 0.0
  nspectral = n_elements(image(0,*))-1
  for i = 0,n_elements(image(0,*))-1 do begin
     thisimage1[*,i] = image[*,nspectral-i]
  endfor
  image=thisimage1
  
  thisimage1 = image * 0.0
  nspectral = n_elements(image(0,*))-1
  for i = 0,n_elements(image(0,*))-1 do begin
     thisimage1[*,i] = invvar[*,nspectral-i]
  endfor
  invvar=thisimage1
           
  return
END
  
