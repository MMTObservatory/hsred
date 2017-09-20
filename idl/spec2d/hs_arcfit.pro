;+
; NAME:
;	IARCFIT
;
; PURPOSE:
;	Generate a two dimensional map between pixel position and
;	wavelength from an arc calibration lamp.
;
; CALLING SEQUENCE:
;       iarcfit, arcfile, wset, xcen, ycen, [datapath=, lamp=, wmapname=], $'
;          [row=, nmed_arc=, xcoeff=, dxcoeff=, func_arc=,
;          norder_arc=, mintol=], $
;          arclambda=arclambda, xdiff=xdiff, gauss=gauss, $
;          write=write, doplot=doplot'
;
; INPUTS:
;	arcfile  - arc lamp FITS file from which to derive the
;                  wavelength map
;
; OPTIONAL INPUTS:
;	datapath   - path to the data and the WRITE directory
;	lamp       - type of comparison lamp (e.g., 'HeAr')
;       wmapname   - name of the output wavelength map
;       row        - row to use for initial wavelength solution
;                    (default to the middle row) 
;	nmed_arc   - number of rows to median filter around ROW for
;                    initial wavelength solution (default 5)
;	xcoeff     - two-element array giving the [starting wavelength
;                    (Angstrom), dispersion (Angstrom/pixel)] initial
;                    guesses (default [3600.0,2.75])
;	dxcoeff    - deviations from XCOEFF (default [50.0,0.1])
;	func_arc   - name of fitting function (default 'legendre')
;	norder_arc - order of the fit of column number versus
;                    wavelength (default 4) 
;	mintol     - minimum tolerance for acceptable arc lines
;                    (default 0.2 Angstrom) 
;	
; KEYWORD PARAMETERS:
;	gauss    - use gaussian profile fitting for final centroid fit
;                  (default to flux-weighted centering)
;	write    - write the wavelength map and the QA plots to disk
;       doplot   - generate QA plots on the screen (will be generated
;                  anyway if WRITE=1)
;
; OUTPUTS:
;	wset    - 2D wavelength map (pixel -> lambda)
;	xcen    - column pixel positions of arc lines [nrows,nlines] 
;	ycen    - row positions of arc lines [nrows,nlines]
;
; OPTIONAL OUTPUTS:
;	arclambda - wavelengths of good lamp lines (Angstrom)
;	xdiff     - residuals of fitted lamp lines to fitted positions
;                   (pixels)
;
; COMMENTS:
;       A 2nd order Legendre function is fitted to the row-dependent
;       arc line position, while an NORDER_ARC Legendre polynomial is
;       fitted in the wavelength dimension.
;
;	Also see FITARCIMAGE in Dave Schlegel's IDLSPEC2D package.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	RD2DSPEC(), READCOL, IARCFIT_GUESS(), DJS_MEDIAN(),
;	TRACE_CRUDE(), FIND_NPEAKS(),  TRACESET2XY(), TRACESET2PIX(),
;	XY2TRACESET, TRACE_FWEIGHT(), TRACE_GWEIGHT(), ICLEANUP,
;	DFPSPLOT, DFPSCLOSE, CWD(), GET_ELEMENT, CLEANPLOT,
;	LINEID_PLOT, DJS_ITERSTAT, ERRPLOT, LEGEND, CMSAVE
;
; DATA FILES:
;	${ISPEC_DIR}/etc/lamphear.dat
;
; OLD MODIFICATION HISTORY:
;	J. Moustakas, 2001 August, U of A
;       jm02nov15uofa - added a quality assurance plot
;       jm03jan12uofa - added support for a HeNeAr lamp
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 April 9, U of A, checked out ISPEC v1.0.0
;       jm03dec7uofa - added DEBUG keyword
;-

pro hs_arcfit, arcfile, wset, xcen, ycen, datapath=datapath, lamp=lamp, $
               wmapname=wmapname, row=row, nmed_arc=nmed_arc, $
               xcoeff=xcoeff, dxcoeff=dxcoeff, func_arc=func_arc, $
               norder_arc=norder_arc, mintol=mintol, arclambda=arclambda, $
               xdiff=xdiff, gauss=gauss, write=write, debug=debug, $
               doplot=doplot

    narc = n_elements(arcfile)
    if narc eq 0L then begin
        print, 'Syntax - iarcfit, arcfile, wset, xcen, ycen, '+ $
          '[datapath=, lamp=, wmapname=], $'
        print, '   row=, nmed_arc=, xcoeff=, dxcoeff=, func_arc=,' + $
          ' norder_arc=, mintol=], $'
        print, '   arclambda=arclambda, xdiff=xdiff, /gauss, $'
        print, '   /write, /debug, /doplot'
        return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(mintol) eq 0L then mintol = 0.5       ; [Angstrom]
    if n_elements(lamp) eq 0L then lamp = 'HeAr'
    if n_elements(xcoeff) eq 0L then xcoeff = [3600.0,2.75]
    if n_elements(dxcoeff) eq 0L then dxcoeff = [50.0,0.1]
    if n_elements(nmed_arc) eq 0L then nmed_arc = 5L
    if n_elements(norder_arc) eq 0L then norder_arc = 4L
    if n_elements(func_arc) eq 0L then $
      func_arc = 'legendre' else func_arc = strlowcase(func_arc)

    if narc gt 1L then begin
        
        if n_elements(wmapname) eq 0L then wmapname = replicate('',narc) else $
          if n_elements(wmapname) ne narc then $
          message, 'Incompatible dimensions in ARCFILE and WMAPNAME.'
        
        for k = 0L, narc-1L do begin
            
            iarcfit, arcfile[k], datapath=datapath, lamp=lamp,$
              wmapname=wmapname[k], $
              row=row, nmed_arc=nmed_arc, xcoeff=xcoeff, $
              dxcoeff=dxcoeff, func_arc=func_arc, $
              norder_arc=norder_arc, mintol=mintol, gauss=gauss, $
              write=write, debug=debug, $
              doplot=doplot

       endfor
       return

    endif

    cube = rd2dspec(arcfile,datapath=datapath)

    if strmatch(strcompress(cube.type,/remove),'*comp*') ne 1B then $
      splog, 'WARNING: '+arcfile+' may not be a comparison lamp.'

    arc = cube.image
    arcinvvar = 1.0/cube.sigmamap^2.0D ; inverse variance
    archead = *cube.header
    ncols = cube.naxis1
    nrows = cube.naxis2

    if n_elements(row) eq 0L then row = nrows/2L else row = long(row)

; read in the lamplist for wavelength calibration into a structure.
; only HEAR and HENEAR are supported

    case strupcase(lamp) of
       'HEAR': lampfile = 'lamphear.dat'
       'HENEAR': lampfile = 'lamphenear.dat'
       else: begin
          splog, 'Unsupported comparison arc lamp.'
          return
       end
    endcase

    lamppath = getenv('HSRED_DIR')+'/etc/'
    print, lamppath+lampfile
    if file_test(lamppath+lampfile,/regular) eq 0L then begin
       splog, 'Lampfile '+lampfile+' not found.'
       return
    endif
    
    splog, 'Reading the lamp file '+lampfile+'.'
    readcol, lamppath+lampfile, lampwave, lampinten, lampquality, $
      lampelement, format='F,F,A,A', /silent
    nlamp = n_elements(lampwave)

    lamp = {element: '', lambda: 0.0D, loglamp: 0.0D, $
            intensity: 0.0D, good: 0B}
    lamp = replicate(lamp,nlamp)
    lamp.element = lampelement
    lamp.lambda = lampwave
    lamp.loglamp = alog10(lampwave)
    lamp.intensity = lampinten
    lamp.good = strupcase(lampquality) eq 'GOOD' or $
      ((lampinten gt 100.0) and (lampquality eq 'BLEND'))

    glamp = where(lamp.good eq 1B,nlamp)
    if nlamp eq 0L then begin
       splog, 'No good arc lines in the lamp list.'
       return
    endif else lamp = lamp[glamp]

; median filter NMED_ARC rows centered on ROW; this median spectrum
; determines the initial wavelength solution
    
    spec = djs_median(arc[*,row-nmed_arc/2L:row+nmed_arc/2L],2) 

; initialized the input wavelength solution based on the dispersion
; and starting wavelength.  the linear transformations below are being
; computed to be compatible with the way the tracesets are stored

    xrange = ncols-1.0 & xmid = 0.5*xrange

    acoeff = fltarr(4)
    acoeff[0] = xcoeff[0] + xmid * xcoeff[1]
    acoeff[1] = xcoeff[1] * xrange / 2.0
    acoeff[2:3] = [5.0,0.0] ; [5.0,-0.8]

    dcoeff = fltarr(4)
    dcoeff[0] = dxcoeff[0] + xmid * dxcoeff[1]
    dcoeff[1] = dxcoeff[1] * xrange / 2.0
    dcoeff[2:3] = [20.0,5.0]
    
    splog, 'Searching for the initial wavelength solution.'
;   nsteps = [1, 80, 40, 1]
    aset = iarcfit_guess(spec,lamp.lambda,lamp.intensity,acoeff=acoeff,$
       dcoeff=dcoeff,nsteps=nsteps,bestcor=bestcor,model=model)
    initcoeff = [aset.coeff[0]-xmid*xcoeff[1],aset.coeff[1]*2.0 / $
                 xrange,aset.coeff[2:3]]
    splog, 'Initial wavelength fit = ', initcoeff, format='(a,99f12.5)'

;   traceset2xy, aset, xdum, rowlambda
;   splog, 'Wavelength range in ROW ', minmax(rowlambda)
;   splog, 'Dispersion in ROW ', (max(rowlambda)-min(rowlambda))/ncols

; trim the lamplist to only those wavelengths within the wavelength
; range.  predict the pixel positions of the lines in our lamp list
; using the initial wavelength solution and attempt to trace those
; lines on the arc lamp image

    xstart = traceset2pix(aset,lamp.lambda,/silent)
    qtrim = (xstart gt 1L) and (xstart lt ncols-2L) and lamp.good
    itrim = where(qtrim,nlamp)
    if nlamp eq 0L then message, 'No arc lines in wavelength range.'
    xstart = xstart[itrim]
    lamp = lamp[itrim]
    
; trace arc lines as a function of row number on the 2D arc image.
; allow for a shift of up to 2 pixels in the initial centers, but only
; 0.5 pixels while tracing, starting at the middle row

    splog, 'Tracing ', nlamp, ' arc lines.'
    xcen = trace_crude(arc,arcinvvar,xstart=xstart,ystart=row,radius=radius,$
      nave=1.0,nmed=1.0,maxerr=maxerr,maxshift0=1.0D,maxshifte=0.2D,$
      yset=ycen)
;   print, total(xcen,1)/nrows

; trace from the bottom of the CCD using the previous trace as the
; initial line centers

    xcen = trace_crude(arc,arcinvvar,xstart=xcen[0,*],ystart=0L,radius=radius,$
      nave=1.0,nmed=1.0,maxerr=maxerr,maxshifte=0.2D,maxshift0=1.0D,xerr=xerr,$
      yset=ycen)
;   print, total(xcen,1)/nrows

    if keyword_set(debug) then begin
    
       djs_plot, findgen(ncols), arc[*,row], xsty=3, ysty=3, ps=10, $
         charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, $
         ytitle='Counts', $
         xtitle='Column Number', xrange=[100,200]
       djs_oplot, findgen(ncols), arc[*,20], color='yellow'
       djs_oplot, findgen(ncols), arc[*,nrows-20], color='red', ps=10
       for i = 0L, nlamp-1L do djs_oplot, [xcen[row,i],xcen[row,i]], $
         !y.crange, line=1, thick=0.5, color='cyan'
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
       
    endif

; fix bad traces (optional)

;   xcen = trace_fix(xcen,minsep=minsep,ngrow=ngrow,ycen=ycen,xerr=xerr)

; for each line, fit a legendre function to the trace to make the
; center a smoothly varying function and to interpolate over bad rows

    ncoeff_init_trace = 7L ; 3L
    xy2traceset, ycen, xcen, crudeset, ncoeff=ncoeff_init_trace, $
      yfit=xcrudefit, /silent

    if keyword_set(debug) then begin
    
       for i = 0L, nlamp-1L do begin
          djs_plot, findgen(nrows), xcen[*,i], ps=4, xsty=3, ysty=3, $
            charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, $
            ytitle='Column Position (pixel)', xtitle='Row Number'
          djs_oplot, findgen(nrows), xcrudefit[*,i], line=0,$
            thick=2.0, color='red'
          legend, 'Line '+string(i,format='(I3.3)')+'/'+$
            string(nlamp,format='(I3)'), $
            /left, /top, box=0, charsize=1.5, charthick=2.0
          cc = get_kbrd(1)
      endfor

  endif
    
; iterate the flux-weighted centers.  in the last iteration, use the
; formal errors in the arc image
  
    xnew = trace_fweight(arc,xcrudefit,ycen)
    xnew = trace_fweight(arc,xnew,ycen)

    if (keyword_set(gauss)) then $
      xcen = trace_gweight(arc,xnew,ycen,sigma=1.0,$
                           invvar=arcinvvar,xerr=xerr) else $
      xcen = trace_fweight(arc,xnew,ycen,radius=2.0,invvar=arcinvvar,xerr=xerr)

; to obtain the final two-dimensional wavelength map fit column
; position (pixels) as a function of line list wavelength (Angstrom).
; use the formal errors in the centroids and do an NORDER_ARC order
; fit to the function between wavelength and column number at each row 

    if (nlamp eq 1L) then $
      ytmp = transpose(lamp.lambda*(dblarr(nrows)+1)) else $
      ytmp = lamp.lambda # (dblarr(nrows)+1)

    xy2traceset, transpose(double(xcen)), ytmp, wset, func=func_arc, $
      ncoeff=(norder_arc+1), maxiter=nlamp, maxrej=1.0, /sticky, xmin=0, $
      xmax=ncols-1, yfit=yfit, /silent

;   plot, xcen[63,*], ytmp[*,63]-yfit[*,63], ps=4

; check for bad or mismatched lines and reject them

    resid2d = ytmp - yfit
    resid1d = total(resid2d,2)/nrows ; mean residuals
    
    repeat begin
      
       dev = max(abs(resid1d),indx)
;      print, 'Bad line - off by '+strn(dev,length=5)+' Angstrom.'

       lmask = bytarr(nlamp)+1B
       lmask[indx] = 0B
       w = where(lmask,nlamp)
       lamp = lamp[w]
       xcen = xcen[*,w]
       
       if (nlamp eq 1L) then $
         ytmp = transpose(lamp.lambda*(dblarr(nrows)+1)) else $
         ytmp = lamp.lambda # (dblarr(nrows)+1)
       
       xy2traceset, transpose(double(xcen)), ytmp, wset, $
         func=func_arc, ncoeff=(norder_arc+1), $
         maxiter=nlamp, maxrej=1.0, /sticky, xmin=0, xmax=ncols-1, $
         yfit=yfit, /silent
       
       resid2d = ytmp - yfit
       resid1d = total(resid2d,2)/nrows ; mean residuals
       
    endrep until total(abs(resid1d) gt mintol) eq 0

    splog, 'Final arcfit complete based on '+string(nlamp,format='(I0)')+$
      ' lines.'
; quality assurance plots

    ifitname = 'arc_'+strmid(arcfile,0,strpos(arcfile,'.fits'))
    

    traceset2xy, wset, xmap, ymap
    xwave = ymap[*,row]
    yflux = arc[*,row]/max(arc[*,row])

    peaks = find_npeaks(yflux,xwave,nfind=10,minsep=5,width=10,$
                        ypeak=ypeak,npeak=npeak)
    get_element, lamp.lambda, peaks, xx
    srt = sort(xx)
    xx = xx[srt] & peaks = peaks[srt]

    crop = uniq(xx)
    xx = xx[crop] & peaks = peaks[crop]
    npeak = n_elements(peaks)
    
; crop PEAKS to actually matched lines (using MINTOL)

    indx = replicate(-1L,npeak)
    for ipeak = 0L, npeak-1L do if (abs(peaks[ipeak]-lamp[xx[ipeak]].lambda)$
                                    lt 2*mintol) then indx[ipeak] = ipeak

    good = where(indx ne -1L,ngood)
    if ngood eq 0L then splog, 'No arc lines have been matched for labeling.' $
    else begin

       xx = xx[indx[good]]
;      print_struct, lamp[xx]
    
       if keyword_set(doplot) or keyword_set(write) then begin

          if keyword_set(doplot) then iter = 1L else iter = 0L

          for i = 0L, iter do begin
             
             if i eq 0L then begin
                !p.font = 0
                dfpsplot, datapath+'qaplot_'+ifitname+'_lineid.ps', $
                  /square, /color, /isolatin1
             endif else window, 2, xs=500, ys=500

             lineid_plot, xwave, yflux, lamp[xx].lambda, lamp[xx].element, $
               string(lamp[xx].lambda,format='(F6.1)'), $
               charsize=1.8, charthick=2.0, xsty=3, ysty=3, ps=10,$
               yrange=[-0.05,1.3], /extend, $
               xtitle='Wavelength ('+angstrom()+')',$
               ytitle='Normalized Flux', $
               lcharsize=1.8, lcharthick=1.9,$
               position=[0.15,0.12,0.95,0.65], title=ifitname

             if i eq 0L then begin
                dfpsclose
                !p.font = -1
             endif
             
          endfor
          cleanplot, /silent

       endif

    endelse

; plot 2
    
    arclambda = lamp.lambda     ; wavelengths of good arc lines
    xarc = transpose(traceset2pix(wset,arclambda,/silent)) ; arc pixel 
                                ;positions [nrows,nlamp]
    xdiff = xcen-xarc           ; residuals between the traced and 
                                ;the predicted positions

; compute offset and stddev for each line center in pixels and in Angstroms

    meanpix = fltarr(nlamp)
    sigmapix = fltarr(nlamp)

    meanlam = fltarr(nlamp)
    sigmalam = fltarr(nlamp)
    
    for k = 0L, nlamp-1L do begin

; pixels
       
       djs_iterstat, xdiff[*,k], mean=mn, sigma=sg
       meanpix[k] = mn
       sigmapix[k] = sg

; wavelength

       djs_iterstat, yfit[k,*]-ytmp[k,*], mean=mn, sigma=sg
       meanlam[k] = mn
       sigmalam[k] = sg

    endfor

; determine the starting and ending wavelengths.  set the plot range

    traceset2xy, wset, transpose(replicate(wset.xmin,nrows)), minwave
    traceset2xy, wset, transpose(replicate(wset.xmax,nrows)), maxwave
   
    waverange = [min([minwave,maxwave]),max([minwave,maxwave])]
    splog, 'Wavelength range (Angstrom): ', $
      strn(waverange[0]), ', ', strn(waverange[1])
    splog, 'Mean disperion (Angstrom/pixel): ', $
      strn((waverange[1]-waverange[0])/ncols)
    
; plot the residuals in pixel position as a function of wavelength

    if keyword_set(doplot) or keyword_set(write) then begin

       if keyword_set(doplot) then iter = 1L else iter = 0L

       for i = 0L, iter do begin
      
          if i eq 0L then begin
             !p.font = 0
             dfpsplot, datapath+'qaplot_'+ifitname+'.ps', $
               /square, /color, /isolatin1
         endif else window, 0, xs=500, ys=500
         
         yminmax = max(abs([max(meanlam+sigmalam),min(meanlam-sigmalam)]))
         
         plot, arclambda, meanlam, xrange=waverange, $
           yrange=yminmax*[-1,1], ps=6, $
           xsty=3, ysty=3, xtit='Wavelength ('+angstrom()+')', $
           ytit='Mean Residuals ('+angstrom()+')', $
           title=ifitname, charsize=1.8, charthick=2.0, xthick=2.0,$
           ythick=2.0, thick=2.0
         oplot, waverange, [0,0]
         errplot, lamp.lambda, meanlam-sigmalam, meanlam+sigmalam
         legend, [textoidl('\Delta')+' = '+ $
                  string(mean(meanlam),format='(F8.5)'),$
                  textoidl('\sigma_{\Delta}')+' = '+ $
                  string(stddev(meanlam),format='(F8.5)')], $
           /right, /top, box=0, charsize=1.7, charthick=2.0
         
         if i eq 0L then begin
             dfpsclose
             !p.font = -1
         endif
         
     endfor
     
 endif
 
 
 if keyword_set(write) then begin
     
     if n_elements(wmapname) eq 0L then wmapname = 'wmap_'+ifitname+'.idlsave'
     if strcompress(wmapname,/remove) eq '' then $
       wmapname = 'wmap_'+ifitname+'.idlsave'
     
     splog, 'Writing '+datapath+wmapname+'.'
     cmsave, filename=datapath+wmapname, wset, arclambda, xcen, archead, $
       names=['wset','arclambda','xarc','header']
     
     openw, lun, datapath+'qalog_'+ifitname+'.log', /get_lun
     printf, lun, '# Two dimensional wavelength map derived from ', arcfile
     printf, lun, '# '+strn(ncols), ' columns and ', strn(nrows), ' rows'
     printf, lun, '# Wavelength range (Angstrom): ', $
       strn(waverange[0]), ', ', strn(waverange[1])
     printf, lun, '# Mean dispersion (Angstrom/pixel): ', $
       strn(mean(wset.coeff[1,*]*2/xrange))
     printf, lun, '# '+strn(nlamp), ' arc lines used in the final solution'
     printf, lun, '# Minimum wavelength tolerance for a good arc line: ',$
       strn(mintol), ' Angstrom'
     printf, lun, '# Mean residuals and RMS of the fit (Angstrom): ', $
       strn(mean(meanlam)), ', ', strn(stddev(meanlam))
     printf, lun, '# Fitting function and order: ', strn(func_arc), ', ', $
       strn(norder_arc)
     printf, lun, '#  '
     printf, lun, '# Nominal wavelength, mean measured wavelength,' + $
       ' pixel position:'
     printf, lun, '#  '
     
     for i = 0L, nlamp-1L do printf, lun, arclambda[i], $
       mean(yfit[i,*]), mean(xarc[*,i])
     
     free_lun, lun
       
    endif
    
    
    icleanup, cube
    
return
end
