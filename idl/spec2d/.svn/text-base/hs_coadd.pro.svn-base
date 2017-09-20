PRO hs_coadd, filelist, outname, outdir=outdir, rerun=rerun, flux=flux, $
              notweak=notweak
  
  ;;This program will take two HS files (presumably spHect files) and coadd
  ;;them  this is really only needed if you have plates that were taken on 
  ;;different nights and want to coadd them. You would think that you could
  ;;simply flux them all together, but they may be affected differently the 
  ;;any problems with the adc
  
  if not keyword_set(outdir) then outdir=cwd()
  nfib = 300
  
  ;;Read in the wavelength for the first file to get some dimensions down
  lam1= mrdfits(filelist(0), 0)
  
  lamout = dblarr(n_elements(filelist), n_elements(lam1[*,0]),$
                  n_elements(lam1[0,*]))
  fluxout = lamout
  ivarout = lamout
  pixmaskout = lamout
  plugmapout = mrdfits(filelist(0), 5)
  
    
  for i = 0, n_elements(filelist) -1 do begin
     
     lamout[i,*,*] = mrdfits(filelist(i), 0)
     fluxout[i,*,*] = mrdfits(filelist(i), 1)
     ivarout[i,*,*] = mrdfits(filelist(i), 2)
     pixmaskout[i,*,*] = mrdfits(filelist(i), 3)
     
     if keyword_set(flux) then begin
        junk = strsplit(filelist(i), '/', /extract)
        if n_elements(junk) gt 1 then begin
           dir = ''
           for j1 = 0, n_elements(junk) -2 do begin
              dir = dir + '/' + junk(j1)
           endfor
        endif else dir= cwd()
        junk1 = strsplit(junk(n_elements(junk)-1), ' ' , length=length)
        root = strsplit(junk(n_elements(junk)-1), '.-', /extract)
        root = root(n_elements(root)-2)
        fcal = findfile(dir + '/hsFluxcalib*' + root + '*fits')
        fcal = fcal(0)
        if file_test(fcal) eq 0 then stop
        
        sfactor = mrdfits(fcal, 1)
        factor = bspline_valu(alog10(lamout[i,*,*]), sfactor)
        
        f= reform(fluxout[i,*,*])
        iv = reform(ivarout[i,*,*])
        divideflat, f, reform(factor), invvar=iv, $
           minval = 0.01*median(factor[*,0])
        fluxout[i,*,*] = f
        ivarout[i,*,*] = iv
        l = lamout[i,*,*]
        plot, l, f
     endif
     
  endfor
  
  
  sn = dblarr(n_elements(filelist))
  for im = 0, n_elements(filelist) -1 do begin
     
     hs_snr,  transpose(transpose(lamout[im,*,*])), $
        transpose(transpose(fluxout[im,*,*])), plugmapout, $
        sn2=sn21, /noplot
     
     sn(im) =sn21
     
  endfor
  
  expid = filelist
  
  sn_exp = sn
  maxval = max(sn_exp/median(sn_exp), iframe_best)
  best_exposure = expid[iframe_best]
  splog, 'Best Exposure is: ' + best_exposure
  
  if not keyword_set(notweak) then begin
     ;; Compute the exposure-to-exposure corrections
     
     lam = lamout
     flux = fluxout
     ivar = ivarout
     mask = pixmaskout
     
     ndim = size(lam, /n_dimensions)
     lamorig = lam
     airtovac, lamorig
     
     if ndim eq 3 then begin
        
        size1 = size(lam)
        nim = size1(1)
        npix = size1(2)
        nfib = size1(3)
        newlam = dblarr(npix, nfib*nim)
        newflux = newlam
        newivar = newlam
        newmask = newlam
        count = 0
        
        for i = 0, nim -1 do begin
           newlam(*,count:count+(nfib-1L))  =  $
              transpose(transpose(lam(i,  *, *)))
           newflux(*,count:count+(nfib-1L)) = $
              transpose(transpose(flux(i, *, *)))
           newivar(*,count:count+(nfib-1L)) =  $
              transpose(transpose(ivar(i, *, *)))
           newmask(*,count:count+(nfib-1L)) =  $
              transpose(transpose(mask(i, *, *)))
           plugmapout.frames = i 
           plugmapout.expid = filelist(i)
           
           if i eq 0 then newplug = plugmapout
           if i ne 0 then newplug = [newplug, plugmapout]
           count = nfib*(i+1)
           
        endfor
        
        lam = newlam
        flux = newflux
        ivar = newivar
        mask = newmask
        plugmap = newplug
        
        junk1 = filelist
        for itest = 0, n_elements(filelist) -1 do begin
           junk = strsplit(filelist(itest), '/', /extract)
           junk1(itest) = junk(n_elements(junk)-1) 
        endfor
        
        corrfiles = 'spFluxcorr-coadd-'+$
           string(indgen(n_elements(filelist)), format='(I2.2)') $
           + '.fits'
        
        hs_frame_flux_tweak, alog10(lam), flux, ivar,$
           best_exposure, plugmap, corrfiles, /diag
        
        
        ;;Read back the corrections and apply
        
        for iexp=0, n_elements(filelist) -1 do begin
           corrfile = corrfiles(iexp)
           corrset = mrdfits(corrfile,1, /silent)
           
           indx = where(plugmap.expid eq expid[iexp])
           traceset2xy, corrset, alog10(lam[*,indx]), corrimg
           
           ;;Don't let the flux correction be more than a factor of 10
           invertcorr = 1.0 / corrimg
           tempflux = flux[*,indx]
           tempivar = ivar[*,indx]
           
           divideflat, tempflux, invvar=tempivar, $
              invertcorr, minval=0.1
           
           tempmask = mask[*, indx] OR (corrimg GE 10) $
              * pixelmask_bits('BADFLUXFACTOR')
           tempmask = mask[*, indx] OR (corrimg LE 0.1) $
              * pixelmask_bits('BADFLUXFACTOR') ;
           
           fluxout[iexp,*,*] = tempflux
           lamout[iexp,*,*] = lam[*,indx]
           ivarout[iexp,*,*] = tempivar
           pixmaskout[iexp,*,*] = tempmask
           
        endfor
        
     endif
     
  endif           
           
  ;;Now combine each of the fibers into a composite spectrum
  newlam = transpose(transpose(lamout[0,*,*]))
  newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
  newfluxivar = transpose(transpose(ivarout[0,*,*])) * 0.0
  newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
  ormask = newmask
  
  for i = 0, 299 do begin
     
     lamin = transpose(lamout(*,*,i))
     
     print, format='("Co-adding ",i4," of ",i4,a1,$)', $
        i, nfib, string(13b)
     
     hs_combine1fiber,(lamin), $
        transpose(fluxout(*,*,i)), transpose(ivarout(*,*,i)), $
        newlam=(newlam(*,i)), $
        newflux=newflux1, newivar=newivar1, $
        finalmask = transpose(pixmaskout(*,*,i)), $
        andmask = newmask1, ormask=ormask1
     newflux(*,i) = newflux1
     newfluxivar(*,i) = newivar1
     newmask(*,i) = newmask1
     ormask(*,i) = ormask1
  endfor
  
  fluxout = newflux
  ivarout = newfluxivar
  lamout = newlam
  pixmaskout = newmask
  ormaskout = ormask
  
                  
  outname = outdir + outname
  sxaddpar, objhdr, 'NAXIS2', n_elements(lamout(*,0)), $
     after='NAXIS'
  sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
                    
  mwrfits, float(lamout), outname,/create, objhdr
  mwrfits, float(fluxout), outname
  mwrfits, float(ivarout), outname
  mwrfits, float(pixmaskout), outname
  mwrfits, float(ormaskout), outname
  mwrfits, plugmapout, outname
  
  
  
END
               
                           
