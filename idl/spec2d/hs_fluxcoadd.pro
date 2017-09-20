;;This program is designed to take the spObs files and flux them and coadd them

PRO hs_fluxcoadd, filelist, rerun=rerun, outname=outname
   
   IF NOT keyword_set(outname) THEN outname = 'spHect-test.fits'
   
   junk = strsplit(outname, ' ', length=ct)
   root = strmid(outname, 7, ct-5-7)
   dfpsplot, 'qaplot-' + root + '.ps', /landscape, /color
   
   
   
   
   IF NOT keyword_set(rerun) THEN rerun = 0
   IF size(rerun, /tname) ne 'STRING' then $
      rerun=string(rerun, format='(i4.4)')

   if not keyword_set(plateid) then begin
      hdr = headfits(filelist(0))
      plateid = long(sxpar(hdr, 'mjd'))
   endif
   plateid = strn(plateid)
   
   lam = mrdfits(filelist(0), 0, objhdr)
   plugmapout = mrdfits(filelist(0), 5) 
   npix = n_elements(lam(*,0))
   nfile = n_elements(filelist)
   
   lamout = dblarr(nfile, npix, 300)
   fluxout = lamout *0.0
   ivarout = lamout * 0.0
   pixmaskout = lamout * 0.0
   
   FOR i = 0, n_elements(filelist) -1 DO BEGIN
      lamout(i,*,*) = mrdfits(filelist(i), 0)
      fluxout(i,*,*) = mrdfits(filelist(i), 1)
      ivarout(i,*,*) = mrdfits(filelist(i), 2)
      pixmaskout(i,*,*) = mrdfits(filelist(i), 3)
   ENDFOR
   
   
   
   
   for im = 0, n_elements(filelist) -1 do begin
      xxtemp = indgen(4608,300)
   for ifib = 0, 299 do begin & xxtemp(*,ifib) = xxtemp(*,0) & endfor
      xy2traceset, xxtemp, transpose(transpose(lamout[im,*,*])), $ 
         tempset, ncoeff=6
      
      flux =  transpose(transpose(fluxout[im,*,*])) 
      ivar =  transpose(transpose(ivarout[im,*,*]))
      
      telluric_factor = hs_telluric_corr(flux, ivar, tempset,$
                                         plugmapout)
      divideflat, flux, telluric_factor, invvar = ivar
      
      fluxout[im,*,*] = flux
      ivarout[im,*,*] = ivar
      
   ENDFOR
   
   
   
   
   fcalfile = strarr(n_elements(filelist))
   scalfile = fcalfile(0)
   
   for i = 0, n_elements(filelist) -1 do begin
      junk = strsplit(filelist(i), ' ', length=ct)
      filelist1 = strmid(filelist(i), 6, ct-6-5)
      
      
      
      fcalfile(i) = 'hsFluxcalib-' + $
         filelist1 + '.fits'
      if i eq 0 then scalfile = 'hsStd-' + filelist1 $
         + '.fits'
   endfor
       
   flux_factor = hs_fluxcorr(lamout, fluxout,$
                             ivarout, plugmapout, $
                             pixmaskout, outname = filelist, $
                             plateid=plateid, $
                             rerun=rerun,/psplot, checkaverage=checkaverage,$
                             stand=stand, header=objhdr, cut=cut, $
                             superfile=superfile, useaverage=useaverage, $
                             fcalfile=fcalfile, scalfile=scalfile) 
   
   fcalfile = strarr(n_elements(filelist))
   sphoto_err = dindgen(n_elements(filelist))
   sn=sphoto_err
    
   for im = 0, n_elements(filelist) -1 do begin
      outname1 = strsplit(filelist(im), '/,.fits', /extract)
      ct = n_elements(outname1)
      outname1= outname1(ct-1)
      
      calibhdr = headfits(fcalfile(im))
      
      sphoto_err(im) = sxpar(calibhdr, 'SPHOROERR')
      flux =  transpose(transpose(fluxout[im,*,*]))
      factor = transpose(transpose(flux_factor[im,*,*]))
      ivar = transpose(transpose(ivarout[im,*,*]))
      
      hs_snr,  transpose(transpose(lamout[im,*,*])), $
         transpose(transpose(fluxout[im,*,*])), plugmapout, $
         sn2=sn21, /noplot
      
      sn(im) =sn21
      
      divideflat, flux, factor,$
         invvar=ivar, minval=0.05*median(factor[*,0])
      
      fluxout[im,*,*] = flux
      ivarout[im,*,*] = ivar
   endfor
    
   expid = filelist
   sn_exp = sn
   sphoto_exp = (1/sphoto_err)
   maxval = max(sn_exp/median(sn_exp)+$
                sphoto_exp/median(sphoto_exp), iframe_best)
   
   if not keyword_set(best_exposure) then $
      best_exposure = expid[iframe_best]
   splog, 'Best Exposure is: ' + best_exposure
   
   ;;--------------------------------------------------------------
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
         newlam(*,count:count+(nfib-1L))  = $
            transpose(transpose(lam(i,  *, *)))
         newflux(*,count:count+(nfib-1L)) = $ 
            transpose(transpose(flux(i, *, *)))
         newivar(*,count:count+(nfib-1L)) = $
            transpose(transpose(ivar(i, *, *)))
         newmask(*,count:count+(nfib-1L)) = $ 
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
      
      
      
      corrfiles = 'spFluxcorr-' + junk1
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
         
         tempmask = mask[*, indx] OR $
            (corrimg GE 10) * pixelmask_bits('BADFLUXFACTOR')
         tempmask = mask[*, indx] OR $
            (corrimg LE 0.1) * pixelmask_bits('BADFLUXFACTOR') ;
         
         fluxout[iexp,*,*] = tempflux
         lamout[iexp,*,*] = lam[*,indx]
         ivarout[iexp,*,*] = tempivar
         pixmaskout[iexp,*,*] = tempmask
         
      endfor
   endif
   
   size1 = size(lamout)
   nim = size1(1)
   npix = size1(2)
   nfib = size1(3)
   
   
   
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
         newlam=(newlam(*,i)), newflux=newflux1,$
         newivar=newivar1, $
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
   
   finalflux=fluxout
   finalivar = ivarout
   finalandmask = pixmaskout
   finalormask = ormaskout
   finalloglam = alog10(newlam)
   
   
   
   
   
   ;; ----------------------------------------------------------------------
   ;; Correct objects for foregound Milky Way dust using SFD Maps
   ;;------------------------------------------------------------------------
   
   
   newplug_struct = create_struct(plugmapout[0],$
                                  'tsobjid', lonarr(5), $
                                  'tsobj_mag', fltarr(5), 'ebv_sfd', 0.0)
   newplug = make_array(val=newplug_struct, dim=nfib)
   struct_assign, plugmapout, newplug
   newplug.tsobj_mag = plugmapout.mag
   plugmapout = newplug
   
   dustdir = getenv('DUST_DIR')
   if keyword_set(dustdir) then $
      dustmaps = $
      file_test(dustdir + '/maps/SFD_dust_4096_ngp.fits') $
   else dustmaps = 0
   
   if dustmaps then begin
      ;; Get values of E(B-V) for each ra/dec from SFD maps
      glactc, plugmapout.ra, plugmapout.dec, 2000, gall, $
         galb, 1, /degre
      ebv_sfd = dust_getval(gall, galb, ipath = dustdir $
                            + '/maps/', /interp)
      plugmapout.ebv_sfd = ebv_sfd 
      
      ;; Deredden the flux and ivar - use O'donnel extinction curve 
      A_v = 3.1 * ebv_sfd
      
      for i = 0, nfib -1 do begin
         a_odonnell = ext_odonnell(10.0^finalloglam[*,i], 3.1)
         e_tau = exp(a_odonnell # A_v / 1.086)
         finalflux[*,i] = finalflux[*,i] * e_tau
         finalivar[*,i] = finalivar[*,i] / e_tau^2
         
      endfor
   endif else begin
      splog, 'WARNING!!!!  Foreground Dust correction NOT applied' 
      ebv_sfd = 0.0
   endelse
   fluxout = finalflux
;   finalivar = ivarout
   ivarout = finalivar
   pixmaskout = finalandmask
   ormaskout = finalormask
   
      
   splog, 'Writing coadded sphect file ' + outname
   
   
   sxaddpar, objhdr, 'NAXIS2', n_elements(lamout(*,0)), after='NAXIS'
   sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
   mwrfits, float(lamout), outname, /create, objhdr
   mwrfits, float(fluxout), outname
   mwrfits, float(ivarout), outname
   mwrfits, float(pixmaskout), outname
   
   
   mwrfits, float(ormaskout), outname
   mwrfits, plugmapout, outname
   
   k = where(plugmapout.objtype EQ  'SPECTROPHOTO_STD')
   
   
   FOR ii = 0, n_elements(k) -1 DO BEGIN
      
      standfile = getenv('HSRED_DIR')+'/etc/standstar.dat'
      readcol, standfile, sra, sdec, u, g, r, i, z, $
         reddening, format='D,D,D,D,D,D,D,D', $
         comment='#', /silent
      
      dist = sqrt( (plugmapout(k(ii)).ra-sra)^2 + $
                   (plugmapout(k(ii)).dec-sdec)^2)*3600.0
      junk = min(dist, j)
  
      maggies = k_project_filters(lamout(*,k(ii)), $
                                  1d-17*fluxout(*,k(ii)))
      mag = -2.5*alog10(maggies)
      smag = [u(j), g(j), r(j), i(j), z(j)]
      leff = k_lambda_eff()
      
      diff = smag - mag
      IF ii EQ 0 THEN out = dblarr(3, n_elements(k))
      out(*,ii) = diff(1:3)
      plotsym, 0, /fill
      color = ['red', 'orange', 'green', 'blue', 'magenta']
      color = [color, color, color, color]
      IF ii EQ 0 THEN $
         djs_plot, leff(1:3), diff(1:3), ps=8, yrange=[-0.2, 0.2], $
         xrange=[4000, 8000], color=color(ii), $
         xtitle='Effective Wavelength', $
         ytitle='De-reddened SDSS Mag - AGES Synth Mag', $
         charthick=2, charsize=1.2
      
      IF ii NE 0 THEN djs_oplot, leff(1:3), diff(1:3), ps=8, color=color(ii)
      
   ENDFOR
   
   djs_iterstat, out, median=median, sigma=sigma
   djs_oplot, [0, 1e9], median*[1,1], linestyle=3, thick=3, color='red'
   djs_oplot, [0, 1e9], median + sigma * [-1, -1], color='blue', $
      linestyle=2, thick=3
   
    djs_oplot, [0, 1e9], median + sigma * [1, 1], color='blue', $
      linestyle=2, thick=3
   
    dfpsclose
end
