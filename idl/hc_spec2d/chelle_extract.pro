PRO chelle_extract, filelist, outroot=outroot, rerun=rerun, $
                    nosky=nosky, noskysub=noskysub
   
   IF NOT keyword_set(rerun) THEN rerun = '0000' 
   IF size(rerun, /tname) NE 'STRING' THEN $
      rerun=string(rerun, format='(i4.4)')
   
   IF NOT keyword_set(arclist) THEN arclist = 'calibration/'+rerun+'/wset.fits'
   nfile = n_elements(filelist)
   IF n_elements(arclist) EQ 1 THEN arclist = replicate(arclist(0), nfile)    
   
   IF direxist('reduction') EQ 0 THEN spawn, 'mkdir reduction'
   IF direxist('reduction/' + rerun) EQ 0 THEN $
      spawn, 'mkdir reduction/' + rerun
   
   IF NOT keyword_set(outroot) THEN outroot='test'
   outname = 'spChelle-' + outroot + '.fits'
   outname = 'reduction/'+ rerun + '/' + outname
   
   
   
   FOR ifile = 0, nfile -1 DO BEGIN
      
      FOR ccdnum = 1, 2 DO BEGIN
         
         ;;read the wavelength solution
         wset = mrdfits(arclist(ifile), ccdnum)
         
         ;;Extract the spectra
         hc_getspec, filelist(ifile), ccdnum, specflux=specflux, $
            specivar=specivar, header=objhdr, rerun=rerun
         
         plugmap = mrdfits(filelist(ifile), 4+ccdnum)
         objtype = strupcase(plugmap.target)
         temp = {objtype:objtype(0), fiberid:0L, xfocal:0.0, yfocal:0.0}
         temp = replicate(temp, 120)
         temp.objtype = strupcase(plugmap.object)
         temp.fiberid = findgen(120)
         temp.xfocal = plugmap.platex
         temp.yfocal = plugmap.platey
         plugmap = struct_addtags(plugmap, temp)
         
        ;;------------------
        ;; Apply heliocentric correction
        ;; Correct LAMBDA, which is used to shift to vacuum wavelengths.
 
        helio=0.0
        ra = sxpar(objhdr,'RA')
        dec = sxpar(objhdr,'DEC')
        
        
        ;;--------------------------------------------------------
        ;; Call standard proc to determine time-stamps
        
        mjd_in = double(sxpar(objhdr,'MJD'))
        jd = 2400000.5D + mjd_in
        
        if ( size(ra, /tname) NE 'INT' $
             AND size(dec, /tname) NE 'INT' $
             AND size(tai_mid, /tname) NE 'INT' ) then begin
           helio = heliocentric(ra, dec, jd=jd)
           splog, 'Heliocentric correction = ', helio, ' km/s'
           sxaddpar, objhdr, 'HELIO_RV', helio, $
              ' Heliocentric correction (added to velocities)'
        endif else begin
           splog, 'WARNING: Header info not present' + $
              ' to compute heliocentric correction'
        endelse
    
        ;;Apply the correction
        delvarx, xx
        traceset2xy, wset, xx, lamtemp
        lamtemp = lamtemp / (1+helio/299792.458)
        xx1 = indgen(n_elements(lamtemp(*,0)), 150)
        for i = 0, 150 -1 do begin
           xx1(*,i) = xx1(*,0)
        endfor
        
        xy2traceset, xx1, lamtemp, wset, ncoeff=6
        airset=wset
        
        ;;Flat field the data
        IF file_test('calibration/' + rerun + '/traceflat.fits') EQ 0L $
           THEN BEGIN
           hs_makesuperflat, rerun=rerun, /chelle
        ENDIF
        
        fflat = mrdfits('calibration/' + rerun + '/traceflat.fits', $
                        ccdnum-1, /silent)
        divideflat, specflux, fflat, invvar=specivar
        
        
        ;;Apply sky flat
        IF file_test('calibration/' + rerun + '/sflat.fits') EQ 0 THEN $
           nosky = 1
        if not keyword_set(nosky) then begin
           ;;The fiberflat.fits is the higher order 
           ;;correction usually made from the skyflat
           if file_test('calibration/'+rerun+'/'+'fiberflat.fits')$
              EQ 0L then begin 
              hs_makefiberflat, rerun=rerun, root=root, /usesky, /chelle
           endif
           
           fflat = mrdfits('calibration/'+rerun+'/'+$
                           'fiberflat.fits',ccdnum-1, /silent)
          
           divideflat, specflux, fflat, invvar=specivar
           
        ENDIF
        
                
        ;;DO the sky subtraction
        
        traceset2xy, airset, xx, lam
        
        IF NOT keyword_set(noskysub) THEN begin
           
           k = where(strupcase(plugmap.object) EQ 'SKY', ct)
           FOR i = 0, ct -1 DO BEGIN
              IF i EQ 0 THEN BEGIN
                 superwave = lam(*,k(0))
                 superflux = specflux(*,k(0))
                 superivar = specivar(*,k(0))
              endif
              
              IF i NE 0 THEN BEGIN
                 superwave = [superwave, lam(*,k(i))]
                 superflux = [superflux, specflux(*,k(i))]
                 superivar = [superivar, specivar(*,k(i))]
              endif
              
           endfor
           
           bkpts=0
           everyn = (2*ct/3) > 1
           sset = bspline_iterfit(superwave, superflux, $
                                  invvar=superivar, $
                                  yfit = skyfit, bkpts=bkpts, $
                              everyn=everyn)
           sky = bspline_valu(lam, sset)
           objsub = specflux - sky
           objsubivar = specivar + 1/sky
           newmask = objsub*0
           
        ENDIF ELSE BEGIN
           objsub = specflux
           objsubivar = specivar
           newmask = objsub*0
        endelse
        

        
           
        ;skystruct = hs_skysubtract(specflux, $
        ;                           specivar, plugmap, $
        ;                           airset, objsub, objsubivar, $
        ;                           newmask=newmask, sset=sset1)
       ; 
        IF ccdnum EQ 1 THEN BEGIN
           lam1 = lam
           flux1 = objsub
           fluxivar1 = objsubivar
           pixmask1 = newmask
           plugmap1 = plugmap
        ENDIF
        
        IF ccdnum EQ 2 THEN BEGIN
           lam2 = lam
           flux2 =objsub
           fluxivar2 = objsubivar
           pixmask2 = newmask
           plugmap2 = plugmap
           
           specflux = [[flux1], [flux2]]
           specfluxivar = [[fluxivar1], [fluxivar2]]
           pixmask = [[pixmask1], [pixmask2]]
           lam = [[lam1],[lam2]]
           plugmapout = [plugmap1, plugmap2]
           
           temp = {frames : 1, snr:0.0}
           temp = replicate(temp, 240)
           plugmapout = struct_addtags(plugmapout, temp)
           
           
        ENDIF
        
     ENDFOR
     
     
     
     IF ifile EQ 0 THEN BEGIN
        lamout = dblarr(nfile, n_elements(flux2(*,0)), 240)
        fluxout = lamout *0.0
        ivarout = lamout *0.0
        pixmaskout = lamout *0.0
     ENDIF
     
     
     lamout[ifile, *,*] = lam
     fluxout[ifile,*,*] = specflux
     ivarout[ifile, *,*] = specfluxivar
     pixmaskout[ifile, *, *] = pixmask
     
     hs_snr,  transpose(transpose(lamout[ifile,*,*])), $
        transpose(transpose(fluxout[ifile,*,*])), plugmapout, $
        snr=snr, /noplot, /nomag, /middle
     plugmapout.snr = snr
          
     junk = strsplit(filelist(ifile), './', /extract)
     outnum = junk(n_elements(junk)-2)
     outfile = outroot + '-' + outnum + '.fits'
     outfile = 'spObs-' + outfile
     outfile = 'reduction/' + rerun + '/'+ outfile
     
     sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
     mwrfits, lam, outfile, objhdr,  /create
     mwrfits, specflux, outfile
     mwrfits, specfluxivar, outfile
     mwrfits, pixmaskout, outfile
     mwrfits, pixmaskout, outfile
     mwrfits, plugmapout, outfile
      
  ENDFOR
  
  
  
  
  
  ;;Now combine each of the fibers into a composite spectrum
  newlam = transpose(transpose(lamout[0,*,*]))
  newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
  newivar = transpose(transpose(ivarout[0,*,*])) * 0.0
  newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
  ormask = newmask
    
  
  nfib = 240
  for i = 0, nfib-1 do begin
     
  
     lamin = transpose(lamout(*,*,i))
     print, format='("Co-adding ",i4," of ",i4,a1,$)', $
        i, nfib, string(13b)
     
     nim = n_elements(fluxout(*,0,0))
     val = fltarr(n_elements(fluxout(*,0,0)))
     
;     FOR ii = 0, nim -1 DO BEGIN
;        val(ii) = djs_median(fluxout(ii,*,i))
;        fluxout(ii,*,i) = fluxout(ii,*,i)/abs(val(ii)/val(0))
;     ENDFOR
     

    
     hs_combine1fiber,(lamin), $
        transpose(fluxout(*,*,i)), transpose(ivarout(*,*,i)), $
        newlam=(newlam(*,i)), $
        newflux=newflux1, newivar=newivar1, $
        finalmask = transpose(pixmaskout(*,*,i)), $
        andmask = newmask1, ormask=ormask1
     newflux(*,i) = newflux1
     newivar(*,i) = newivar1
     newmask(*,i) = newmask1
     ormask(*,i) = ormask1
  endfor
  
  
  
   fluxout = newflux
   ivarout = newivar
   lamout = newlam
   pixmaskout = newmask
   ormaskout = ormask
   
   
   hs_snr,  lamout, fluxout, plugmapout, $                    
      snr=snr, /noplot, /nomag, /middle     
   plugmapout.snr = snr
   
   if NOT keyword_set(outname) then $
      outname = plugfile + '-' + rerun + '.fits'
   
   splog, 'Writing coadded sphect file ' + outname
   
   
   sxaddpar, objhdr, 'NAXIS2', n_elements(lamout(*,0)), after='NAXIS'
   sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
   mwrfits, float(lamout), outname, /create, objdr
   mwrfits, float(fluxout), outname
   mwrfits, float(ivarout), outname
   mwrfits, float(pixmaskout), outname
   mwrfits, float(ormaskout), outname
   mwrfits, plugmapout, outname
       
END




  
