;+
; NAME:
;   hs_readspec
;
; PURPOSE:
;  Main reduction workhorse
;
; CALLING SEQUENCE:
;  hs_readspec, filelist,  skysubtract=skysubtract, $
;     fiberflat=fiberflat,  tweaksky=tweaksky, $
;     plugcat=plugcat, combine=combine, $
;     lamout=lamout, fluxout=fluxout, fluxcorr=fluxcorr, $
;     tellcorr=tellcorr, uberextract=uberextract, $
;     ivarout=ivarout, rerun=rerun,  outname=outname,  $
;     plateid=plateid, $
;     qaplot=qaplot, quicklook=quicklook, quickplot=quickplot,$
;     stand=stand,  checkaverage=checkaverage, $
;     docosmic=docosmic, superfile=superfile, dostand=dostand
;
; INPUTS:
;  filelist - list of files to be reduced
;
; OPTIONAL KEYWORDS:
;   skysubtract - turn on sky subtraction
;   fiberflat   - turn on flat fielding
;   tweaksky    - use skylines to tweak the intial wavelength 
;                 solution
;   plugcat     - turn this on if you are using specially formated
;                 plugcat files for your reductions. If you are not sure, 
;                 then you are probably not using these files
;   plugfile    - If plugcat is set, then plugfile must be set to the plugcat 
;                 file for the reduction
;   combine     - Set this to combine the multiple observations. The spectra
;                 are extracted and then combined
;   fluxcorr    - Set this to flux the data with F star standards
;                 NOTE: Fluxing must have at least two observations
;   tellcorr    - Set this to telluric correct the data.  You need 
;                 standard stars in order for this to work
;   uberextract - Set this to turn on tweaksky, fiberflat, skysubtract
;                 tellcorr, and fluxcorr. Obivously, you don't want this 
;                 if you don't have standards
;   rerun       - This should be set to the rerun number you would like to 
;                 for the reduction. If not set rerun=0000
;   outname     - Root outname you would like to use for the data.  If
;                 this is not set, then outname='test'
;   qaplot      - Setting this results in the quality assurance plots to be
;                 be created
;   superfile   - If this is set, the file $HSRED_DIR/etc/superstand.fits
;                 instead of typing the stars
;   docosmic    - If docosmic is set then the cosmic ray rejection is turned 
;                 on.
;   checkaverage- This will check the current flux soln against an archival 
;                 solution and will truncate the data if the counts of the 
;                 standard stars are fewer than 100 
;   dostand    - This will provide the standard flags that you need
;                 (skyssubtract, fiberflat, tweaksky, combine,
;                 writeout, writeall)
; 
; OUTPUTS:
;      Creates a variety of files (depending on settings) in the 
;      /reduction/YYYY/ direcory where YYYY is the rerun
;
; OPTIONAL OUTPUTS:
;   lamout       - the final wavelength map
;   fluxout      - the final flux
;   ivarout      - the final inverse variance map
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;    The flat field is generated after tracing.  A proper flat field would
;    alsop have a 2d component.  There are issues with this.  
;   
;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA         
;                
;-
;------------------------------------------------------------------------------
PRO hs_readspec, filelist, plugfile=plugfile, skysubtract=skysubtract, $
                fiberflat=fiberflat,  tweaksky=tweaksky, $
                plugcat=plugcat, combine=combine, $
                lamout=lamout, fluxout=fluxout, fluxcorr=fluxcorr,$
                tellcorr=tellcorr, uberextract=uberextract, $
                ivarout=ivarout, rerun=rerun, outname=outname, $
                plateid=plateid, $
                qaplot=qaplot, quicklook=quicklook, quickplot=quickplot, $
                stand=stand,  checkaverage=checkaverage, $
                docosmic=docosmic, doplot=doplot, $
                superfile=superfile, useaverage=useaverage, $
                noftweak=noftweak, nobias=nobias, dodark=dodark, $
                nosky=nosky, debug=debug, doskyscale=doskyscale, $
                skyscale=skyscale, dostand=dostand, writeall=writeall, $
                writeout=writeout, oldmap=oldmap
  
  
  
  
  if NOT keyword_set(plugcat) then begin
     plugfile = strsplit(filelist(0), '/', /extract)
     plugfile = plugfile(n_elements(plugfile)-1)
     junk  = strsplit(plugfile, ' ', length=ct)
     plugfile = strmid(plugfile,0,ct-5)
  endif
  
  if keyword_set(quicklook) then begin   
     hs_readspec, filelist, plugfile=plugfile, /tweaksky, $
        /fiberflat, /skysubtract, plugcat=plugcat, /writeall, $
        rerun=rerun, stand=stand, checkaverage=checkaverage, $
        docosmic=docosmic, /writeout, outname=outname, qaplot=qaplot, $
        superfile=superfile, useaverage=useaverage, nobias=nobias, $
        dodark=dodark, nosky=nosky, debug=debug, $
        doskyscale=doskyscale, skyscale=skyscale, oldmap=oldmap
  endif
     
     
        
  if keyword_set(uberextract) then begin
     hs_readspec,  filelist, plugfile=plugfile, /skysubtract, $
        /fiberflat,  /tweaksky, $
        plugcat=plugcat, /combine, /writeout,$
        lamout=lamout, rerun=rerun,  fluxout=fluxout, $
        /fluxcorr, /tellcorr, ivarout=ivarout, $
        outname=outname, /writeall, $
        plateid=plateid, qaplot=qaplot, quickplot=quickplot, $
        stand=stand, checkaverage=checkaverage, docosmic=docosmic, $
        superfile=superfile, useaverage=useaverage, noftweak=noftweak, $
        nobias=nobias, dodark=dodark, nosky=nosky, debug=debug, $
        doskyscale=doskyscale, skyscale=skyscale, oldmap=oldmap
     return
  endif
  
  if keyword_set(dostand) then begin
     hs_readspec, filelist, plugfile=plugfile, /skysubtract, $
        /fiberflat, /tweaksky, plugcat=plugcat, /combine, $
        /writeout, /writeall, fluxout=fluxout, lamout=lamout, rerun=rerun, $
        ivarout=ivarout, outname=outname, plateid=plateid, $
        qaplot=qaplot, stand=stand, docosmic=docosmic, superfile=superfile, $
        noftweak=noftweak, oldmap=oldmap
     return
  endif
  
  ;;These are thing you will always want - they are now hardwired
  writeout = 1
  writeall = 1
  
  if not keyword_set(rerun) then rerun= string(0,format='(i4.4)')
  if size(rerun, /tname) ne 'STRING' then rerun=string(rerun, format='(i4.4)')
  if direxist('reduction') eq 0L then spawn, 'mkdir reduction'
  if direxist('reduction/' + rerun) eq 0L then spawn, 'mkdir reduction/'+rerun
  outdir = 'reduction/' + rerun + '/'
  
  if not keyword_set(plateid) then begin
     hdr = headfits(filelist(0), exten=1)
     plateid = long(sxpar(hdr, 'mjd'))
  endif
  
  if keyword_set(outname) then begin
     outname = strsplit(outname, '/', /extract)
     outname = outname(n_elements(outname)-1)
     if strmid(outname, 0, 7) ne 'spHect-' then $
        outname = 'spHect-'+outname
  endif else begin
     outname = 'spHect-test.fits'
     splog, 'No out name given - using spHect-test.fits'
  endelse
  
  outname = outdir + outname   
  
  if keyword_set(qaplot) then begin
     
     ;;Start the qaplot - 
     
     outplot = strsplit(outname, '/', /extract)
     junk = strsplit(outplot(n_elements(outplot)-1), ' ', length=ct)
     
     outplot = strmid(outplot(n_elements(outplot)-1), 7, ct-7)
     junk = strsplit(outplot, ' ', length=ct)
     outplot = strmid(outplot, 0, ct-5)
     
     dfpsplot, 'reduction/' + rerun + '/qaplot_' + outplot + '.ps', $
        /color, /isolatin1
  endif
  
  if file_test('calibration/'+rerun+'/traceset.fits') eq 0L then begin
     ;;Make the traceset
     hs_maketraceset, /write, rerun=rerun, mjd=root
  endif
  
  skysingle = ''
  skymulti = ''
  
  for j = 0, n_elements(filelist) -1 do begin
     for ccdnum = 1, 2 do begin
        


         hs_proc, filelist(j), ccdnum, spec, specivar, pixelmask=pixmask, $
           rerun=rerun, root=root, header=objhdr, docosmic=docosmic,$
           nobias=nobias
        
        if keyword_set(dodark) then begin
           darkfile = lookforgzip('calibration/' + rerun + '/dark.fits')
           if file_test(darkfile) eq 0 then begin
              splog, 'NO DARKS FOUND!!'
              return
           endif
           
           hs_proc, darkfile,ccdnum, dark,$
              darkivar, rerun=rerun, root=root, header=darkheader, $
              nobias=nobias
           
           dtime = sxpar(darkheader, 'DARKTIME')
           time  =sxpar(objhdr, 'DARKTIME')
           scale = time/dtime
           dark = dark*scale
           splog, 'Applying Dark Correction with scale = ' + string(scale)
           spec = spec-dark
          delvarx, darkcase
        endif
        
        sigma= 2.0
        proftype = 2.0
        highrej = 30
        lowrej = 30
        npoly = 4
        wfixed = [1, 1, 1 ]
        
        
        hs_readtraceset, ccdnum , xsol=xsol, rerun=rerun
        
        hs_extract_image, spec, specivar, xsol, sigma,$
           specflux, specfluxivar, $
           proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
           npoly=npoly, relative=1, ymodel=ymodel
                
        delvarx, spec
        delvarx, specivar
        
        if file_test('calibration/'+rerun+'/'+'wset.fits') EQ 0L then begin
           hs_findwave, rerun=rerun, root=root
        endif
        
        wset = mrdfits('calibration/'+rerun+'/'+'wset.fits', ccdnum, /silent)
        split = strsplit(filelist(j), '/,.fits', /extract)
			 
        if keyword_set(plugcat) then begin
           hs_cattoplug, plugfile, plugmap, ccdnum=ccdnum, oldmap=oldmap
           if ccdnum eq 1 then other = 2
           if ccdnum eq 2 then other = 1
           hs_cattoplug, plugfile, plugcheck, ccdnum= other, oldmap=oldmap
        endif else begin
           if ccdnum eq 1 then other = 2
           if ccdnum eq 2 then other = 1
           hs_maptoplug, plugfile + '_map', plugmap, ccdnum=ccdnum
           hs_maptoplug, plugfile + '_map', plugcheck ,ccdnum=other
        endelse
                

        
        ratio = hs_scalefix( plugmap)
        scale = specflux*0.0
        for i = 0, n_elements(specflux(0,*))-1 do $
           scale(*,i) = replicate(ratio(i), n_elements(specflux(*,0)))
        divideflat, specflux, scale, invvar=specivar
        ;;Apply the fiberflat now if specified
        if keyword_set(fiberflat) then begin
           ;;Traceflat.fits is the flat made with the dome flat (basic
           ;;correction plug fringing)
           if file_test('calibration/'+rerun+'/'+'traceflat.fits')$
              EQ 0L then begin 
              hs_makesuperflat, rerun=rerun, root=root
           endif
           fflat = mrdfits('calibration/'+rerun+'/'+ $
                           'traceflat.fits',ccdnum-1, /silent)
           print, ccdnum
           divideflat, specflux, fflat, invvar=specivar

           if file_test('calibration/' + rerun + '/sflat.fits') EQ 0 THEN $
             nosky = 1


           if not keyword_set(nosky) then begin
              ;;The fiberflat.fits is the higher order 
              ;;correction usually made from the skyflat
              if file_test('calibration/'+rerun+'/'+'fiberflat.fits')$
                 EQ 0L then begin 
                 hs_makefiberflat, rerun=rerun, root=root, /usesky
              endif
              
              fflat = mrdfits('calibration/'+rerun+'/'+$
                              'fiberflat.fits',ccdnum-1, /silent)
              divideflat, specflux, fflat, invvar=specivar
           endif
           
	 if keyword_set(qaplot) then begin
              traceset2xy, wset, xx, temp
              temp1 = fltarr(n_elements(temp(*,0)), n_elements(temp(0,*)))
              for ifib = 0, n_elements(temp1(0,*)) -1 do begin
                 temp1(*,ifib) = findgen(n_elements(temp1(*,0)))
              endfor
              temp = alog10(temp)
              xy2traceset, temp1, temp, tempset, ncoeff=6
               
              hs_qaplot_fflat, fflat, tempset
              DELVARX, FLAT
              delvarx, fflat
              delvarx, tempset
           endif
           
        endif
        if keyword_set(doskyscale) then begin
           sscale = mrdfits(skyscale, ccdnum-1)
           divideflat, specflux, sscale, invvar=specfluxivar
        endif
        
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
        
        if keyword_set(tweaksky) then begin
           iskies = where(strmatch(plugmap.objtype, 'SK*'))
           
           hs_tweaksky, specflux, specfluxivar, wset, ccdnum=ccdnum, $
              airset=airset, xsky=xsky, rerun=rerun, $
              filename=filelist(j), skywaves=xwaves, $
              arcshift=arcshift
           
           delvarx, xx
           
           if keyword_set(qaplot) then begin
              traceset2xy, wset, xx, temp
              temp1 = fltarr(n_elements(temp(*,0)), n_elements(temp(0,*)))
              
              for ifib = 0, n_elements(temp1(0,*)) -1 do begin
                 temp1(*,ifib) = findgen(n_elements(temp1(*,0)))
              endfor
              temp = alog10(temp)
              xy2traceset, temp1, temp, tempset, ncoeff=6   
              hs_qaplot_skyshift, tempset, xsky, xwaves, arcshift
              delvarx, xsky
              delvarx, xwaves
           endif
           
        endif else begin
           airset = wset
           
        endelse
        
   
        
        nccd = n_elements(where(strmid(plugmap.objtype,0,3) eq 'SKY'))
        nother = n_elements(where(strmid(plugcheck.objtype,0,3) eq 'SKY'))
        traceset2xy, wset, xx, lam11
       
        k = where(strmatch(plugmap.objtype, 'SKY*'))
        
        if k(0) ne -1 then begin
           
           ;;Check at 4000 Angstroms
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1 do begin
              kcheck(ifib) = where(abs(lam11[*,ifib]-4000.0) $
                                   eq min(abs(lam11[*,ifib]-4000.0)))
           endfor         
           
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 2*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT'
           
           
           ;;Check at 5577 (this will keep the bspline from 
           ;;dying if there is a zero fiber 
           
           kcheck = indgen(n_elements(k))
           
           for ifib = 0,  n_elements(k)-1 do begin
              kcheck(ifib) = where(abs(lam11[*,ifib]-5577.0) $
                                   eq min(abs(lam11[*,ifib]-5577.0)))
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 2*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT'
           
           ;;Check at 5800 for continuum sources
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1  do begin
              kcheck(ifib) = where(abs(lam11[*,ifib]-5800.0) $
                                   eq min(abs(lam11[*,ifib]-5800.0)))
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 2*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT'
           
           ;;Check at 8725 and 8975 for the Red monster
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1 do begin
              kcheck(ifib) = where(abs(lam11[*,ifib]-8725.0) $
                                   eq min(abs(lam11[*,ifib]-8725.0)))
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 2*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT'
           
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1  do begin
              kcheck(ifib) = where(abs(lam11[*,ifib]-9150.0) $
                                   eq min(abs(lam11[*,ifib]-9150.0)))
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 2*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT'
           
        endif
        
        junk = strsplit(filelist(j), ' ', length=filecount)
        junk = strsplit(cwd(), ' ', length=pathcount)
        
        if keyword_set(qaplot) or keyword_set(doplot)  then begin
           k = where(strmatch(plugmap.objtype, 'SKY*'))
           if k(0) ne -1 then PLOT, LAM11[*,K],specFLUX[*,K], $
              xrange=[5000,6000]
           k = where(strmatch(plugmap.objtype, 'RE*'))
           if k(0) ne -1 then djs_oPLOT, LAM11[*,K],specFLUX[*,K], color='red'
           k = where(strmatch(plugmap.objtype, 'SKY*'))
           if k(0) ne -1 then PLOT, LAM11[*,K],specFLUX[*,K], $
              xrange=[6000,9200]
           k = where(strmatch(plugmap.objtype, 'RE*'))
           if k(0) ne -1 then djs_oPLOT, LAM11[*,K],specFLUX[*,K], color='red'
        endif
        
        
         if keyword_set(skysubtract) and nother gt 3 and nccd gt 3 $
            then begin
            
            traceset2xy, airset, xx,lam
           ; hs_skyscale, lam, specflux, plugmap, scale=scale, debug=debug
           ; hs_scale, lam, specflux,  scale=scale
            scale_sky, lam, specflux, specfluxivar, scale=scale
           
            divideflat, specflux, scale, invvar=specfluxivar, minval=0.01
            
            k = where(strmatch(plugmap.objtype, 'SKY*'))
            
            traceset2xy, airset, xx,lam
            xx = indgen(4608,n_elements(k))
            for i = 0,n_elements(k)-1 do begin & xx(*,i) = xx(*,0) & endfor
            xy2traceset, xx, lam[*,k], tempset, ncoeff=6
            delvarx, xx
            
            skystuct = hs_skysubtract(specflux[*,k],$
                                    specfluxivar[*,k], plugmap[k], $
                                    tempset,   objsub1, objsubivar1, $
                                    newmask=newmask, sset=sset1)
            skyval = findgen(n_elements(k))
		

            for i = 0, n_elements(k) - 1 do begin
               
               l = where(abs(lam[*,k(i)]-5000) eq min(abs(lam[*,k(i)]-5000)))
               m = where(abs(lam[*,k(i)]-6000) eq min(abs(lam[*,k(i)]-6000)))
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            
            bad = where(skyval gt mean+2.5*sigma)
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT-56'
            skyval = findgen(n_elements(k))
            
            for i = 0, n_elements(k) - 1 do begin
               l = where(abs(lam[*,k(i)]-6000) eq min(abs(lam[*,k(i)]-6000)))
               m = where(abs(lam[*,k(i)]-7000) eq min(abs(lam[*,k(i)]-7000)))
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            bad = where(skyval gt mean+2.5*sigma)
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT67'
            skyval = findgen(n_elements(k))
            
            for i = 0, n_elements(k) - 1 do begin
               l = where(abs(lam[*,k(i)]-7000) eq min(abs(lam[*,k(i)]-7000)))
               m = where(abs(lam[*,k(i)]-8000) eq min(abs(lam[*,k(i)]-8000)))
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            bad = where(skyval gt mean+2.5*sigma)
         
                        
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT67'
            skystuct = hs_skysubtract(specflux, specfluxivar, plugmap, $
                                      airset,   objsub1, objsubivar1, $
                                      newmask=newmask1, sset=sset)
            
          
            skymulti='true'
            if ccdnum eq 1 then sset1 = sset
            if ccdnum eq 1 then scale1 = scale
            if ccdnum eq 2 then sset2 = sset
            if ccdnum eq 2 then scale2 = scale
            
            pixmask = newmask1
            
            objsub = objsub1
            objsubivar = objsubivar1
          
            if keyword_set(qaplot) then begin
             split = strsplit(filelist(j), '/,.fits', /extract)
             
             traceset2xy,airset, xx, temp
             temp1 = temp
             for ifib = 0, n_elements(temp(0,*)) -1 do begin
                temp1(*,ifib) = findgen(n_elements(temp(*,0)))
             endfor
             temp = alog10(temp)
             xy2traceset, temp1, temp, tempset, ncoeff=5
             iskies = where(strmatch(plugmap.objtype, 'SKY'))
             
             hs_qaplot_skysub, specflux, specfluxivar, $
                objsub, objsubivar, tempset, iskies
             traceset2xy, airset, xx, lam
             
             k = where(strmatch(plugmap.objtype, 'SKY*'))
             
             for i = 0, n_elements(k) -1 do begin
                if i eq 0 then plot, lam[*,k(i)],$
                   median(objsub[*,k(i)],100), /xstyle
                if i gt 0 then oplot, lam[*,k(i)], median(objsub[*,k(i)],100)
                
                if ccdnum eq 1 then factor=0
                if ccdnum eq 2 then factor=150
                djs_xyouts, max(smooth(objsub[*,k(i)],100))*0.90, 4300, $
                   'Aperture ' + string(k(i)+factor, $
                                        format='(i3.3)'), color='red'
             endfor
             
          endif
          
       endif else begin
          objsub = specflux
          objsubivar = specfluxivar
          pixmask=0.0*objsub
       endelse
       delvarx, specflux
       delvarx, specfluxivar
       if ccdnum eq 1 then begin
          
          flux1 = objsub
          fluxivar1 = objsubivar
          airset1 = airset    
          pixmask1 = pixmask
          plugmap1=plugmap
          
          
       endif else if ccdnum eq 2 then begin
          
          flux2 = objsub
          fluxivar2 = objsubivar
          airset2 = airset
          pixmask2 = pixmask
          plugmap2=plugmap
          
          ncol = n_elements(flux2(*,0))
          specflux = dblarr(ncol, 300)
          specfluxivar = dblarr(ncol, 300)
          specflux(*,0:149) = flux1
          specflux(*,150:299) = flux2
          specfluxivar(*,0:149) = fluxivar1
          specfluxivar(*,150:299) = fluxivar2
          
          pixmask = dblarr(ncol,300)
          pixmask(*,0:149) = pixmask1
          pixmask(*,150:299) = pixmask2
          
          coord1 = airset1.coeff
          coord2 = airset2.coeff
          coord = dblarr(n_elements(coord1[*,0]), 300)
          coord[*,0:149] = coord1
          coord[*,150:299] = coord2
          
          airset = create_struct('func', airset1.func, $
                                 'xmin', airset1.xmin, $
                                 'xmax', airset1.xmax, $
                                 'coeff', coord)
          
          
          temp = {expid:filelist(j(0))}
          temp = replicate(temp, 300)
          temp1 = [plugmap1, plugmap2]
          plugmapout = struct_addtags(temp, temp1)

                    
          delvarx, plugmap1
          delvarx, plugmap2
          delvarx, flux1
          delvarx, flux2
          delvarx, fluxivar1
          delvarx, fluxivar2
          delvarx, pixmask1
          delvarx, pixmask2
          
          
          if keyword_set(skysubtract) and (nother le 3 or nccd le 3) $
             then begin
             
             delvarx, xx 
             traceset2xy, airset, xx,lam
             hs_scale, lam, specflux, sigma=3.0, scale=scale
             
             divideflat, specflux, scale, invvar=specfluxivar, minval=0.1
             k = where(strmatch(plugmapout.objtype, 'SKY*'))
             
             
             traceset2xy, airset, xx,lam
             xx = indgen(4608,n_elements(k))
          for i = 0,n_elements(k)-1 do begin & xx(*,i) = xx(*,0) & endfor
             xy2traceset, xx, lam[*,k], tempset, ncoeff=6
             delvarx, xx
             
             skystuct = hs_skysubtract(specflux[*,k],$
                                       specfluxivar[*,k], plugmapout[k], $
                                       tempset,   objsub1, objsubivar1, $
                                       newmask=newmask)
             
             skyval = findgen(n_elements(k))
             
             for i = 0, n_elements(k) - 1 do begin
                l = where(abs(lam[*,k(i)]-5000) eq min(abs(lam[*,k(i)]-5000)))
                m = where(abs(lam[*,k(i)]-6000) eq min(abs(lam[*,k(i)]-6000)))
                djs_iterstat, objsub1[l:m,i], median=median
                skyval(i) = median
             endfor
             
             djs_iterstat, skyval, mean = mean, sigma=sigma
             bad = where(skyval gt mean+2.5*sigma)
             if bad(0) ne -1 then plugmapout[k(bad)].objtype='REJECT-OBJECT-56'
             
             skyval = findgen(n_elements(k))
             for i = 0, n_elements(k) - 1 do begin
                l = where(abs(lam[*,k(i)]-6000) eq min(abs(lam[*,k(i)]-6000)))
                m = where(abs(lam[*,k(i)]-7000) eq min(abs(lam[*,k(i)]-7000)))
                djs_iterstat, objsub1[l:m,i], median=median
                skyval(i) = median
             endfor
             
             djs_iterstat, skyval, mean = mean, sigma=sigma
             bad = where(skyval gt mean+2.5*sigma)
             If bad(0) ne -1 then plugmapout[k(bad)].objtype='REJECT-OBJECT67'
             skyval = findgen(n_elements(k))
             for i = 0, n_elements(k) - 1 do begin
                l = where(abs(lam[*,k(i)]-7000) eq min(abs(lam[*,k(i)]-7000)))
                m = where(abs(lam[*,k(i)]-8000) eq min(abs(lam[*,k(i)]-8000)))
                djs_iterstat, objsub1[l:m,i], median=median
                skyval(i) = median
             endfor
             
             djs_iterstat, skyval, mean = mean, sigma=sigma
             bad = where(skyval gt mean+2.5*sigma)
             if bad(0) ne -1 then plugmapout[k(bad)].objtype='REJECT-OBJECT67'
             
             splog, 'Mean =', mean
             splog, 'Sigma = ', sigma
             splog, 'Doing the full CCD skysubtract'
             
             skystuct = hs_skysubtract(specflux, specfluxivar, $
                                       plugmapout, airset,$
                                       objsub1, objsubivar1, $
                                       newmask=newmask1, sset=sset) 
             skysingle = 'true'
             pixmask = newmask1
             
             if keyword_set(qaplot) then begin
                split = strsplit(filelist(j), '/,.fits', /extract)
                traceset2xy,airset, xx, temp
                temp1 = temp
                for ifib = 0, n_elements(temp(0,*)) -1 do begin
                   temp1(*,ifib) = findgen(n_elements(temp(*,0)))
                endfor
                temp = alog10(temp)
                xy2traceset, temp1, temp, tempset, ncoeff=5
                iskies = where(strmatch(plugmapout.objtype, 'SKY'))
                
                hs_qaplot_skysub, specflux, specfluxivar, objsub, $
                   objsubivar, tempset, iskies
                
                
             endif
             
             specflux = objsub1
             specfluxivar = objsubivar1
          endif
       endif
       delvarx, xx
       
       traceset2xy, airset1, xx, lam1
       
       xx = dblarr(n_elements(lam1(*,0)), 2*n_elements(lam1(0,*)))
       for i = 0, n_elements(xx(0,*)) -1 do begin
          xx(*,i) = lindgen(n_elements(lam1(*,0)))
       endfor
       
       traceset2xy, airset, xx, lam
       
       if j eq 0 and ccdnum eq 1 then begin
          ;;Setup the output arrays that we need
          lamout = dblarr( n_elements(filelist), $
                           n_elements(lam(*,0)), 300)
          fluxout = dblarr( n_elements(filelist), $
                            n_elements(lam(*,0)), 300)
          ivarout=  dblarr( n_elements(filelist), $
                            n_elements(lam(*,0)), 300)
          pixmaskout=  dblarr( n_elements(filelist), $
                               n_elements(lam(*,0)), 300)
          headout = strarr(n_elements(filelist),$
                           n_elements(objhdr))
          fullfitout=  dblarr( n_elements(filelist), $
                               n_elements(lam(*,0)), 300)
       endif
    endfor
    lamout[j,*, *] = lam
    fluxout[j,*, *] = specflux
    ivarout[j,*, *] = specfluxivar
    pixmaskout[j,*,*] = pixmask
    headout[j,*] = objhdr
    
    if keyword_set(writeall) then begin
       
       file1 = strsplit(filelist(j),'/',/extract)
       file1=file1(n_elements(file1)-1)
       file1=cwd()+file1
       junk = strsplit(file1, ' ', length=count)
       junk = strsplit(cwd(), ' ', length=pathcount)
       outfile = outdir + strmid(file1, pathcount, $
                                 count-pathcount-5) $
          + '-'+ rerun + '.fits'
       
       ;;if the outname keyword is set, this will make a new outname
       if keyword_set(outname) then begin
          junk = strsplit(outname, '/.', /extract)
          if n_elements(junk) gt 1 then begin
             outroot = junk(n_elements(junk)-2)
          endif
          
          junk = strsplit(filelist(j), '.', /extract)
          outnum = junk(n_elements(junk)-2)
          outfile = outroot + '-' + outnum + '.fits'
          
          ;;Now check to see if 'spHect-' has been appended
          if strmatch(strmid(outfile,0,7), 'spHect-') then begin
             junk = strsplit(outfile, ' ', length=outlength)
             outfile = strmid(outfile, 7, outlength-7)
             
          endif
          outfile = 'spObs-' + outfile
          outfile = outdir + outfile
          
       endif
       
       splog, 'Writing single obs file ' + outfile
       
       sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
       
       mwrfits, float(transpose(transpose(lamout[j,*,*]))),$
          outfile, objhdr, /create
       mwrfits, float(transpose(transpose(fluxout[j,*,*]))), outfile
       mwrfits, float(transpose(transpose(ivarout[j,*,*]))), outfile
       mwrfits, float(transpose(transpose(pixmaskout[j,*,*]))), outfile
       mwrfits, float(transpose(transpose(pixmaskout[j,*,*]))), outfile
       mwrfits, plugmapout, outfile
       
       
       if keyword_set(skysubtract) then begin
          if skymulti eq 'true' then begin
             mwrfits, sset1, outfile
             mwrfits, scale1, outfile
             mwrfits, sset2, outfile
             mwrfits, scale2, outfile
          endif
          if skysingle eq 'true' then begin
             mwrfits, sset, outfile
             mwrfits, scale, outfile
             mwrfits,  0.0 , outfile
             mwrfits,  0.0 , outfile
             
          endif
       endif
    endif
    
    if keyword_set(qaplot) then begin
       
       test = strsplit(outname, '/', /extract)
       test = test(n_elements(test)-1)
       test= test(0)
       junk = strsplit(test, ' ', length=ct)
       test = strmid(test, 0, ct-5)
       
       if strmid(test,0,7) eq 'spHect-' then test = strmid(test, 7, 3)
       num = j
       hs_qaplot_snr, flux=specflux, lam=lam,$
          plugmap=plugmapout, rerun=rerun, plateid=test, $
          num = num, filename=filelist(j)
       
    endif
    
    if keyword_set(quickplot) then begin
       hs_quickplot, flux=specflux, lam=lam, $
          plugmap=plugmapout, rerun=rerun
    endif
    delvarx, lam
    delvarx, specflux
    delvarx, specfluxivar
    delvarx, pixmask
 endfor
 if keyword_set(tellcorr) then begin
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
       
    endfor
    
 endif
 
 if keyword_set(fluxcorr) then begin
    
    flux_factor = hs_fluxcorr(lamout, fluxout,$
                              ivarout, plugmapout, $
                              pixmaskout, outname = filelist, $
                              plateid=plateid, $
                              rerun=rerun,/psplot, checkaverage=checkaverage,$
                              stand=stand, header=objhdr, cut=cut, $
                              superfile=superfile, useaverage=useaverage) 
    
    if keyword_set(cut)  then begin
       npix = n_elements(lamout(0,*,0))-1
       lamout = lamout[*,cut:npix,*]
       fluxout = fluxout[*,cut:npix,*]
       ivarout = ivarout[*,cut:npix,*]
       pixmaskout = pixmaskout[*,cut:npix,*]
    endif
            
    fcalfile = strarr(n_elements(filelist))
    sphoto_err = dindgen(n_elements(filelist))
    sn=sphoto_err
    
    for im = 0, n_elements(filelist) -1 do begin
       outname1 = strsplit(filelist(im), '/,.fits', /extract)
       ct = n_elements(outname1)
       outname1= outname1(ct-1)
       
       fcalfile(im) = 'reduction/' + rerun +'/hsFluxcalib-' + $
          outname1 + '-' + plateid + '-' + $
          rerun +'.fits'
       
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
       
       corrfiles = 'reduction/' + rerun +'/spFluxcorr-' + junk1
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
    
    if keyword_set(combine) AND n_elements(filelist) GT 1 then begin
       
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
       
    endif
  
    
      
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
          galb, 1, /degree
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
;    finalivar = ivarout
    ivarout = finalivar
    pixmaskout = finalandmask
    ormaskout = finalormask
    
 endif
 
 ;;This should handle the combining of non-fluxed data 
 ;;in a way that will not kill the data
 
 if not keyword_set(fluxcorr) and keyword_set(combine) then begin
    
    if not keyword_set(noftweak) then begin
       
       sn = fltarr(n_elements(filelist))
       for im = 0, n_elements(filelist) -1 do begin
          hs_snr,  transpose(transpose(lamout[im,*,*])), $
             transpose(transpose(fluxout[im,*,*])), plugmapout, $
             sn2=sn21, /noplot, /nomag
          
          sn(im) =sn21  
       endfor
       
       expid = filelist
       sn_exp = sn
       maxval = max(sn_exp/median(sn_exp), iframe_best)
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
          
          corrfiles = 'reduction/' + rerun +'/spFluxcorr-' + junk1
          if not keyword_set(noftweak) then begin
             
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
                   (corrimg LE 0.1) * pixelmask_bits('BADFLUXFACTOR');
                
                fluxout[iexp,*,*] = tempflux
                lamout[iexp,*,*] = lam[*,indx]
                ivarout[iexp,*,*] = tempivar
                pixmaskout[iexp,*,*] = tempmask
                
             endfor
             
          endif else begin
             for iexp=0, n_elements(filelist) -1 do begin
                
                indx = where(plugmap.expid eq expid[iexp])
                tempflux = flux[*,indx]
                tempivar = ivar[*,indx]
                tempmask = mask[*, indx]
                fluxout[iexp,*,*] = tempflux
                lamout[iexp,*,*] = lam[*,indx]
                ivarout[iexp,*,*] = tempivar
                pixmaskout[iexp,*,*] = tempmask
                   
             endfor
          endelse
       endif   
    ENDIF
    
    
    ;;Now combine each of the fibers into a composite spectrum
    newlam = transpose(transpose(lamout[0,*,*]))
    newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
    newfluxivar = transpose(transpose(ivarout[0,*,*])) * 0.0
    newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
    ormask = newmask
    
    
    
    
    for i = 0, 299 do begin
       
       nfib = 300
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
 ENDIF
 
 
 if keyword_set(writeout) AND n_elements(filelist) GT 1 then begin
    
    if NOT keyword_set(outname) then $
       outname = plugfile + '-' + rerun + '.fits'
        
    splog, 'Writing coadded sphect file ' + outname
    
    
    sxaddpar, objhdr, 'NAXIS2', n_elements(lamout(*,0)), after='NAXIS'
    sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
    mwrfits, float(lamout), outname, /create, objhdr
    mwrfits, float(fluxout), outname
    mwrfits, float(ivarout), outname
    mwrfits, float(pixmaskout), outname
    
    
    if keyword_set(combine) then mwrfits, float(ormaskout), outname
    if NOT keyword_set(combine) then mwrfits, float(pixmaskout), outname
    mwrfits, plugmapout, outname
    if keyword_set(skysubtract) then begin
       if skymulti eq 'true' then begin
          mwrfits, sset1, outname
          mwrfits, scale1, outname
          mwrfits, sset2, outname
          mwrfits, scale2, outname
       endif
       if skysingle eq 'true' then begin
          mwrfits, sset, outname
          mwrfits, scale, outname
          mwrfits,  0.0 , outname
          mwrfits,  0.0 , outname
          
       endif
    endif
    
 endif
 
 
 if keyword_set(qaplot) then dfpsclose
 
 scale = 0.0
 scale1 = 0.0
 scale2 = 0.0
 skystuct = 0.0
 skystuct1 = 0.0
 skystuct2 = 0.0
 
 
 
 return
END



