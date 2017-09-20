;+
; NAME:
;   hs_extract
;
; PURPOSE:
;  Main reduction workhorse
;
; CALLING SEQUENCE:
;       hs_extract, filelist, plugfile=plugfile, skysubtract=skysubtract, $
;                fiberflat=fiberflat,  tweaksky=tweaksky, $
;                plugcat=plugcat, combine=combine, $
;                lamout=lamout, fluxout=fluxout, fluxcorr=fluxcorr,$
;                tellcorr=tellcorr, uberextract=uberextract, $
;                ivarout=ivarout, rerun=rerun, outname=outname, $
;                plateid=plateid, sky2d=sky2d, $
;                qaplot=qaplot, quicklook=quicklook, quickplot=quickplot, $
;                stand=stand,  checkaverage=checkaverage, $
;                docosmic=docosmic, doplot=doplot, customskies=customskies, $
;                superfile=superfile, useaverage=useaverage, $
;                noftweak=noftweak, nobias=nobias, dodark=dodark, $
;                nosky=nosky, debug=debug, doskyscale=doskyscale, $
;                mmtcam=mmtcam, skyscale=skyscale, dostand=dostand, $
;                writeall=writeall, doredleak=doredleak,writeout=writeout,$
;                oldmap=oldmap, dummy=dummy, dofringes=dofringes
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
;   tellcorr    - Set this to telluric correct the data. Works for
;                 fluxed or unfluxed data, but so far only on 270-line
;                 data (supercedes old code that only worked on fluxed 
;                 data, but may be problematic if telluric correction
;                 is needed for fluxed, 600-line data)
;   uberextract - Set this to turn on tweaksky, fiberflat, skysubtract
;                 tellcorr, doredleak, and fluxcorr. Obivously,
;                 you don't want this if you don't have standards
;   rerun       - This should be set to the rerun number you would like to 
;                 for the reduction. If not set rerun=0000
;   outname     - Root outname you would like to use for the data.  If
;                 this is not set, then outname='test'
;   qaplot      - Setting this results in the quality assurance plots to be
;                 be created
;   superfile   - If this is set, the file $HSRED_DIR/etc/superstand.fits
;                 instead of typing the stars
;   docosmic    - If docosmic is set then the cosmic ray rejection is turned 
;                 on. docosmic=2 uses HCOSMIC if > 1 image present
;                 (You want this)
;   checkaverage- This will check the current flux soln against an archival 
;                 solution and will truncate the data if the counts of the 
;                 standard stars are fewer than 100 
;   dostand    - This will provide the standard flags that you need
;                 (skyssubtract, fiberflat, tweaksky, combine, noftweak
;                 tellcor, doredleak, writeout, writeall)
;   dofringes - Estimates and subtracts red fringing contribution,
;               before applying flat field. Testing has shown this 
;               makes very little difference, but option is retained
;   noftweak - enables new, simpler coadding of exposures, without any
;              scaling of fluxes between exposures. Cannot be set with /fluxcorr
;   doredleak - set to remove contaminating red flux from robot
;               positioner LEDs, longward of ~8500A. 
;   mmtcam   - placeholder to eventually remove the additional red
;              leak introduced while the mmtcam camera is turned
;              on. Not yet functional
;   dummy    - useful when trying to force the pipeline to treat a
;              calibration or other non-standard file as if it were an
;              actual target observation
;   sky2d    - turn on enhanced sky subtraction where the 'supersky'
;              model is allowed to vary smoothly as a function of 
;              fiber number. If there are too few sky fibers, in most
;              cases the code with revert automatically to normal 1D sky
;              subtraction; unset this keyword to force that behavior
;   customskies - Set when configuration contains custom-placed sky fibers
;        Normal behavior is to treat all fibers named 'SKY' or 'sky*' as a
;        sky fiber. With this flag set, also includes skies named '*sky' in
;        the map files. Note that it is the name that matters, not the 0/1/2
;        designation of fiber type in the map files
;
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
;    New telluric correction code supercedes old code that only worked
;    on fluxed data, but may be problematic if telluric correction
;    is needed for fluxed, 600-line data
;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA         
;   June  2014 - Many updates for HSRED v1.2, S. Moran, SAO TDC
;                
;-
;------------------------------------------------------------------------------
PRO hs_extract, filelist, plugfile=plugfile, skysubtract=skysubtract, $
                fiberflat=fiberflat,  tweaksky=tweaksky, $
                plugcat=plugcat, combine=combine, $
                lamout=lamout, fluxout=fluxout, fluxcorr=fluxcorr,$
                tellcorr=tellcorr, uberextract=uberextract, $
                ivarout=ivarout, rerun=rerun, outname=outname, $
                plateid=plateid, sky2d=sky2d, $
                qaplot=qaplot, quicklook=quicklook, quickplot=quickplot, $
                stand=stand,  checkaverage=checkaverage, $
                docosmic=docosmic, doplot=doplot, customskies=customskies, $
                superfile=superfile, useaverage=useaverage, $
                noftweak=noftweak, nobias=nobias, dodark=dodark, $
                nosky=nosky, debug=debug, doskyscale=doskyscale, $
                mmtcam=mmtcam, skyscale=skyscale, dostand=dostand, $
                writeall=writeall, doredleak=doredleak,writeout=writeout,$
                oldmap=oldmap, dummy=dummy, dofringes=dofringes
  
  
  
  if NOT keyword_set(plugcat) then begin
     plugfile = strsplit(filelist(0), '/', /extract)
     plugfile = plugfile(n_elements(plugfile)-1)
     junk  = strsplit(plugfile, '.',  /extract)
     junk1  = strsplit(plugfile, ' ', length=ct)
     if junk[n_elements(junk)-1] eq 'gz' then begin
        plugfile = strmid(plugfile, 0, ct-8)
     endif else begin
        plugfile = strmid(plugfile,0,ct-5)
     end
     
  endif
  
  if keyword_set(quicklook) then begin   
     hs_extract, filelist, plugfile=plugfile, /tweaksky, $
        /fiberflat, /skysubtract, plugcat=plugcat, /writeall, $
        rerun=rerun, stand=stand, checkaverage=checkaverage, $
        docosmic=docosmic, /writeout, outname=outname, qaplot=qaplot, $
        superfile=superfile, useaverage=useaverage, nobias=nobias, $
        dodark=dodark, nosky=nosky, debug=debug, sky2d=sky2d, customskies=customskies, $
        doskyscale=doskyscale, skyscale=skyscale, oldmap=oldmap, dummy=dummy, dofringes=dofringes, doredleak=doredleak
  endif
     
     
        
  if keyword_set(uberextract) then begin
     hs_extract,  filelist, plugfile=plugfile, /skysubtract, $
        /fiberflat,  /tweaksky, $
        plugcat=plugcat, /combine, /writeout,$
        lamout=lamout, rerun=rerun,  fluxout=fluxout, $
        /fluxcorr, /tellcorr, ivarout=ivarout, $
        outname=outname, /writeall, /sky2d, customskies=customskies, $
        plateid=plateid, qaplot=qaplot, quickplot=quickplot, $
        stand=stand, checkaverage=checkaverage, docosmic=2, $
        superfile=superfile, useaverage=useaverage, noftweak=noftweak, $
        nobias=nobias, dodark=dodark, nosky=nosky, debug=debug, $
        doskyscale=doskyscale, skyscale=skyscale, oldmap=oldmap, $
        dummy=dummy, dofringes=dofringes, /doredleak, mmtcam=mmtcam
     return
  endif
  
  if keyword_set(dostand) then begin
     hs_extract, filelist, plugfile=plugfile, /skysubtract, $
        /fiberflat, /tweaksky, plugcat=plugcat, /combine, $
        /writeout, /writeall, fluxout=fluxout, lamout=lamout, rerun=rerun, $
        ivarout=ivarout, outname=outname, plateid=plateid, /sky2d, $
        qaplot=qaplot, stand=stand, docosmic=2, superfile=superfile, $
        customskies=customskies, /tellcorr, mmtcam=mmtcam, $
        /noftweak, oldmap=oldmap, dummy=dummy, dofringes=dofringes, /doredleak
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

;if we are doing CR_rejection, need to call hs_proc to load all the images
;up front and pass them all to HCOSMIC
if (docosmic eq 2) and (n_elements(filelist) gt 1) then begin
  exptimes=fltarr(n_elements(filelist))
  for p=0, n_elements(filelist)-1 do begin
    if p eq 0 then crfiles1='reduction/'+rerun+'/cr1-'+strn(p)+'.fits' else crfiles1=[crfiles1,'reduction/'+rerun+'/cr1-'+strn(p)+'.fits']
    if p eq 0 then crfiles2='reduction/'+rerun+'/cr2-'+strn(p)+'.fits' else crfiles2=[crfiles2,'reduction/'+rerun+'/cr2-'+strn(p)+'.fits']
    hs_proc, filelist[p], 1, spec1, specivar1, pixelmask=pixmask, $
           rerun=rerun, root=root, header=objhdr1, docosmic=0,$
           nobias=nobias
    exptimes[p]=sxpar(objhdr1, 'EXPTIME')
    hs_proc, filelist[p], 2, spec2, specivar2, pixelmask=pixmask, $
           rerun=rerun, root=root, header=objhdr2, docosmic=0,$
           nobias=nobias
    if p eq 0 then begin 
       dimens=size(spec1)
       nx=dimens[1]
       ny=dimens[2]
       nim=n_elements(filelist)
       images1=fltarr(nx,ny,nim)
       ivar1=fltarr(nx,ny,nim)
       images2=fltarr(nx,ny,nim)
       ivar2=fltarr(nx,ny,nim)
     endif     
     images1[*,*,p]=spec1
     ivar1[*,*,p]=specivar1
     images2[*,*,p]=spec2
     ivar2[*,*,p]=specivar2
  endfor
  hcosmic, images=images1, ivar=ivar1
  hcosmic, images=images2, ivar=ivar2
endif

;now proceed with the main loop
  for j = 0, n_elements(filelist) -1 do begin
     for ccdnum = 1, 2 do begin
        hs_proc, filelist(j), ccdnum, spec, specivar, pixelmask=pixmask, $
           rerun=rerun, root=root, header=objhdr, nobias=nobias, $
           docosmic=((docosmic eq 1) or ((n_elements(filelist) eq 1) and (docosmic eq 2)))
          mjd_in = double(sxpar(objhdr,'MJD'))
          jd = 2400000.5D + mjd_in
        if (docosmic eq 2) and (n_elements(filelist) gt 1) then begin
          if ccdnum eq 1 then begin
            spec=images1[*,*,j]
            specivar=ivar1[*,*,j]
         endif else begin
            spec=images2[*,*,j]
            specivar=ivar2[*,*,j]
          endelse
          sxaddpar, objhdr, 'COSMIC', 'CR-rejected using HCOSMIC'
        endif
        ;trim blanks from fits header
        nonblank=where(strtrim(objhdr,2) ne '', blcnt)
        if blcnt ge 1 then objhdr=objhdr[nonblank]
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
           sxaddpar, objhdr, 'DARKCOR', 'Dark correction applied'
        endif
        
;        sigma= 2.0
        sigma=1.5
        proftype = 2.0
        highrej = 30
        lowrej = 30
        npoly = 4
        wfixed = [1, 1, 1 ]
        
        
        hs_readtraceset, ccdnum , xsol=xsol, rerun=rerun
        
ansimage=1
ymodel=1
        extractmask=fix(pixmask*0+1)
        ivarzero=where(specivar eq 0.0, ivarcnt)       
        if ivarcnt ge 1 then extractmask[ivarzero]=0
        hs_extract_image, spec, specivar, xsol, sigma,$
           specflux, specfluxivar, mask=extractmask, $
           proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
           npoly=npoly, relative=1, ymodel=ymodel, ansimage=ansimage

        bad=where(specfluxivar eq 0.0, badcnt)
        if badcnt ge 1 then begin
             pmask=(specfluxivar eq 0.0)
             specflux=djs_maskinterp(specflux,pmask,iaxis=0, /spline)
        endif
        sxaddpar, objhdr, 'EXTRACT', 'Spectra traced and extracted using hs_extract_image' 
                
        delvarx, spec
        delvarx, specivar
        
        if file_test('calibration/'+rerun+'/'+'wset.fits') EQ 0L then begin
           hs_findwave, rerun=rerun, root=root
        endif

        if file_test('calibration/'+rerun+'/'+'dispset.fits') EQ 0L then begin
           hs_findarcdisper, rerun=rerun, /write
        endif
        
        wset = mrdfits('calibration/'+rerun+'/'+'wset.fits', ccdnum, /silent)
        dispset = mrdfits('calibration/'+rerun+'/'+'dispset.fits', ccdnum, /silent)
        split = strsplit(filelist(j), '/,.fits', /extract)
        
        if keyword_set(plugcat) then begin
           hs_cattoplug, plugfile, plugmap, ccdnum=ccdnum, oldmap=oldmap
           if ccdnum eq 1 then other = 2
           if ccdnum eq 2 then other = 1
           hs_cattoplug, plugfile, plugcheck, ccdnum= other, oldmap=oldmap
        endif else begin
           if ccdnum eq 1 then other = 2
           if ccdnum eq 2 then other = 1
           hs_maptoplug, plugfile + '_map', plugmap, ccdnum=ccdnum, dummy=dummy
           hs_maptoplug, plugfile + '_map', plugcheck ,ccdnum=other, dummy=dummy
           ;fix for old data where map files do not have ra, dec for sky/unused
           empty=where((plugmap.ra eq 0.0) and (plugmap.dec eq 0.0),empty_cnt)
           if empty_cnt ge 1 then begin
                plugmap[empty].ra=hms2dec(sxpar(objhdr,'RA'))*15.
                plugmap[empty].dec=hms2dec(sxpar(objhdr,'DEC'))
           endif
           empty=where((plugcheck.ra eq 0.0) and (plugcheck.dec eq 0.0),empty_cnt)
           if empty_cnt ge 1 then begin
                plugcheck[empty].ra=hms2dec(sxpar(objhdr,'RA'))*15.
                plugcheck[empty].dec=hms2dec(sxpar(objhdr,'DEC'))
           endif
           if keyword_set(customskies) then begin
                custom= where(strmatch(plugmap.objtype, '*SKY*'), cu_cnt)
                if cu_cnt ge 1 then plugmap[custom].objtype='SKY'
                custom = where(strmatch(plugcheck.objtype, '*SKY*'), cu_cnt)
                if cu_cnt ge 1 then plugcheck[custom].objtype='SKY'
            endif else begin
                custom= where(strmatch(plugmap.objtype, 'SKY*'), cu_cnt)
                if cu_cnt ge 1 then plugmap[custom].objtype='SKY'
                custom = where(strmatch(plugcheck.objtype, 'SKY*'), cu_cnt)
                if cu_cnt ge 1 then plugcheck[custom].objtype='SKY'
            endelse
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
           ;first do a fringing correction if desired
           if keyword_set(dofringes) then begin
             if file_test('calibration/'+rerun+'/'+'fringes_pattern.fits')$
                EQ 0L then begin 
                hs_fringing, rerun
             endif
             fringes = mrdfits('calibration/'+rerun+'/'+ $
                           'fringes_pattern.fits',ccdnum-1, /silent)
            
             tmpsmooth=specflux*0.0
             nsmooth=n_elements(specflux[*,0])/50
             for q=0, 149 do tmpsmooth[*,q]=djs_median(specflux[*,q],width=nsmooth, boundary='reflect') 
             specflux=specflux-fringes*tmpsmooth
             for q=0, 149 do tmpsmooth[*,q]=djs_median(fflat[*,q],width=nsmooth, boundary='reflect') 
             fflat=fflat-fringes*tmpsmooth ;*fflat             
              sxaddpar, objhdr, 'FRINGES', 'Fringe pattern removed from flat field and science spectra'
           endif
           divideflat, specflux, fflat, invvar=specivar
           sxaddpar, objhdr, 'DOMEFLAT', 'Dome flat applied'

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
              sxaddpar, objhdr, 'FIBERFLAT', 'Fiber flat applied'
              endif else begin
              ;;make the fiberflats from domeflat instead of skyflat
              if file_test('calibration/'+rerun+'/'+'fiberflat.fits')$
                 EQ 0L then begin 
                 hs_makefiberflat, rerun=rerun, root=root, usesky=0
              endif
              
              fflat = mrdfits('calibration/'+rerun+'/'+$
                              'fiberflat.fits',ccdnum-1, /silent)
              if (size(fflat))[0] eq 2 then divideflat, specflux, fflat, invvar=specivar
              sxaddpar, objhdr, 'FIBERFLAT', 'Fiber flat applied, based on dome not sky'
           endelse
           
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
        

;;---------------------------------
;; tweak the wavelength solution based on sky line positions
;;        
        if keyword_set(tweaksky) then begin
           iskies = where(strmatch(plugmap.objtype, 'SK*'))
           
           hs_tweaksky, specflux, specfluxivar, wset, ccdnum=ccdnum, $
              airset=airset, xsky=xsky, rerun=rerun, $
              filename=filelist(j), skywaves=xwaves, $
              arcshift=arcshift
           sxaddpar, objhdr, 'SKYTWEAK', 'Wavelengths tweaked based on sky lines'

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
              delvarx, xx
           endif
           
        endif else begin
           airset = wset
           
        endelse

 ;;-------------------------------------
 ;;  correct for the red leak here, before sky subtraction, and before rejecting bad skies
 ;; 
        if keyword_set(doredleak) then begin
           specflux_orig=specflux
           specfluxivar_orig=specfluxivar
           hs_remove_redleak, airset, plugmap, sxpar(objhdr,'EXPTIME'), specflux=specflux, invvar=specfluxivar, mmtcam=mmtcam
           sxaddpar, objhdr, 'REDLEAK', 'Red leak removed using template for intensity vs fiber position'
        endif

;;--------------------------------
;;begin sky subtraction part of code
;;
        ;do a first pass at rejecting contaminated sky fibers
        nccd = n_elements(where(strmid(plugmap.objtype,0,3) eq 'SKY'))
        nother = n_elements(where(strmid(plugcheck.objtype,0,3) eq 'SKY'))
        delvarx, xx
        traceset2xy, airset, xx, lam11
        
        k = where(strmatch(plugmap.objtype, 'SKY*'))
        
        if k(0) ne -1 then begin
           
           ;;Check at 4000 Angstroms
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1 do begin
              kcheck(ifib) = (where(abs(lam11[*,ifib]-4000.0) $
                                   eq min(abs(lam11[*,ifib]-4000.0))))[0]
           endfor         
           
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 3*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT4'
           
           
           ;;Check at 5577 (this will keep the bspline from 
           ;;dying if there is a zero fiber 
           
           kcheck = indgen(n_elements(k))
           
           for ifib = 0,  n_elements(k)-1 do begin
              kcheck(ifib) = (where(abs(lam11[*,ifib]-5577.0) $
                                   eq min(abs(lam11[*,ifib]-5577.0))))[0]
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 3*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT5'
           
           ;;Check at 5800 for continuum sources
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1  do begin
              kcheck(ifib) = (where(abs(lam11[*,ifib]-5800.0) $
                                   eq min(abs(lam11[*,ifib]-5800.0))))[0]
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where(specflux[kcheck,k] gt mean + 3*sigma)
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT6'
           
           ;;Check at 8725 and 8975 for the Red monster-have to check
           ;;above and below the mean now, since we have already applied our
           ;;red leak correction
           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1 do begin
              kcheck(ifib) = (where(abs(lam11[*,ifib]-8725.0) $
                                   eq min(abs(lam11[*,ifib]-8725.0))))[0]
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where((specflux[kcheck,k] gt mean + 3*sigma) or (specflux[kcheck,k] lt mean - 3*sigma))
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT8'

;           w = where((specflux[kcheck,k] gt mean + 1.5*sigma) or (specflux[kcheck,k] lt mean - 1.5*sigma))
;           if w(0) ne -1 then begin
;                   for zz=0, n_elements(w)-1 do begin
;                      vv=where(lam11[*,k[w[zz]]] ge 8750.)
;                     specfluxivar[zz,k[w[zz]]]=0.0 
;                   endfor
;           endif

           kcheck = indgen(n_elements(k))
           for ifib = 0, n_elements(k)-1  do begin
              kcheck(ifib) = (where(abs(lam11[*,ifib]-9150.0) eq min(abs(lam11[*,ifib]-9150.0))))[0]
           endfor     
           djs_iterstat, specflux[kcheck,k], mean=mean, sigma=sigma
           w = where((specflux[kcheck,k] gt mean + 3*sigma) or (specflux[kcheck,k] lt mean - 3*sigma))
           if w(0) ne -1 then plugmap[k(w)].objtype = 'REJECT9'
;           w = where((specflux[kcheck,k] gt mean + 1.5*sigma) or (specflux[kcheck,k] lt mean - 1.5*sigma))
;           if w(0) ne -1 then begin
;                  for zz=0, n_elements(w)-1 do begin
;                      vv=where(lam11[*,k[w[zz]]] ge 8750.)
;                     specfluxivar[zz,k[w[zz]]]=0.0 
;                   endfor
;           endif
        endif

        nccd = n_elements(where(strmid(plugmap.objtype,0,3) eq 'SKY'))

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
            
            ;;
            ;;Scale fiber transmission by fitting sky lines in each fiber to the average supersky
            ;;--------------------------------------------
            delvarx, xx
            traceset2xy, airset, xx,lam
            k = where(strmatch(plugmap.objtype, 'SKY*'))
            xy2traceset, xx[*,k], lam[*,k], tempset, ncoeff=6
            if keyword_set(dummy) then jd=0
            skystuct = hs_skysubtract(specflux[*,k],$
                                    specfluxivar[*,k], plugmap[k], $
                                    tempset,   objsub1, objsubivar1, fullfit=fullfit,$
                                    newmask=newmask, jd=jd, sset=sset)


            hs_skyfit, lam, specflux, scale=scale, skyset=sset
            divideflat, specflux, scale, invvar=specfluxivar, minval=0.1
            ;but if we have removed red leak, lets
            ;apply fiber scaling, and then re-do red leak  
            ;DOESNT SEEM TO IMPROVE THINGS
            ;if keyword_set(doredleak) then begin
            ;    specflux=specflux_orig
            ;    specfluxivar=specfluxivar_orig
            ;    divideflat, specflux, scale, invvar=specfluxivar, minval=0.1
            ;    hs_remove_redleak, airset, plugmap, sxpar(objhdr,'EXPTIME'), specflux=specflux, invvar=specfluxivar, mmtcam=mmtcam
            ;endif



            k = where(strmatch(plugmap.objtype, 'SKY*'))
            low=indgen(75) ;for top and bottom halves of a ccd
            hi=indgen(75)+75
            traceset2xy, airset, xx,lam
            xy2traceset, xx[*,low], lam[*,low], airset_low, ncoeff=6
            xy2traceset, xx[*,hi], lam[*,hi], airset_hi, ncoeff=6
            traceset2xy, dispset, xx,disp
            xy2traceset, xx[*,k], lam[*,k], tempset, ncoeff=6
            xy2traceset, xx[*,k], disp[*,k], tempdisp, ncoeff=3
            tempdisp = struct_addtags(tempdisp $
                                  , create_struct('MEDIAN', dispset.median[k]))
            delvarx, xx
            if keyword_set(dummy) then jd=0
            skystuct = hs_skysubtract(specflux[*,k],$
                                    specfluxivar[*,k], plugmap[k], $
                                    tempset,   objsub1, objsubivar1, fullfit=fullfit,$
                                    newmask=newmask, jd=jd, sset=sset1) ;, $
                                    ;dispset=tempdisp)

            skyval = findgen(n_elements(k))
            
            for i = 0, n_elements(k) - 1 do begin
               
               l = where(abs(lam[*,k(i)]-5000) eq min(abs(lam[*,k(i)]-5000)))
               m = where(abs(lam[*,k(i)]-6000) eq min(abs(lam[*,k(i)]-6000)))
               m=m[0]
               l=l[0]
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            
            bad = where((skyval gt mean+2.5*sigma) or (skyval lt mean-2.5*sigma))
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT-56'

            skyval = findgen(n_elements(k))
            for i = 0, n_elements(k) - 1 do begin
               
               l = where(abs(lam[*,k(i)]-6500) eq min(abs(lam[*,k(i)]-6500)))
               m = where(abs(lam[*,k(i)]-6600) eq min(abs(lam[*,k(i)]-6600)))
               m=m[0]
               l=l[0]
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            
            bad = where((skyval gt mean+2.5*sigma) or (skyval lt mean-2.5*sigma))
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT-HA'

            skyval = findgen(n_elements(k))
            
            for i = 0, n_elements(k) - 1 do begin
               l = where(abs(lam[*,k(i)]-6000) eq min(abs(lam[*,k(i)]-6000)))
               m = where(abs(lam[*,k(i)]-7000) eq min(abs(lam[*,k(i)]-7000)))
               m=m[0]
               l=l[0]
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            bad = where((skyval gt mean+2.5*sigma) or (skyval lt mean-2.5*sigma))
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT67'
            skyval = findgen(n_elements(k))
            
            for i = 0, n_elements(k) - 1 do begin
               l = where(abs(lam[*,k(i)]-7000) eq min(abs(lam[*,k(i)]-7000)))
               m = where(abs(lam[*,k(i)]-8000) eq min(abs(lam[*,k(i)]-8000)))
               m=m[0]
               l=l[0]
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor
            
            djs_iterstat, skyval, mean = mean, sigma=sigma
            bad = where((skyval gt mean+2.5*sigma) or (skyval lt mean-2.5*sigma))
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT78'

            for i = 0, n_elements(k) - 1 do begin
               l = where(abs(lam[*,k(i)]-8000) eq min(abs(lam[*,k(i)]-8000)))
               m = where(abs(lam[*,k(i)]-9000) eq min(abs(lam[*,k(i)]-9000)))
               m=m[0]
               l=l[0]
               djs_iterstat, objsub1[l:m,i], median=median
               skyval(i) = median
            endfor

            djs_iterstat, skyval, mean = mean, sigma=sigma
            bad = where((skyval gt mean+2.5*sigma) or (skyval lt mean-2.5*sigma))
            if bad(0) ne -1 then plugmap[k(bad)].objtype='REJECT-OBJECT89'
            bad = where((skyval gt mean+1.5*sigma) or (skyval lt mean-1.5*sigma))
            if bad[0] ne -1 then begin
                 for zz=0, n_elements(bad)-1 do begin
                     vv=where(lam[*,k[bad[zz]]] ge 8750.)
                     specfluxivar[zz,k[bad[zz]]]=0.0 
                 endfor
            endif

            nccd = n_elements(where(strmid(plugmap.objtype,0,3) eq 'SKY'))
            if ccdnum eq 1 then nccd1=nccd
            if ccdnum eq 1 then skysubq1=0
            if ccdnum eq 1 then skysubq2=0
            if ccdnum eq 2 then sky2d=sky2d_old
            if ccdnum eq 1 then  sky2d_old=sky2d
            if keyword_set(sky2d) and (nccd lt 6) then sky2d=0            


           if keyword_set(sky2d) then begin
              skystuct = hs_skysubtract(specflux, $
                                      specfluxivar, plugmap, $
                                      airset,   objsub1, objsubivar1, upper=5, lower=5, $
                                    newmask=newmask1, sset=sset, jd=jd, nord=3, npoly=2, dispset=dispset) 
              skyflux=specflux-objsub1
              sxaddpar, objhdr, 'SKYSUB', 'Sky spline-fit and subtracted, including 2nd dimen on fiber num'

              ;lets also revert back to 1d sky if one *half* of the ccd contains 3 or fewer
             ;this fixes an edge case that comes up from time to time, where the sky
             ;distribution is lopsided
             nccda = n_elements(where(strmid(plugmap[low].objtype,0,3) eq 'SKY'))
             nccdb = n_elements(where(strmid(plugmap[hi].objtype,0,3) eq 'SKY'))

             if (nccda le 3) then begin            
              k = where(strmatch(plugmap.objtype, 'SKY*'), kcnt)
              low_sky=indgen((k[5]>74)+1) ;make sure we have the  nearest 6 skies, at least
              delvarx, xx
              traceset2xy, airset, xx,lam
              xy2traceset, xx[*,low_sky], lam[*,low_sky], airset_low, ncoeff=6

              skystucta = hs_skysubtract(specflux[*,low_sky], $
                                      specfluxivar[*,low_sky], plugmap[low_sky], $
                                      airset_low,   objsub1a, objsubivar1a, $
                                    newmask=newmask1a, sset=sseta, jd=jd, nord=3, upper=4, lower=4) 
              skyflux[*,low]=specflux[*,low]-objsub1a[*,low]
              objsub1[*,low]=objsub1a[*,low]
              objsubivar1[*,low]=objsubivar1a[*,low]
              newmask1[*,low]=newmask1a[*,low]

              if ccdnum eq 1 then skysubq1=1
              if (ccdnum eq 2) and (skysubq1 eq 1) then sxaddpar, objhdr, 'SKYSUBQ1', 'Sky fit reverted to 1D for fibers 1-75'
              if ccdnum eq 2 then sxaddpar, objhdr, 'SKYSUBQ3', 'Sky fit reverted to 1D for fibers 151-225'
             endif else if (nccdb le 3) then begin

              k = where(strmatch(plugmap.objtype, 'SKY*'), kcnt)
              hi_sky=indgen(150)
              hi_sky=hi_sky[(k[(kcnt-6)]<75):*] ;make sure we have the  nearest 6 skies, at least
              delvarx, xx
              traceset2xy, airset, xx,lam
              xy2traceset, xx[*,hi_sky], lam[*,hi_sky], airset_hi, ncoeff=6
              skystuctb = hs_skysubtract(specflux[*,hi_sky], $
                                      specfluxivar[*,hi_sky], plugmap[hi_sky], $
                                      airset_hi,   objsub1b, objsubivar1b, $
                                    newmask=newmask1b, sset=ssetb, jd=jd, nord=3, upper=4, lower=4) 
              skyflux[*,hi]=specflux[*,hi]-objsub1b[*,(low+(n_elements(hi_sky)-n_elements(hi)))]
              objsub1[*,hi]=objsub1b[*,(low+(n_elements(hi_sky)-n_elements(hi)))]
              objsubivar1[*,hi]=objsubivar1b[*,(low+(n_elements(hi_sky)-n_elements(hi)))]
              newmask1[*,hi]=newmask1b[*,(low+(n_elements(hi_sky)-n_elements(hi)))]

              if ccdnum eq 1 then skysubq2=1
              if (ccdnum eq 2) and (skysubq2 eq 1) then sxaddpar, objhdr, 'SKYSUBQ2', 'Sky fit reverted to 1D for fibers 76-150'
              if ccdnum eq 2 then sxaddpar, objhdr, 'SKYSUBQ4', 'Sky fit reverted to 1D for fibers 226-300'

             endif

              if ccdnum eq 2 then begin
                 sxaddpar, objhdr, 'NSKIES1', 'Used '+strn(nccd1)+' good skies on CCD1'
                 sxaddpar, objhdr, 'NSKIES2', 'Used '+strn(nccd)+' good skies on CCD2'
              endif
           endif  else begin
              skystuct = hs_skysubtract(specflux, $
                                      specfluxivar, plugmap, $
                                      airset,   objsub1, objsubivar1, $
                                    newmask=newmask1, sset=sset, jd=jd, upper=4, lower=4) 
              skyflux=specflux-objsub1
              sxaddpar, objhdr, 'SKYSUB', 'Sky spline-fit and subtracted'
              if ccdnum eq 2 then begin
                sxaddpar, objhdr, 'NSKIES1', 'Used '+strn(nccd1)+' good skies on CCD1'
                sxaddpar, objhdr, 'NSKIES2', 'Used '+strn(nccd)+' good skies on CCD2'
              endif
           endelse
            if ccdnum eq 2 then sky2d=sky2d_old

            pixmask = newmask1            
            objsub = objsub1
            objsubivar = objsubivar1

            skymulti='true'
            if ccdnum eq 1 then sset1 = sset1
            if ccdnum eq 1 then sky1 = skyflux
            if ccdnum eq 1 then scale1 = scale
            if ccdnum eq 2 then sky2 = skyflux
            if ccdnum eq 2 then sset2 = sset1
            if ccdnum eq 2 then scale2 = scale
            

            if keyword_set(qaplot) then begin
             split = strsplit(filelist(j), '/,.fits', /extract)
             ;airset=wset
             delvarx, xx
             traceset2xy,airset, xx, temp
             temp1 = temp
             for ifib = 0, n_elements(temp(0,*)) -1 do begin
                temp1(*,ifib) = findgen(n_elements(temp(*,0)))
             endfor
             temp = alog10(temp)
             xy2traceset, temp1, temp, tempset, ncoeff=5
             iskies = where(strmatch(plugmap.objtype, 'SKY*'))
             
             hs_qaplot_skysub, specflux, specfluxivar, $
                objsub, objsubivar, tempset, iskies
             traceset2xy, airset, xx, lam
             
             k = where(strmatch(plugmap.objtype, 'SKY*'))
             
             for i = 0, n_elements(k) -1 do begin
                if i eq 0 then plot, lam[*,k[i]],$
                   median(objsub[*,k[i]],100), /xstyle
                if i gt 0 then oplot, lam[*,k[i]], median(objsub[*,k[i]],100)
                
                if ccdnum eq 1 then factor=0
                if ccdnum eq 2 then factor=150
                djs_xyouts, max(smooth(objsub[*,k(i)],100))*0.90, 4300, $
                   'Aperture ' + string(k(i)+factor, $
                                        format='(i3.3)'), color='red'
             endfor
             
          endif
          
         endif else begin
          delvarx, xx 
          traceset2xy, airset, xx,lam
          objsub = specflux
          objsubivar = specfluxivar
          pixmask=0.0*objsub
          sky1=pixmask
          sky2=pixmask
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
          specsky = dblarr(ncol, 300)
          specfluxivar = dblarr(ncol, 300)
          specflux(*,0:149) = flux1
          specflux(*,150:299) = flux2
          if keyword_set(skysubtract) then begin
             specsky(*,0:149) = sky1
             specsky(*,150:299) = sky2
          endif
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
          delvarx, sky1
          delvarx, sky2
          delvarx, fluxivar1
          delvarx, fluxivar2
          delvarx, pixmask1
          delvarx, pixmask2
          

          if keyword_set(skysubtract) and (nccd le 3) or (nother le 3) $
             then begin
             
             delvarx, xx 
             traceset2xy, airset, xx,lam
            k = where(strmatch(plugmap.objtype, 'SKY*'))
            xy2traceset, xx[*,k], lam[*,k], tempset, ncoeff=6
            if keyword_set(dummy) then jd=0
            skystuct = hs_skysubtract(specflux[*,k],$
                                    specfluxivar[*,k], plugmap[k], $
                                    tempset,   objsub1, objsubivar1, fullfit=fullfit,$
                                    newmask=newmask, jd=jd, sset=sset)

;             hs_scale, lam, specflux, sigma=3.0, scale=scale
             hs_skyfit, lam, specflux ,scale=scale, skyset=sset
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
             sky1=specflux-objsub1
             skysingle = 'true'
             pixmask = newmask1
             sxaddpar, objhdr, 'SKYCRUDE', 'Sky subtracted with low number of sky fibers' 
             
             if keyword_set(qaplot) then begin
                split = strsplit(filelist(j), '/,.fits', /extract)
                traceset2xy,airset, xx, temp
                temp1 = temp
                for ifib = 0, n_elements(temp(0,*)) -1 do begin
                   temp1(*,ifib) = findgen(n_elements(temp(*,0)))
                endfor
                temp = alog10(temp)
                xy2traceset, temp1, temp, tempset, ncoeff=5
                iskies = where(strmatch(plugmapout.objtype, 'SKY*'))
                
                hs_qaplot_skysub, specflux, specfluxivar, objsub, $
                   objsubivar, tempset, iskies
                
                
             endif
             
             specflux = objsub1
             specfluxivar = objsubivar1
             specsky=sky1
          endif

         ;do the telluric correction on the whole frame at once
         if keyword_set(tellcorr) then begin
           delvarx, xx
           traceset2xy, airset, xx, lambda
           airmass=tai2airmass(plugmapout.ra, plugmapout.dec, jd=jd, altitude=2606., longitude=(360.-(110.+53./60.+4.4/3600.)),latitude=(31.+41/60.+19.6/3600.))
           hs_remove_telluric, lambda, flux=specflux, ivar=specfluxivar, airmass=airmass
           sxaddpar, objhdr, 'TELLCOR', 'Telluric correction for A and B bands applied'
         endif

         ;and then do the heliocentric correction the the lambda values
         ;last, for the full frame
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
             helio = heliocentric(hms2dec(ra)*15., hms2dec(dec), jd=jd,altitude=2606., longitude=(360.-(110.+53./60.+4.4/3600.)),latitude=(31.+41/60.+19.6/3600.))
             splog, 'Heliocentric correction = ', helio, ' km/s'
             sxaddpar, objhdr, 'HELIO_RV', helio, $
              ' Heliocentric correction (km/s, added)'
          endif else begin
             splog, 'WARNING: Header info not present' + $
                ' to compute heliocentric correction'
          endelse
        
        
          ;;Apply the correction
          delvarx, xx
          traceset2xy, airset, xx, lamtemp
          lamtemp = lamtemp / (1+helio/299792.458)
          xy2traceset, xx, lamtemp, airset, ncoeff=6

         ;take the just-corrected lambdas to put in the final array
          lam=lamtemp

       endif

       

       if j eq 0 and ccdnum eq 1 then begin
          ;;Setup the output arrays that we need
          lamout = dblarr( n_elements(filelist), $
                           n_elements(lam(*,0)), 300)
          fluxout = dblarr( n_elements(filelist), $
                            n_elements(lam(*,0)), 300)
          skyout = dblarr( n_elements(filelist), $
                            n_elements(lam(*,0)), 300)
          ivarout=  dblarr( n_elements(filelist), $
                            n_elements(lam(*,0)), 300)
          pixmaskout=  dblarr( n_elements(filelist), $
                               n_elements(lam(*,0)), 300)
          exptimes=fltarr(n_elements(filelist))
;          headout = strarr(n_elements(filelist),$
;                           n_elements(objhdr))
          fullfitout=  dblarr( n_elements(filelist), $
                               n_elements(lam(*,0)), 300)
       endif
    endfor
    lamout[j,*, *] = lam
    fluxout[j,*, *] = specflux
    skyout[j,*, *] = specsky
    ivarout[j,*, *] = specfluxivar
    pixmaskout[j,*,*] = pixmask
;    headout[j,*] = objhdr
    
    if keyword_set(writeall) then begin
       
       file1 = strsplit(filelist(j),'/',/extract)
       file1=file1(n_elements(file1)-1)
       
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

       if j eq 0 then texptime=sxpar(objhdr, 'EXPTIME') else texptime=texptime+sxpar(objhdr, 'EXPTIME')       
       exptimes[j]=sxpar(objhdr, 'EXPTIME')
       splog, 'Writing single obs file ' + outfile
       
       sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
       sxaddpar, objhdr, 'COMMENT', 'HDU0: wavelengths (Angstroms)'
if not keyword_set(noftweak) then begin
       sxaddpar, objhdr, 'COMMENT', 'HDU1: sky-sub, weighted coadded spectra (total counts)'
       sxaddpar, objhdr, 'COMMENT', 'HDU2: inverse variance (counts)'
       sxaddpar, objhdr, 'COMMENT', 'HDU3: AND bad pixel mask'
       sxaddpar, objhdr, 'COMMENT', 'HDU4: OR bad pixel mask (coadded) or sky spectra (single-exposure)'
       sxaddpar, objhdr, 'COMMENT', 'HDU5: Plugmap structure (fiber info)'
       sxaddpar, objhdr, 'COMMENT', 'HDU6: Combined sky spectra'
       sxaddpar, objhdr, 'COMMENT', 'HDU7: Summed (unweighted) spectra'
    endif else begin
       sxaddpar, objhdr, 'COMMENT', 'HDU1: sky-sub, fluxed spectra'
       sxaddpar, objhdr, 'COMMENT', 'HDU2: inverse variance'
       sxaddpar, objhdr, 'COMMENT', 'HDU3: AND bad pixel mask'
       sxaddpar, objhdr, 'COMMENT', 'HDU4: OR bad pixel mask (coadded) or sky spectra (single-exposure)'
       sxaddpar, objhdr, 'COMMENT', 'HDU5: Plugmap structure (fiber info)'
endelse
       
       mwrfits, float(transpose(transpose(lamout[j,*,*]))),$
          outfile, objhdr, /create
       mwrfits, float(transpose(transpose(fluxout[j,*,*]))), outfile
       mwrfits, float(transpose(transpose(ivarout[j,*,*]))), outfile
       mwrfits, float(transpose(transpose(pixmaskout[j,*,*]))), outfile
       if keyword_set(skysubtract) then mwrfits, specsky, outfile else  mwrfits, float(transpose(transpose(pixmaskout[j,*,*]))), outfile
       mwrfits, plugmapout, outfile
       
       
       if keyword_set(skysubtract) then begin
          ;writing out these redundant extensions in case of an
          ;observation with a single exposure,to match rough format of spHect files
          mwrfits, specsky, outfile
          mwrfits, float(transpose(transpose(fluxout[j,*,*]))), outfile
          specsky_out=specsky
          ;if skymulti eq 'true' then begin
          ;   mwrfits, sset1, outfile
          ;   mwrfits, scale1, outfile
          ;   mwrfits, sset2, outfile
          ;   mwrfits, scale2, outfile
          ;endif
          ;if skysingle eq 'true' then begin
          ;   mwrfits, sset, outfile
          ;   mwrfits, scale, outfile
          ;   mwrfits,  0.0 , outfile
          ;   mwrfits,  0.0 , outfile
          ;   
          ;endif
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
    delvarx, specsky
    delvarx, specfluxivar
    delvarx, pixmask
 endfor


 if keyword_set(fluxcorr) then begin
    
    flux_factor = hs_fluxcorr(lamout, fluxout,$
                              ivarout, plugmapout, $
                              pixmaskout, outname = filelist, $
                              plateid=plateid, $
                              rerun=rerun,psplot=qaplot, checkaverage=checkaverage,$
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
    
    if keyword_set(combine) then begin
       
       ;;Now combine each of the fibers into a composite spectrum
       newlam = transpose(transpose(lamout[0,*,*]))
       newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
       newfluxivar = transpose(transpose(ivarout[0,*,*])) * 0.0
       newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
       ormask = newmask
       
       for i = 0, nfib-1 do begin
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
 
       sxaddpar, objhdr, 'FLUXCOR','Exposures flux-calibrated and averaged'
       
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
;    finalivar = ivarout
    ivarout = finalivar
    pixmaskout = finalandmask
    ormaskout = finalormask
    
 endif
 
 ;;This should handle the combining of non-fluxed data 
 ;;in a way that will not kill the data
 size1 = size(lamout)
 nim = size1[1]
 npix = size1[2]
 nfib = size1[3]
 if not keyword_set(fluxcorr) and keyword_set(combine) then begin   
    if keyword_set(noftweak) then begin
       ;do a very simple coadd here, and skip the rest
      if nim eq 1 then begin
       fluxout = transpose(transpose(fluxout))
       fluxout_boxcar = fluxout
       if keyword_set(skysubtract) then skyout=specsky_out else skyout=fluxout*0.0
       ivarout = transpose(transpose(ivarout))
       pixmaskout = transpose(transpose(pixmaskout))
       ormaskout = pixmaskout
       lamout=transpose(transpose(lamout))
       sxaddpar, objhdr, 'SINGLE_EXP', 'Combined-format file, but only a single exposure was available for the coadd'
       sxaddpar, objhdr, 'COMMENT', 'Reduced '+systime()
       sxaddpar, objhdr, 'COMMENT', 'using SAO/TDC fork of IDL Hectospec pipeline created by R. Cool'
      endif else begin
       newflux=fltarr(npix,nfib)
       newfluxivar=newflux
       newflux_alt=newflux
       newfluxivar_alt=newflux
       newsky = newflux
       newlam=newflux
       if nim eq 3 then newlam[*,*] = lamout[1,*,*] else if nim eq 5 then newlam[*,*] = lamout[2,*,*] else newlam[*,*] = lamout[0,*,*] 
       newmask=newflux
       ormask=newflux
       tempflux = fluxout 
       tempsky=skyout
       tempivar = ivarout 
       for ifib=0, nfib-1 do begin                
          print, format='("Co-adding ",i4," of ",i4,a1,$)', ifib+1, nfib, string(13b)
         ;if there are 3 or 5 exposures, minimize shifts by using wavelengths from the middle one
         if nim eq 3 then templam = lamout[1,*,ifib] else if nim eq 5 then templam = lamout[2,*,ifib] else templam = lamout[0,*,ifib] 
         for iexp=0, nim-1 do begin                
            if ((nim eq 3) and (iexp eq 1)) or ((nim eq 5) and (iexp eq 2)) or ((nim ne 3) and (nim ne 5) and (iexp eq 0)) then begin
              tempflux[iexp,*,ifib]=fluxout[iexp,*,ifib]
              tempsky[iexp,*,ifib]=skyout[iexp,*,ifib]
              tempivar[iexp,*,ifib]=ivarout[iexp,*,ifib]
              newmask1=pixmaskout[iexp,*,ifib]
            endif else begin
              ;use hs_combine1fiber to just resample spectrum onto the reference
              ;wavelength grid
              hs_combine1fiber,transpose(lamout[iexp,*,ifib]), transpose(fluxout[iexp,*,ifib]), transpose(ivarout[iexp,*,ifib]), skyflux=transpose(skyout[iexp,*,ifib]), $
                newlam=transpose(templam), newflux=newflux1,$
                newivar=newivar1, newsky=newsky1, finalmask = transpose(pixmaskout[iexp,*,ifib]), $
                andmask = newmask1
              tempflux[iexp,*,ifib]=newflux1
              tempsky[iexp,*,ifib]=newsky1
              tempivar[iexp,*,ifib]=newivar1
            endelse
            if iexp eq 0 then newmask[*,ifib]=newmask1 else newmask[*,ifib]=(newmask[*,ifib] AND newmask1)
            if iexp eq 0 then ormask[*,ifib]=newmask1 else ormask[*,ifib]=(newmask[*,ifib] OR newmask1)
         endfor
         newflux[*,ifib]=total(tempflux[*,*,ifib], 1)
         newsky[*,ifib]=total(tempsky[*,*,ifib], 1)
         newfluxivar[*,ifib]=1./total(1./tempivar[*,*,ifib],1)
         ;alternative: weighted coadd
         ;first switch to counts per second
         for z=0, nim-1 do begin
            tempflux[z,*,ifib]=tempflux[z,*,ifib]/exptimes[z]
            tempivar[z,*,ifib]=tempivar[z,*,ifib]*exptimes[z]^2.
         endfor
         newflux_alt[*,ifib]=total(tempflux[*,*,ifib]*tempivar[*,*,ifib], 1)/total(tempivar[*,*,ifib], 1)*texptime
;         newfluxivar_alt[*,ifib]=1./total(1./tempivar[*,*,ifib],1)/(texptime/nim)^2.
         newfluxivar_alt[*,ifib]=total(tempivar[*,*,ifib],1)/texptime^2.

         bad=where(finite(newfluxivar[*,ifib],/nan), badcnt)
         if badcnt ge 1 then newfluxivar[bad,ifib]=0.0
         bad=where(finite(newfluxivar[*,ifib],/inf), badcnt)
         if badcnt ge 1 then newfluxivar[bad,ifib]=0.0
       endfor
       fluxout = newflux_alt
       fluxout_boxcar = newflux
       skyout=newsky
       ivarout = newfluxivar_alt
;stop
       pixmaskout = newmask
       ormaskout = ormask
       lamout=newlam
       sxaddpar, objhdr, 'COMBINED', 'Weighted sum of'+string(nim, format='(I2)')+' exposures, tot exptime: '+string(texptime,format='(I8.1)')
       sxaddpar, objhdr, 'COMMENT', 'Reduced '+systime()
       sxaddpar, objhdr, 'COMMENT', 'using SAO/TDC fork of IDL Hectospec pipeline created by R. Cool'
       sxaddpar, objhdr, 'EXPTIME', texptime
     endelse
    endif else begin   
       ;otherwise, do a crude flux scaling here and
       ;proceed with the old way of coadding
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
     endif   
     ;;Now combine each of the fibers into a composite spectrum
     newlam = transpose(transpose(lamout[0,*,*]))
     newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
     newfluxivar = transpose(transpose(ivarout[0,*,*])) * 0.0
     newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
     ormask = newmask
     size1 = size(lamout)
     nim = size1(1)
     npix = size1(2)
     nfib = size1(3)
     for i = 0, nfib-1 do begin
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
   endelse    
 endif

 
 if keyword_set(writeout) then begin
    
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
    if keyword_set(skysubtract) and keyword_set(noftweak) and keyword_set(combine) and not keyword_set(fluxcorr) then begin
       mwrfits, skyout, outname
       mwrfits, fluxout_boxcar, outname
       ;if skymulti eq 'true' then begin
       ;   mwrfits, sset1, outname
       ;   mwrfits, scale1, outname
       ;   mwrfits, sset2, outname
       ;   mwrfits, scale2, outname
       ;endif
       ;if skysingle eq 'true' then begin
       ;   mwrfits, sset, outname
       ;   mwrfits, scale, outname
       ;   mwrfits,  0.0 , outname
       ;   mwrfits,  0.0 , outname
       ;   
       ;endif
    endif ;else begin
   ;      mwrfits, float(pixmaskout), outname
   ;      mwrfits, float(pixmaskout), outname
 ;   endelse
    
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



