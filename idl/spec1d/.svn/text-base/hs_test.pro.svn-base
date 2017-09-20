;+
; NAME:
;   hs_reduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   hs_reduce1d, platefile, wavemin=wavemin, wavemax=wavemax, rerun=rerun, 
;                pseudoflux=pseudoflux
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platefile  - Plate file(s) from spectro-2D; default to all files
;                matching 'spPlate*.fits'
;   wavemin    - Minimum wave to use in the redshift
;   wavemax    - maximum wave to use in the redshift algorithm
;   rerun      - reduction rerun
;   pseudoflux - apply an averaged vector if the data is not fluxed
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Names of output files are derived from PLATEFILE.
;   For example, if PLATEFILE='spHect-0306-51690.fits', then
;     ZALLFILE = 'spZall-0306-51690.fits'
;     ZBESTFILE = 'spZbest-0306-51690.fits'
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   cpbackup
;   dfpsclose
;   dfpsplot
;   elodie_best()
;   filter_thru()
;   mrdfits()
;   mwrfits
;   qaplot_fcalibvec
;   splog
;   skymask()
;   speclinefit
;   struct_addtags()
;   sxaddpar
;   sxdelpar
;   sxpar()
;   synthspec()
;   vdispfit
;   zfind()
;   zrefind()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro hs_test, platefile, doplot=doplot, debug=debug, fiberid=fiberid, $
                 wavemin=wavemin, wavemax=wavemax, rerun=rerun, pass=pass, $
                 field=field, pseudoflux=pseudoflux,airmass=airmass
  
  
  if NOT keyword_set(wavemin) then wavemin=4100
  if NOT keyword_set(wavemax) then wavemax=8500
  
   if (NOT keyword_set(platefile)) then begin
      splog,'Please specify an input file'
   endif else begin
      if (size(platefile,/tname) NE 'STRING') then $
       message, 'PLATEFILE must be a file name'
      if (keyword_set(platefile)) then nplate = n_elements(platefile) $
       else nplate = 0
   endelse
   if (keyword_set(debug)) then doplot = 1

   ;----------
   ; If multiple plate files exist, then call this script recursively
   ; for each such plate file.

   if (nplate EQ 0) then begin
      splog, 'No plate files specified or found'
      return
   endif else if (nplate EQ 1) then begin
      platefile = platefile[0]
   endif else begin
      for i=0, nplate-1 do begin
         hs_reduce1d, platefile[i], fiberid=fiberid, doplot=doplot, debug=debug
      endfor
      return
   endelse

   ;----------
   ; Determine names of output files
   
 
   if keyword_set(pass) and keyword_set(field) then begin
      platemjd = string(pass, format='(i1.1)') + string(field, format='(i2.2)')
      outkey = platemjd + '-' + rerun
   endif else begin
      platemjd= '0000'
      junk = strsplit(platefile, '/', /extract)
      junk1 = strsplit(junk(n_elements(junk)-1), ' ', length=ct)
      outkey = strmid(junk(n_elements(junk)-1), 0, ct-5)
      if strmid(outkey, 0, 7) eq 'spHect-' THEN $
         outkey = strmid(outkey, 7, ct-5)
 
   endelse
   
   if not keyword_set(rerun) then rerun = string(0, format='(i4.4)')
   
   zallfile = 'spZall-' + outkey+'.fits'
   zbestfile = 'spZbest-' + outkey +'.fits'
   zlinefile = 'spZline-' + outkey +'.fits'
   if (NOT keyword_set(logfile)) then $
    logfile = 'spDiag1d-' + outkey +'.log'
   plotfile = 'spDiag1d-' + outkey +'.ps'

   if (keyword_set(doplot) AND NOT keyword_set(debug)) then begin
      debugfile = 'spDiagDebug1d-' + platemjd +'-' + rerun+ '.ps'
      cpbackup, debugfile
      dfpsplot, debugfile, /color
   endif

   stime0 = systime(1)

   if (keyword_set(logfile)) then begin
      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + systime()
   endif
   if (keyword_set(plotfile)) then $
    splog, 'Plot file ' + plotfile
   if (keyword_set(debugfile)) then $
    splog, 'Debug plot file ' + debugfile
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   ;----------
   ; Read the 2D output file
   
   wave = mrdfits(platefile,0)
   airtovac, wave
   objflux = mrdfits(platefile,1,hdr)
   if (NOT keyword_set(hdr)) then $
    message, 'Plate file not valid: ' + platefile
   npixobj = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
   objivar = mrdfits(platefile,2)
   andmask = mrdfits(platefile,3)
   ormask = mrdfits(platefile,4)
   plugmap = mrdfits(platefile,5)
   
   
   
   npix = n_elements(wave(*,0))
   
   ;;July 23 04 - if the wavemin and wavemax are 
   ;;outside the data range, go to the data range
   
   wavemin = max(wave[2,*]) >  wavemin
   wavemax = min(wave[npix-2, *]) < wavemax   
   
   print, wavemin, wavemax
   
   
   if Keyword_set(pseudoflux) then begin
      
      fcalfile = getenv('HSRED_DIR') + '/etc/average_flux.fits'
      tcalfile = getenv('HSRED_DIR') + '/etc/average_tell.fits'
      
      fcal = mrdfits(fcalfile, 1)
      tcal = mrdfits(tcalfile, 1)
      
      flux_factor = bspline_valu(alog10(wave), fcal)
      tell_factor = bspline_valu((wave), tcal)
      
      divideflat, objflux, flux_factor, invvar=objivar
      divideflat, objflux, tell_factor, invvar=objivar
      
      
      ;;Changed on 9-24 to mask out region of telluric
      
      nfiber = n_elements(wave(0,*))
      for ifib = 0, nfiber-1 do begin
         k = where(abs(tell_factor(*,ifiber)-1.0) gt 0.01)
         objivar(k,ifiber) = 0.0 ;; Kill the inverse variance
         andmask(k,ifiber) = andmask(k,ifiber)+134217728L ;; Set the mask bit 
         ormask(k, ifiber) = ormask(k,ifiber) +134217728L ;; Set the mask bit
      endfor
      
      
   endif
   
   
   
   if keyword_set(wavemin) then begin
      k = where( abs(wave(*,0) - wavemin) eq min(abs(wave(*,0) - wavemin)))
      
      npix = n_elements(wave(*,0))
      wave = wave(k(0):npix-1,*)
      objflux = objflux(k(0):npix-1,*)
      objivar = objivar(k(0):npix-1,*)
      andmask = andmask(k(0):npix-1,*)
      ormask= ormask(k(0):npix-1,*)
   endif
     if keyword_set(wavemax) then begin
      k = where( abs(wave(*,0) - wavemax) eq min(abs(wave(*,0) - wavemax)))
      
      npix = n_elements(wave(*,0))
      wave = wave(0:k(0),*)
      objflux = objflux(0:k(0),*)
      objivar = objivar(0:k(0),*)
      andmask = andmask(0:k(0),*)
      ormask= ormask(0:k(0),*)
   endif
   
   
   
   
   npixobj = n_elements(wave(*,0))
   sxaddpar, hdr, 'NAXIS1', npixobj

   
  
   
   anyandmask = transpose(andmask[0,*])
   anyormask = transpose(ormask[0,*])
   for ipix=1, npixobj-1 do $
    anyandmask = anyandmask OR transpose(andmask[ipix,*])
   for ipix=1, npixobj-1 do $
      anyormask = anyormask OR transpose(ormask[ipix,*])
   objivar = skymask(objivar, andmask, ormask)+objivar
   
   
   ;----------
   ; Trim to specified fibers if FIBERID is set

   if (keyword_set(fiberid)) then begin
      if (min(fiberid) LE 0 OR max(fiberid) GT nobj) then $
         message, 'Invalid value for FIBERID: must be between 0 and '$
         +string(nobj)
      objflux = objflux[*,fiberid-1]
      objivar = objivar[*,fiberid-1]
      wave = wave[*,fiberid-1]
      anyandmask = anyandmask[fiberid-1]
      anyormask = anyormask[fiberid-1]
      plugmap = plugmap[fiberid-1]
      nobj = n_elements(fiberid)
   endif else begin
      fiberid = lindgen(nobj) + 1
   endelse
   splog, 'Number of fibers = ', nobj

   ;----------
   ; Look for where the S/N is unreasonably large
   ; or where flux is unphysically negative.

   for iobj=0L, nobj-1 do begin
      junk = where(abs(objflux[*,iobj]) * sqrt(objivar[*,iobj]) GT 200., ct)
      if (ct GT 0) then $
       splog, 'WARNING: Fiber #', fiberid[iobj], $
        ' has ', ct, ' pixels with S/N > 200'

      junk = where(objflux[*,iobj] * sqrt(objivar[*,iobj]) LE -10., ct)
      if (ct GT 0) then $
       splog, 'WARNING: Fiber #', fiberid[iobj], $
        ' has ', ct, ' pixels with Flux < -10*Noise'
   endfor

   ;----------
   ; Mask out points that are unphysically negative (10-sigma negatives),
   ; and mask the neighboring 2 pixels in each direction.
   
   
   for iobj=0L, nobj-1 do begin
      thismask = objflux[*,iobj] * sqrt(objivar[*,iobj]) LE -10.
      thismask = smooth(float(thismask),5) GT 0
      objivar[*,iobj] = objivar[*,iobj] * (1 - thismask)
   endfor
   
      
   ;;Now transform the spectra to the correct scaling
                                ;(using the galaxy eigentemplate)
   npix = n_elements(wave(*,0))
   coeff0 = alog10(wave(0,0))
   coeff1 = 1.0e-4
   naxis = (alog10(wave(npix-1,0)) - coeff0) / coeff1
   
   
   sxaddpar, hdr,'COEFF0', coeff0
   sxaddpar, hdr,'COEFF1', coeff1
   sxaddpar, hdr,'NAXIS1', naxis
   
   npixobj=naxis
   
   lam = 10D^(coeff0+coeff1*findgen(naxis)) ; logarithmic wavelength vector
   newloglam1 = (lam)
   
   logwave = alog10(wave)
   newloglam = dblarr(n_elements(newloglam1), n_elements(logwave(0,*)))
   for i = 0, n_elements(logwave(0,*)) -1 do begin
      newloglam(*,i) = newloglam1
   endfor
   
   
   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')
   
   
   newflux = newloglam * 0.0
   newivar = newloglam * 0.0
   newand = newloglam * 0.0
   newor = newloglam * 0.0
   
   
   for i = 0, n_elements(logwave(0,*)) -1 do begin
      hs_combine1fiber, wave(*,i), objflux(*,i), objivar(*,i), $
         newlam=newloglam1, finalmask = andmask(*,i), $
         newflux=newflux1, newivar=newivar1, andmask=andmask1, $
         ormask=ormask1
      
      newflux(*,i) = newflux1
      newivar(*,i) = newivar1
      newand(*,i) = andmask1
      newor(*,i) = ormask1
      
   endfor
   

   
   
   objflux = newflux
   objivar = newivar
   andmask = newand
   ormask =  newor
   


   
   ;----------
   ; Find GALAXY redshifts

   npoly = 5
   zmin = 0.3 ; -3000 km/sec
   zmax = 0.7 ; 
   pspace = 10
   nfind = 5
   plottitle = 'Galaxy Redshift'

   eigenfile = 'spEigenGal-*.fits'

   splog, 'Compute GALAXY redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_gal = zfind(objflux, objivar, hdr=hdr, $
    eigenfile=eigenfile, npoly=npoly, zmin=zmin, zmax=zmax, pspace=pspace, $
    nfind=nfind, width=5*pspace, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to compute GALAXY redshifts = ', systime(1)-t0

   splog, 'Locally re-fitting GALAXY redshifts'
   t0 = systime(1)
   res_gal = zrefind(objflux, objivar, hdr=hdr, $
    pwidth=5, pspace=1, width=5, zold=res_gal, $
    plottitle=plottitle, doplot=doplot, debug=debug)
   splog, 'CPU time to re-fit GALAXY redshifts = ', systime(1)-t0

   
   res_gal.class = 'GALAXY'
   res_gal.subclass = ' '

   res_all = res_gal ; Append results
   nper = (size(res_all,/dimens))[0]

   
   ;----------
   ; Sort results for each object by ascending order in chi^2/DOF,
   ; but putting any results with zero degrees-of-freedom at the end.

   minvdiff = 1000.0 ; km/s
   cspeed = 2.99792458e5

   for iobj=0, nobj-1 do begin
      res1 = res_all[*,iobj]

      rchi2 = res1.rchi2

      isort = sort(rchi2 + (res1.dof EQ 0)*max(rchi2))
      for ii=0, nper-1 do begin
         res_all[ii,iobj] = res1[isort[ii]]
      endfor

      ; Find the difference in reduced chi^2 between each result and the next
      res1 = res_all[*,iobj]
      rchi2 = res1.rchi2
      for ii=0, nper-2 do begin
         inext = (where( $
          abs(res1[ii+1:nper-1].z - res1[ii].z) GT minvdiff/cspeed $
          AND res1[ii+1:nper-1].dof GT 0))[0]
         if (inext NE -1) then $
          res_all[ii,iobj].rchi2diff = rchi2[ii+1+inext] - rchi2[ii]
      endfor
   endfor

   ;----------
   ; Generate the synthetic spectra, and count the fraction of points
   ; that deviate more than N sigma (where N goes from 1 to NFSIG).

   t0 = systime(1)
   nfsig = 10
   chi68p = fltarr(nper,nobj)
   fracnsigma = fltarr(nfsig,nper,nobj)
   fracnsighi = fltarr(nfsig,nper,nobj)
   fracnsiglo = fltarr(nfsig,nper,nobj)
   counts_spectro = fltarr(5,nper,nobj)
   counts_synth = fltarr(5,nper,nobj)

   objloglam = objloglam0 + lindgen(npixobj) * objdloglam
   wavevec = 10d^objloglam
   flambda2fnu = wavevec^2 / 2.99792e18

   fthru = filter_thru(objflux * rebin(flambda2fnu,npixobj,nobj), $
    waveimg=wavevec, mask=(objivar EQ 0))
   counts_spectro[*,0,*] = transpose(fthru) * 10^((48.6 - 2.5*17.)/2.5)

   ; Loop in reverse order, so that we look at the best-fit spectra last,
   ; and keep those spectra around for later.

; Save time for now and only look at best fit, since SYNTHSPEC and
; FILTER_THRU are so slow ???
;   for iper=nper-1, 0, -1 do begin
   for iper=0, 0, -1 do begin
      ; Copy this for all fits, since the measured magnitudes are the same
      counts_spectro[*,iper,*] = counts_spectro[*,0,*]

      synflux = synthspec(res_all[iper,*], loglam=objloglam)

      for iobj=0, nobj-1 do begin
         igood = where(objivar[*,iobj] GT 0, ngood)
         if (ngood GT 0) then begin
            chivec = (objflux[igood,iobj] - synflux[igood,iobj]) $
             * sqrt(objivar[igood,iobj])
            abschivec = abs(chivec)
            chi68p[iper,iobj] = (abschivec[sort(abschivec)])[floor(0.68*ngood)]
            for isig=0, nfsig-1 do begin
               fracnsigma[isig,iper,iobj] = total(abschivec GT isig+1) / ngood
               fracnsighi[isig,iper,iobj] = total(chivec GT isig+1) / ngood
               fracnsiglo[isig,iper,iobj] = total(chivec LT (-isig-1)) / ngood
            endfor
         endif
      endfor

      fthru = filter_thru(synflux * rebin(flambda2fnu,npixobj,nobj), $
       waveimg=wavevec)
      counts_synth[*,iper,*] = transpose(fthru) * 10^((48.6 - 2.5*17.)/2.5)
   endfor

   splog, 'CPU time to generate chi^2 statistics = ', systime(1)-t0

   ;----------
   ; Zero-out the dispersion template if the best-fit was not a galaxy.

   for iobj=0, nobj-1 do begin
      if (strtrim(res_gal[iobj].class,2) NE 'GALAXY') then $
       dispflux[*,iobj] = 0
   endfor

   ;----------
   ; Add other fields to the output structure

   splog, 'Adding other fields to output structure'
   res1 = { plate:    platemjd, $ ;;long(sxpar(hdr, 'PLATEID')), $
           ;; tile:     long(sxpar(hdr, 'TILEID')), $
            mjd:      long(sxpar(hdr, 'MJD')), $
            fiberid:  0L        , $
            objid:    lindgen(5), $
            objtype:  ' '       , $
            plug_ra:  0.0d      , $
            plug_dec: 0.0d      }
   res_prepend = make_array(value=res1, dimension=size(res_all,/dimens))
   res_all = struct_addtags(res_prepend, res_all)

   for iobj=0, nobj-1 do begin
      res_all[*,iobj].fiberid = fiberid[iobj]
      res_all[*,iobj].objid = plugmap[iobj].icode
      res_all[*,iobj].objtype = plugmap[iobj].objtype
      res_all[*,iobj].plug_ra = plugmap[iobj].ra
      res_all[*,iobj].plug_dec = plugmap[iobj].dec
   endfor

   res1 = { wavemin:   0.0, $
            wavemax:   0.0, $
            wcoverage: 0.0, $
            zwarning:  0L, $
            sn_median: 0.0, $
            chi68p: 0.0, $
            fracnsigma: fltarr(nfsig), $
            fracnsighi: fltarr(nfsig), $
            fracnsiglo: fltarr(nfsig), $
            counts_spectro: fltarr(5), $
            counts_synth: fltarr(5), $
            counts_sky: fltarr(5), $
            anyandmask: 0L, $
            anyormask:  0L}
            ;;spec1_g:   sxpar(hdr, 'SPEC1_G'), $
            ;;spec1_r:   sxpar(hdr, 'SPEC1_R'), $
            ;;spec1_i:   sxpar(hdr, 'SPEC1_I'), $
            ;;spec2_g:   sxpar(hdr, 'SPEC2_G'), $
            ;;spec2_r:   sxpar(hdr, 'SPEC2_R'), $
            ;;spec2_i:   sxpar(hdr, 'SPEC2_I') }
   res_append = make_array(value=res1, dimension=size(res_all,/dimens))
   res_all = struct_addtags(res_all, res_append)

   for iobj=0, nobj-1 do begin
      igood = where(objivar[*,iobj] NE 0, ngood)
      res_all[*,iobj].wavemin = $
       10^(objloglam0 + (igood[0]>0)*objdloglam) * (ngood NE 0)
      res_all[*,iobj].wavemax = $
       10^(objloglam0 + (igood[(ngood-1)>0])*objdloglam) * (ngood NE 0)
      res_all[*,iobj].wcoverage = ngood * objdloglam
      res_all[*,iobj].anyandmask = anyandmask[iobj]
      res_all[*,iobj].anyormask = anyormask[iobj]
      if (ngood GT 0) then $
       res_all[*,iobj].sn_median = $
        median( sqrt(objivar[igood,iobj]) * abs(objflux[igood,iobj]))
   endfor

   res_all.chi68p = chi68p
   res_all.fracnsigma = fracnsigma
   res_all.fracnsighi = fracnsighi
   res_all.fracnsiglo = fracnsiglo
   res_all.counts_spectro = counts_spectro
   res_all.counts_synth = counts_synth

   ;----------
   ; Generate output headers for spZbest, spZall, spZline files.

   sxaddpar, hdr, 'NAXIS', 0
   sxdelpar, hdr, 'NAXIS1'
   sxdelpar, hdr, 'NAXIS2'
   sxaddpar, hdr, 'EXTEND', 'T', after='NAXIS'
   sxaddpar, hdr, 'VERS1D', idlspec2d_version(), $
    ' Version of idlspec2d for 1D reduction', after='VERSCOMB'
   spawn, 'uname -n', uname
   sxaddpar, hdr, 'UNAME', uname[0]
   ww = strsplit(uname[0], '.', /extract)
   if (ww[1<(n_elements(ww)-1)] EQ 'fnal') then return

   ;----------
   ; Call the line-fitting code for this plate

   splog, 'Call line-fitting code'

   
   ;----------
   ; Set ZWARNING flags.

   splog, 'Setting flags'
   zwarning = lonarr(nper,nobj)

   ; Warning: Sky fiber.
   for iobj=0, nobj-1 do begin
      if (strtrim(plugmap[iobj].objtype,2) EQ 'SKY') then $
       zwarning[*,iobj] = zwarning[*,iobj] OR sdss_flagval('ZWARNING', 'SKY')
   endfor

   
   ;;July 23 04 - removed wavelength coverage flag code as I am
   ;;not sure how to get it to know about the wavemin and wavemax
   ;;this is not a big needed flag!
   
   
   
   ; Warning: too little wavelength coverage.
   ;qflag = res_all.wcoverage LT 0.18
   ;zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'LITTLE_COVERAGE')

   ; Warning: delta-chi^2 is too small as compared to the next best ID.
   minrchi2diff = 0.01
   qflag = res_all.rchi2diff LT minrchi2diff $
    OR res_all.rchi2diff LT minrchi2diff * res_all.rchi2
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'SMALL_DELTA_CHI2')

   ; Warning: synthetic spectrum is negative (for STAR only).
   qflag = (strtrim(res_all.class) EQ 'STAR' $
    AND strtrim(res_all.subclass) NE 'CV' $
    AND res_all.theta[0] LE 0)
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'NEGATIVE_MODEL')

   ; Warning: Fraction of points above 5 sigma is too large (> 5%),
   ; except for QSO's where we just look at the fraction of high outliers
   ; since we expect absorption lines that could give many low outliers.
   qflag = (strtrim(res_all.class) NE 'QSO' AND fracnsigma[4,*,*] GT 0.05) $
    OR (strtrim(res_all.class) EQ 'QSO' AND fracnsighi[4,*,*] GT 0.05)
   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'MANY_OUTLIERS')

   ; Warning: Redshift-error warning flag set to -1, which means that
   ; the chi^2 minimum was at the edge of the redshift-fitting range.
;   qflag = res_all.z_err EQ -1
;   zwarning = zwarning OR qflag * sdss_flagval('ZWARNING', 'Z_FITLIMIT')

   ; Warning: For QSOs, if C_IV, CIII], Mg_II, H_beta or H_alpha are negative
   ; and have at least a few pixels on each side of the fit (LINENPIXLEFT >= 4,
   ; LINENPIXRIGHT >= 4, and DOF >= 4).  Must be at least 3-sigma negative.
   for iobj=0, nobj-1 do begin
      if (strtrim(res_all[0,iobj].class,2) EQ 'QSO') then begin
         indx = where(zline.fiberid EQ res_all[0,iobj].fiberid $
          AND (strmatch(zline.linename, 'C_IV 1549*') $
            OR strmatch(zline.linename, 'C_III] 1908*') $
            OR strmatch(zline.linename, 'Mg_II 2799*') $
            OR strmatch(zline.linename, 'H_beta*') $
            OR strmatch(zline.linename, 'H_alpha*')) )
         if (indx[0] NE -1) then begin
            qflag = total( $
             zline[indx].linearea + 3*zline[indx].linearea_err LT 0 $
             AND zline[indx].linearea_err GT 0 $
             AND zline[indx].linenpixleft GE 4 $
             AND zline[indx].linenpixright GE 4 $
             AND zline[indx].linedof GE 4) NE 0
            zwarning[0,iobj] = zwarning[0,iobj] $
             OR qflag * sdss_flagval('ZWARNING', 'NEGATIVE_EMISSION')
         endif
      endif
   endfor

   res_all.zwarning = zwarning

   ;----------
   ; Write the output files

   splog, 'Writing output files'
   sxaddpar, hdr, 'NAXIS', 0
   sxdelpar, hdr, 'NAXIS1'
   sxdelpar, hdr, 'NAXIS2'
   sxaddpar, hdr, 'EXTEND', 'T', after='NAXIS'
   sxaddpar, hdr, 'VERS1D', idlspec2d_version(), $
    'Version of idlspec2d for 1D reduction', after='VERSCOMB'
   spawn, 'uname -n', uname
   sxaddpar, hdr, 'UNAME', uname[0]
   ww = strsplit(uname[0], '.', /extract)
   if (ww[1<(n_elements(ww)-1)] EQ 'fnal') then return

;   zans = struct_addtags((res_all[0,*])[*], res_elodie)
   zans = (res_all[0,*])[*]
   mwrfits, 0, zbestfile, hdr, /create ; Retain the original header in first HDU
   mwrfits, zans, zbestfile
   mwrfits, synflux, zbestfile
   mwrfits, dispflux, zbestfile

   sxaddpar, hdr, 'DIMS0', nper, ' Number of fits per objects'
   sxaddpar, hdr, 'DIMS1', nobj, ' Number of objects'
   mwrfits, 0, zallfile, hdr, /create ; Retain the original header in first HDU
   mwrfits, res_all, zallfile

   if (keyword_set(debugfile)) then dfpsclose

;   ;----------
;   ; Generate final QA plots
;
;   splog, 'Generating QA plots';

;   if (keyword_set(plotfile)) then begin
;      cpbackup, plotfile
;      dfpsplot, plotfile, /color
;   endif;

;   plottitle = string(zans[0].plate, zans[0].mjd, $
;    format='("Flux-Calibration Errors Plate=", i4, " MJD=", i5)')
;   qaplot_fcalibvec, objloglam, objflux, objivar, synflux, plugmap, zans, $
;    plottitle=plottitle
;
;   if (keyword_set(plotfile)) then dfpsclose

   ;----------
   ; Close log file

   splog, 'Total time for SPREDUCE1D = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE1D at ' + systime()
   if (keyword_set(logfile)) then splog, /close

   return
end
;------------------------------------------------------------------------------
