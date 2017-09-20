;+
; NAME:
;   hs_batch_extract
;
; PURPOSE:
;   Reduce all of the files as layd out in the cal-YYYY.list
;   file in the list directory
;
; CALLING SEQUENCE:
; hs_batch_extract, path=path, root=root,  skysubtract=skysubtract, $
;                      fiberflat=fiberflat,  tweaksky=tweaksky, $
;                      plugcat=plugcat, combine=combine, writeout=writeout,$
;                      lamout=lamout, fluxout=fluxout, fluxcorr=fluxcorr, $
;                      tellcorr=tellcorr, uberextract=uberextract, $
;                      ivarout=ivarout, rerun=rerun,  outname=outname,  $
;                      writeall=writeall,  plateid=plateid, sky2d=sky2d, $
;                      qaplot=qaplot, quicklook=quicklook, quickplot=quickplot,$
;                      stand=stand,  checkaverage=checkaverage, $
;                      customskies=customskies, docosmic=docosmic, $
;                      superfile=superfile,doredleak=doredleak,$
;                      noftweak=noftweak, nobias=nobias, dodark=dodark, $
;                      mmtcam=mmtcam, dostand=dostand, oldmap=oldmap, $
;                      dummy=dummy, dofringes=dofringes
;
; INPUTS:
;   lists/cal-YYYY.list - this file is automatically used
;                          as a reduction plan
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
;   
;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA         
;   Jan 2013 - smoran; fixed bug by replacing reference to 'oldcat' with 'oldmap'             
;   June 2014 - Updated lists of keywords to pass the right things on
;               to hs_extract.pro 
;-
;------------------------------------------------------------------------------
PRO hs_batch_extract, path=path, root=root,  skysubtract=skysubtract, $
                      fiberflat=fiberflat,  tweaksky=tweaksky, $
                      plugcat=plugcat, combine=combine, writeout=writeout,$
                      lamout=lamout, fluxout=fluxout, fluxcorr=fluxcorr, $
                      tellcorr=tellcorr, uberextract=uberextract, $
                      ivarout=ivarout, rerun=rerun,  outname=outname,  $
                      writeall=writeall,  plateid=plateid, sky2d=sky2d, $
                      qaplot=qaplot, quicklook=quicklook, quickplot=quickplot,$
                      stand=stand,  checkaverage=checkaverage, $
                      customskies=customskies, docosmic=docosmic, $
                      superfile=superfile,doredleak=doredleak,$
                      noftweak=noftweak, nobias=nobias, dodark=dodark, $
                      mmtcam=mmtcam, dostand=dostand, oldmap=oldmap, $
                      dummy=dummy, dofringes=dofringes
  
  if NOT keyword_set(path) then path=cwd()
  if NOT keyword_set(listpath) then listpath = path+'lists/'
  if NOT keyword_set(rerun) then rerun = string(0, format='(i4.4)')
  if NOT keyword_set(calpath) then calpath = path + $
     'calibration/' + rerun + '/'
  
  if size(rerun, /tname) ne 'string' then $
     rerun = string(rerun, format='(i4.4)')
  
  
  ;;Check to see if the calibration directory exists

  if direxist(path+'calibration') eq 0L then $
     spawn, 'mkdir ' + path+'calibration'
  if direxist(path+'calibration/' + rerun) EQ 0L then $
     spawn, 'mkdir ' + path+'calibration/' + rerun

  if keyword_set(outname) then noname =1
  if not keyword_set(outname) then noname =0
  

;;Create the biasfiles
  
  biasfile = findfile(listpath+'bias*')
  junk = strsplit(listpath, ' ', length = pathcount)
  mjd = (strmid(biasfile, pathcount+5, 5))[0]
  root = mjd
  
  readcol, listpath + 'cal.list', filelist,$
     plugmap1, format='A,A', comment='#', $
     delimiter=' '
  

  
   ifile = filelist
   plugmap1 = 'plugdir/' + plugmap1
   for i = 0, n_elements(ifile) -1 do begin
    filelist = strsplit(ifile(i), ',', /extract)
    file = strsplit(filelist(0), '/', /extract)
    file = file(n_elements(file)-1)
    junk = strsplit(file, ' ', length=ct)
    if direxist('reduction') eq 0L then spawn, 'mkdir reduction'
    if direxist('reduction/' + rerun) eq 0L then $
       spawn, 'mkdir reduction/' + rerun
    junk = strsplit(file, '.', /extract)
    if junk[n_elements(junk)-1] eq 'gz' then begin
       root = strmid(file, 0, ct-8)
    endif else begin
       root = strmid(file, 0, ct-5)
    endelse


    
    if noname eq 0 then outname = 'spHect-' + root + $
       '-' + rerun + '.fits'
        
    hs_extract, filelist, plugfile=plugmap1(i), skysubtract=skysubtract, $
       fiberflat=fiberflat,  tweaksky=tweaksky, sky2d=sky2d, $
       plugcat=plugcat, combine=combine, writeout=writeout, lamout=lamout, $
       fluxout=fluxout, fluxcorr=fluxcorr, tellcorr=tellcorr, $
       uberextract=uberextract, ivarout=ivarout, rerun=rerun, customskies=customskies, $
       outname=outname,  writeall=writeall, plateid=plateid, mmtcam=mmtcam, $
       qaplot=qaplot, quicklook=quicklook, quickplot=quickplot, stand=stand, $
       checkaverage=checkaverage, docosmic=docosmic, superfile=superfile, $
       noftweak=noftweak, nobias=nobias, dodark=dodark, dostand=dostand, $
       oldmap=oldmap, dummy=dummy, dofringes=dofringes,doredleak=doredleak
 endfor
 
  
END
