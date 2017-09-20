;+
; NAME:
;   hs_pipeline_wrap
;
; PURPOSE:
;  Wrapper script for HSRED v1.2 pipeline. Runs hs_preroc to generate
;  file lists to run on, then runs hs_calibproc on biases, flats, and
;  arcs, then hs_batch_extract to extract and reduce the science frames
;
; CALLING SEQUENCE:
;  hs_pipeline_wrap, dostand=dostand, uberextract=uberextract, $
;     do600=do600, rerun=rerun, docosmic=docosmic, sky2d=sky2d, $
;     doredleak=doredleak, customskies=customskies, noftweak=noftweak
;
; REQUIRED KEYWORDS:
;  rerun
;
; OPTIONAL KEYWORDS:
;   rerun       - a 4-digit number or string to tag the current run of
;                 the pipeline. If not specified, defaults to '0100'
;   dostand     - Set to perform a standard reduction with red leak
;                 removal, CR-rejection using HCOSMIC, new enhanced
;                 sky-subtraction (/sky2d), and simple weighted
;                 coaddition of exposures. For 270-line data this also
;                 includes telluric absorption correction, but this
;                 option is not yet available for 600-line reductions
;   uberextract - Set to perform a reduction including flux
;                 calibration, if you have included F-stars among your
;                 targets. The file $HSRED_DIR/etc/standardstars.dat
;                 must first be updated with the relevant info for your
;                 standard stars. /uberextract also turns on all the
;                 same options as /dostand, except that combining
;                 exposures is done using the old method, and the
;                 resulting outputs are in flux-calibrated units, not
;                 coadded counts.
;   do600       - Set this keyword if you are reducing
;                 data taken with the 600-line grating 
;   qaplot      - Setting this results in the quality assurance plots to be
;                 be created
;   docosmic    - If docosmic is set then the cosmic ray rejection is turned 
;                 on. docosmic=2 uses HCOSMIC if > 1 image present
;                 (You want this). docosmic=1 uses old method
;   noftweak - enables new, simpler coadding of exposures, without any
;              scaling of fluxes between exposures. Cannot be set with 
;              /uberextract.
;   doredleak - set to remove contaminating red flux from robot
;               positioner LEDs, longward of ~8500A. 
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
;      ./lists/, ./reduction/YYYY/ and ./calibration/YYYY 
;       directories where YYYY is the rerun
;
; REVISION HISTORY:
;   April 2014 - Written by S. Moran         
;                
;-
;------------------------------------------------------------------------------
pro hs_pipeline_wrap, rerun=rerun, sky2d=sky2d, docosmic=docosmic, qaplot=qaplot, dodark=dodark, doredleak=doredleak, dostand=dostand, customskies=customskies, nonlinear=nonlinear, do600=do600, noftweak=noftweak, uberextract=uberextract

if keyword_set(dostand) then begin
  if not keyword_set(do600) then $
hs_pipeline_wrap, qaplot=qaplot, /sky2d,/doredleak, docosmic=2, rerun=rerun, customskies=customskies, nonlinear=nonlinear, /noftweak else hs_pipeline_wrap, qaplot=qaplot, /sky2d,doredleak=doredleak, docosmic=2, rerun=rerun, customskies=customskies, nonlinear=nonlinear, /do600, /noftweak
return
endif

if keyword_set(uberextract) then begin
  if not keyword_set(do600) then $
hs_pipeline_wrap, qaplot=qaplot, /sky2d,/doredleak, docosmic=2, rerun=rerun, customskies=customskies, nonlinear=nonlinear, noftweak=-1 else hs_pipeline_wrap, qaplot=qaplot, /sky2d,doredleak=doredleak, docosmic=2, rerun=rerun, customskies=customskies, nonlinear=nonlinear, /do600, noftweak=-1
return
endif

if not keyword_set(rerun) then rerun=100 ;change this ID number to save multiple reductions for the same data
if not keyword_set(sky2d) then sky2d=0 
if not keyword_set(do600) then do600=0 
if not keyword_set(noftweak) then noftweak=1
if noftweak eq -1 then noftweak=0
if not keyword_set(doredleak) then doredleak=0 
if not keyword_set(qaplot) then qaplot=0 
if not keyword_set(docosmic) then docosmic=0 
tweaksky=1
if not keyword_set(dodark) then dodark=0
if not keyword_set(customskies) then customskies=0
mmtcam=0 
if qaplot eq 1 then set_plot, 'ps' else set_plot, 'x'
if qaplot eq 1 then begin
  device, file='waveplot.ps', xsize=8.5, ysize=11, /inches
  !P.multi=[2,1,2]
endif
;now make lists of files to run on
hs_preproc ;, /overwrite ;can add ,/overwrite

;process calibrations for this night
hs_calibproc, /doall, qaplot=qaplot, rerun=rerun ;/dodark ;can change rerun number for multiple versions of the reduction
if qaplot eq 1 then device, /close


if keyword_set(do600) then begin
hs_batch_extract, /fiberflat, /combine, /skysubtract, tweaksky=tweaksky, rerun=rerun, docosmic=docosmic, qaplot=qaplot, dummy=0, dofringes=0, tellcorr=0, noftweak=noftweak, sky2d=sky2d, dodark=dodark,doredleak=doredleak, mmtcam=mmtcam, customskies=customskies
endif else begin
hs_batch_extract, /fiberflat, /combine, /skysubtract, tweaksky=tweaksky, rerun=rerun, docosmic=docosmic, qaplot=qaplot, dummy=0, dofringes=0, /tellcorr, noftweak=noftweak, sky2d=sky2d, dodark=dodark,doredleak=doredleak, mmtcam=mmtcam, customskies=customskies
endelse

if size(rerun, /tname) ne 'string' then rerun = string(rerun, format='(i4.4)')
cd, 'reduction/'+rerun
hectfiles=findfile('spHect*fits')
for j=0, n_elements(hectfiles)-1 do begin
obj=strmid(hectfiles[j],7,(strlen(hectfiles[j])-7-15))
outfile=obj+'_skysub.ms.fits'
if keyword_set(do600) then hs_toiraf, hectfiles[j], outfile, dx=0.55, /multi, /hsplit else hs_toiraf, hectfiles[j], outfile, wavemin=3700., wavemax=9150., dx=1.1999119, /multi, /hsplit
if keyword_set(nonlinear) then begin
    if keyword_set(do600) then hs_toiraf, hectfiles[j], outfile, dx=0.55, /multi, hsplit=2 else hs_toiraf, hectfiles[j], outfile, wavemin=3700., wavemax=9150., dx=1.1999119, /multi, hsplit=2
    spawn, 'mkdir skysub_nonlin_'+obj+'/skydir'
    spawn, 'mv skysub_nonlin_'+obj+'/*sky*fits skysub_nonlin_'+obj+'/skydir/'
    spawn, 'rm -rf skysub_nonlin_'+obj+'/*unused*fits'
    spawn, 'rm -rf skysub_nonlin_'+obj+'/*rejec*fits'
endif
endfor
cd, '../..'
;optionally, can run the 1D pipeline for measuring redshifts
;files = findfile('reduction/'+rerun+'/spHect*fits')
;hs_reduce1d, files, /pseudoflux
!P.multi=0
set_plot, 'x'

end
