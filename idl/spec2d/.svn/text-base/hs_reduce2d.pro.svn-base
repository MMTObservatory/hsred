;+
; NAME:
;   hs_reduce2d
;
; PURPOSE:
;   Run the 2d reduction pipeline
;
; CALLING SEQUENCE:
;  hs_reduce2d, plugplan=, write=, checkaverage=, docosmic=, superfil=, 
;
; INPUTS:
;  Uses the cal.list file to do the reductions
;
; OPTIONAL KEYWORDS:
;  plugplan - translation file for the reductions
;  checkaverage - same as in hs_extract
;  docosmic - same as in hs_extract
;  superfile - same as in hs_extract
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA

;------------------------------------------------------------------------------
PRO hs_reduce2d, path=path, rerun=rerun, outdir=outdir, $
                 plugplan=plugplan, write=write, $
                 checkaverage=checkaverage, docosmic=docosmic,$
                 superfile=superfile
  
  
  if NOT keyword_set(path) then path = cwd()
  if NOT keyword_set(rerun) then rerun = 000
  
  rerun = string(rerun, format='(i4.4)')
  
  ;;Check to see if the calibration directory exists
  
  if direxist(path+'calibration') eq 0L then spawn,$
     'mkdir ' + path+'calibration'
  if direxist(path+'calibration/' + rerun) EQ 0L then $
     spawn, 'mkdir ' + path+'calibration/' + rerun
  
  calpath = path+'calibration/' + rerun+'/'
  if direxist(calpath) eq 0 then spawn, 'mkdir ' + calpath
  listpath = path + 'lists/'
  
  
  hs_calibproc, listpath=listpath, calpath=calpath, $
     mjd=mjd, rerun=rerun, /doall
  
  ;;Now run the extractor on each set of spectra in the cal file.
  
  
  if direxist(path + 'reduction') EQ 0 then spawn, $
     'mkdir ' + path + 'reduction'
  if direxist(path + 'reduction/'+rerun) EQ 0 then $
     spawn, 'mkdir ' + path + 'reduction/' + rerun
  
  redpath =   path+'reduction/' + rerun+'/'
  
  readcol, listpath + 'cal.list', $
     filelist, plugmap1, format='A,A', comment='#', $
     delimiter=' '
  
  ifile = filelist
  for i = 0, n_elements(ifile) -1 do begin
     
     filelist = strsplit(ifile(i), ',', /extract)
     
     plugmap = path + 'plugdir/' + plugmap1(i)
     
     if keyword_set(plugplan) then begin
        readcol, plugplan, dir, plug, field, pass, format='A,A,I,I'
        
        for k = 0, n_elements(plug)  -1 do begin
           if strmatch(plugmap, '*'+plug(k)) eq 1l then match = k
        endfor
     endif 
     
     out = 'spHect-' + string(pass(match),format='(i1.1)') $
        + string(field(match),format='(i2.2)') + $
        '-' + rerun + '.fits'
     
     print, out
     
     hs_extract, filelist, plugfile=plugmap, /uberextract, $
        /writeout, /writeall,$
        /combine, /plugcat,  outname=out , $
        rerun=rerun, plateid=mjd, $
        /qaplot, checkaverage=checkaverage, docosmic=docosmic, $
        superfile=superfile
  endfor
  
END

 
 
