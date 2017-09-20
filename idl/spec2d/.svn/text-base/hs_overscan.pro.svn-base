;+
; NAME:
;   hs_overscan
;
; PURPOSE:
;   Overscan correct a hectospec frame
;
; CALLING SEQUENCE:
;  hs_overscan, image, head, overscan=overscan, 
;              checkoverscan=checkoverscan, $
;              norder_overscan=norder_overscan
;
; INPUTS:
;   image - input image
;   head  - header for the input file
;   
;
; OPTIONAL KEYWORDS:
;   overscan - overscan region for the data
;   
;
; OUTPUTS:
;   
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
;   
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in ISPEC
;   
;-
;------------------------------------------------------------------------------
PRO hs_overscan, image, head, overscan=overscan, checkoverscan=checkoverscan, $
                 norder_overscan=norder_overscan
  
  origimage = image             ;keep the original image
  
  imsize = size(origimage, /dimension)
  ncols = imsize[0]
  nrows = imsize[1]
  
  ;;pad the defaul overscan region by two columns
  if keyword_set(overscan) then overscan = long(overscan) else begin
     
     overscan = long(strsplit(sxpar(head, 'BIASSEC'), $
                              '[:,]', /EXTRACT))-1l
     
     overscan[0] = overscan[0] + 2L
     overscan[1] = overscan[1] - 2L
     
  endelse
  
  x1 = overscan[0] > 0L
  x2 = overscan[1] < (ncols-1L)
  
  y1 = overscan[2] > 0L
  y2 = overscan[3] < (imsize[1] - 1L)
  
  colaxis = findgen(ncols)
  rowaxis = findgen(nrows)
  
  Imtype = strlowcase(strcompress((sxpar(head,'imagetyp')),/remove))
  object = strlowcase(strcompress((sxpar(head,'object')),/remove))
  if n_elements(norder_overscan) eq 0L then norder_overscan = 2
  
  
;;Median of the overscan along perpendicular to readout. 
  
  ny = n_elements(rowaxis[y1:y2])
  nx = n_elements(colaxis[x1:x2])
  
  scanvector = djs_median(origimage[x1:x2,y1:y2],1)
  
  scanerror = scanvector*0.0
  for i = 0L, ny-1L do scanerror[i] = stddev(origimage[x1:x2,y1+i])
  
  fnorder = norder_overscan
  coef = func_fit(rowaxis[y1:y2],scanvector,fnorder+1, $
                  function_name='flegendre',yfit=fit)
  
  yrange = [min(scanvector-scanerror),max(scanvector+scanerror)]$
     *[0.999,1.001]   
  
  oscanstr = '['+strjoin(strcompress(overscan[0:1],/remove),',')+':'+$
     strjoin(strcompress(overscan[2:3],/remove),',')+']'
  
  if keyword_set(checkoverscan) and (not keyword_set(qflag)) then begin
     
     
     if !d.window ne 3L then window, 3, xs=600, ys=600
     plot, rowaxis[y1:y2], scanvector, xsty=3, ysty=3, ps=10, charsize=1.5, $
        charthick=2.0, xtitle='Row Number', ytitle='Counts', yrange=yrange, $
        title=strupcase(object)+' ('+strupcase(imtype)+')'
     oplot, rowaxis[y1:y2], fit, thick=3.0
     oplot, rowaxis[y1:y2], scanvector+scanerror, line=2, thick=3.0
     oplot, rowaxis[y1:y2], scanvector-scanerror, line=2, thick=3.0
     cc = get_kbrd(1)
     if strupcase(cc) eq 'Q' then qflag = 1L
     
  endif
  
  residuals = 100.0*(scanvector-fit)/fit
  
  ;; subtract the overscan fit
  image = origimage*0.0
  ;; subtract the overscan fit
  
  image = origimage*0.0
  for k = 0L, ncols-1L do image[k,y1:y2] = origimage[k,y1:y2] - fit
      
  ;; update the head, remembering that IRAF indices start at one
  
  sxdelpar, head, 'BIASSEC'
  sxaddpar, head, 'OVERSCAN', im_today()+' Overscan section '+$
     '['+strn(x1+1)+':'+strn(x2+1)+','+strn(y1+1)+':'+$
     strn(y2+1)+'] with mean=+ ' + $
     strn(mean(scanvector),format='(F10.3)'), before='HISTORY'
  
  
  return
END

  
  
