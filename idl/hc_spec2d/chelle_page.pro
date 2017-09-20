;+
; NAME:
;   hs_page_file
;
; PURPOSE:
;   Examine reduced hs data
;
; CALLING SEQUENCE:
;  hs_page_file, filename, nstart=, plugmap=, ivar=, nsmooth=, xrange=, 
;
; INPUTS:
;  filename 6   - filename for for the Hectospec data
;
; OPTIONAL KEYWORDS:
;   nstart - starting aperture
;   nsmooth - number of pixels to smooth over
;   xrange - plotting xrange
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
PRO chelle_page, filename, nstart=nstart, plugmap=plugmap, ivar=ivar, $
                  zmap=zmap, autoline=autoline, smin = smin, synth=synth, $
                  nsmooth=nsmooth, psym=psym,title=title, xrange=xrange, $
                  clipfact=clipfact, dosnr=dosnr, pseudoflux=pseudoflux, $
                  synthmag=synthmag
  
   if not keyword_set(psym) then psym=10
   if not keyword_set(smin) then smin = 3650
   if NOT keyword_set(nstart) then nstart = 0
   if NOT keyword_set(clipfact) then clipfact =0.6
  
   ;;Read the file
   lam = mrdfits(filename, 0)
   object = mrdfits(filename, 1)
   plugmap = mrdfits(filename, 5)
     

     
   test = ''
   splog, 'Commands'
   splog, 'N - next'
   splog, 'P - previous'
   splog, 'Stop - stop'
   splog, 'C - clear lines'
  
   i = nstart
   stop = 'go'
   if NOT keyword_set(nsmooth) then nsmooth=0
   while stop ne 'stop' do begin
      if keyword_set(title) then title1=title(i)
      if not keyword_set(title) then title1=''
      ymin = min(smooth(object[*,i], nsmooth))
      ymax = max(smooth(object[*,i], nsmooth))        
      xrange= [min(lam[*,i]), max(lam[*,i])]
      
      djs_plot, lam[*,i], smooth(object[*,i],nsmooth),$
         /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
         xrange=xrange,/nstyle, /ystyle
      djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
         ymin + (ymax-ymin)*0.95, 'Aperture  ' + strn([i]), $
         color='yellow', thick=2, charsize=2, /isolatin1
      
      if keyword_set(plugmap) then begin
         
         djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
            ymin + (ymax-ymin)*0.90  ,  plugmap[i].objtype, $
            color='yellow', thick=2, charsize=2, /isolatin1
         djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.05,$
            ymin + (ymax-ymin)*0.85, $
            'snr = ' + strn(plugmap(i).snr), $
            color='yellow', thick=2, charsize=2, /isolatin1
      endif
     
     loop = 0 
     while loop eq 0 do begin
        read, prompt='Command: ', test
        
        if strupcase(test) eq 'N' then begin
           i = i+ 1
           loop = 1
        endif
        
        if strupcase(test) eq 'P' then begin
           i = i - 1
           loop = 1
        endif
        
        IF strupcase(test) EQ 'A' THEN BEGIN
           
           lines=[3933.66, 3968.469]
           splog, 'Enter the Velocity you would like plotted'
           read, prompt='Velocity : ', v
           c = astroconst()
           c = c.c / 1e5
           FOR ii = 0, n_elements(lines) -1 DO BEGIN
              djs_oplot, lines(ii)*(1+v/c)*[1,1], [-1e4, 1e8], $
                 color='red', thick=2
           ENDFOR
           
        ENDIF
        
           
           

        
        if strupcase(test) eq 'C' then begin
           i = i 
           loop = 1
        endif
        
        
        if strupcase(test) eq 'STOP' then begin
           loop = 1
           stop = 'stop'
        endif
        
     endwhile
     
  endwhile
  
END

