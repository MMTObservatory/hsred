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

PRO hs_page_file, filename, nstart=nstart, plugmap=plugmap, ivar=ivar, $
                  zmap=zmap, autoline=autoline, smin = smin, synth=synth, $
                  nsmooth=nsmooth, psym=psym,title=title, xrange=xrange, $
                  clipfact=clipfact, dosnr=dosnr, pseudoflux=pseudoflux, $
                  synthmag=synthmag
  
  if not keyword_set(xrange) then xrange=[3800,8500]
  if not keyword_set(psym) then psym=0
  if not keyword_set(smin) then smin = 3650
  if NOT keyword_set(nstart) then nstart = 0
  if NOT keyword_set(clipfact) then clipfact =0.6
  
  ;;Read the file
  lam = mrdfits(filename, 0)
  object = mrdfits(filename, 1)
  ivar = mrdfits(filename, 2)
  mask = mrdfits(filename, 4)*0.0
  plugmap = mrdfits(filename, 5)
     
     
  if keyword_set(pseudoflux) then begin
     fcalfile = getenv('HSRED_DIR') + '/etc/average_flux.fits'
     tcalfile = getenv('HSRED_DIR') + '/etc/average_tell.fits'
     
     fcal = mrdfits(fcalfile, 1)
     tcal = mrdfits(tcalfile, 1)
     flux_factor = bspline_valu(alog10(lam), fcal)
     tell_factor = bspline_valu((lam), tcal)
     divideflat, object, flux_factor, invvar=ivar
     divideflat, object, tell_factor, invvar=ivar
  endif
  
  
  test = ''
  splog, 'Commands'
  splog, 'N - next'
  splog, 'P - previous'
  splog, 'Stop - stop'
  splog, 'A - plot  absorption lines'
  splog, 'E - plot  emission lines'
  splog, 'C - clear lines'
  
  i = nstart
  stop = 'go'
  
  if keyword_set(dosnr) then hs_snr, lam, object, plugmap, snr=snr
  if not keyword_set(dosnr) then snr = fltarr(n_elements(object(0,*)))
  
  if NOT keyword_set(nsmooth) then nsmooth=0
  while stop ne 'stop' do begin
     
     ;;Mask pixels were the ivar is lt 0.6 the smoothed value
     temp = ivar[*,i] / smooth(ivar[*,i], 21)
     newmask = temp lt clipfact
     
     junk = min(abs(lam[*,i]-5588),pix)
     mask[(pix-10):(pix+10),i]=mask[(pix-10):(pix+10),i]+$
        newmask[(pix-10):(pix+10)]
     
     junk = min(abs(lam[*,i]-5577), pix)
     
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     junk = min(abs(lam[*,i]-6300), pix)
     
     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     
     junk = min(abs(lam[*,i]-6363), pix)

     mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
     
     pixels = where(mask[*,i] eq 0)  
     
   

     l1 =  min(abs(lam[pixels,i]-xrange(0)),l)
     m1 =  min(abs(lam[pixels,i]-xrange(1)),m)
     
     l = l(0)
     m = m(0)
     
     ymax = max(smooth(object[pixels(l):pixels(m),i],nsmooth)) *1.1
     ymin = min(smooth(object[pixels(l):pixels(m),i],nsmooth))
     diff = ymax-ymin    
     
     if keyword_set(title) then title1=title(i)
     if not keyword_set(title) then title1=''
     djs_plot, lam[pixels,i], smooth(object[pixels,i],nsmooth),$
        /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
        xrange=xrange,/nstyle, /ystyle, thick=2
     djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
        ymin + (ymax-ymin)*0.95, 'Aperture  ' + strn([i]), $
        color='yellow', thick=2, charsize=2, /isolatin1
     
     if keyword_set(plugmap) then begin
        
        djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
           ymin + (ymax-ymin)*0.90  ,  plugmap[i].objtype, $
           color='yellow', thick=2, charsize=2, /isolatin1
        djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.05,$
           ymin + (ymax-ymin)*0.85, $
           'snr = ' + strn(snr(i)), $
           color='yellow', thick=2, charsize=2, /isolatin1
        djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
           ymin + (ymax-ymin)*0.95, $
           'icode = ' + strn(plugmap[i].icode), $
           color='yellow', thick=2, charsize=2, /isolatin1
        djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
           ymin + (ymax-ymin)*0.9, $
           'ap mag = ' + strn(plugmap[i].rapmag), $
           color='yellow', thick=2, charsize=2, /isolatin1
        
        ;;Here, I am going to do something a bit strange and go ahead and
        ;;integrate the spectrum at the redshift observed.  This only works 
        ;;if the /synthmag is set
        
        if keyword_set(synthmag) then begin
           
           k_load_filters, 'sdss_r0.par', nlam, lam1, trans
           integrand = object[*,i]*1e-17
           intertrans = interpol(trans, lam1, lam[*,i])
           k = where(intertrans lt 0)
           intertrans(k) = 0
           rmag = int_tabulated(lam[*,i], integrand*intertrans)
           djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.45,$
              ymin + (ymax-ymin)*0.95, $
              'synth mag = ' + strn(-2.5*alog10(rmag)), $
              color='yellow', thick=2, charsize=2, /isolatin1
        endif
        
        
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
        
        if strupcase(test) eq 'C' then begin
           i = i 
           loop = 1
        endif
        
        
        if strupcase(test) eq 'STOP' then begin
           loop = 1
           stop = 'stop'
        endif
        
        if strupcase(test) eq 'A' then begin
           splog, 'Enter the Redshift you would like plotted'
           read, prompt='Redshift : ', z
           
           lines= [3797.90,3834.00,3889.00,3933.70,3968.50,4101.70,$
                   4304.40,4340.48,4861.34,5175.36,5268.98,5892.50,$
                   4172.00,4226.74]
           labels = ['H8', 'H7' , 'H6' ,  'K',  'H',  'H \delta',  $
                     'G',  'H \gamma',  'H \beta',  'Mg',  'CaFe', $
                     'Na',  'CaI',  'CaI']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='magenta', linestyle=2
              djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
                 labels(j), charsize=2, color='yellow'
           endfor
        endif
        
        if strupcase(test) eq 'E' then begin
           splog, 'Enter the Redshift you would like plotted'
           read, prompt='Redshift : ', z
           
           lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
                    1908.734, 2326.0, 2798.745, 3727.08, 4341.68, $
                    4862.68, 4960.3,5008.24, 6549.86,$
                    6564.61, 6585.27, 6718.29 ,6732.67]
           labels = ['Ly \beta', 'Ly \alpha', 'NV', 'Si V', 'C IV', $
                     'He II', 'C III', 'C II','MgII',  'O II', 'H \gamma', $
                     'H \beta', 'O III', 'O III', 'N II', 'H \alpha', $
                     'N II', 'S II', 'S II']
           for j = 0, n_elements(lines) -1 do begin
              djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
                 color='yellow', linestyle=2
              djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
                 labels(j), charsize=2, color='yellow'
           endfor 
        endif
        
     endwhile
     
  endwhile
  
END

