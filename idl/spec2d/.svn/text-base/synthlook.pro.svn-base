PRO synthlook, filename
  
   dfpsplot, '/tmp/spectra.ps', /landscape, /color
  !P.multi=[0,1,4]
  
  for ifile=0, n_elements(filename)-1 do begin
     print, filename(ifile)
     
;;Read the file
  lam = mrdfits(filename(ifile), 0)
  object = mrdfits(filename(ifile), 1)
  ivar = mrdfits(filename(ifile), 2)
  mask = mrdfits(filename(ifile), 4)*0.0
  plugmap = mrdfits(filename(ifile), 5)
  
  k = [42,44]
  
  ;;Some parameters
  xrange=[3800,8500]
  clipfact = 0.6
  nsmooth=11
  

    
  for j = 0, n_elements(k) -1 do begin
     i = k(j)
     
     ;;Mask pixels were the ivar is lt 0.6 the smoothed value
        temp = ivar[*,i] / smooth(ivar[*,i], 21)
        newmask = temp lt clipfact
     
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
        djs_plot, lam[pixels,i], smooth(object[pixels,i],nsmooth),$
       /xstyle , yrange=[ymin, ymax],  ps=psym, title=filename(ifile), $
        xrange=xrange,/nstyle, /ystyle, thick=2
        djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
       ymin + (ymax-ymin)*0.90, 'Aperture  ' + strn([i]), $
        color='red', thick=2, charsize=1, /isolatin1
        djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
          ymin + (ymax-ymin)*0.80  ,  plugmap[i].objtype, $
           color='red', thick=2, charsize=1, /isolatin1
     endfor
  endfor
  
     dfpsclose
     

   
END


