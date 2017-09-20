;;x_specobjcr_grow
;;shamelessly copied from the XIDL routine x_specobjcr by X. Prochaska
;;et al

pro x_specobjcr_grow, grow, crmap=crmap

  ;; Announce
  print, 'Growing CR hits by ', grow

  ;;
  sz = size(crmap, /dimens)
  cridx = where(crmap LT 0., ncr)
  if ncr EQ 0 then return
  xcr = cridx mod sz[0]
  ycr = cridx / sz[0]

  ximg = lindgen(sz[0]) # replicate(1L,sz[1])
  yimg = replicate(1L,sz[1]) # lindgen(sz[1])
  x0 = (xcr - grow) > 0
  x1 = (xcr + grow) < (sz[0]-1)
  y0 = (ycr - grow) > 0
  y1 = (ycr + grow) < (sz[1]-1)

  ;; Loop
  for qq=0L,ncr-1 do begin
      crmap[x0[qq]:x1[qq],y0[qq]:y1[qq]] = -99.
  endfor
return
end
 
;;
;PROCEDURE: HCOSMIC
;Routine to remove cosmic rays modelled after hcosmic.pl/.cl by N
;Caldwell and J. Mink
;
; DESCRIPTION - pass a list of files, or NX x NY x NIMAGES data
;               cube, to detect and remove cosmic rays from a set of
;               several un-extracted spectra, with the dispersion
;               direction running along the y-axis
;
pro hcosmic, files=files, images=images, ivar=ivar, minflux=minflux, limit=limit, grow=grow, mask=mask

if not keyword_set(limit) then limit=20.
if not keyword_set(minflux) then minflux=20.
if not keyword_set(grow) then grow=1
if (not keyword_set(files)) and (not keyword_set(images)) then begin
   print, 'HCOSMIC: No images passed. Must set IMAGES keyword with an array of [NX,NY,NIM], or else a string array containing file names'
   return
endif
;If we've been passed a string array, check and load all the files
if keyword_set(files) then begin
  for i=0, n_elements(files)-1 do begin
     exist=file_test(files[i])
     if exist then begin
        tmpim=mrdfits(files[i],0,hdr) 
        tmpivar=mrdfits(files[i],1) 
     endif else begin
        print, 'File: '+files[i]+' not found; Please check path'
        return
     endelse
     if i eq 0 then begin 
       dimens=size(tmpim)
       nx=dimens[1]
       ny=dimens[2]
       nim=n_elements(files)
       images_new=fltarr(nx,ny,nim)
       ivar_new=fltarr(nx,ny,nim)
     endif     
     images_new[*,*,i]=tmpim
     ivar_new[*,*,i]=tmpivar
  endfor
  images=images_new
  ivar=ivar_new
endif else begin
  dimens=size(images)
  nx=dimens[1]
  ny=dimens[2]
  nim=dimens[3]
endelse
;construct the median and the statistics image derived from it
if nim lt 2 then begin
   print, 'Only one valid image; use LACosmic instead'
   return
endif

medim=((sqrt(abs(median(images,dimen=3)))*1.4) > 5.)
cr_cnt=fltarr(nim)

;now loop through one 3x20pixel box at a time
for ix=0, nx-1, 3 do begin
     ix2=(ix+2)<(nx-1)
     for iy=0,ny-1,20 do begin
       iy2=(iy+19)<(ny-1)
       imsec1=images[ix:ix2,iy:iy2,0]
       imsec2=images[ix:ix2,iy:iy2,1]
       ivarsec1=ivar[ix:ix2,iy:iy2,0]
       ivarsec2=ivar[ix:ix2,iy:iy2,1]
       medsec=medim[ix:ix2,iy:iy2] 
       nkeep=n_elements(imsec1)-5
       stat1=total(imsec1[(sort(imsec1))[0:nkeep]])/float(nkeep)
       stat2=total(imsec2[(sort(imsec2))[0:nkeep]])/float(nkeep)
       ;only calculate ratio between images if counts are above minflux
       if (stat1 gt minflux) or (stat2 gt minflux) then ratio=stat1/stat2 else ratio=1.
       ;now find cosmics according to the formula pix1-pix2*ratio>median*limit
       cosmic1=where((imsec1-imsec2*ratio) gt (limit*medsec),cr1_cnt)
       cosmic2=where((imsec2-imsec1/ratio) gt (limit*medsec),cr2_cnt)
       cr_cnt[0]=cr_cnt[0]+cr1_cnt
       cr_cnt[1]=cr_cnt[1]+cr2_cnt
       ;put off growing and interpolating over CRs until whole image has had CRs detected
        if cr1_cnt ge 1 then begin
          imsec1[cosmic1]=0.0
          ivarsec1[cosmic1]=-99.0
          images[ix:ix2,iy:iy2,0]=imsec1
          ivar[ix:ix2,iy:iy2,0]=ivarsec1
        endif
        if cr2_cnt ge 1 then begin
          imsec2[cosmic2]=0.0
          ivarsec2[cosmic2]=-99.0
          images[ix:ix2,iy:iy2,1]=imsec2
          ivar[ix:ix2,iy:iy2,1]=ivarsec2
        endif
        if nim gt 2 then begin
           ;temporarily fill in bad parts of imsec1 with corresponding pixels from imsec2
           if cr1_cnt ge 1 then $
              imsec1[cosmic1]=imsec2[cosmic1]*ratio
           stat1=total(imsec1[(sort(imsec1))[0:nkeep]])/float(nkeep)
           ;then loop through the rest of the images
           for i=2, nim-1 do begin
              imsec2=images[ix:ix2,iy:iy2,i]
              ivarsec2=ivar[ix:ix2,iy:iy2,i]
              stat2=total(imsec2[(sort(imsec2))[0:nkeep]])/float(nkeep)
              if (stat1 gt minflux) or (stat2 gt minflux) then ratio=stat1/stat2 else ratio=1.
              cosmic2=where((imsec2-imsec1/ratio) gt (limit*medsec),cr2_cnt)
              cr_cnt[i]=cr_cnt[i]+cr2_cnt
              if cr2_cnt ge 1 then begin
                imsec2[cosmic2]=0.0
                ivarsec2[cosmic2]=-99.0  
                images[ix:ix2,iy:iy2,i]=imsec2
                ivar[ix:ix2,iy:iy2,i]=ivarsec2
              endif
           endfor
        endif
     endfor    
endfor
;now grow CR mask and linearly interpolate over CR pixels
for i=0, nim-1 do begin
  print, 'HCOSMIC: identified '+strn(cr_cnt[i])+' cosmic rays in image #'+strn(i+1)
  ivar1=reform(ivar[*,*,i])
  x_specobjcr_grow, grow, crmap=ivar1
  ivar[*,*,i]=ivar1
endfor
mask=(ivar lt -98.)
;stop
;images=djs_maskinterp(images,mask,iaxis=1, /spline)
ivar[where(mask)]=0.0
;stop
return
end


