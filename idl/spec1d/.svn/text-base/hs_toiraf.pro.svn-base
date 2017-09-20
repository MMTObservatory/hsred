;;+
; NAME:
;   hs_toiraf
;
; PURPOSE:
;  Convert HSRED-format multi-extension fits files into linearized
;  spectra of various formats more readable by iraf or ds9
;
; CALLING SEQUENCE:
;       hs_toiraf, sphectfile, outname, wavemin=wavemin, dx=dx, 
;       wavemax=wavemax, simple=simple, multispec=multispec, hsplit=hsplit
;
; INPUTS:
;   sphectfile - file name of input spHect (or spObs-) fits file
;   outname    - name for output 2D fits file
; OPTIONAL KEYWORDS:
;   simple    - If set, writes out only the simpled 2D file with
;               NLAMBDAxNFIBERS dimensions
;   multispec - If set, instead writes out an NLAMBDAxNFIBERSx4
;               multispec data-cube, where the 4 spectra in the 3rd
;               dimension are 1) weighted-coadded spectrum 2) summed
;               spectrum, 3) sky spectrum, 4) sigma spectrum
;   hsplit    - If set with multispec, further writes out one file per
;               target spectrum, all stored in the skysub_<config>
;               subdirectory, in the same NLAMBDAx1x4 data cube format
;               If hsplit=2, will instead write everything in the iraf
;               .ec. format, where the spectra are not linearized, and
;               the wavelengths are stored as legandre coefficients in
;               the fits header.
;   wavemin,wavemax   - Min and max wavelengths for the output
;                       spectra. Defaults to minmax(flux)
;   dx        -       - pixel size in angstroms of the output,
;               defaults to 1.2, good for 270-line data, but must be
;               changed (to ~0.55) for 600-line data

Pro hs_toiraf,  sphectfile, outname, wavemin=wavemin, dx=dx, wavemax=wavemax, simple=simple, multispec=multispec, hsplit=hsplit ;hsplit=2 means write non-linear iraf-style .ec.fits files
  
  
  lam = mrdfits(sphectfile, 0, header)
  flux = mrdfits(sphectfile, 1, hdr1)
  ivar = mrdfits(sphectfile, 2)
  nfib = n_elements(lam[0,*])
  plugmap = mrdfits(sphectfile,5)
if keyword_set(multispec) then begin
   FITS_INFO, sphectfile,N_ext=n_ext, /silent
   if n_ext ge 6 then sky = mrdfits(sphectfile,6) else sky=flux*0.0
   if n_ext ge 7 then flux_box = mrdfits(sphectfile,7) else flux_box=flux
endif


  if not keyword_set(wavemin) then wavemin=max(lam[0,*])
  
  if wavemin lt  min(lam) then wavemin = min(lam)
  
  if not keyword_set(wavemax) then wavemax=(min(lam[n_elements(lam[*,0])-1,*]))
  if not keyword_set(dx) then dx=1.2
  
  npix = long( (wavemax-wavemin)/dx) + 1
  newlam = wavemin + dx * lindgen(npix)
  newflux = fltarr(npix, nfib)
  newivar = fltarr(npix, nfib)
  newsky=fltarr(npix, nfib)
  newflux_box=fltarr(npix, nfib)
  
  for i = 0, nfib-1 do begin
     splog, 'Resampling fiber' + string(i)
     
     if keyword_set(multispec) then begin
        goodpix=where(ivar[*,i] ne 0.0, ngood)
        if ngood gt 1 then begin
          hs_combine1fiber, lam[goodpix,i], flux[goodpix,i], ivar[goodpix,i], skyflux=sky[goodpix,i], $
          newlam=newlam, newflux=newflux1, newivar=newivar1, newsky=newsky1, nord=4
          hs_combine1fiber, lam[goodpix,i], flux_box[goodpix,i], ivar[goodpix,i], newlam=newlam,newflux=newflux2, nord=4
          newsky[*,i]=newsky1
          newflux_box[*,i]=newflux2
        endif
     endif else begin         
         goodpix=where(ivar[*,i] ne 0.0, ngood)
         if ngood gt 1 then begin
           hs_combine1fiber, lam[goodpix,i], flux[goodpix,i], ivar[goodpix,i], $
            newlam=newlam, newflux=newflux1, newivar=newivar1, nord=4
         endif
      endelse
     newflux[*,i] = newflux1
     newivar[*,i] = newivar1
  endfor
  
  ;;Add IRAF keywords
  if keyword_set(multispec) then sxaddpar, header, 'NAXIS', 3 else sxaddpar, header, 'NAXIS', 2
  sxaddpar, header, 'NAXIS1', npix
  sxaddpar, header, 'NAXIS2', nfib
  if keyword_set(multispec) then   sxaddpar, header, 'NAXIS3', 4
  sxaddpar, header, 'SIMPLE', "T", before='NAXIS'
if keyword_set(simple) then  sxaddpar, header, 'EXTEND', 'F', after='NAXIS2' else sxaddpar, header, 'EXTEND', 'T', after='NAXIS2'
  sxaddpar, header, 'CTYPE1', 'LINEAR'
  sxaddpar, header, 'CRVAL1', wavemin
  sxaddpar, header, 'CRPIX1', 1
  sxaddpar, header, 'CDELT1', dx
  sxaddpar, header, 'CDELT2', 1
  sxaddpar, header, 'CD1_1', dx
  sxaddpar, header, 'CD2_2', 1
  sxaddpar, header, 'LTM1_1', 1
  sxaddpar, header, 'LTM2_2', 1
  sxaddpar, header, 'DISPAXIS', 1
  sxaddpar, header, 'DC-FLAG', 0
  sxaddpar, header, 'WAT0_001', 'system=equispec'
  sxaddpar, header, 'WAT1_001', 'wtype=linear label=wavelength units=angstroms'
  sxaddpar, header, 'WAT2_001', 'wtype=linear'
 sxaddpar, header, 'BANDID1', 'spectrum - sky, weighted, IDL pipeline'
if keyword_set(multispec) then begin
 sxaddpar, header, 'BANDID2', 'spectrum - sky, unweighted, IDL pipeline'
 sxaddpar, header, 'BANDID3', 'sky spectrum, IDL pipeline'
 sxaddpar, header, 'BANDID4', 'sigma, IDL pipeline'
 sxaddpar, header, 'CD3_3', 1
 sxaddpar, header, 'LTM3_3', 1
  sxaddpar, header, 'WAT3_001', 'wtype=linear'
  sxaddpar, header, 'WCSDIM', 3
endif
;sxaddpar, header, 'IMAGEID', '1'
sxdelpar, header, 'IMAGEID'
sxdelpar, header, 'LTM1_1'
sxdelpar, header, 'LTM2_2'
sxdelpar, header, 'LTV1'
sxdelpar, header, 'LTV2'
sxdelpar, header, 'DTM1_1'
sxdelpar, header, 'DTM2_2'
sxdelpar, header, 'DTV1'
sxdelpar, header, 'DTV2'
sxdelpar, header, 'ATM1_1'
sxdelpar, header, 'ATM2_2'
sxdelpar, header, 'ATV1'
sxdelpar, header, 'ATV2'
sxdelpar, header, 'EXTNAME'
sxdelpar, header, 'EXTVER'
sxdelpar, header, 'CCDSEC'
sxdelpar, header, 'DETSIZE'
sxdelpar, header, 'DISPAXIS'
sxdelpar, header, 'DETSEC'
sxdelpar, header, 'DATASEC'
sxdelpar, header, 'XTENSION'
sxdelpar, header, 'INHERIT'
sxdelpar, header, 'CCDSUM'
blanks=where(strtrim(header,2) eq '', blank_cnt, complement=nonblank)
if blank_cnt ge 1 then   header=header[nonblank]
header_root=header

  for i = 0, nfib -1 do begin
     id=' 2 '
     if strtrim(strlowcase(plugmap[i].objtype),2) eq 'sky' then id=' 0 ' else if strtrim(strlowcase(plugmap[i].objtype),2) eq 'target' then id=' 1 '
     apnumstring=strn(i+1)+id+strn(i)+' '+strn(i+1)
     sxaddpar, header, 'APNUM'+strn(i+1), apnumstring 
  endfor
  for i = 0, nfib -1 do begin  
    id=' '+strn(fix(plugmap[i].icode))+' '
     if strtrim(strlowcase(plugmap[i].objtype),2) eq 'sky' then id=' -1 ' else if (strtrim(strlowcase(plugmap[i].objtype),2) ne 'target') and strtrim(strlowcase(plugmap[i].objtype),2) ne 'opptar' then id=' 0 '  
    apidstring=strtrim(strlowcase(plugmap[i].objtype),2)+' '+hs_dec2hms(plugmap[i].ra/15., /double)+' '+hs_dec2hms(plugmap[i].dec, /double)+id+string(fix(plugmap[i].rcode), format='(I3)') +' '+string(plugmap[i].xfocal, format='(f8.3)')+' '+string(plugmap[i].yfocal, format='(f8.3)')
    sxaddpar, header, 'APID'+strn(i+1), apidstring 
  endfor

if NOT keyword_set (multispec) then begin
  mwrfits, newflux,outname, header, /create
  if NOT keyword_set(simple)  then begin
    mwrfits, newivar,outname, header
    mwrfits, plugmap, outname
  endif
endif else begin
  multi=fltarr(npix, nfib, 4)
  multi[*,*,0]=newflux
  multi[*,*,1]=newflux_box
  multi[*,*,2]=newsky
  sigma=sqrt(1./newivar)
  bad=where(finite(sigma,/nan), badcnt)
  if badcnt ge 1 then sigma[bad]=abs(newflux[bad])
  bad=where(finite(sigma,/inf), badcnt)
  if badcnt ge 1 then sigma[bad]=abs(newflux[bad])
  multi[*,*,3]=sigma
  mwrfits, multi, outname, header, /create
  if keyword_set(hsplit) then begin
    if n_ext gt 5 then ntrim=15 else ntrim=8
    obj=strmid(outname,0, strlen(outname)-ntrim)
    if hsplit eq 2 then begin
       if n_ext gt 5 then if direxist('skysub_nonlin_'+obj) eq 0 then spawn, 'mkdir skysub_nonlin_'+obj 
       if n_ext le 5 then if direxist('1d.nonlin_'+obj) eq 0 then spawn, 'mkdir 1d.nonlin_'+obj 
    endif
    if hsplit eq 1 then begin 
       if n_ext gt 5 then if direxist('skysub_'+obj) eq 0 then spawn, 'mkdir skysub_'+obj
       if n_ext le 5 then if direxist('1d.'+obj) eq 0 then spawn, 'mkdir 1d.'+obj 
    endif
    hdr1d=header_root
    sxaddpar, hdr1d, 'NAXIS2', 1
    for i=0, 299 do begin   
       id=' 2 '
       if strtrim(strlowcase(plugmap[i].objtype),2) eq 'sky' then id=' 0 ' else if strtrim(strlowcase(plugmap[i].objtype),2) eq 'target' then id=' 1 '
       apnumstring=strn(1)+id+strn(0)+' '+strn(1)
       sxaddpar, hdr1d, 'RA', hs_dec2hms(plugmap[i].ra/15.,sep=':',/double)
       sxaddpar, hdr1d, 'DEC', hs_dec2hms(plugmap[i].dec,sep=':',/double)
       sxaddpar, hdr1d, 'APNUM1', apnumstring 
       id=' '+strn(fix(plugmap[i].icode))+' '
       if strtrim(strlowcase(plugmap[i].objtype),2) eq 'sky' then id=' -1 ' else if (strtrim(strlowcase(plugmap[i].objtype),2) ne 'target') and strtrim(strlowcase(plugmap[i].objtype),2) ne 'opptar' then id=' 0 '  
       apidstring=strtrim(strlowcase(plugmap[i].objtype),2)+' '+hs_dec2hms(plugmap[i].ra/15., /double)+' '+hs_dec2hms(plugmap[i].dec, /double)+id+string(fix(plugmap[i].rcode), format='(I3)') +' '+string(plugmap[i].xfocal, format='(f8.3)')+' '+string(plugmap[i].yfocal, format='(f8.3)')
       sxaddpar, hdr1d, 'APID'+strn(1), apidstring 
       if strtrim(strlowcase(plugmap[i].objtype),2) eq 'target' then begin
          if strtrim(strlowcase(plugmap[i].bcode),2) eq 'target' then newname=strtrim(sxpar(hdr1d,'CAT-ID'),2)+'_'+strn(fix(plugmap[i].icode)) else newname=strtrim(plugmap[i].bcode,2)
      endif else newname=strtrim(strlowcase(plugmap[i].objtype),2)
      sxaddpar, hdr1d, 'OBJECT', newname
      if hsplit eq 2 then begin
        sxaddpar, hdr1d, 'WAT0_001', 'system=multispec'
        sxaddpar, hdr1d, 'WAT1_001', 'wtype=multispec label=wavelength units=angstroms'
       sxaddpar,hdr1d, 'WAT3_001', 'wtype=linear'
        sxaddpar, hdr1d, 'CTYPE1', 'MULTISPE'
        sxaddpar, hdr1d, 'CTYPE2', 'MULTISPE'
        sxaddpar, hdr1d, 'CTYPE3', 'LINEAR'
       sxaddpar, hdr1d, 'CDELT1', 1.
        sxaddpar, hdr1d, 'CD1_1', 1.
       sxaddpar, hdr1d, 'WCSDIM', 3
;        sxdelpar, hdr1d, 'CDELT1'
;        sxdelpar, hdr1d, 'CD1_1'
        sxdelpar, hdr1d, 'CRVAL1'
        sxdelpar, hdr1d, 'CRPIX1'
        sxaddpar, hdr1d, 'LTM1_1', 1.0
        sxaddpar, hdr1d, 'LTM2_2', 1.0
        lam1d=lam[*,i]
        xx=findgen(n_elements(lam[*,i]))
        flux1d=flux[*,i]
        flux1d_box=flux_box[*,i]
        sky1d=sky[*,i]
        xy2traceset, xx, lam1d, airset, ncoeff=6
        npix_nl=n_elements(flux1d)
        sxaddpar, hdr1d, 'NAXIS2', npix_nl
        badpix=where(ivar[*,i] eq 0.0, nbad)
        sigma1d=sqrt(1./ivar[*,i])
;        if nbad ge 1 then for q=0, nbad-1 do sigma1d[badpix[q]]=(sigma1d[((badpix[q]-1)>0)]+sigma1d[((badpix[q]+1)<(npix_nl-1))])/2.
        if nbad ge 1 then for q=0, nbad-1 do sigma1d[badpix[q]]=abs(flux1d[badpix[q]])
        nonlin_cube=fltarr(npix_nl, 1, 4)
        nonlin_cube[*,0,0]=flux1d
        nonlin_cube[*,0,1]=flux1d_box
        nonlin_cube[*,0,2]=sky1d
        nonlin_cube[*,0,3]=sigma1d
        funcstr='"1 1 2 '+string(min(lam1d),format='(f10.5)')+' '+string(dx,format='(f9.7)')+' '+strn(npix_nl)+' 0. 0.5 1.5 1. 0. 2 6 1 '+strn(npix_nl)+' '+strjoin(string(airset.coeff, format='(f13.7)'))+'"'
;        funcstr='"1 1 2 1. '+string(dx,format='(f9.7)')+' '+strn(npix_nl)+' 0. 0.5 1.5 1. 0. 2 6 1 '+strn(npix_nl)+' '+strjoin(string(airset.coeff, format='(f13.7)'))+'"'
        lamStr='wtype=multispec spec1 = '+funcstr
;        lamStr=lamStr+' spec2 = '+funcstr
;        lamStr=lamStr+' spec3 = '+funcstr
;        lamStr=lamStr+' spec4 = '+funcstr
        for k=0, long(strlen(lamStr)/(80-12)) do begin
		name = string(format='("WAT2_", i03)', k+1)
		sxaddpar, hdr1d, name, strmid(lamStr, k*(80-12), 80-12)
	endfor
        if n_ext gt 5 then outfile1d='skysub_nonlin_'+obj+'/'+string((i+1),format='(I3.3)')+'.'+newname+'.ec.fits'
        if n_ext le 5 then outfile1d='1d.nonlin_'+obj+'/'+string((i+1),format='(I3.3)')+'.'+newname+'.ec.fits'
;        mwrfits, multi[*,i,*], outfile1d, hdr1d, /create
        mwrfits, nonlin_cube, outfile1d, hdr1d, /create
      endif else begin
        if n_ext gt 5 then outfile1d='skysub_'+obj+'/'+string((i+1),format='(I3.3)')+'.'+newname+'.ms.fits' else outfile1d='1d.'+obj+'/'+string((i+1),format='(I3.3)')+'.'+newname+'.ms.fits'
        mwrfits, multi[*,i,*], outfile1d, hdr1d, /create
      endelse
    endfor
  endif
endelse

END
