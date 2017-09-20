;;This program makes the traceset needed for a hectospec dataset and writes 
;;it out.

PRO hs_findarcdisper, arcfile=arcfile, $
                      write=write, xcen=xcen, $
                      tset=tset, xsol=xsol, rerun=rerun
  
  if NOT keyword_set(arcfile) then begin
;     file = findfile('*fit*')
;     N = n_elements(file)/2
;     header = headfits(file(N), exten=1,  /silent)
;     mjd = sxpar(header, 'MJD')
     arcfile = 'calibration/'+rerun+'/arc.fits'
  endif
  
  filelist = findfile('calibration/'+rerun+'/*fit*')
  k = where(strmatch(filelist, arcfile)) 
  if k(0) eq -1 then begin
     splog, 'This arcfile does not exist'
     return
  endif
  
  for ccdnum = 1, 2 do begin
     
     hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun
     sigma= 2.0
     proftype = 1.0
     highrej = 15
     lowrej = 15
     npoly = 1
     wfixed = [1,1]
     
     hs_proc, arcfile, ccdnum, arcim, arcivar, rerun=rerun
     
     extract_image, arcim, arcivar, xsol, sigma, flux, fluxivar, $
        proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
        npoly=npoly, relative=1
     
     if file_test('calibration/'+rerun+'/lineloc.fits') eq 0L then begin
        splog, 'The wavelength calibrator has not been run as'
        splog, 'lineloc.fits is not present'
        return
     endif
     
     xcen = mrdfits('calibration/'+rerun+'/lineloc.fits',ccdnum-1)
     n = n_elements(xcen(0,*))
     indn=n/3
     ind = indgen(indn)
     for i = 0, indn -1 do begin
        ind(i) = 3*i
     endfor
     
     xcen=xcen(*,ind+1)
     dispset = hs_fitdispersion(flux,fluxivar, xcen) 
     if keyword_set(write) then begin
        if ccdnum eq 1 then begin
           mwrfits, dispset, 'calibration/'+rerun+'/dispset.fits', /create
        endif
        
        if ccdnum eq 2 then begin
           mwrfits, dispset, 'calibration/'+rerun+'/dispset.fits'
        endif
     endif
  endfor
  
  return
  
END

  
  
  
