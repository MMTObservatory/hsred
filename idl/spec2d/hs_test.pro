PRO hs_test.pro
  
  hs_proc, 'dflat-53110.fits', 1, flatim, flativar
  
   xsol = trace150crude(flatim, flativar, yset=ycen)
  xy2traceset, ycen, xsol, tset, ncoeff=5
  traceset2xy, tset, ycen, xsol
  
  hs_proc, 'arc-53110.fits', 1, arcimg, arcivar
 
  
  
  sigma= 2.0
  proftype = 3
  highrej = 15
  lowrej = 15
  npoly = 1   
  wfixed = [1,1]
             
  extract_image, arcimg, arcivar, xsol, sigma, flux, fluxivar, $
     proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
     npoly=npoly, relative=1
  

                     
  
     hs_wavecal, flux, fluxivar, xpeak, $
        ypeak, wset, ncoeff=arccoeff, $
        bestcorr=corr, row=37, parity='odd', ccdnum=1
     


      hs_proc, 'kfield1_1.0617.fits', 1, spec, specivar
      
      extract_image, spec, specivar, xsol, sigma, specflux, specfluxivar, $
         proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
                                npoly=npoly, relative=1
      
      
      
      
      
     
     dims = size(flatim, /dimens)
nx = dims[0]
ny = dims[1]

waveimg = fltarr(nx,ny)
traceset2xy, wset, xx, lam
 plot, 10^loglam[*,10]

xy2traceset, transpose(xsol), transpose(loglam), tmpset, $
   func='legendre', ncoeff=5, xmin=0, xmax=nx-1, upper=2.0, lower=2.0

xtmp=0
      traceset2xy, tmpset, xtemp, waveimg
djs_plot, 10^(loglam[*,100])/10, flux[*,100]

lampdefault = filepath('lamphenear.dat', $
                      root_dir=getenv('IDLSPEC2D_DIR'), $
                      subdirectory='etc')
      lampfilename = lampdefault
      
      splog, 'Reading lamp file ', lampfilename
   readcol, lampfilename, lampwave, lampinten, lampquality, format='D,F,A'
   lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0.0d0}
   lamps = replicate(lamps, N_elements(lampwave))
   lamps.lambda = lampwave
   lamps.loglam = alog10(lampwave)
   lamps.intensity = lampinten
   lamps.good = strupcase(lampquality) EQ 'GOOD' AND lampinten GT 0     
   djs_oplot, lamps.lambda/10, lamps.intensity/max(lamps.intensity)*$
      max(flux[*,100]), color='red'
  
      
END


