PRO hc_getspec, file, ccdnum, docosmic=docosmic, header=objhdr, $
                specflux=specflux, specivar=specivar, $
                pixelmask=pixelmask, $
                rerun=rerun, nobias=nobias
   
   ;;This program will do the hecto-chelle extraction and return the 
   ;;extracted spectra.   I do this here rather than in the standard 
   ;;hc_extract as there are some problems with junk fibers that must be
   ;;discarded properly
   
   
   ;;First read the science frame -- note the /chelle is 
   ;;used to tell this to use the hectochelle bpms - they are
   ;;currently blank - this will LIKELY BREAK IF YOU BIN -- as the 
   ;;BPMS assume 1x1 binning.
   IF NOT keyword_set(rerun) THEN rerun = '0000'
  
   hs_proc, file, ccdnum, spec, specimivar, pixelmask=pixelmask, $
      rerun=rerun, header=objhdr, docosmic=docosmic, nobias=nobias, $
      /chelle
   
   ;;Setup the extraction parameters - check the hs_extract_image
   ;;docs to see what these mean
   
   sigma= 3.0
   proftype = 2.0
   highrej = 100
   lowrej = 100
   npoly = 1
   wfixed = [1,0,1]
   
   hs_readtraceset, ccdnum , xsol=xsol, rerun=rerun
   
   hs_extract_image, spec, specimivar, xsol, sigma,$
      specflux, specivar, $
      proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
      npoly=npoly, relative=1, ymodel=ymodel
   
   mwrfits, spec, 'test.fits', /create
   mwrfits, ymodel, 'test.fits'
   
   
   
   ;;Now chuck the spectra that are either bad or are not fibers that are 
   ;;completely matched on the mapfile
   IF ccdnum EQ 1 THEN BEGIN
      specflux = specflux(*,3:122)
      specivar = specivar(*,3:122)
   ENDIF
   IF ccdnum EQ 2 THEN BEGIN
      specflux = specflux(*,0:119)
      specivar = specivar(*,0:119)
   ENDIF
   
   
END

