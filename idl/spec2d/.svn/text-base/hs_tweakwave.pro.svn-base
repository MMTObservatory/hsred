;;This program will be used to create a tweak to the wavelegnth solution
PRO hs_tweakwave, arcguess, rerun=rerun, root=root
  
  
  
 ;;Read the first ccd and extract the spectra
  
  arcname = 'calibration/' + rerun + '/' + 'arc-' + root + '.fits'
  
  ccdnum=1
  
   hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun
         
     
     sigma= 2.0
     proftype = 2.0
     highrej = 15
     lowrej = 15
     npoly = 1
     wfixed = [1,1,1]

     hs_proc, arcname, ccdnum, arcimg, arcivar, rerun=rerun, root=root
 
     hs_extract_image, arcimg, arcivar, xsol, sigma, flux, fluxivar, $
        proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
        npoly=npoly, relative=1
     
     
     ;;Now get the solution for the first chip first fiber
     
     
     stop = ''
     
     
     splog, 'WELCOME TO THE INTERACTIVE WAVELENGTH TWEAKER'
   
 
     
     test = fltarr(n_elements(flux[*,0]), 6)
     for i = 0, 5 do begin
        test(*,i) = flux(*,0)
     endfor
     
     
     mwrfits, test, 'tmp-wavetest.fits', /create
          
     
     while stop eq '' do begin
        
        splog, 'Please enter an estimate for the MAXIMUM lambda:'

        
        read, prompt='Min Wave : ', minwave
        splog, 'Please enter your guess at the mean dispersion (1.21 A):'
        read, prompt='Disp : ', disp
        
     
        iarcfit, 'tmp-wavetest.fits', /doplot, xcoeff=[minwave, -1*disp], $
           dxcoeff=[100,0.5], wset
        
        ent = ''
        splog, 'Does this solution seem like a good first guess?'
        splog, '1 - Yes'
        splog, '2 - No'
        read, prompt = 'Good? : ', ent
        
        if strmatch(ent,'1') then stop = 'stop'
        
     endwhile
     
     wset10 = wset
     
       
     test = fltarr(n_elements(flux[*,1]), 6)
     for i = 0, 5 do begin
        test(*,i) = flux(*,1)
     endfor
     
     
     mwrfits, test, 'tmp-wavetest.fits', /create
          
     
     splog, 'Now fitting CCD #1 Fiber #1'
     stop = ''
     
     while stop eq '' do begin
        
        
        splog, 'Please enter an estimate for the MAXIMUM lambda:'

        read, prompt='Min Wave : ', minwave
        splog, 'Please enter your guess at the mean dispersion (1.20 A):'
        read, prompt='Disp : ', disp
     
        iarcfit, 'tmp-wavetest.fits', /doplot, xcoeff=[minwave, -1*disp], $
           dxcoeff=[100,0.5], wset
        
        ent = ''
        splog, 'Does this solution seem like a good first guess?'
        splog, '1 - Yes'
        splog, '2 - No'
        read, prompt = 'Good? : ', ent
        
        if strmatch(ent,'1') then stop = 'stop'
        
     endwhile
     
     wset11 = wset
     
     
     
     ccdnum=2
  
   hs_readtraceset, ccdnum, xsol=xsol, rerun=rerun
         
     
     sigma= 2.0
     proftype = 3
     highrej = 15
     lowrej = 15
     npoly = 1
     wfixed = [1,1]

     hs_proc, arcname, ccdnum, arcimg, arcivar, rerun=rerun, root=root
 
     hs_extract_image, arcimg, arcivar, xsol, sigma, flux, fluxivar, $
        proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
        npoly=npoly, relative=1
     
     
     ;;Now get the solution for the first chip first fiber
     
     
     stop = ''
     
     
     splog, 'WELCOME TO THE INTERACTIVE WAVELENGTH TWEAKER'
   
 
     
     test = fltarr(n_elements(flux[*,0]), 6)
     for i = 0, 5 do begin
        test(*,i) = flux(*,0)
     endfor
     
     
     mwrfits, test, 'tmp-wavetest.fits', /create
          
     
     splog, 'Now fitting CCD#2 Fiber #0'
     
     while stop eq '' do begin
        
        splog, 'Please enter an estimate for the MAXIMUM lambda:'

        
        read, prompt='Min Wave : ', minwave
        splog, 'Please enter your guess at the mean dispersion (1.21 A):'
        read, prompt='Disp : ', disp
     
        iarcfit, 'tmp-wavetest.fits', /doplot, xcoeff=[minwave, -1*disp], $
           dxcoeff=[100,0.5], wset
        
        ent = ''
        splog, 'Does this solution seem like a good first guess?'
        splog, '1 - Yes'
        splog, '2 - No'
        read, prompt = 'Good? : ', ent
        
        if strmatch(ent,'1') then stop = 'stop'
        
     endwhile
     
     wset20 = wset
     
       
     test = fltarr(n_elements(flux[*,1]), 6)
     for i = 0, 5 do begin
        test(*,i) = flux(*,1)
     endfor
     
     
     mwrfits, test, 'tmp-wavetest.fits', /create
          
     
     splog, 'Now fitting CCD #2 Fiber #1'
     stop = ''
     
     while stop eq '' do begin
        
        
        splog, 'Please enter an estimate for the MAXIMUM lambda:'

        read, prompt='Min Wave : ', minwave
        splog, 'Please enter your guess at the mean dispersion (1.20 A):'
        read, prompt='Disp : ', disp
     
        iarcfit, 'tmp-wavetest.fits', /doplot, xcoeff=[minwave, -1*disp], $
           dxcoeff=[100,0.5], wset
        
        ent = ''
        splog, 'Does this solution seem like a good first guess?'
        splog, '1 - Yes'
        splog, '2 - No'
        read, prompt = 'Good? : ', ent
        
        if strmatch(ent,'1') then stop = 'stop'
        
     endwhile
     
     wset21 = wset
     
     
     arcguess = fltarr(4,4)
     arcguess(0,*) = (wset10.coeff)[0:3,3]
     arcguess(1,*) = (wset11.coeff)[0:3,3]
     arcguess(2,*) = (wset20.coeff)[0:3,3]
     arcguess(3,*) = (wset21.coeff)[0:3,3]
     
     
     
     
     
     
  END
  
  
  
        
