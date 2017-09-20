Pro hs_compare, filelist
  
  ;;This program will take an input list of spZbest files and compare the 
  ;;the results to the OSU redshifts.
  nfile = n_elements(filelist)
  readcol, '/home/rcool/catalog.cat', ra, dec, type, rank, mag, apmag, $
     ncode, icode, rcode, etime, redshift, inversnr, $
     format='D,D,A,D,D,D,D,D,D,D,D', /silent
  
  readcol, '/home/rcool/catalog.zOSU', junk, junk, zosu, junk, junk, junk, $
     format='D,D,D,D,D,D', /silent
  
  
  for ifile = 0, nfile -1 do begin
     
     spec = mrdfits(filelist(ifile), 1)
     
     spherematch, ra, dec, spec.plug_ra/15., spec.plug_dec, 1.0/3600., $
        match1, match2, dist
     
     spec = spec[match2]
     
     if n_elements(output) ne 0 then begin
        newoutput = create_struct(output[0])
       
        output = make_array(val=newoutput, dim=n_elements(match2))
     endif
     
    if n_elements(output) eq 0 then $
        create_struct, output, 'output', ['file', 'ra', 'dec', 'zBs', 'zosu', $
                                       'mag', 'apmag', 'aperture', 'zwarning',$
                                       'rchi2diff', 'class'], $
        'A,D,D,D,D,D,D,D,D,D,A', dimen = n_elements(match2)
     
     file = strarr(n_elements(match2))
     file(*) = filelist(ifile)
     output.file = file 
     output.ra = spec.plug_ra
     output.dec = spec.plug_dec
     output.zbs = spec.z
     output.zosu = zosu(match1)
     output.mag = mag(match1)
     output.apmag = apmag(match1)
     output.aperture = spec.fiberid
     output.zwarning = spec.zwarning
     output.rchi2diff = spec.rchi2diff
     output.class = spec.class
     
     if n_elements(output1) ne 0 then begin
        outputtemp = create_struct(output[0])
        output2 = make_array(val=outputtemp, $
                             dim = (n_elements(output)+n_elements(output1)))
        output2[0:n_elements(output)-1] = output
        output2[n_elements(output):n_elements(output2)-1] = output1
        
        output1 = output2
     endif
     
     if n_elements(output1) eq 0 then output1 = output
     
     
  endfor
  
  
  mwrfits, output1, '/home/rcool/check.fits', /create
  
END


    
