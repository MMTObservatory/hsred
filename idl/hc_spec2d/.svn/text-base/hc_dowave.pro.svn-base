PRO hc_dowave, file, outroot=outroot, rerun=rerun, debug=debug, xrange=xrange
   
   IF NOT keyword_set(rerun) THEN rerun = '0000'
   IF size(rerun, /tname) NE 'STRING' THEN $
      rerun=string(rerun, format='(i4.4)')
   IF NOT keyword_set(outroot) THEN outroot = 'wset'
   outfile = 'calibration/' + rerun + '/' + outroot +'.fits'
   
   IF file_test(outfile) THEN rmfile, outfile
   
   IF file_test('calibration/' + rerun + '/traceset.fits') EQ 0l THEN $
      hs_maketraceset, /chelle, rerun=rerun
   
   rmfile, 'calibration/0000/wave.fits'
   
   FOR ccdnum = 1, 2 DO BEGIN
    
      hc_getspec, file, ccdnum, specflux=specflux, rerun=rerun
;      hc_initialwave, specflux, rerun=rerun
      
      IF ccdnum EQ 1 THEN idfile = 'reduction/0000/arc/database/id000'
      IF ccdnum EQ 2 THEN idfile = 'reduction/0000/arc/database/id120'
      
      readcol, idfile, pix, wave
      test = [[pix], [wave]]
      mwrfits, test, '/tmp/hsred_arcguess.fits', /create
      
      
      hc_initfit, specflux, oddcent=oddcent, evencent=evencent, $
         oddwave=oddwave, evenwave=evenwave, debug=debug, xrange=xrange
      
      delvarx, wset
      wset = hc_construct_wset(specflux, oddcent=oddcent, $
                               oddwave=oddwave, evencent=evencent, $
                               evenwave=evenwave)

      IF ccdnum EQ 1 THEN mwrfits, wset, outfile, /create
      IF ccdnum EQ 2 THEN mwrfits, wset, outfile
      
   ENDFOR
   
END

