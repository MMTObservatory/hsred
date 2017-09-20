PRO hs_correctplugmap, filelist
   
   FOR ifile = 0, n_elements(filelist) -1 DO BEGIN
      
      plugmap = mrdfits(filelist(ifile), 5)
      
      
      standfile = getenv('HSRED_DIR')+'/etc/standstar.dat'
      readcol, standfile, sra, sdec, u, g, r, i, z, $
         reddening, format='D,D,D,D,D,D,D,D', $
         comment='#', /silent
      
      spherematch, plugmap.ra, plugmap.dec, sra, sdec, 1.0/3600., $
         match1, match2, dist
     
      reddening2 = reddening
      reddening0 = reddening2*1.874
      reddening1 = reddening2*1.379
      reddening3 = reddening2*0.758
      reddening4 = reddening2*0.538
      
      mag = dblarr(5,300)*0.0
      mag(0, match1) = u(match2) + reddening0(match2)
      mag(1, match1) = g(match2) + reddening1(match2)
      mag(2, match1) = r(match2) + reddening2(match2)
      mag(3, match1) = i(match2) + reddening3(match2)
      mag(4, match1) = z(match2) + reddening4(match2)
      
      plugmap.mag = mag
      
      rmfile, '/tmp/hsred_tmp.fits'
      FOR ii = 0, 4 DO BEGIN
         test = mrdfits(filelist(ifile), ii)
         mwrfits, test, '/tmp/hsred_tmp.fits'
      endfor
      
      mwrfits, plugmap, '/tmp/hsred_tmp.fits'
      
      spawn, 'mv /tmp/hsred_tmp.fits ' + filelist(ifile)
      
   ENDFOR
   
END

         
